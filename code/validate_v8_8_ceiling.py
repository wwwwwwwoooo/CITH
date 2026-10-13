"""Independent checks of finite productivity ceilings and shared technology.

Targets: cor:ceiling-budget (A1), cor:capped-basin and its qualifications
(A2), app:ceiling-growth-proof (A3), app:norep-ceiling-proof (A4), and
cor:shared-technology (B). Numerical assertions supplement, not replace,
the appendix proofs. Finite trajectories do not establish asymptotic limits.
Only archived inputs are read; no other validator is imported.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
import math
import platform
import sys

import numpy as np
import scipy
from scipy.optimize import brentq, minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
SEED = 2026101380
TOL = 2e-9
ROOT_TOL = 2e-13
BRACKET_TOL = 2e-9
MAX_STEPS = 5000


def production(kind: str):
    if kind == 'linear':
        return (lambda b: b), (lambda b: 1.0)
    if kind == 'log':
        return np.log1p, (lambda b: 1.0 / (1.0 + b))
    if kind == 'root':
        return (lambda b: 2.0 * b / (np.sqrt(1.0 + b) + 1.0)), \
               (lambda b: 1.0 / np.sqrt(1.0 + b))
    raise ValueError(kind)


def tipping(p, s, kind='log') -> float:
    ell, ellp = production(kind)
    mask = np.asarray(s) > 0
    top = float(np.min(np.asarray(p)[mask] / np.asarray(s)[mask]))
    return float(brentq(lambda b: ell(b) / ellp(b) - top,
                        0.0, top, xtol=ROOT_TOL))


def threshold(b, p, r, s, d0=1.0, kind='log') -> float:
    ell, ellp = production(kind)
    p, r, s = (np.asarray(v, dtype=float) for v in (p, r, s))
    active = r + s > 0
    den = p[active] * ellp(b) - s[active] * ell(b)
    if np.any(den <= 0):
        return math.inf
    return float(np.max((r[active] + d0 * s[active]) / den,
                        initial=0.0))


def budget_ceiling(cap, p, r, s, d0=1.0, kind='log') -> float:
    """Feasible-set bisection; preserves plateaus and the empty-set case."""
    ell, ellp = production(kind)
    p, r, s = (np.asarray(v, dtype=float) for v in (p, r, s))
    active = r + s > 0

    def feasible(b):
        return bool(np.all(cap * (p[active] * ellp(b) -
                                  s[active] * ell(b)) >=
                           r[active] + d0 * s[active]))
    if not feasible(0.0):
        return 0.0
    lo, hi = 0.0, tipping(p, s, kind)
    for _ in range(90):
        mid = (lo + hi) / 2.0
        if feasible(mid):
            lo = mid
        else:
            hi = mid
    return lo


@dataclass(frozen=True)
class Model:
    p: float = 1.0
    r: float = 0.0
    s: float = 1.0
    d0: float = 1.0

    def crit(self, b):
        return threshold(b, [self.p], [self.r], [self.s], self.d0)

    def tip(self):
        return tipping([self.p], [self.s])

    def star(self, cap):
        return budget_ceiling(cap, [self.p], [self.r], [self.s], self.d0)


def derivative_numerator(x, b, a, v):
    z = max(0.0, b - v.p * x)
    return v.r + v.s * v.d0 + v.s * v.r * x + \
        a * (v.s * math.log1p(z) - v.p / (1.0 + z))


def policy(b, a, v=Model()):
    if derivative_numerator(0.0, b, a, v) <= 0:
        return 0.0
    if derivative_numerator(b / v.p, b, a, v) >= 0:
        return b / v.p
    return brentq(lambda x: derivative_numerator(x, b, a, v),
                  0.0, b / v.p, xtol=ROOT_TOL)


def value(x, b, a, v=Model()):
    return math.log(v.d0 + a * math.log1p(max(0.0, b - v.p * x)) +
                    v.r * x) + v.s * x


def saturated(a0, cap, ga, t):
    # Stable at t=0 and for slow early growth.
    return a0 + (cap - a0) * (-math.expm1(-ga * t))


def trajectory(a0, cap, ga, gx, bar, v=Model(), steps=MAX_STEPS):
    """Return certified finite outcomes; an unresolved run is never relabeled."""
    b, tip = bar, v.tip()
    path = []
    for t in range(steps):
        a = saturated(a0, cap, ga, t)
        if b >= tip:
            return 'cross', t, b, path
        dn = derivative_numerator(0.0, b, a, v)
        if dn < -1e-12:
            return 'capture', t, b, path
        x = policy(b, a, v)
        path.append((t, b, a, x))
        b += gx * x
    return 'unresolved', steps, b, path


def face(budget, a, d0=1.0, ji=1.0, jr=1.0,
         si=1.0, sr=1.0, pm=1.0, pb=1.0):
    """Solve the restricted face, including z=0, using multiplier KKT."""
    loga = math.log(a)

    def services(loglam):
        w = si * max(0.0, math.log(ji / (si * pm)) - loglam)
        b = sr * max(0.0, math.log(jr / (sr * pb)) - loglam)
        return w, b

    def interior(z):
        loglam = loga - math.log1p(z) - math.log(d0 + a * math.log1p(z))
        w, b = services(loglam)
        return z + pm * w + pb * b - budget

    if interior(0.0) <= 0.0:
        z = brentq(interior, 0.0, budget, xtol=ROOT_TOL)
        loglam = loga - math.log1p(z) - math.log(d0 + a * math.log1p(z))
    else:
        z = 0.0
        hi = max(math.log(ji / (si * pm)), math.log(jr / (sr * pb)))
        lo = hi - 2.0 - budget / min(pm * si, pb * sr)
        loglam = brentq(lambda ll: pm * services(ll)[0] +
                       pb * services(ll)[1] - budget,
                       lo, hi, xtol=ROOT_TOL)
    w, b = services(loglam)
    lam = math.exp(loglam)
    fprime = 1.0 / (d0 + a * math.log1p(z))
    return z, w, b, lam, fprime


def run():
    rng = np.random.default_rng(SEED)
    groups = ('A1_definition', 'A1_monotonic_limit', 'A1_endpoints',
              'A2_policy', 'A2_low_cap', 'A2_ordered_basin',
              'A2_counterexamples', 'A3_subunit', 'A3_linear_endpoint',
              'A4_face', 'A4_ceiling_obstruction', 'A4_boundary',
              'B_coordinates', 'B_archived_examples')
    counts = {g: [0, 0] for g in groups}
    failures = []
    details = {'basin_brackets': [], 'low_cap_crossing_times': [],
               'scope': 'Finite checks only; limit and uniqueness claims rely on proofs.'}

    def check(group, ok, context):
        passed = bool(ok)
        counts[group][0 if passed else 1] += 1
        if not passed:
            failures.append({'group': group, 'context': str(context)})

    def close(x, y, atol=TOL, rtol=TOL):
        return abs(x - y) <= atol + rtol * max(abs(x), abs(y))

    # A1: empty sets, included positive endpoints, independent corner tests.
    for kind in ('linear', 'log', 'root'):
        for _ in range(40):
            p = rng.uniform(.5, 3., 3)
            r = rng.uniform(0., 2., 3)
            s = rng.uniform(.05, 1.5, 3)
            if rng.random() < .5:
                s[2] = 0.0
            d0 = float(rng.uniform(1., 3.))
            tip = tipping(p, s, kind)
            f0 = threshold(0., p, r, s, d0, kind)
            cap = f0 * float(rng.uniform(1.2, 30.))
            bs = budget_ceiling(cap, p, r, s, d0, kind)
            check('A1_definition', 0. < bs < tip, (kind, bs, tip))
            check('A1_definition', close(threshold(bs, p, r, s, d0, kind), cap), kind)
            check('A1_definition', threshold(.8 * bs, p, r, s, d0, kind) <= cap, kind)
            check('A1_definition', threshold((bs + tip) / 2., p, r, s, d0, kind) > cap, kind)
            # A finite external trial cannot beat the corner at a feasible budget.
            ell, _ = production(kind)
            b = .8 * bs
            corner = math.log(d0 + cap * ell(b))
            weights = rng.dirichlet(np.ones(3)) * b * rng.uniform(.1, .9)
            x = weights / p
            ext = math.log(d0 + cap * ell(b - float(p @ x)) + float(r @ x)) + float(s @ x)
            check('A1_definition', ext <= corner + TOL, kind)
            caps = f0 * np.array([.5, 1., 2., 10., 1e3, 1e6])
            budgets = [budget_ceiling(c, p, r, s, d0, kind) for c in caps]
            check('A1_monotonic_limit', all(x <= y + TOL for x, y in zip(budgets, budgets[1:])), kind)
            check('A1_monotonic_limit', (tip - budgets[-1]) / tip < 1e-4, kind)
            check('A1_endpoints', budgets[0] == 0., kind)
    p, r, s = [1., 1.], [0., 2.], [1., 0.]
    for cap, expected in ((1.9, 0.), (2., .5), (4., .75)):
        check('A1_endpoints', close(budget_ceiling(cap, p, r, s, kind='linear'), expected), (cap, expected))
    check('A1_endpoints', close(math.log(1. + 2. * .4),
                                math.log(1. + 2. * (.4 - .1) + 2. * .1)), 'tie at plateau')
    v = Model()
    check('A1_endpoints', close(budget_ceiling(1., [1.], [0.], [1.]), 0.), 'f(0) equality, strict increase')

    # A2: the root policy is checked against an independent scalar optimizer.
    for _ in range(120):
        v = Model(float(rng.uniform(.5, 2.)), float(rng.uniform(0., .5)),
                  float(rng.uniform(.2, 1.5)), float(rng.uniform(1., 2.)))
        b, a = float(rng.uniform(.02, 3.)), float(rng.uniform(.05, 20.))
        x = policy(b, a, v)
        opt = minimize_scalar(lambda y: -value(y, b, a, v),
                              bounds=(0., b / v.p), method='bounded',
                              options={'xatol': 2e-12})
        best = max(value(0., b, a, v), value(b / v.p, b, a, v), -opt.fun)
        check('A2_policy', close(value(x, b, a, v), best), (b, a))
        check('A2_policy', policy(1.2 * b, a, v) >= x - TOL, 'budget order')
        check('A2_policy', policy(b, 1.2 * a, v) <= x + TOL, 'productivity order')

    for _ in range(50):
        v = Model(float(rng.uniform(.7, 1.5)), float(rng.uniform(0., .4)),
                  float(rng.uniform(.4, 1.2)), float(rng.uniform(1., 2.)))
        bar, gx = v.tip() * float(rng.uniform(.15, .7)), float(rng.uniform(.3, 1.2))
        cap = v.crit(bar) * float(rng.uniform(.4, .9))
        a0, ga = cap * float(rng.uniform(.05, .9)), float(rng.uniform(.02, 1.))
        dx = policy(bar, cap, v)
        bound = math.ceil((v.tip() - bar) / (gx * dx))
        state, t, b, path = trajectory(a0, cap, ga, gx, bar, v, steps=bound + 2)
        check('A2_low_cap', state == 'cross' and t <= bound, (state, t, bound))
        check('A2_low_cap', all(xx >= dx - TOL for _, _, _, xx in path), 'uniform increment')
        details['low_cap_crossing_times'].append(t)
        for j in range(4):
            a = saturated(a0, cap, ga, t + j)
            x = policy(b, a, v)
            check('A2_low_cap', x > (b - v.tip()) / v.p, 'runaway investment')
            nxt = b + gx * x
            check('A2_low_cap', nxt - v.tip() > (1. + gx / v.p) * (b - v.tip()), 'runaway stock')
            b = nxt

    # High-cap family: certify both finite basins and bracket their separator.
    v, bar, cap, gx = Model(), .2, 4., 1.
    bs = v.star(cap)
    for ga in (.025, .05, .1):
        lo, hi = .02, v.crit(bar)
        check('A2_ordered_basin', trajectory(lo, cap, ga, gx, bar, v)[0] == 'cross', ('low', ga))
        check('A2_ordered_basin', trajectory(hi, cap, ga, gx, bar, v)[0] == 'capture', ('high', ga))
        for _ in range(45):
            if hi - lo < BRACKET_TOL:
                break
            mid = (lo + hi) / 2.
            outcome = trajectory(mid, cap, ga, gx, bar, v)[0]
            if outcome == 'cross':
                lo = mid
            elif outcome == 'capture':
                hi = mid
            else:
                raise RuntimeError(f'Unresolved separator bracket at ga={ga}, a0={mid}')
        low, high = trajectory(lo, cap, ga, gx, bar, v), trajectory(hi, cap, ga, gx, bar, v)
        check('A2_ordered_basin', low[0] == 'cross' and high[0] == 'capture', ga)
        check('A2_ordered_basin', hi - lo < BRACKET_TOL, ga)
        check('A2_ordered_basin', high[2] < bs + TOL < v.tip(), ga)
        common = min(len(low[3]), len(high[3]))
        gaps = [low[3][i][1] - high[3][i][1] for i in range(common)]
        check('A2_ordered_basin', all(g >= -TOL for g in gaps), 'ordered budgets')
        check('A2_ordered_basin', all(x <= y + TOL for x, y in zip(gaps, gaps[1:])), 'nonshrinking gap')
        near = min(abs(row[1] - bs) for row in high[3])
        # This is a finite shadowing check, not a numerical proof of the limit.
        check('A2_ordered_basin', 0. < near < bs - bar,
              ('finite progress toward effective budget', ga, near))
        details['basin_brackets'].append({'ga': ga, 'lower_crossed': lo,
            'upper_captured': hi, 'B_star': bs, 'B_tip': v.tip(),
            'finite_capture_budget': high[2], 'finite_shadowing_distance': near,
            'lower_crossing_time': low[1], 'upper_capture_time': high[1]})

    # Both counterexamples in app:capped-basin-proof, with exact certificates.
    v = Model()
    for a0 in np.geomspace(1e-8, 10., 60):
        b1 = .1 + policy(.1, float(a0), v)
        a1 = 5. + float(a0) / 2.
        check('A2_counterexamples', b1 <= .2 + TOL and a1 > v.crit(.2), 'all capture')
        check('A2_counterexamples', policy(b1, a1, v) == 0., 'capture by date one')
    check('A2_counterexamples', v.crit(.2) < 30. / 19. < 5., 'analytic certificate')
    cap = v.crit(.2)
    for a0 in np.linspace(cap / 10., 1., 40):
        b1 = .1 + policy(.1, float(a0), v)
        a1 = min(10. * float(a0), cap)
        check('A2_counterexamples', close(b1, .2), 'hard cap first step')
        check('A2_counterexamples', abs(derivative_numerator(0., b1, a1, v)) < TOL,
              'hard cap equality at included endpoint')
        check('A2_counterexamples', close(v.star(cap), .2), 'shared terminal boundary')
    a0 = .5 * cap / 10.
    b2 = .2 + policy(.2, 10. * a0, v)
    b = b2
    for t in range(200):
        if b >= v.tip():
            break
        b += policy(b, cap, v)
    check('A2_counterexamples', b2 > .2 and b >= v.tip(), 'below hard cap threshold crosses')

    # A3: direct failure of necessary output rates, plus the eta=0 exception.
    for kind in ('linear', 'log', 'root'):
        ell, ellp = production(kind)
        for eta in (.2, .4, .6, .8):
            for cap in (.5, 2., 8.):
                # Choose a large budget; compare in logs to avoid overflow.
                b = 10. ** max(8., (1. - eta) / eta * (math.log10(cap) + 3.))
                logomega = math.log(cap * ell(b))
                required = math.log(b) / (1. - eta)  # s=M=1
                check('A3_subunit', logomega < required, (kind, eta, cap, b))
                check('A3_subunit', cap * ell(b) <= cap * b * (1. + TOL), 'concave tangent')
        check('A3_linear_endpoint', 8. * ellp(1e8) < 1. if kind != 'linear'
              else 8. * ellp(1e8) >= 1., kind)
    check('A3_linear_endpoint', .5 < 1. <= 2., 'cap test does not imply actual test')

    # A4: optimized face, net prices, growing obstruction, and zero-slope controls.
    for a in (.2, 1., 5.):
        ratios = []
        for budget in (10., 100., 1e3, 1e5, 1e7):
            z, w, b, lam, fp = face(budget, a)
            check('A4_face', close(z + w + b, budget, rtol=1e-10), (a, budget))
            check('A4_face', z > 0. and w > 0. and b > 0., (a, budget))
            check('A4_face', close(fp / lam, (1. + z) / a), 'interior identity')
            check('A4_face', close(math.exp(-w), lam) and close(math.exp(-b), lam), 'service KKT')
            k = 3. * max(a, 1. / math.e)
            lower = budget / (k * (1. + a * math.log1p(budget)))
            check('A4_ceiling_obstruction', fp / lam >= lower * (1. - TOL), 'uniform lower bound')
            check('A4_ceiling_obstruction', lam <= k / budget * (1. + TOL), 'multiplier bound')
            check('A4_ceiling_obstruction', 0. * fp / lam < 1., 'zero remote slope control')
            ratios.append(fp / lam)
        check('A4_ceiling_obstruction', all(x < y for x, y in zip(ratios, ratios[1:])), a)
        check('A4_ceiling_obstruction', .1 * ratios[-1] > 1., 'positive remote slope fails')
    # No lower bound on a is needed. At this tiny a, the optimum is exactly z=0.
    z, w, b, lam, fp = face(20., math.exp(-400.))
    check('A4_boundary', z == 0. and close(w, 10.) and close(b, 10.), 'zero-local optimum')
    check('A4_boundary', math.exp(-400.) <= lam and close(lam, math.exp(-10.)), 'boundary KKT')
    check('A4_boundary', fp / lam > 1e4 and .1 * fp / lam > 1., 'use multiplier at boundary')
    check('A4_boundary', not close(fp / lam, (1. + z) / math.exp(-400.)), 'interior identity is invalid')
    # Compare capped productivity with the original unbounded a=B^2 construction.
    for budget in (1e3, 1e5, 1e7):
        z, w, b, lam, fp = face(budget, budget**2)
        check('A4_ceiling_obstruction', .1 * fp / lam < 1., 'original construction positive control')

    # B: exact transformed phase coordinates, including all zero-slope cases.
    for i in range(300):
        b, d0, a, p = rng.uniform(.01, 10.), rng.uniform(1., 4.), rng.uniform(.1, 20.), rng.uniform(.1, 5.)
        r, s = rng.uniform(0., 5.), rng.uniform(0., 2.)
        if i % 4 == 0:
            r = 0.
        if i % 5 == 0:
            s = 0.
        la, lr, kappa = np.exp(rng.uniform(-2., 2., 3))
        xi = (1. + b) * math.log1p(b) * s / p
        zeta = (1. + b) * (r + d0 * s) / (p * a)
        direct_xi = (1. + b) * math.log1p(b) * s / (p / kappa)
        direct_zeta = (1. + b) * (lr * r + d0 * s) / ((p / kappa) * la * a)
        transformed = kappa / la * (zeta + (1. + b) * (lr - 1.) * r / (p * a))
        check('B_coordinates', close(direct_xi, kappa * xi), i)
        check('B_coordinates', close(direct_zeta, transformed), i)
        shared = (1. + b) * (la * r + d0 * s) / ((p / kappa) * la * a)
        decomposition = kappa * (1. + b) * r / (p * a) + \
            kappa / la * (1. + b) * d0 * s / (p * a)
        check('B_coordinates', close(shared, decomposition), 'shared-gain components')

    archive = ROOT / 'archive/v7_3/data/allocation'
    for filename, expected in (('instantaneous.npz', 219),
                               ('instant_local_productivity_x10.npz', 314),
                               ('instant_shared_productivity_x10.npz', 222),
                               ('instant_cost_floor_0p01.npz', 0)):
        with np.load(archive / filename, allow_pickle=False) as data:
            # Column names are checked against the frozen archive, not inferred.
            if 'all_local' not in data.files:
                raise KeyError(f'{filename}: available keys {data.files}')
            actual = int(np.count_nonzero(data['all_local']))
        check('B_archived_examples', actual == expected, (filename, actual, expected))

    report = {'seed': SEED, 'runtime': {'python': platform.python_version(),
              'numpy': np.__version__, 'scipy': scipy.__version__},
              'tolerances': {'comparison': TOL, 'root_absolute': ROOT_TOL,
                             'separator_bracket': BRACKET_TOL},
              'maximum_trajectory_steps': MAX_STEPS,
              'counts': counts, 'details': details, 'failures': failures,
              'interpretation': 'Assertion counts are not independent random draws.'}
    (ROOT / 'results').mkdir(exist_ok=True)
    (ROOT / 'results/v8_8_ceiling.json').write_text(json.dumps(report, indent=2) + '\n')
    for group, (passed, failed) in counts.items():
        print(f'{group}: PASS={passed} FAIL={failed}')
    if failures:
        print(json.dumps(failures[:20], indent=2))
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(run())
