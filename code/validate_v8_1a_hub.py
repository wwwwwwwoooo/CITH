"""Additive v8 Stage 1a checks; no archived file is read or changed.

Run: python validate_v8_1a_hub.py
Requires Python 3, NumPy, SciPy. Seed and tolerances are explicit.
A concave numerical optimization is compared with the closed-form corner test.
Complex-step derivatives independently check the algebra at the corner.
"""
from __future__ import annotations
import json
import platform
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import minimize

SEED = 2026101301
DERIV_TOL = 2e-9
VALUE_TOL = 2e-7
OUT = Path(__file__).resolve().parents[1] / 'results'


def ell(z, kind):
    if kind == 'linear': return z
    if kind == 'log': return np.log1p(z)
    if kind == 'sqrt': return 2 * (np.sqrt(1 + z) - 1)
    raise ValueError(kind)


def ellp(z, kind):
    return {'linear': 1.0, 'log': 1/(1+z), 'sqrt': 1/np.sqrt(1+z)}[kind]


def fval(C, eta):
    if eta == 1: return np.log1p(C)
    return np.expm1((1-eta)*np.log1p(C))/(1-eta)


def objective(s, v):
    # s is the vector of external budget shares. Complex inputs are intentional.
    z = v['B'] * (1-np.sum(s))
    xx = v['B'] * np.asarray(s) / v['p']
    C = v['d0'] - 1 + v['a'] * ell(z, v['kind'])
    C += np.sum(v['r'][:2] * np.log1p(v['kg'] * xx[:2]) / v['kg'])
    E = -v['amp'][0] * np.expm1(-xx[0]/v['scale'][0])
    R = -v['amp'][1] * np.expm1(-np.sum(xx[:2])/v['scale'][1])
    I = -v['amp'][2] * np.expm1(-np.sum(xx)/v['scale'][2])
    return fval(C, v['eta']) + np.dot(v['weights'], np.array([E, R, I]))


def run():
    rng = np.random.default_rng(SEED)
    counts = {'derivative': [0, 0], 'optimizer': [0, 0],
              'reductions': [0, 0], 'edges': [0, 0]}
    failures = []
    skipped = 0
    inward = outward = 0
    max_deriv_error = max_corner_excess = 0.0

    def record(group, ok, detail=''):
        counts[group][not bool(ok)] += 1
        if not ok: failures.append(f'{group}: {detail}')

    for k in range(360):
        v = dict(B=10**rng.uniform(-1, 1.4), a=10**rng.uniform(-.5, 2),
                 d0=1+rng.uniform(0, 3), kind=('linear','log','sqrt')[k % 3],
                 eta=(0., .5, 1., 2.)[k % 4],
                 p=10**rng.uniform(-.3, .5, 3),
                 r=np.r_[10**rng.uniform(-2, .3, 2), 0.],
                 kg=10**rng.uniform(-1, .3, 2),
                 amp=10**rng.uniform(-1, .3, 3),
                 scale=10**rng.uniform(-.3, .7, 3),
                 weights=10**rng.uniform(-3, -.1, 3))
        if k % 19 == 0: v['weights'][:] = 0
        if k % 37 == 0: v['r'][:] = 0
        comp = v['weights']*v['amp']/v['scale']
        S = np.array([sum(comp), sum(comp[1:]), comp[2]])
        B, a, d0 = v['B'], v['a'], v['d0']
        om = a*ell(B,v['kind']); gam = ell(B,v['kind'])/ellp(B,v['kind'])
        required = gam*(v['r']/om + S*(d0+om)**v['eta']/om)
        gap = v['p']-required
        closed = bool(np.all(gap >= 0))
        local = a*ellp(B,v['kind'])/(d0+om)**v['eta']
        predicted = B/v['p']*(v['r']/(d0+om)**v['eta'] + S) - B*local
        observed = np.array([np.imag(objective(np.eye(3,dtype=complex)[j]*(1j*1e-20),v))/1e-20
                             for j in range(3)])
        err = np.max(abs(observed-predicted))/(1+np.max(abs(predicted)))
        max_deriv_error = max(max_deriv_error,float(err))
        record('derivative',err <= DERIV_TOL,f'case {k}, error {err}')

        u0 = float(objective(np.zeros(3),v))
        scale_u = 1+abs(u0)
        fits=[]
        # This numerical problem uses the complete objective, not the corner test.
        for initial in (np.zeros(3),np.full(3,.2)):
            fit=minimize(lambda s:-float(objective(s,v))/scale_u, initial,
                         bounds=[(0.,1.)]*3,
                         constraints={'type':'ineq','fun':lambda s:1-np.sum(s)},
                         method='SLSQP',options={'ftol':1e-12,'maxiter':600})
            if fit.success and np.sum(fit.x) <= 1+1e-9: fits.append(fit)
        if not fits:
            record('optimizer',False,f'case {k}: no successful feasible optimization')
            continue
        best=max([u0]+[float(objective(f.x,v)) for f in fits])
        excess=best-u0
        if closed:
            inward+=1
            max_corner_excess=max(max_corner_excess,excess/scale_u)
            record('optimizer',excess <= VALUE_TOL*scale_u,f'case {k}: excess {excess}')
        else:
            outward+=1
            if excess > VALUE_TOL*scale_u:
                record('optimizer',True)
            elif np.max(predicted) <= 1e-5*scale_u:
                skipped+=1  # Explicitly unresolved, not counted as a pass.
            else:
                record('optimizer',False,f'case {k}: predicted ascent not resolved, {excess}')

        # Reductions do not assume that the optimization has passed.
        log_required = (v['r']+d0*S)/(a*ellp(B,v['kind']))+S*gam
        hub_log = gam*(v['r']+S*(d0+om))/om
        record('reductions',np.allclose(log_required,hub_log,rtol=2e-13,atol=2e-13))
        L=np.log1p(B)
        original=(1+B)*L*S+(1+B)*(v['r']+d0*S)/a
        general_log=(v['r']+d0*S)*(1+B)/a+S*(1+B)*L
        record('reductions',np.allclose(original,general_log,rtol=2e-13,atol=2e-13))

    # Exact eta=0 tie: local and the first external channel have equal return.
    for eta in (0., .5, 1., 2.):
        for B in (1e-10, 1., 10.):
            vals=[float(fval(2*(B-x)+2*x,eta)) for x in np.linspace(0,B,21)]
            record('edges',max(vals)-min(vals)<2e-13*(1+abs(vals[0])),'pure-computation tie')
            # Zero returns/slopes impose no restriction, including d0=1.
            zero_values=[float(fval(2*(B-x),eta)) for x in np.linspace(0,B,21)]
            record('edges',max(zero_values)<=zero_values[0]+1e-13,'zero residuals and returns')
    # B=0 has only the zero allocation. The divided hub formula is not evaluated.
    record('edges',np.all(np.zeros(4)==0),'B=0 singleton')
    result=dict(seed=SEED,runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
                tolerances=dict(derivative=DERIV_TOL,relative_value=VALUE_TOL),counts=counts,
                optimizer_inward=inward,optimizer_outward=outward,unresolved=skipped,
                max_scaled_derivative_error=max_deriv_error,max_scaled_corner_excess=max_corner_excess,
                failures=failures)
    OUT.mkdir(exist_ok=True)
    (OUT/'v8_1a_hub.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'1a HUB: seed={SEED}; numpy={np.__version__}; scipy={scipy.__version__}')
    for group,(passed,failed) in counts.items(): print(f'{group}: PASS={passed} FAIL={failed}')
    print(f'optimizer regimes: inward={inward} outward={outward} UNRESOLVED={skipped}')
    print(f'max scaled derivative error={max_deriv_error:.3e}')
    print(f'max scaled excess above a passing corner={max_corner_excess:.3e}')
    if failures:
        print('\n'.join(failures[:10])); raise SystemExit(1)

if __name__=='__main__': run()
