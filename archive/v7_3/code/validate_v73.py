"""Recompute only the new v7.3 examples (F3 and O1).

The v7.2 data, result tables, and Figures 1, 2, and 4 are not written.
Analytical endpoint/derivative solutions are checked with scalar and joint
optimizers. All values quoted in the new prose are emitted here.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
from pathlib import Path
import json
import platform
import numpy as np
import scipy
from scipy.optimize import brentq, minimize, minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'

@dataclass(frozen=True)
class Parameters:
    B1: float = 10.0
    p1: float = 1.0
    p2: float = 1.0
    r1: float = 1.0
    r2: float = 1.0
    a1: float = 5.0
    a2: float = 100.0
    C01: float = 1.0
    C02: float = 1.0
    B20: float = 0.0
    gL: float = 0.1
    gX: float = 0.0
    beta: float = 0.9


def stage2(B: float, p: Parameters) -> tuple[float, float, float, float]:
    if B < -1e-12:
        raise ValueError('Negative continuation budget')
    B = max(0.0, B)
    cap = max(p.a2 * p.p2 / p.r2 - 1.0, 0.0)
    z = min(B, cap)
    x = (B - z) / p.p2
    denom = 1.0 + p.C02 + p.a2 * np.log1p(z) + p.r2 * x
    marginal = max(p.a2 / (1.0 + B), p.r2 / p.p2) / denom
    return float(np.log(denom)), float(marginal), z, x


def objective_derivative(x: float, p: Parameters) -> tuple[float, float, float, float]:
    z = p.B1 - p.p1 * x
    d1 = 1.0 + p.C01 + p.a1 * np.log1p(z) + p.r1 * x
    B2 = p.B20 + p.gL * z + p.gX * x
    V2, V2prime, _, _ = stage2(B2, p)
    static = (p.r1 - p.p1 * p.a1 / (1.0 + z)) / d1
    future = p.beta * V2prime * (p.gX - p.p1 * p.gL)
    return float(np.log(d1) + p.beta * V2), float(static + future), float(static), float(future)


def check_example(name: str, p: Parameters) -> dict:
    xmax = p.B1 / p.p1
    z_static = min(p.B1, max(p.a1 * p.p1 / p.r1 - 1.0, 0.0))
    x_static = (p.B1 - z_static) / p.p1
    d0 = objective_derivative(0.0, p)[1]
    dlast = objective_derivative(xmax, p)[1]
    if d0 <= 0.0:
        xstar = 0.0
    elif dlast >= 0.0:
        xstar = xmax
    else:
        xstar = float(brentq(lambda x: objective_derivative(x, p)[1], 0.0, xmax, xtol=1e-13))
    value = objective_derivative(xstar, p)[0]
    scalar = minimize_scalar(lambda x: -objective_derivative(x, p)[0],
                             bounds=(0.0, xmax), method='bounded',
                             options={'xatol': 1e-12})
    if not scalar.success:
        raise RuntimeError(scalar.message)
    candidates = [0.0, float(scalar.x), xmax]
    scalar_x = max(candidates, key=lambda x: objective_derivative(x, p)[0])
    scalar_error = abs(value - objective_derivative(scalar_x, p)[0])
    B2star = p.B20 + p.gL * (p.B1 - p.p1 * xstar) + p.gX * xstar
    _, _, z2star, x2star = stage2(B2star, p)

    # Independent joint problem in (x1, z2, x2), without stage2().
    def joint_neg(q: np.ndarray) -> float:
        x1, z2, x2 = q
        c1 = 1.0 + p.C01 + p.a1 * np.log1p(p.B1 - p.p1 * x1) + p.r1 * x1
        c2 = 1.0 + p.C02 + p.a2 * np.log1p(z2) + p.r2 * x2
        return -float(np.log(c1) + p.beta * np.log(c2))

    def joint_jac(q: np.ndarray) -> np.ndarray:
        x1, z2, x2 = q
        z1 = p.B1 - p.p1 * x1
        c1 = 1.0 + p.C01 + p.a1 * np.log1p(z1) + p.r1 * x1
        c2 = 1.0 + p.C02 + p.a2 * np.log1p(z2) + p.r2 * x2
        return np.array([-(p.r1-p.p1*p.a1/(1+z1))/c1,
                         -p.beta*p.a2/((1+z2)*c2), -p.beta*p.r2/c2])

    def budget(q: np.ndarray) -> float:
        x1, z2, x2 = q
        return p.B20 + p.gL * (p.B1 - p.p1*x1) + p.gX*x1 - z2 - p.p2*x2

    joint_results = []
    for x_start in (0.0, xmax/2.0, xmax):
        b_start = p.B20 + p.gL*(p.B1-p.p1*x_start) + p.gX*x_start
        initial = np.array([x_start, b_start/2, b_start/(2*p.p2)])
        fit = minimize(joint_neg, initial, jac=joint_jac, method='SLSQP',
                       bounds=((0.0,xmax),(0.0,None),(0.0,None)),
                       constraints={'type':'ineq','fun':budget,
                                    'jac':lambda q: np.array([p.gX-p.p1*p.gL,-1.0,-p.p2])},
                       options={'ftol':1e-12,'maxiter':1000})
        if not fit.success:
            raise RuntimeError(f'{name}: {fit.message}')
        if budget(fit.x) < -1e-9:
            raise AssertionError('Infeasible joint result')
        joint_results.append({'x1':float(fit.x[0]),'z2':float(fit.x[1]),'x2':float(fit.x[2]),
                              'value':float(-fit.fun),'budget_slack':float(budget(fit.x))})
    joint_error = max(abs(value-r['value']) for r in joint_results)
    grid = np.linspace(0.0,xmax,2001)
    grid_values = np.array([objective_derivative(float(x),p)[0] for x in grid])
    grid_derivatives = np.array([objective_derivative(float(x),p)[1] for x in grid])
    assert scalar_error < 1e-9 and joint_error < 1e-9
    assert np.max(grid_values) <= value + 1e-10
    assert np.max(np.diff(grid_derivatives)) <= 1e-10
    np.savez_compressed(OUT/f'v73_{name}_grid.npz',x=grid,value=grid_values,derivative=grid_derivatives)
    _, total, static, future = objective_derivative(0.0,p)
    return {'parameters':asdict(p),'static_term_at_zero':static,'future_term_at_zero':future,
            'total_derivative_at_zero':total,'static_xstar':x_static,'dynamic_xstar':xstar,
            'dynamic_z1star':p.B1-p.p1*xstar,'B2star':B2star,'z2star':z2star,'x2star':x2star,
            'dynamic_value':value,'scalar_xstar':scalar_x,'scalar_objective_error':scalar_error,
            'joint_results':joint_results,'max_joint_objective_error':joint_error,
            'grid_points':len(grid),'max_grid_excess':float(np.max(grid_values)-value)}


def check_production_bound() -> list[dict]:
    # sqrt(z) is differentiable at B>0, not at zero; the proof needs only B.
    functions = [
        ('sqrt(z)',lambda B:np.sqrt(B),lambda B:0.5/np.sqrt(B),lambda B:2*B),
        ('1-exp(-z)',lambda B:-np.expm1(-B),lambda B:np.exp(-B),lambda B:np.expm1(B)),
        ('z',lambda B:B,lambda B:1.0,lambda B:B),
    ]
    rows=[]
    S=0.001
    for name,ell,derivative,identity in functions:
        for B in (1.0,2.0,10.0):
            ratio=float(ell(B)/derivative(B))
            floor=S*ratio
            lower=S*B
            tangent_slack=float(ell(B)-B*derivative(B))
            assert ratio >= B-1e-12 and tangent_slack >= -1e-12
            assert np.isclose(ratio,identity(B),rtol=2e-14,atol=1e-14)
            rows.append({'ell':name,'B':B,'S':S,'ell_B':float(ell(B)),
                         'ell_prime_B':float(derivative(B)),'ell_over_derivative':ratio,
                         'residual_floor':floor,'linear_lower_bound':lower,
                         'tangent_slack':tangent_slack})
    return rows


def main() -> None:
    OUT.mkdir(exist_ok=True)
    forward = check_example('static_inward_dynamic_outward',Parameters(a1=100.0,B20=1.0,gL=0.2,gX=2.0))
    reverse = check_example('static_outward_dynamic_inward',Parameters())
    assert forward['static_term_at_zero'] < 0 < forward['total_derivative_at_zero']
    assert forward['static_xstar'] == 0 and forward['dynamic_xstar'] > 0
    assert reverse['static_term_at_zero'] > 0 > reverse['total_derivative_at_zero']
    assert reverse['static_xstar'] == 6 and reverse['dynamic_xstar'] == 0
    result={'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
            'two_period':{'static_inward_dynamic_outward':forward,'static_outward_dynamic_inward':reverse},
            'production_lower_bound':check_production_bound()}
    # Directly generated LaTeX prevents transcription errors in new numbers.
    macros=(
        '% Generated by code/validate_v73.py; original result_macros.tex is unchanged.\n'
        +f'\\newcommand{{\\ReverseStaticTerm}}{{{reverse["static_term_at_zero"]:+.5f}}}\n'
        +f'\\newcommand{{\\ReverseFutureTerm}}{{{reverse["future_term_at_zero"]:+.5f}}}\n'
        +f'\\newcommand{{\\ReverseTotalTerm}}{{{reverse["total_derivative_at_zero"]:+.5f}}}\n'
        +f'\\newcommand{{\\ReverseStaticX}}{{{reverse["static_xstar"]:.0f}}}\n'
        +f'\\newcommand{{\\ReverseDynamicX}}{{{reverse["dynamic_xstar"]:.0f}}}\n'
    )
    (ROOT/'manuscript/v73_result_macros.tex').write_text(macros)
    data=json.dumps(result,indent=2,allow_nan=False)+'\n'
    (OUT/'v73_checks.json').write_text(data)
    print(data,end='')

if __name__=='__main__':
    main()
