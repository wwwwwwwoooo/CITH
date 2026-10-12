"""Executable checks for added CITH growth/curvature/production statements.

All examples are normalized mathematical constructions, not empirical estimates.
Outputs go to ../results/additions_checks.json. No manuscript numbers are typed
into the output independently of the calculations.
"""
from __future__ import annotations
import json
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import brentq, minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'
OUT.mkdir(exist_ok=True)


def critical_price(a: float, B: float, r: float, S: float,
                   d0: float, eta: float) -> float:
    if a <= 0 or B <= 0 or d0 < 1 or S < 0 or r < 0 or eta < 0:
        raise ValueError('Outside the concave model domain')
    return (1 + B) * (r + S * (d0 + a*np.log1p(B))**eta) / a


def main() -> None:
    B, r, S, d0, p = 10., .1, .001, 2., 1.
    L = np.log1p(B)
    H = (1+B)*L*S
    A = (1+B)*(r+d0*S)
    ans: dict = {
        'environment': {'python': platform.python_version(),
                        'numpy': np.__version__, 'scipy': scipy.__version__},
        'parameters': dict(B=B, r=r, S=S, d0=d0, p=p),
    }
    ans['growth'] = {}
    rows = []
    for budget in [10., 100., 1000.]:
        floor = (1+budget)*np.log1p(budget)*S
        price_at_a100 = critical_price(100., budget, r, S, d0, 1.)
        # Construct growing-budget paths satisfying the exact condition.
        # Price growth: p(t)=2H(t); choose a(t)=B(t)+A(t)/H(t).
        productivity = budget + (r+d0*S)/(S*np.log1p(budget))
        growing_price = 2*floor
        need = critical_price(productivity,budget,r,S,d0,1.)
        assert need < growing_price
        # Residual decay: S(t)=1/[2(1+B)log(1+B)], fixed p=1.
        decaying_slope = 1/(2*(1+budget)*np.log1p(budget))
        decay_A = (1+budget)*(r+d0*decaying_slope)
        decay_a = budget + 2*decay_A
        decay_need = critical_price(decay_a,budget,r,decaying_slope,d0,1.)
        assert decay_need < 1
        rows.append(dict(B=budget,H=floor,critical_price_at_a100=price_at_a100,
                         growing_price=growing_price,price_path_a=productivity,
                         price_path_required=need,decaying_slope=decaying_slope,
                         slope_path_a=decay_a,slope_path_required=decay_need))
    floor_crossing=brentq(lambda b:(1+b)*np.log1p(b)*S-p,10,1000,xtol=1e-12)
    ans['growth']={'rows':rows,'fixed_price_floor_crossing_B':floor_crossing,
                   'H_at_crossing':(1+floor_crossing)*np.log1p(floor_crossing)*S}

    curvature=[]
    for eta in [0., .5, 1., 2.]:
        fun=lambda a:critical_price(a,B,r,S,d0,eta)-p
        if eta<1:
            high=1.
            while fun(high)>0:
                high*=2
            root=brentq(fun,1e-12,high,xtol=1e-13)
            sample=[root/2,root,root*2,100.]
            assert fun(root/2)>0 and fun(root*2)<0
            roots=[root]
        elif eta==1:
            root=A/(p-H)
            assert abs(fun(root))<1e-12
            sample=[root/2,root,root*2,100.]
            roots=[root]
        else:
            # S L^2 a^2 + (2 S d0 L - p/(1+B)) a + r+S d0^2 = 0.
            coeff=[S*L**2,2*S*d0*L-p/(1+B),r+S*d0**2]
            roots=sorted(float(x) for x in np.roots(coeff))
            sample=[1.,2.,10.,20.,100.]
            for root in roots:
                assert abs(fun(root))<1e-12
            assert fun(1)>0 and fun(2)<0 and fun(10)<0 and fun(20)>0
        curvature.append({'eta':eta,'roots':roots,
                          'checks':[{'a':a,'required_price':critical_price(a,B,r,S,d0,eta),
                                     'inward':bool(fun(a)<=1e-12)} for a in sample]})
    ans['curvature']={'H':H,'A':A,'cases':curvature}

    # A concave one-channel global objective with the same residual slope at zero.
    # E is absent; residual utility is S*(1-exp(-x)); external production is r*x.
    global_checks=[]
    grid=np.linspace(0,B/p,20001)
    for eta,productivity in [(0.5,2.),(1.,2.),(2.,1.),(2.,2.),(2.,10.),(2.,20.)]:
        def value(x):
            c=d0-1+productivity*np.log1p(B-p*x)+r*x
            f=np.log1p(c) if eta==1 else np.expm1((1-eta)*np.log1p(c))/(1-eta)
            return f+S*(-np.expm1(-x))
        opt=minimize_scalar(lambda x:-value(x),bounds=(0,B/p),method='bounded',
                            options={'xatol':1e-12})
        candidates=np.array([0.,float(opt.x),B/p])
        v=np.array([value(x) for x in candidates]); idx=int(np.argmax(v))
        xstar=float(candidates[idx]); difference=float(v[idx]-value(0))
        predicts=critical_price(productivity,B,r,S,d0,eta)<=p
        if predicts:
            assert difference<1e-11
        else:
            assert difference>1e-9
        assert float(np.max(value(grid)))-float(v[idx])<1e-10
        global_checks.append({'eta':eta,'a':productivity,'corner_predicted':bool(predicts),
                              'x_star':xstar,'improvement_over_corner':difference})
    ans['curvature']['global_optimization_checks']=global_checks

    # Concave external output changes interior allocations but not the corner test.
    # Linear g=x vs concave g=log(1+x), both g'(0)=1; pure computation.
    production=[]
    for budget in [8.,10.,100.]:
        linear_x=max(0.,budget-9.)
        concave_x=max(0.,(budget-9.)/11.)
        corner_derivative=1-10/(1+budget)
        for label,expected in [('linear',linear_x),('concave',concave_x)]:
            def cfun(x):
                ext=x if label=='linear' else np.log1p(x)
                return 10*np.log1p(budget-x)+ext
            opt=minimize_scalar(lambda x:-cfun(x),bounds=(0,budget),method='bounded',
                                options={'xatol':1e-12})
            best=max(cfun(0),cfun(budget),cfun(float(opt.x)))
            assert abs(cfun(expected)-best)<1e-10
            production.append({'B':budget,'external_production':label,
                               'x_star':expected,'z_star':budget-expected,
                               'corner_derivative':corner_derivative,
                               'all_inward':bool(expected==0)})
    # Direct check of the generalized local-production floor S*ell/ell'.
    B2,a2,r2,S2,d02,p2=7.,5.,.2,.003,2.,1.
    local_checks=[]
    for name,ell,ellp in [('logarithmic',np.log1p(B2),1/(1+B2)),
                         ('linear',B2,1.),('square_root_shifted',np.sqrt(1+B2)-1,
                                         .5/np.sqrt(1+B2))]:
        threshold=(r2+d02*S2)/(a2*ellp)+S2*ell/ellp
        direct=(r2/(d02+a2*ell)+S2)/(a2*ellp/(d02+a2*ell))
        assert abs(threshold-direct)<1e-12
        local_checks.append({'local_form':name,'price_threshold':threshold,
                             'direct_marginal_ratio':direct,'residual_floor':S2*ell/ellp})
    ans['production']={'allocation_checks':production,'general_local_form_checks':local_checks,
                       'general_local_form_parameters':dict(B=B2,a=a2,r=r2,S=S2,d0=d02,p=p2)}

    # Exact recurrence and conditional workload/value comparison.
    b,d,c=3,4,.7
    w=(b**(d+1)-1)//(b-1); wn=(b**(d+2)-1)//(b-1)
    assert wn==b*w+1
    ans['recursion']={'b':b,'d':d,'c':c,'W_d':w,'W_next':wn,
                      'cost_ratio':wn/w,'identity_cost_ratio':b+1/w,
                      'log_value_increment':float(np.log1p(c*wn)-np.log1p(c*w)),
                      'limiting_log_value_increment':float(np.log(b))}

    # Extend the matched-efficiency support construction beyond the archived horizon.
    support=[]
    for t in [0.,10.,20.,30.,40.,50.,60.]:
        depth=2+int(t//10); w=2**(depth+1)-1
        ratios=np.array([w/(32*np.exp(.12*t)),w/256,w/(24*np.exp(.12*t))])
        support.append({'t':t,'depth':depth,'worlds':w,'ratios':ratios.tolist(),
                        'max_ratio':float(ratios.max()),'feasible':bool(ratios.max()<=1)})
    assert all(row['feasible'] for row in support[:-1]) and not support[-1]['feasible']
    ans['support_extension']=support
    Svec=np.array([.012,.002,.001]); rvec=np.array([.1,.1,0.])
    ans['joint_witness']=[{'worlds':w,'a':float(w/np.log(11)),
                           'marginal_ratios':(11*np.log(11)*(Svec+(rvec+Svec)/w)).tolist()}
                          for w in [7,127]]
    (OUT/'additions_checks.json').write_text(json.dumps(ans,indent=2,allow_nan=False)+'\n')
    print(json.dumps(ans,indent=2,allow_nan=False))

if __name__=='__main__':
    main()
