"""Self-contained Stage 1d global checks, including smooth counterexamples.

The grid-plus-local-refinement global search does not assume endpoint optimality.
The endpoint theorem is tested ONLY with affine external production and linear
residual utility. Two counterexamples test extensions outside those assumptions.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar

SEED=2026101304
TOL=2e-8
OUT=Path(__file__).resolve().parents[1]/'results'

def grid_max(fun, lo, hi):
    grid=np.linspace(lo,hi,2001);vv=fun(grid)
    candidates=[(float(grid[0]),float(vv[0])),(float(grid[-1]),float(vv[-1]))]
    idx=np.flatnonzero((vv[1:-1]>vv[:-2]) & (vv[1:-1]>vv[2:]))+1
    for k in idx:
        fit=minimize_scalar(lambda x:-float(fun(x)),bounds=(grid[k-1],grid[k+1]),
                            method='bounded',options={'xatol':1e-12})
        if not fit.success:raise RuntimeError(fit.message)
        candidates.append((float(fit.x),float(-fit.fun)))
    k=int(np.argmax(vv));candidates.append((float(grid[k]),float(vv[k])))
    return max(candidates,key=lambda row:row[1])

def run():
    rng=np.random.default_rng(SEED)
    counts={'endpoint_equivalence':[0,0],'Gamma_identities':[0,0],'edges':[0,0],
            'smooth_counterexamples':[0,0]}
    def ck(g,ok):counts[g][not bool(ok)]+=1
    foc_pass=foc_pass_global_fail=unresolved=0
    for i in range(2400):
        B=10**rng.uniform(-1,1);p=10**rng.uniform(-.5,.5);psi=10**rng.uniform(-1,.3)
        a=10**rng.uniform(-3,1.7);d0=10**rng.uniform(0,1)
        r=0. if i%11==0 else 10**rng.uniform(-2,.5)
        S=0. if i%13==0 else 10**rng.uniform(-3,.7)
        def U(x):return np.log(d0+a*np.expm1(psi*(B-p*x))/psi+r*x)+S*x
        D0=d0+a*np.expm1(psi*B)/psi
        f0=(r-a*p*np.exp(psi*B))/D0+S
        u0=float(U(0));ue=float(U(B/p))
        xm,um=grid_max(U,0.,B/p)
        if f0<=0:
            foc_pass+=1
            if um>u0+TOL*(1+abs(u0)):foc_pass_global_fail+=1
        predicted=(f0<=0 and u0>=ue)
        observed=um<=u0+TOL*(1+abs(u0))
        if predicted:
            ck('endpoint_equivalence',observed)
        elif not observed:
            ck('endpoint_equivalence',True)
        elif abs(f0)<1e-6 or abs(ue-u0)<TOL*(1+abs(u0)):
            unresolved+=1
        else:ck('endpoint_equivalence',False)
        gam=-np.expm1(-psi*B)/psi
        raw=(np.expm1(psi*B)/psi)/np.exp(psi*B)
        ck('Gamma_identities',np.isclose(gam,raw,rtol=2e-13,atol=2e-13) and gam<1/psi)
    for power in (1.,2.,5.):
        for B in (.01,1.,100.):
            gam=B**power/(power*B**(power-1))
            ck('Gamma_identities',np.isclose(gam,B/power) and gam<=B*(1+1e-12))
    # Every allocation ties: a=psi*d0, S=psi*p, r=0.
    xx=np.linspace(0,3,301)
    vals=np.log(1+np.expm1(3-xx))+xx
    ck('edges',np.max(abs(vals-3))<1e-12)
    # A strict FOC does not suffice even with a linear residual.
    B=5.;a=2.;r=.3
    D=1+a*np.expm1(B)
    f0=(r-a*np.exp(B))/D+1
    u0=np.log(D);ue=np.log(1+r*B)+B
    ck('edges',f0<0 and ue>u0)
    linear_example=dict(derivative=float(f0),U0=float(u0),Uend=float(ue))
    # Smooth, increasing, concave residual; r=0.
    def R(x):return .76/20*(np.logaddexp(0,140)-np.logaddexp(0,20*(7-x)))
    def UR(x):return np.log(1+.01*np.expm1(10-x))+R(x)
    f0=-.01*np.exp(10)/(1+.01*np.expm1(10))+.76/(1+np.exp(-140))
    xm,um=grid_max(UR,0.,10.)
    residual_example=dict(derivative=float(f0),U0=float(UR(0)),U7=float(UR(7)),
                          Uend=float(UR(10)),xstar=xm,Ustar=um)
    ck('smooth_counterexamples',f0<0 and UR(10)<UR(0)<UR(7) and um>UR(0))
    # Smooth, increasing, concave external production; S=0.
    def G(x):return 9/50*(np.logaddexp(0,100)-np.logaddexp(0,50*(2-x)))
    def UG(x):return np.log(np.exp(3-x)+G(x))
    f0=-1+9/np.exp(3)
    xm,um=grid_max(UG,0.,3.)
    production_example=dict(derivative=float(f0),U0=float(UG(0)),U2=float(UG(2)),
                            Uend=float(UG(3)),xstar=xm,Ustar=um)
    ck('smooth_counterexamples',f0<0 and UG(3)<UG(0)<UG(2) and um>UG(0))
    result=dict(seed=SEED,tolerance=TOL,counts=counts,FOC_pass=foc_pass,
                FOC_pass_but_not_global=foc_pass_global_fail,unresolved=unresolved,
                linear_example=linear_example,residual_counterexample=residual_example,
                production_counterexample=production_example)
    OUT.mkdir(exist_ok=True)
    (OUT/'v8_1d_exponential.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'1d EXPONENTIAL: seed={SEED}')
    for g,(passed,failed) in counts.items():print(f'{g}: PASS={passed} FAIL={failed}')
    print(f'FOC passes={foc_pass}; of these not globally inward={foc_pass_global_fail}; UNRESOLVED={unresolved}')
    for name,v in [('linear residual, FOC-only',linear_example),('concave residual',residual_example),
                   ('concave external production',production_example)]:
        print(name+': '+', '.join(f'{k}={value:.9f}' for k,value in v.items()))
    if any(v[1] for v in counts.values()):raise SystemExit(1)

if __name__=='__main__':run()
