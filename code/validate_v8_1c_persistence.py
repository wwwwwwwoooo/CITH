"""Self-contained additive checks of Stage 1c and bounded-budget edge cases."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar

SEED=2026101303
TOL=3e-8
OUT=Path(__file__).resolve().parents[1]/'results'

def run():
    rng=np.random.default_rng(SEED)
    counts={'minimum':[0,0],'unbounded_obstruction':[0,0],
            'sublog_scaling':[0,0],'necessary_rate':[0,0],
            'eta_zero':[0,0],'counterexamples':[0,0]}
    details=[]
    def ck(g,ok):counts[g][not bool(ok)]+=1
    for _ in range(240):
        eta=1+10**rng.uniform(-1.5,.6); d0=10**rng.uniform(0,1.2)
        omstar=d0/(eta-1)
        logmin=(eta-1)*np.log(d0)+eta*np.log(eta)-(eta-1)*np.log(eta-1)
        fit=minimize_scalar(lambda v:eta*np.logaddexp(np.log(d0),v)-v,
                            bounds=(np.log(omstar)-15,np.log(omstar)+15),
                            method='bounded',options={'xatol':1e-12})
        ck('minimum',fit.success and abs(fit.fun-logmin)<TOL)
    for eta in (1.,1.2,2.,4.):
        d0=1.;s=.02;M=1.
        m=1 if eta==1 else d0**(eta-1)*eta**eta/(eta-1)**(eta-1)
        for B in (1e3,1e5,1e8):
            for L,Lp in ((B,1.),(np.log1p(B),1/(1+B)),
                         (2*(np.sqrt(1+B)-1),1/np.sqrt(1+B))):
                gam=L/Lp
                ck('unbounded_obstruction',s*gam*m>M)
                for om in np.logspace(-4,10,15):
                    req=s*gam*(d0+om)**eta/om
                    ck('unbounded_obstruction',req>M)
    p=np.array([1.,1.4,.8]);r=np.array([.05,.03,0.]);S=np.array([.002,.001,.0005])
    d0=1.
    for eta in (.25,.5,.75):
        com=max(1.,2*max(S/p)**(1/(1-eta)))
        for B in (1e2,1e4,1e8):
            for name,L,Lp in (('linear',B,1.),('log',np.log1p(B),1/(1+B)),
                              ('sqrt',2*(np.sqrt(1+B)-1),1/np.sqrt(1+B))):
                gam=L/Lp; om=com*gam**(1/(1-eta)); a=om/L
                req=gam*(r/om+S*np.exp(eta*np.log(d0+om))/om)
                ck('sublog_scaling',np.all(req<p))
                ck('necessary_rate',np.all(om+1e-12 >= (S*gam/p)**(1/(1-eta))))
                if eta==.5 and name=='linear':details.append(dict(B=B,a=float(a),max_price_ratio=float(max(req/p))))
    # Original proposed linear-production growth rate at eta=1/2 through B=10^8.
    for B in np.logspace(1,8,29):
        a=.01*B; om=a*B
        req=(.05+.002*np.sqrt(1+om))/a
        ck('sublog_scaling',req<1)
    for _ in range(200):
        B=10**rng.uniform(-2,4); eta=rng.uniform(0,.95)
        L=np.log1p(B);Lp=1/(1+B);gam=L/Lp
        om=10**rng.uniform(-2,12);p0=10**rng.uniform(-1,1);s=10**rng.uniform(-4,0)
        req=s*gam*(1+om)**eta/om
        if req<=p0: ck('necessary_rate',om>=(s*gam/p0)**(1/(1-eta))*(1-1e-12))
    for B in (.01,1.,100.):
        for r0,s,p0 in ((0.,0.,1.),(.2,0.,1.),(0.,.2,1.),(.2,.3,1.)):
            for kind in ('linear','log'):
                L=B if kind=='linear' else np.log1p(B)
                Lp=1 if kind=='linear' else 1/(1+B)
                acr=(r0+s)/(p0*Lp)
                a=max(acr,1e-4)
                ck('eta_zero',r0+s<=p0*a*Lp+1e-12)
    # Counterexample to sufficiency with only an upper price bound.
    for B in (10.,30.,100.):
        a=10*B;p0=np.exp(-B)
        ck('counterexamples',np.sqrt(1+a*B)>p0*a)
    # B bounded above is not enough for the eta>1 high-a exclusion.
    for a in (4.,10.,100.,1e6):
        B=a**-2
        ck('counterexamples',(1+a*B)**2<=a)
    # A zero residual slope removes the high-a obstruction at every eta.
    for eta in (0.,1.,2.,5.):ck('counterexamples',.1<=100/(1+10))
    OUT.mkdir(exist_ok=True)
    (OUT/'v8_1c_persistence.json').write_text(json.dumps({'seed':SEED,'tolerance':TOL,
            'counts':counts,'eta_half_examples':details},indent=2)+'\n')
    print(f'1c PERSISTENCE: seed={SEED}')
    for group,(passed,failed) in counts.items():print(f'{group}: PASS={passed} FAIL={failed}')
    print('eta=0.5, ell(z)=z, a=0.01 B: checked through B=1e8')
    print('counterexamples: vanishing prices; B=a^-2 at eta=2; zero residual slopes')
    if any(v[1] for v in counts.values()):raise SystemExit(1)

if __name__=='__main__':run()
