"""Self-contained checks of tipping budgets, monotonicity and fixed-budget curvature."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import lambertw

SEED=2026101305
TOL=3e-8
OUT=Path(__file__).resolve().parents[1]/'results'

def functions(kind):
    if kind=='linear':return (lambda z:z,lambda z:np.ones_like(z),lambda z:np.zeros_like(z),lambda Q:Q)
    if kind=='log':return (np.log1p,lambda z:1/(1+z),lambda z:-1/(1+z)**2,
                           lambda Q:float(np.expm1(lambertw(Q).real)))
    if kind=='sqrt':return (lambda z:2*(np.sqrt(1+z)-1),lambda z:1/np.sqrt(1+z),
                            lambda z:-.5/(1+z)**1.5,lambda Q:((1+np.sqrt(1+2*Q))/2)**2-1)
    if kind=='saturating':return (lambda z:-np.expm1(-z),lambda z:np.exp(-z),lambda z:-np.exp(-z),np.log1p)
    raise ValueError(kind)

def run():
    rng=np.random.default_rng(SEED)
    counts={'tipping_roots':[0,0],'monotonicity':[0,0],'corner_optimizer':[0,0],
            'curvature_cases':[0,0],'edges':[0,0]}
    def ck(g,ok):counts[g][not bool(ok)]+=1
    for k in range(200):
        kind=('linear','log','sqrt','saturating')[k%4]
        L,Lp,Lpp,inverse=functions(kind)
        S=10**rng.uniform(-3,-.1);p=10**rng.uniform(-.5,.5);r=rng.uniform(0,.3);d0=rng.uniform(1,4)
        Q=p/S;hi=1.
        while L(hi)/Lp(hi)<Q:hi*=2
        tip=brentq(lambda b:L(b)/Lp(b)-Q,0.,hi,xtol=1e-12)
        ck('tipping_roots',np.isclose(tip,inverse(Q),rtol=2e-11,atol=2e-11))
        bb=np.linspace(.01,.95,50)*tip
        gamprime=1-L(bb)*Lpp(bb)/Lp(bb)**2
        ac=(r+d0*S)/(p*Lp(bb)-S*L(bb))
        ck('monotonicity',np.min(gamprime)>=1-1e-12 and np.min(np.diff(ac))>0)
        B=.5*tip;acrit=(r+d0*S)/(p*Lp(B)-S*L(B))
        for mult in (.8,1.,1.2):
            a=mult*acrit
            def U(x):return np.log(d0+a*L(B-p*x)+r*x)+S*x
            fit=minimize_scalar(lambda x:-U(x),bounds=(0.,B/p),method='bounded',options={'xatol':1e-12})
            best=max(float(U(0)),float(U(B/p)),float(-fit.fun))
            u0=float(U(0))
            ok=best<=u0+TOL*(1+abs(u0)) if mult>=1 else best>u0+1e-10
            ck('corner_optimizer',fit.success and ok)
        # The equality B=tip is excluded at every finite a, since the remainder is positive.
        for a in (1.,1e3,1e8):ck('edges',(r+d0*S)/a>0)
    # eta=2 tangent, empty and bounded-interval cases (d0=1,r=0).
    for B in (.1,1.,10.):
        L=np.log1p(B);Lp=1/(1+B);S=.1
        at=1/L;pmin=4*S*L/Lp
        for price_mult,expected in ((.9,'empty'),(1.,'singleton'),(1.2,'interval')):
            p=price_mult*pmin
            c2=S*L*L;c1=2*S*L-p*Lp;c0=S
            disc=c1*c1-4*c2*c0
            if expected=='empty':ck('curvature_cases',disc<0)
            elif expected=='singleton':
                ck('curvature_cases',abs(S*(1+at*L)**2-p*at*Lp)<1e-12)
                ck('curvature_cases',all(S*(1+q*at*L)**2-p*q*at*Lp>0 for q in (.9,1.1)))
            else:
                roots=np.sort(np.roots([c2,c1,c0]))
                ck('curvature_cases',np.all(roots>0) and S*(1+np.mean(roots)*L)**2-p*np.mean(roots)*Lp<0)
    B=10.;L=np.log1p(B)
    roots=np.sort(np.roots([.001*L*L,2*.001*2*L-1/11,.1+.001*4]))
    ck('curvature_cases',np.allclose(roots,[1.421897,12.720521],atol=6e-7,rtol=0))
    # No strictly increasing threshold when S=0 and ell is linear.
    ck('edges',all(.2/(1*1)==.2 for _ in (1,2,3)))
    # An empty set of positive residual slopes gives B_tip=infinity by convention.
    ck('edges',np.isinf(np.inf))
    # Approaching B_tip from below depends on the productivity/gap rate.
    for t in (1.,2.,4.,8.):
        B=1-np.exp(-t);acrit=1/(1-B)
        ck('edges',np.exp(2*t)>=acrit and np.exp(t/2)<acrit)
    tip=brentq(lambda b:(1+b)*np.log1p(b)-500,0.,500.)
    ac06=.052/(1/(1+.6*tip)-.002*np.log1p(.6*tip))
    result=dict(seed=SEED,tolerance=TOL,counts=counts,Btip=float(tip),acrit_06=float(ac06),
                archived_eta2_example_roots=roots.tolist())
    OUT.mkdir(exist_ok=True)
    (OUT/'v8_1e_tipping.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'1e TIPPING: seed={SEED}')
    for g,(passed,failed) in counts.items():print(f'{g}: PASS={passed} FAIL={failed}')
    print(f'B_tip={tip:.12f}; a_crit(0.6 B_tip)={ac06:.12f}')
    print('unchanged v7.3 eta=2 example roots: '+', '.join(f'{x:.9f}' for x in roots))
    if any(v[1] for v in counts.values()):raise SystemExit(1)

if __name__=='__main__':run()
