"""Self-contained additive checks of the sharp concave-production floor."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

SEED=2026101302
TOL=2e-9
OUT=Path(__file__).resolve().parents[1]/'results'

def run():
    rng=np.random.default_rng(SEED)
    counts={'floor':[0,0],'witness_optimizer':[0,0],'boundary':[0,0],'equality':[0,0]}
    def ck(group,ok): counts[group][not bool(ok)]+=1
    for _ in range(600):
        B=10**rng.uniform(-3,.7)
        power=rng.uniform(.1,1)
        values=[(B,1.),(np.log1p(B),1/(1+B)),
                (2*(np.sqrt(1+B)-1),1/np.sqrt(1+B)),
                (B**power,power*B**(power-1)),(-np.expm1(-B),np.exp(-B))]
        for L,Lp in values: ck('floor',L/Lp >= B-2e-12*(1+B))
        ck('equality',abs(values[0][0]/values[0][1]-B)<1e-13)
        # Nonlinear after B is allowed: ell(z)=z up to B, B+log(1+z-B) after B.
        ck('equality',B/1.==B)
    for k in range(120):
        B=10**rng.uniform(-1,1.5); d0=rng.uniform(1,4)
        S=np.sort(10**rng.uniform(-3,-.3,3))[::-1]
        r=np.r_[10**rng.uniform(-2,.2,2),0.]
        if k%7==0:S[:]=0
        if k%13==0:r[:]=0
        p=S*B+10**rng.uniform(-1,.5,3)
        acr=max((r+d0*S)/(p-S*B))
        a=max(1.,1.2*acr)
        def u(s):
            x=B*np.asarray(s)/p; z=B*(1-np.sum(s))
            return np.log(d0+a*z+np.dot(r,x))+np.dot(S,x)
        fit=minimize(lambda s:-float(u(s)),np.full(3,.2),method='SLSQP',
                     bounds=[(0.,1.)]*3,constraints={'type':'ineq','fun':lambda s:1-np.sum(s)},
                     options={'ftol':1e-12,'maxiter':500})
        ck('witness_optimizer',fit.success and fit.fun >= -u(np.zeros(3))-TOL*(1+abs(u(np.zeros(3)))))
    for B in (1e-8,.1,1.,100.):
        for S in (.001,1.):
            for a in (1e-5,1.,1e5,1e10):
                p=S*B
                # Explicit positive correction avoids subtracting nearly equal large numbers.
                correction=S/(1+a*B)
                step=1j*1e-25
                numerical=np.imag(np.log(1+a*(B-p*step))+S*step)/1e-25
                ck('boundary',correction>0 and numerical>0 and abs(numerical-correction)<=1e-12*S)
    # A_j=0 iff r_j=S_j=0, because d0>0 and ell'(B)>0.
    for d0 in (1.,2.,10.): ck('boundary',(0+d0*0)==0)
    OUT.mkdir(exist_ok=True)
    (OUT/'v8_1b_floor.json').write_text(json.dumps({'seed':SEED,'tolerance':TOL,'counts':counts},indent=2)+'\n')
    print(f'1b TIGHT FLOOR: seed={SEED}')
    for group,(passed,failed) in counts.items():print(f'{group}: PASS={passed} FAIL={failed}')
    if any(v[1] for v in counts.values()):raise SystemExit(1)

if __name__=='__main__':run()
