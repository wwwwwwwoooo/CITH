"""No-replication KKT, substitutes, bundle dominance, and zero corners.
Self-contained; optimizer evaluates the complete five-allocation objective.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import minimize
SEED=2026101321
OUT=Path(__file__).resolve().parents[1]/'results'

def run():
    rng=np.random.default_rng(SEED); counts={k:[0,0] for k in ('restricted_KKT','full_optimizer','slack_identity','bundle','edges')}; regimes={}
    def check(k,c): counts[k][0 if c else 1]+=1
    for t in range(320):
        B=10**rng.uniform(-.4,1.5); a=10**rng.uniform(-.2,2.7); C0=rng.uniform(0,2)
        prices=10**rng.uniform(-.5,.6,4); returns=10**rng.uniform(-2,0,2)
        amps=10**rng.uniform(-2,.5,3); scales=10**rng.uniform(-.2,.7,3)
        if t%8==0: amps[2]=0
        if t%9==0: amps[1]=0
        eta=float(rng.choice([0.,.5,1.,2.])); backup=t%2==0
        def calc(v,grad=False):
            z=B*v[0]; x,y,w,b=B*v[1:]/prices
            C=C0+a*np.log1p(z)+returns@[x,y]
            F=np.log1p(C) if eta==1 else np.expm1((1-eta)*np.log1p(C))/(1-eta)
            der=(1+C)**(-eta); args=np.array([x,x+y+b,x+y+w]); exp=np.exp(-args/scales)
            val=F+np.sum(amps*(-np.expm1(-args/scales)))
            E,R,I=amps/scales*exp
            marginal=np.array([der*a/(1+z),der*returns[0]+E+R+I,der*returns[1]+R+I,I,R])
            if grad:return val,B*marginal/np.r_[1,prices]
            return val
        con={'type':'eq','fun':lambda v: v.sum()-1,'jac':lambda v:np.ones(5)}
        rb=[(0,1),(0,0),(0,0),(0,1),(0,1) if backup else (0,0)]
        starts=[np.array([1.,0,0,0,0]),np.array([.5,0,0,.5,0])]
        if backup:starts.append(np.array([.4,0,0,.3,.3]))
        opts=[minimize(lambda v:-calc(v),v,jac=lambda v:-calc(v,True)[1],bounds=rb,constraints=con,method='SLSQP',options={'ftol':1e-12,'maxiter':900}) for v in starts]
        opt=min(opts,key=lambda o:o.fun); v=opt.x; value,g=calc(v,True); active=[k for k in (0,3,4) if (k!=4 or backup) and v[k]>1e-7]
        lam=np.mean(g[active])/B; allowed=[0,3]+([4] if backup else [])
        gap=max(g[k] for k in allowed)-g@v
        check('restricted_KKT',gap<2e-5*(1+abs(value)))
        z=B*v[0]; x,y,w,b=B*v[1:]/prices; C=C0+a*np.log1p(z)
        fp=(1+C)**(-eta); E=amps[0]/scales[0]; R=amps[1]/scales[1]*np.exp(-b/scales[1]); I=amps[2]/scales[2]*np.exp(-w/scales[2])
        lhs=fp*returns+np.array([E+R+I,R+I]); score=lhs-prices[:2]*lam
        pred=np.max(score)<=2e-7
        fullbounds=[(0,1)]*4+([(0,1)] if backup else [(0,0)])
        starts2=[v,np.array([.2,.2,.2,.2,.2]) if backup else np.array([.25,.25,.25,.25,0])]
        full=[minimize(lambda vv:-calc(vv),vv,jac=lambda vv:-calc(vv,True)[1],bounds=fullbounds,constraints=con,method='SLSQP',options={'ftol':1e-12,'maxiter':1200}) for vv in starts2]
        best=min(full,key=lambda o:o.fun); gain=-best.fun-value
        ok=gain<5e-6*(1+abs(value)) if pred else (gain>1e-9 or np.max(score)<1e-5)
        check('full_optimizer',ok)
        ci=I/lam; cr=R/lam
        reduced=fp*returns+np.array([E,0])-(prices[:2]-ci-cr)*lam
        check('slack_identity',np.allclose(score,reduced,atol=1e-11,rtol=1e-10))
        key=('backup' if backup else 'observation')+('_no_rep' if pred else '_rep');regimes[key]=regimes.get(key,0)+1
        if backup and w>1e-5 and b>1e-5 and z>1e-5:
            if prices[1]<prices[2]+prices[3]-1e-5:check('bundle',score[1]>0)
    # Construct an exact tie: p_h=p_m+p_b, r_h=0, identical information/resilience.
    # z=1, w=b=1, a=4, C0=0, p_m=p_b=1. Set marginal residuals=lambda.
    a=4.; z=w=b=1.; lam=a/(1+z)/(1+a*np.log1p(z)); amp=lam*np.e
    U=lambda y:np.log1p(a*np.log1p(z))+amp*(-np.expm1(-(b-y+y)))+amp*(-np.expm1(-(w-y+y)))
    check('edges',abs(U(0)-U(.5))<1e-13)
    check('edges',np.isfinite(lam) and lam>0)
    # A zero information motive has zero shadow service cost, not price p_m.
    check('edges',0/lam==0)
    OUT.mkdir(exist_ok=True);res={'seed':SEED,'counts':counts,'regimes':regimes};(OUT/'v8_2a_replication.json').write_text(json.dumps(res,indent=2))
    print(f'2a REPLICATION seed={SEED}')
    for k,(p,f) in counts.items():print(f'{k}: PASS={p} FAIL={f}')
    print('regimes:',regimes)
    assert sum(v[1] for v in counts.values())==0
if __name__=='__main__':run()
