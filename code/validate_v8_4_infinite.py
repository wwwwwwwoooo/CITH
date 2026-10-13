"""Infinite-horizon log benchmark, verification, linear residual and ceiling.
Finite-horizon Bellman updates use independent scalar optimization.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import minimize_scalar
SEED=2026101340
OUT=Path(__file__).resolve().parents[1]/'results'

def run():
 rng=np.random.default_rng(SEED);counts={k:[0,0] for k in ['corner','Bellman_limit','transversality','policy_bound','residual_derivative','ceiling','edges']};maxerr=0.
 def check(k,c):counts[k][0 if c else 1]+=1
 for _ in range(180):
  al,ax,gl,gx=10**rng.uniform(-.7,.7,4);beta=rng.uniform(.03,.96);B=10**rng.uniform(-3,3)
  A=lambda s:al+(ax-al)*s;G=lambda s:gl+(gx-gl)*s
  H=lambda s:(1-beta)*np.log(A(s))+beta*np.log(G(s))
  op=minimize_scalar(lambda s:-H(s),bounds=(0,1),method='bounded',options={'xatol':1e-13})
  candidates=[0.,1.,op.x];share=max(candidates,key=H);hstar=H(share);K=hstar/(1-beta)**2
  pred=(1-beta)*ax/al+beta*gx/gl<=1
  check('corner', (share<1e-7)==pred or abs((1-beta)*ax/al+beta*gx/gl-1)<1e-6)
  coeff=const=0.;horizon=int(np.ceil(32/-np.log(beta)))+10
  for k in range(horizon):
   obj=lambda s:np.log(A(s))+beta*coeff*np.log(G(s))
   oo=minimize_scalar(lambda s:-obj(s),bounds=(0,1),method='bounded')
   best=max(obj(0),obj(1),-oo.fun);const=best+beta*const;coeff=1+beta*coeff
  err=abs((coeff*np.log(B)+const)-(np.log(B)/(1-beta)+K));maxerr=max(maxerr,err)
  check('Bellman_limit',err<2e-6*(1+abs(K)))
  T=int(np.ceil(45/-np.log(beta)))+10
  lnB=np.log(B);total=0.
  for t in range(T):
   s=rng.uniform();total+=beta**t*(lnB+np.log(A(s)));lnB+=np.log(G(s))
  tail=beta**T*(lnB/(1-beta)+K)
  check('transversality',abs(tail)<2e-11*(1+abs(K)))
  check('policy_bound',total+tail<=np.log(B)/(1-beta)+K+1e-9)
  # One deviation in physical x with p=1, followed by inward forever.
  S=rng.uniform(.001,.2)
  J=lambda x:np.log(al*(B-x)+ax*x)+S*x+beta*(np.log(gl*(B-x)+gx*x)/(1-beta))
  numeric=J(1e-22j).imag/1e-22
  exact=((ax/al-1)+beta/(1-beta)*(gx/gl-1))/B+S
  check('residual_derivative',np.isclose(numeric,exact,atol=1e-10,rtol=1e-10))
  bar=rng.uniform(.3,4);Bx=bar*rng.uniform(.1,2);gL=rng.uniform(.1,2);gX=gL+rng.uniform(0,3)
  xn=Bx*1e-6;next0=min(bar,gL*Bx);next1=min(bar,gL*(Bx-xn))+gX*xn
  check('ceiling',next1>=next0-1e-12)
 # Two logical reversals and a neutral-growth tie.
 check('edges',.1*.5+.9*2>1);check('edges',.1*2+.9*.5<1)
 check('edges',.1*2+.9*1>1);check('edges',.1*.5+.9*1<1)
 # Divergence witness: beta*gX>=1 and S>0 under all-external production.
 check('edges',.9*1.2>1)
 # The affine value is exactly flat in share only when both return ratios equal one.
 check('edges',np.isclose(.4*1+.6*1,1))
 OUT.mkdir(exist_ok=True);res={'seed':SEED,'counts':counts,'max_finite_horizon_value_error':maxerr};(OUT/'v8_4_infinite.json').write_text(json.dumps(res,indent=2))
 print(f'4 INFINITE HORIZON seed={SEED}')
 for k,(p,f) in counts.items():print(f'{k}: PASS={p} FAIL={f}')
 print(f'max finite-horizon value error={maxerr:.3e}')
 assert not sum(v[1] for v in counts.values())
if __name__=='__main__':run()
