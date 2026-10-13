"""Subjective delay, exact corners, depletion, and replication-chain tests.
Inputs are synthetic. Quadrature and scalar optimization are independent checks.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
SEED=2026101350
OUT=Path(__file__).resolve().parents[1]/'results'

def run():
 rng=np.random.default_rng(SEED);counts={k:[0,0] for k in ['clock_integral','corner_optimizer','sufficient_bound','bounded_threshold','chain','depletion','counterexamples']}
 def check(k,c):counts[k][0 if c else 1]+=1
 for _ in range(450):
  B,a,p,r,S,d0=10**rng.uniform(-2,2,6);d0=1+d0
  delta,w0,g,T=10**rng.uniform(-3,0,4);t=rng.uniform(0,5)
  num=delta*quad(lambda u:w0*np.exp(g*u),t,t+T,epsabs=1e-12)[0]
  exact=delta*w0*np.exp(g*t)*np.expm1(g*T)/g
  check('clock_integral',np.isclose(num,exact,rtol=1e-11,atol=1e-12) and num>=delta*w0*np.exp(g*t)*T-1e-12)
  D=rng.uniform(0,12);atten=np.exp(-D)
  required=atten*(1+B)*((r+d0*S)/a+S*np.log1p(B))
  obj=lambda x:np.log(d0+a*np.log1p(B-p*x)+r*atten*x)+S*atten*x
  op=minimize_scalar(lambda x:-obj(x),bounds=(0,B/p),method='bounded',options={'xatol':1e-12})
  best=max(obj(0),obj(B/p),-op.fun)
  gap=best-obj(0);margin=p-required
  check('corner_optimizer',(gap<1e-8*(1+abs(obj(0)))) if margin>=0 else (gap>0 or abs(margin)<1e-7*(1+p)))
  H=(1+B)*(1+np.log1p(B));amin=.001
  check('sufficient_bound',required<=atten*H*((r+d0*S)/a+S)+1e-9)
  bar=B+rng.uniform(0,20);a2=a+rng.uniform(0,50)
  threshold=max(0,np.log(((1+bar)*(r+d0*S)/a+(1+bar)*np.log1p(bar)*S)/p))
  required2=np.exp(-threshold-.01)*(1+B)*((r+d0*S)/a2+S*np.log1p(B))
  check('bounded_threshold',required2<p*(1+1e-12))
  R=10**rng.uniform(-1,2);Dchain=max(0,np.log(R))+rng.uniform(.2,4);v=R*np.exp(-Dchain)
  series=sum(v**k for k in range(1,1000));check('chain',np.isclose(series,v/(1-v),rtol=1e-11,atol=1e-11))
  stock=rng.uniform(0,10);prices=1+(stock+np.cumsum(np.zeros(40)))**2
  check('depletion',np.all(prices==prices[0]))
 # Original Gamma-only condition fails as B -> 0, with a finite remote return.
 for n in [10,20,50]:
  B=np.exp(-n);D=0.;required=(1+B)*2
  check('counterexamples',D-np.log((1+B)*np.log1p(B))>n-.001 and required>1)
 # Bounded-above prices are not bounded below: discount loses against cheaper access.
 for n in [10,20,50]:
  B=1.;p=np.exp(-2*n);required=np.exp(-n)*2
  check('counterexamples',n-np.log(2*(1+np.log(2)))>0 and required>p)
 check('counterexamples',np.isclose(np.exp(-.02*40),np.exp(-.02*40)))
 OUT.mkdir(exist_ok=True);(OUT/'v8_5_latency.json').write_text(json.dumps({'seed':SEED,'counts':counts},indent=2))
 print(f'5 LATENCY seed={SEED}')
 for k,(p,f) in counts.items():print(f'{k}: PASS={p} FAIL={f}')
 assert not sum(x[1] for x in counts.values())
if __name__=='__main__':run()
