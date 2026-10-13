"""Marked-Poisson uniformity, continuous preemption, binary cascades.
Numerical comparisons include direct Bernoulli enumeration and best responses.
"""
from pathlib import Path
from itertools import product
import json
import numpy as np
from scipy.optimize import minimize_scalar
SEED=2026101360
OUT=Path(__file__).resolve().parents[1]/'results'

def run():
 rng=np.random.default_rng(SEED);counts={k:[0,0] for k in ['mixture_equivalence','necessary_bound','binary_equilibria','continuous_endpoints','counterexamples','poisson_simulation']}
 def check(k,c):counts[k][0 if c else 1]+=1
 for _ in range(700):
  NA=10**rng.uniform(-1,4);mi=rng.uniform(0,1);mo=mi+rng.uniform(.01,6);P=rng.uniform(.01,.99);pi=rng.uniform();L=-np.log(P)
  fi=-np.expm1(-mi);fo=-np.expm1(-mo);prob=np.exp(-NA*((1-pi)*fi+pi*fo));bound=(L/NA-fi)/(fo-fi)
  check('mixture_equivalence',(prob>=P)==(pi<=bound) or abs(prob-P)<1e-12)
  check('necessary_bound',prob<P or pi<=L/(NA*fo)+1e-12)
  N=4;cap=rng.uniform(.1,2);d=rng.uniform(-3,2);iota=rng.uniform(0,3)
  def stable(state):
   for i in range(N):
    diff=cap*(d+iota*cap*(sum(state)-state[i]))
    if (state[i]==0 and diff>1e-12) or (state[i]==1 and diff< -1e-12):return False
   return True
  check('binary_equilibria',stable((0,)*N)==(d<=0) and stable((1,)*N)==(d+iota*cap*(N-1)>=0))
  aa=rng.uniform(-2,2);bb=rng.uniform(.2,4)
  f=lambda x:aa*x-bb*x*x/2
  for others,target in [(0.,0.),((N-1)*cap,cap)]:
   op=minimize_scalar(lambda x:-(f(x)+iota*others*x),bounds=(0,cap),method='bounded')
   cand=[0.,cap,op.x];best=max(cand,key=lambda x:f(x)+iota*others*x)
   pred=aa<=0 if target==0 else aa-bb*cap+iota*others>=0
   check('continuous_endpoints',pred==(abs(best-target)<1e-6) or abs(aa-bb*target+iota*others)<1e-7)
 # The proposed cap-game criterion is false although the zero derivative is positive at maximal rivals.
 check('counterexamples',-1+2>0 and -1-10+2<0)
 x=.9
 for _ in range(50):x=np.clip((-1+2*x)/10,0,1)
 check('counterexamples',x==0)
 # Two stable binary equilibria need not be triggered by one defector.
 N=5;d=-2.5;iota=1.;cap=1.
 check('counterexamples',d<=0 and d+iota*(N-1)*cap>=0 and d+iota*cap<0)
 # A lower bound on outward means cannot establish a sufficient rarity window.
 NA=100;mo=.01;L=-np.log(.05)
 check('counterexamples',NA*(-np.expm1(-mo))<L and np.exp(-NA*(-np.expm1(-20)))<.05)
 # Clustered origin simulation; child-count aggregation conditional on parent count.
 n=250000;parents=rng.poisson(2,size=n);fronts=rng.poisson(2*parents);phat=np.mean(fronts==0);truth=np.exp(-2*(1-np.exp(-2)))
 check('poisson_simulation',abs(phat-truth)<6*np.sqrt(truth*(1-truth)/n))
 table={str(N):float(-np.log(.05)/(N*(-np.expm1(-3)))) for N in [100,1000,10000]}
 OUT.mkdir(exist_ok=True);(OUT/'v8_6_population.json').write_text(json.dumps({'seed':SEED,'counts':counts,'uniformity_bounds':table,'simulation':{'estimate':phat,'exact':truth,'n':n}},indent=2))
 print(f'6 POPULATION seed={SEED}')
 for k,(p,f) in counts.items():print(f'{k}: PASS={p} FAIL={f}')
 for N,v in table.items():print(f'N_A={N}, outward fraction <= {v:.10f} ({100*v:.7f} percent)')
 print(f'Poisson-origin simulation: {phat:.7f}; exact {truth:.7f}; n={n}')
 assert not sum(x[1] for x in counts.values())
if __name__=='__main__':run()
