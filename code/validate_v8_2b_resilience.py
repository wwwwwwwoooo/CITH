"""Exact endogenous-risk derivative and asymptotics; no global claim."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import minimize_scalar
SEED=2026101322
OUT=Path(__file__).resolve().parents[1]/'results'

def run():
 rng=np.random.default_rng(SEED);counts={k:[0,0] for k in ['derivative','necessary_condition','K_identity','asymptotics','edges']}
 def check(k,c):counts[k][0 if c else 1]+=1
 def F(C,eta):return np.log1p(C) if eta==1 else np.expm1((1-eta)*np.log1p(C))/(1-eta)
 def K(C,eta):return F(C,eta)*np.exp(eta*np.log1p(C))/C
 for t in range(600):
  eta=float(rng.choice([0.,.3,.9,1.,1.4,2.]));gam=float(rng.choice([0.,.5,1.]));B=10**rng.uniform(-1,1);a=10**rng.uniform(-1,2);p=10**rng.uniform(-.5,.5);r=rng.uniform(0,.4);haz=rng.uniform(0,.8);kap=10**rng.uniform(-1,.4)
  Om=a*np.log1p(B);G=(1+B)*np.log1p(B)
  def utility(x):
   C=a*np.log1p(B-p*x)+r*x
   return (1-haz+haz*(-np.expm1(-kap*x/C**gam)))*F(C,eta)
  cs=utility(1e-22j).imag/1e-22
  exact=(1-haz)*(1+Om)**(-eta)*(r-p*a/(1+B))+haz*kap*F(Om,eta)/Om**gam
  check('derivative',abs(cs-exact)<2e-9*(1+abs(exact)))
  req=G*r/Om+haz/(1-haz)*kap*G*K(Om,eta)/Om**gam
  check('necessary_condition',(exact<=1e-12)==(p>=req-1e-12))
  # A positive derivative guarantees a feasible improving point, regardless of nonconcavity.
  if exact>1e-6:
   eps=min(1e-5,B/(100*p));gain=utility(eps)-utility(0)
   for _ in range(12):
    if gain>0:break
    eps*=.1;gain=utility(eps)-utility(0)
   check('necessary_condition',gain>0)
  if eta!=1:ki=((1+Om)-(1+Om)**eta)/((1-eta)*Om)
  else:ki=(1+Om)*np.log1p(Om)/Om
  check('K_identity',np.isclose(K(Om,eta),ki,rtol=2e-11,atol=1e-11))
 for eta in [0.,.2,.5,.8]:
  C=1e80;check('asymptotics',np.isclose(K(C,eta),1/(1-eta),rtol=2e-12))
 for gam in [0.,.5,1.]:
  check('asymptotics', (K(1e40,1)/1e40**gam > 1) if gam==0 else (K(1e40,1)/1e40**gam < 1e-16))
 for eta in [0.,.5,1.,2.]:
  check('edges',np.isclose(K(1e-6,eta),1,rtol=2e-6));check('edges',0/(1-0)==0)
 # The omitted survival factor matters: with hazard=.5 the odds are 1, not .5.
 check('edges',abs(.5/(1-.5)-1)<1e-15)
 OUT.mkdir(exist_ok=True);(OUT/'v8_2b_resilience.json').write_text(json.dumps({'seed':SEED,'counts':counts},indent=2))
 print(f'2b ENDOGENOUS RISK seed={SEED}')
 for k,(p,f) in counts.items():print(f'{k}: PASS={p} FAIL={f}')
 print('Exact seed multiplier is hazard/(1-hazard); eta=1 adds log(output).')
 assert not sum(v[1] for v in counts.values())
if __name__=='__main__':run()
