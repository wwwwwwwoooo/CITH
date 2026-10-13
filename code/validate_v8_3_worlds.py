"""Free world-number optimization, logarithmic case, integer correction."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import minimize_scalar
SEED=2026101330
OUT=Path(__file__).resolve().parents[1]/'results'
def run():
 rng=np.random.default_rng(SEED);counts={k:[0,0] for k in ['power_optimizer','shifted_normalization','log_optimizer','integer_worlds','edges']}
 def check(k,c):counts[k][0 if c else 1]+=1
 for t in range(700):
  th=rng.uniform(0,1);al=rng.uniform(.01,1);u=10**rng.uniform(-2,1);C=u*10**rng.uniform(0,5)
  v=lambda logc:th*np.log(C)+(al-th)*logc
  opt=minimize_scalar(lambda x:-v(x),bounds=(np.log(u),np.log(C)),method='bounded')
  best=max(v(np.log(u)),v(np.log(C)),-opt.fun)
  exact=th*np.log(C)+(al-th)*np.log(u) if th>al else al*np.log(C)
  check('power_optimizer',abs(best-exact)<1e-9)
  nu=max(th,al);Cs=C/u-1;baseline=u**al;scaled=np.exp(exact)/baseline
  check('shifted_normalization',np.isclose((scaled-1)/nu,np.expm1(nu*np.log1p(Cs))/nu,rtol=1e-11,atol=1e-10))
  # Optimize the log-welfare objective in log(c/u), including both endpoints.
  th=0 if t%10==0 else th; L=np.log(C/u)
  val=lambda v:np.exp(th*(L-v))*v
  o=minimize_scalar(lambda v:-val(v),bounds=(0,L),method='bounded')
  expected=val(min(L,1/th)) if th>0 else L
  check('log_optimizer',abs(max(val(0),val(L),-o.fun)-expected)<1e-6*(1+abs(expected)))
  Ci=rng.uniform(1,90)*u;n=int(Ci/u);Ws=np.arange(1,n+1);brute=np.max(Ws**th*(Ci/Ws)**al)
  exacti=Ci**al*n**max(th-al,0)
  check('integer_worlds',np.isclose(brute,exacti,rtol=1e-12))
 for th,al in [(0,.5),(1,1),(.5,.5),(1,.1)]:
  check('edges',max(th,al)>0 and 0<=1-max(th,al)<1)
  check('edges',np.isclose(1**th*2**al,2**al))
 OUT.mkdir(exist_ok=True);(OUT/'v8_3_worlds.json').write_text(json.dumps({'seed':SEED,'counts':counts},indent=2))
 print(f'3 WORLD NUMBER seed={SEED}')
 for k,(p,f) in counts.items():print(f'{k}: PASS={p} FAIL={f}')
 print('theta=1, log welfare: exactly linear only for C>=e*u.')
 assert not sum(v[1] for v in counts.values())
if __name__=='__main__':run()
