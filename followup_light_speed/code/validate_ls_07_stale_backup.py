"""Stale geographical backups: v8's endogenous-risk floor with a recency factor psi.
The product objective is differentiated numerically and compared with the closed form;
the recency factor exp(-g*omega*r/c) is checked against an explicit staleness model.
"""
from pathlib import Path
import json
import numpy as np
SEED=2026101407
OUT=Path(__file__).resolve().parents[1]/'results'

def run():
 rng=np.random.default_rng(SEED);counts={k:[0,0] for k in ['derivative','floor_equivalence','hazard_ceiling','staleness_model','monotone_distance']}
 def check(k,c):counts[k][0 if c else 1]+=1
 for _ in range(600):
  B,a,p,r=10**rng.uniform(-1,2),10**rng.uniform(-1,3),10**rng.uniform(-1,1),10**rng.uniform(-3,0)
  eta=rng.uniform(0,1.5);gam=rng.choice([0.,1.]);th=rng.uniform(0,.9);phip=10**rng.uniform(-2,1);psi=rng.uniform(0,1)
  F=(lambda C:np.log1p(C)) if abs(eta-1)<1e-9 else (lambda C:((1+C)**(1-eta)-1)/(1-eta))
  phi=lambda u:1-np.exp(-phip*u)
  C=lambda x:a*np.log1p(B-p*x)+r*x
  U=lambda x:(1-th+th*psi*phi(x/C(x)**gam))*F(C(x))
  h=1e-30;num=U(1j*h).imag/h  # complex-step derivative at the corner
  Om=a*np.log1p(B);Gam=(1+B)*np.log1p(B);K=F(Om)*(1+Om)**eta/Om
  exact=th*psi*phip*Om**-gam*F(Om)+(1-th)*(1+Om)**-eta*(r-a*p/(1+B))
  check('derivative',np.isclose(num,exact,rtol=1e-11,atol=1e-14*(1+abs(exact))))
  floor=Gam*r/Om+th/(1-th)*psi*phip*Gam*Om**-gam*K
  check('floor_equivalence',(exact<=0)==(p>=floor) or abs(p-floor)<1e-9*(1+p))
  if eta<1:
   ceil=p*(1-eta)/(p*(1-eta)+psi*phip*Gam)
   # Limit of the floor as the remote-computation term vanishes and K -> 1/(1-eta), seed backup.
   lim=th/(1-th)*psi*phip*Gam/(1-eta)
   check('hazard_ceiling',(p>=lim)==(th<=ceil) or abs(p-lim)<1e-9*(1+p))
 # Explicit staleness: internal value grows at subjective rate g; backup lags by r/c objective.
 for _ in range(300):
  g,om,r_ly=10**rng.uniform(-9,-3),10**rng.uniform(0,7),10**rng.uniform(-3,1.5)
  V=lambda s:np.exp(g*s);t_sub=rng.uniform(0,1e3);lag=om*r_ly  # c=1 ly/yr, so r/c=r_ly years
  check('staleness_model',np.isclose(V(t_sub-lag)/V(t_sub),np.exp(-g*om*r_ly),rtol=1e-12))
  r2=r_ly*rng.uniform(1.01,10);check('monotone_distance',-g*om*r2<-g*om*r_ly and np.exp(-g*om*r2)<=np.exp(-g*om*r_ly))
 ex={'omega':1e6,'r_ly':4.24,'subjective_lag_years':1e6*4.24}
 OUT.mkdir(exist_ok=True);(OUT/'ls_07_stale_backup.json').write_text(json.dumps({'seed':SEED,'counts':counts,'example':ex},indent=2))
 print(f'LS7 STALE BACKUP seed={SEED}')
 for k,(p_,f_) in counts.items():print(f'{k}: PASS={p_} FAIL={f_}')
 print(f"subjective lag of a 4.24 ly backup at omega=1e6: {ex['subjective_lag_years']:.4g} subjective years")
 assert not sum(x[1] for x in counts.values())
if __name__=='__main__':run()
