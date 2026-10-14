"""Causal acquisition envelope, finite exponential phases and the truncated v8 runaway.
Inputs are synthetic. Site sets, acquisition schedules and growth laws are drawn
independently of the formulas being checked.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import brentq
SEED=2026101401
OUT=Path(__file__).resolve().parents[1]/'results'

def last_root(f,lo,hi,n=4000):
 """Largest t in [lo,hi] with f(t)>=0 for an f that is negative beyond it."""
 ts=np.linspace(lo,hi,n);v=np.array([f(t) for t in ts]);idx=np.nonzero(v>=0)[0]
 if idx.size==0:return None
 i=idx[-1]
 if i==n-1:return hi
 return brentq(f,ts[i],ts[i+1],xtol=1e-12,rtol=1e-12)

def run():
 rng=np.random.default_rng(SEED);counts={k:[0,0] for k in ['influence_envelope','import_envelope','exponential_phase','explicit_bound','runaway_truncation','edges']}
 def check(k,c):counts[k][0 if c else 1]+=1
 c=1.
 for _ in range(300):
  # Random sites outside a home ball; resources and acquisition times obey causality only.
  rh=rng.uniform(0,2);n=rng.integers(5,400);rad=rh+rng.uniform(0,50,n)**rng.uniform(.5,1.5)
  res=10**rng.uniform(-2,1,n);Bh=rng.uniform(0,10)
  M=lambda r:res[rad<=r].sum()
  t_inf=(rad-rh)/c*(1+rng.exponential(1.,n));t_imp=2*(rad-rh)/c*(1+rng.exponential(1.,n))
  for t in np.linspace(0,120,13):
   check('influence_envelope',Bh+res[t_inf<=t].sum()<=Bh+M(rh+c*t)+1e-12)
   check('import_envelope',Bh+res[t_imp<=t].sum()<=Bh+M(rh+c*t/2)+1e-12)
 for _ in range(400):
  B0,g,K,rh,k=10**rng.uniform(-2,2),10**rng.uniform(-2.5,0),10**rng.uniform(-3,2),rng.uniform(0,10),rng.uniform(.5,4)
  Bh=B0*rng.uniform(0,2);cc=10**rng.uniform(-1,1)
  env=lambda t:Bh+K*(rh+cc*t)**k;f=lambda t:env(t)-B0*np.exp(min(g*t,700.))
  T1=max(rh/cc,(Bh/K)**(1/k)/(2*cc)) if Bh>0 else rh/cc
  T2=2*cc*(2*K/B0)**(1/k)*(2*k/(np.e*g))**2;Tb=max(T1,T2)
  Tc=last_root(f,0,1.05*Tb+1,8000)
  if Tc is None:check('exponential_phase',f(0)<0);check('explicit_bound',True);continue
  later=np.linspace(Tc+1e-9*(1+Tc),Tc+5*(1+Tb),400)[1:]
  check('exponential_phase',np.all(np.array([f(t) for t in later])<0))
  check('explicit_bound',Tc<=Tb*(1+1e-10)+1e-10)
 # v8 basin runaway: once B_t>=B_tip, the excess grows at least geometrically.
 p,S,r,d0,ga,gX=1.,.002,.05,1.,.05,.5
 Btip=brentq(lambda B:(1+B)*np.log1p(B)-p/S,1,1e4,xtol=1e-14);Bbar=.6*Btip
 def xstar(B,a):
  T=lambda x:r+S*d0+S*r*x+a*(S*np.log1p(B-p*x)-p/(1+B-p*x))
  if T(0)<=0:return 0.
  if T(B/p)>=0:return B/p
  return brentq(T,0,B/p,xtol=1e-13,rtol=1e-13)
 for a0 in [.5,1.,2.,4.]:
  X=0.;Bs=[]
  for t in range(60):
   B=Bbar+gX*X;Bs.append(B);X+=xstar(B,a0*np.exp(ga*t))
  Bs=np.array(Bs);t1=int(np.argmax(Bs>=Btip))
  geo=all(Bs[t+1]-Btip>(1+gX/p)*(Bs[t]-Btip)-1e-9 for t in range(t1,59))
  check('runaway_truncation',Bs[t1]>=Btip and geo)
  # A cubic causal envelope with unit-period light radius is overtaken in finite time.
  env=Bbar+1e3*(1+np.arange(60))**3;viol=np.nonzero(Bs>env)[0]
  check('runaway_truncation',viol.size>0 and viol[0]<60)
 check('edges',np.isclose(np.exp(3.)*1.,np.exp(3.)) and np.e*2.5<=np.exp(2.5))
 for y in np.linspace(0,50,101):
  for m in [.5,1,2,6]:check('edges',np.exp(y)>=(np.e*y/m)**m*(1-1e-12))
 OUT.mkdir(exist_ok=True);(OUT/'ls_01_envelope.json').write_text(json.dumps({'seed':SEED,'counts':counts,'Btip':Btip},indent=2))
 print(f'LS1 ENVELOPE seed={SEED}')
 for k,(p_,f_) in counts.items():print(f'{k}: PASS={p_} FAIL={f_}')
 assert not sum(x[1] for x in counts.values())
if __name__=='__main__':run()
