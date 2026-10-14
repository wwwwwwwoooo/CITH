"""Causally capped inward basin: domination, monotonicity, capture and the v8 example.
The static policy is solved from the full one-channel objective of v8 eq. (basin-policy);
capture and crossing are certified, not inferred from a finite horizon alone.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import brentq,minimize_scalar
SEED=2026101402
OUT=Path(__file__).resolve().parents[1]/'results'

def model(p,S,r,d0):
 Btip=brentq(lambda B:(1+B)*np.log1p(B)-p/S,1e-12,1e12,xtol=1e-14,rtol=1e-15)
 def acrit(B):return (1+B)*(r+d0*S)/(p-(1+B)*np.log1p(B)*S)
 def xstar(B,a):
  T=lambda x:r+S*d0+S*r*x+a*(S*np.log1p(B-p*x)-p/(1+B-p*x))
  if T(0)<=0:return 0.
  if T(B/p)>=0:return B/p
  return brentq(T,0,B/p,xtol=1e-14,rtol=1e-14)
 return Btip,acrit,xstar

def path(a0,ga,gX,Bbar,E,M,Esup=np.inf,tmax=3000):
 """Returns outcome, period and budget history. E(t)=None means uncapped.
 Capture is certified when x=0 and productivity exceeds the threshold at the largest
 budget reachable without further investment, Bbar+min(gX*X,sup E)."""
 Btip,acrit,xstar=M;X=0.;Bs=[]
 for t in range(min(tmax,int(650/ga))):
  a=a0*np.exp(ga*t);acq=gX*X if E is None else min(gX*X,E(t));B=Bbar+acq;Bs.append(B)
  if B>=Btip:return 'cross',t,np.array(Bs)
  x=xstar(B,a);Bev=Bbar+min(gX*X,Esup)
  if x==0 and Bev<Btip and a>=acrit(Bev):return 'capture',t,np.array(Bs)
  X+=x
 return 'open',tmax,np.array(Bs)

def separator(ga,gX,Bbar,E,M,hi):
 lo=1e-9
 if path(lo,ga,gX,Bbar,E,M)[0]!='cross':return 0.
 for _ in range(55):
  m=.5*(lo+hi);o=path(m,ga,gX,Bbar,E,M)[0]
  if o=='cross':lo=m
  else:hi=m
 return .5*(lo+hi)

def run():
 rng=np.random.default_rng(SEED);counts={k:[0,0] for k in ['policy','domination','monotone_time','monotone_a0','capture_inclusion','no_early_crossing','bounded_envelope','lock_in','v8_example']}
 def check(k,c):counts[k][0 if c else 1]+=1
 for _ in range(160):
  p,S,r=10**rng.uniform(-.5,.5),10**rng.uniform(-3.5,-2),10**rng.uniform(-2,0);d0=1+10**rng.uniform(-2,.5)
  M=model(p,S,r,d0);Btip,acrit,xstar=M;Bbar=Btip*rng.uniform(.2,.95)
  ga,gX=10**rng.uniform(-2,-.5),10**rng.uniform(-1,.5);a0=acrit(Bbar)*10**rng.uniform(-1.5,.3)
  # Static policy agrees with direct scalar maximization of the full objective.
  B=Bbar*rng.uniform(.5,1.2);a=a0*rng.uniform(.5,3)
  obj=lambda x:np.log(d0+a*np.log1p(B-p*x)+r*x)+S*x
  op=minimize_scalar(lambda x:-obj(x),bounds=(0,B/p),method='bounded',options={'xatol':1e-12})
  best=max([obj(0),obj(B/p),-op.fun]);check('policy',obj(xstar(B,a))>=best-1e-9*(1+abs(best)))
  k=rng.uniform(1,3);tE=rng.uniform(2,60);e=(Btip-Bbar)/tE**k;E=lambda t,e=e,k=k:e*t**k
  ou,tu,Bu=path(a0,ga,gX,Bbar,None,M);oc,tc,Bc=path(a0,ga,gX,Bbar,E,M)
  n=min(len(Bu),len(Bc));check('domination',np.all(Bc[:n]<=Bu[:n]+1e-10))
  check('monotone_time',np.all(np.diff(Bc)>=-1e-12) and np.all(np.diff(Bu)>=-1e-12))
  a1=a0*rng.uniform(1.01,2);o1,t1,B1=path(a1,ga,gX,Bbar,E,M);n1=min(len(B1),len(Bc))
  check('monotone_a0',np.all(B1[:n1]<=Bc[:n1]+1e-10) and not (oc=='capture' and o1=='cross'))
  check('capture_inclusion',ou!='capture' or oc=='capture')
  check('no_early_crossing',oc!='cross' or tc>=np.ceil(tE)-1e-9)
  # A causal allowance that never reaches the tipping gap captures every path.
  cap=.9*(Btip-Bbar);Eb=lambda t,e=e,k=k,cap=cap:min(e*t**k,cap);ob=path(a0*rng.uniform(.01,1),ga,gX,Bbar,Eb,M,cap)[0]
  check('bounded_envelope',ob=='capture')
  if oc=='cross':
   # Lock-in: after crossing, the static optimizer stays positive at every finite productivity.
   Bx=Bc[-1];ok=all(xstar(Bx*(1+s),a0*np.exp(ga*(tc+s*40))*10**s)>0 for s in np.linspace(0,1,6))
   check('lock_in',ok)
 # The v8 numerical example and its causally capped counterparts.
 p,S,r,d0,ga,gX=1.,.002,.05,1.,.05,.5;M=model(p,S,r,d0);Btip,acrit,xstar=M;Bbar=.6*Btip
 au=separator(ga,gX,Bbar,None,M,acrit(Bbar))
 check('v8_example',abs(Btip-106.001289)<1e-6 and abs(acrit(Bbar)-7.279612)<1e-6 and abs(au-6.63216)<1e-5)
 check('v8_example',path(.5,ga,gX,Bbar,None,M)[:2]==('cross',2) and path(50.,ga,gX,Bbar,None,M)[0]=='capture')
 table=[]
 for e3 in [1e-5,1e-4,1e-3,3e-3,1e-2,3e-2,1e-1,1.]:
  E=lambda t,e3=e3:e3*t**3;ac=separator(ga,gX,Bbar,E,M,acrit(Bbar))
  table.append({'e3':e3,'t_E':((Btip-Bbar)/e3)**(1/3),'a_sep_capped':ac,'outcome_a0_0.5':path(.5,ga,gX,Bbar,E,M)[0]})
 seps=[row['a_sep_capped'] for row in table]
 check('v8_example',all(s<=au+1e-9 for s in seps) and all(np.diff(seps)>=-1e-9))
 OUT.mkdir(exist_ok=True)
 (OUT/'ls_02_capped_basin.json').write_text(json.dumps({'seed':SEED,'counts':counts,'Btip':Btip,'acrit_Bbar':acrit(Bbar),'a_sep_uncapped':au,'capped_table':table},indent=2))
 print(f'LS2 CAPPED BASIN seed={SEED}')
 for k,(p_,f_) in counts.items():print(f'{k}: PASS={p_} FAIL={f_}')
 print(f'uncapped a_sep={au:.6f}')
 for row in table:print(f"e3={row['e3']:g}: t_E={row['t_E']:.4f} capped a_sep={row['a_sep_capped']:.6f} a0=0.5 -> {row['outcome_a0_0.5']}")
 assert not sum(x[1] for x in counts.values())
if __name__=='__main__':run()
