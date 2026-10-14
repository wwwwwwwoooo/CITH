"""Centrally evaluated frontier expansion stops at a finite radius; stopping-radius asymptotics.
The general bound uses only stated envelopes; the one-channel benchmark is solved exactly
and checked by direct scalar maximization of v8's one-channel objective with discounting.
Lengths are in a fixed unit (pc in the illustration); all parameters are scenarios.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import brentq,minimize_scalar
SEED=2026101406
OUT=Path(__file__).resolve().parents[1]/'results'

def last_crossing(h,lo,hi,n=6000):
 """Largest R in [lo,hi] with h(R)>=0, for h negative beyond it (log grid, then refine)."""
 R=np.geomspace(lo,hi,n);v=np.array([h(x) for x in R]);idx=np.nonzero(v>=0)[0]
 if idx.size==0:return 0.
 i=idx[-1]
 if i==n-1:return np.inf
 return brentq(h,R[i],R[i+1],xtol=1e-13*R[i],rtol=1e-14)

def run():
 rng=np.random.default_rng(SEED);counts={k:[0,0] for k in ['general_bound','benchmark_exact','asymptotic','no_latency_obstruction','realized_radius']}
 def check(k,c):counts[k][0 if c else 1]+=1
 for _ in range(300):
  kind=rng.choice(['log','lin','pow']);th=rng.uniform(.2,.9)
  if kind=='log':ell,dell,m,KG=(lambda z:np.log1p(z)),(lambda z:1/(1+z)),2.,1/np.e
  elif kind=='lin':ell,dell,m,KG=(lambda z:z),(lambda z:1.+0*z),1.,1.
  else:ell,dell,m,KG=(lambda z,t=th:(1+z)**t-1),(lambda z,t=th:t*(1+z)**(t-1)),1.,1/th
  eta=rng.uniform(0,1);d0=1+10**rng.uniform(-2,1);rbar,Sbar=10**rng.uniform(-2,0),10**rng.uniform(-3.5,-1.5)
  amin,pmin,B0,K,k=10**rng.uniform(-1,1.5),10**rng.uniform(-1,.5),10**rng.uniform(-.5,1.5),10**rng.uniform(-2,1),rng.uniform(1,3)
  kap,ls=rng.uniform(1,3),10**rng.uniform(-1,2)
  K1=(rbar+d0*Sbar)/(amin*ell(B0))+Sbar
  lnQ=lambda R:np.log(K1*KG/pmin)+m*np.log1p(B0+K*(1+R)**k)
  Rbar=last_crossing(lambda R:lnQ(R)-kap*R/ls,1e-6,1e9)
  C0=max(0.,np.log(K1*KG/pmin)+m*np.log1p(B0+K))
  Rhat=last_crossing(lambda R:C0+k*m*np.log1p(R)-kap*R/ls,1e-6,1e9)
  check('general_bound',np.isfinite(Rbar) and Rbar<=Rhat*(1+1e-12)+1e-12)
  for R in Rbar+(1+Rbar)*10**rng.uniform(-6,1,8):
   B=B0+rng.uniform()**3*K*(1+R)**k;a=amin*10**rng.uniform(0,2);rj,Sj=rbar*rng.uniform(),Sbar*rng.uniform()
   Om=a*ell(B);Preq=(rj+Sj*(d0+Om)**eta)/(a*dell(B));D=kap*R/ls*(1+rng.exponential(.2));p=pmin*10**rng.uniform(0,1)
   check('general_bound',np.exp(-D)*Preq<p)
 # One-channel benchmark, log production and value, B(R)=Bc+K R^k, delay exponent kappa R/ls.
 def benchmark(p,S,r,d0,a,Bc,K,k,kap,ls):
  phi=lambda R:np.log((1+Bc+K*R**k)*((r+d0*S)/a+S*np.log1p(Bc+K*R**k))/p)
  return last_crossing(lambda R:phi(R)-kap*R/ls,1e-9,1e22,12000),phi
 for _ in range(60):
  p,S,r,a=10**rng.uniform(-.5,.5),10**rng.uniform(-3.5,-1.5),10**rng.uniform(-2,0),10**rng.uniform(-1,1);d0=1+10**rng.uniform(-2,0)
  Bc,K,k,kap,ls=10**rng.uniform(0,2),10**rng.uniform(-1,1),rng.uniform(1,3),rng.uniform(1,3),10**rng.uniform(0,3)
  Rs,phi=benchmark(p,S,r,d0,a,Bc,K,k,kap,ls)
  if Rs==0:check('benchmark_exact',all(phi(R)<kap*R/ls for R in np.geomspace(1e-9,1e6,200)));continue
  for R,want in [(Rs*(1-1e-7),True),(Rs*(1+1e-7),False),(Rs*2,False),(Rs*10,False)]:
   B=Bc+K*R**k;D=kap*R/ls;f=np.exp(-D)
   obj=lambda x:np.log(d0+a*np.log1p(B-p*x)+r*f*x)+S*f*x
   op=minimize_scalar(lambda x:-obj(x),bounds=(0,B/p),method='bounded',options={'xatol':1e-14*B})
   der=(r*f-a*p/(1+B))/(d0+a*np.log1p(B))+S*f
   check('benchmark_exact',(der>0)==want and ((-op.fun>obj(0)) if want else (-op.fun<=obj(0)+1e-12*(1+abs(obj(0))))))
 # Asymptotics R_stop/(l_s ln l_s) -> k/kappa as l_s grows.
 p,S,r,d0,a,Bc,K,k,kap=1.,.002,.05,1.,5.,.6*106.00128942808259,4.1,3.,2.
 asym=[]
 for e in range(1,15):
  ls=10.**e;Rs,_=benchmark(p,S,r,d0,a,Bc,K,k,kap,ls);asym.append({'l_s':ls,'R_stop':Rs,'ratio':Rs/(ls*np.log(ls))})
 dev=np.array([x['ratio']-k/kap for x in asym]);lsv=np.array([x['l_s'] for x in asym])
 scaled=dev*np.log(lsv)/np.log(np.log(lsv))
 check('asymptotic',np.all(dev>0) and np.all(np.diff(dev[2:])<0))
 # The deviation is O(ln ln l_s / ln l_s), as in the proof.
 check('asymptotic',np.all((scaled[2:]>.5)&(scaled[2:]<3)))
 for x,sc_ in zip(asym,scaled):x['scaled_deviation']=float(sc_)
 # Without latency (D=0), the corner fails at all large radii: v8 growth obstruction.
 for R in 10.**np.arange(1,8):
  B=Bc+K*R**k;check('no_latency_obstruction',(1+B)*((r+d0*S)/a+S*np.log1p(B))>p)
 # Illustration: nearest site at 1.30 pc; realized radius from the home versus R_stop.
 r1=1.30;rows=[]
 _,phi=benchmark(p,S,r,d0,a,Bc,K,k,kap,1.);ls_first=kap*r1/phi(r1)
 for ls in np.geomspace(.03,300,81):
  Rs,phi=benchmark(p,S,r,d0,a,Bc,K,k,kap,ls)
  fails=lambda R:phi(R)>kap*R/ls
  if not fails(r1):Rreal=0.
  else:
   grid=np.geomspace(r1,max(10*Rs,10*r1),20000);ok=np.array([not fails(x) for x in grid]);j=int(np.argmax(ok))
   Rreal=brentq(lambda R:kap*R/ls-phi(R),grid[j-1],grid[j],xtol=1e-14*grid[j],rtol=1e-14) if ok.any() else np.inf
  check('realized_radius',Rreal<=Rs*(1+1e-9)+1e-12)
  check('realized_radius',(Rreal>0)==(ls>ls_first))
  rows.append({'l_s_pc':float(ls),'R_stop_pc':float(Rs),'R_stop_beyond_r1_pc':float(Rs) if Rs>r1 else None,'R_realized_pc':Rreal,'asymptote_pc':float(k/kap*ls*np.log(ls)) if ls>1 else None})
 OUT.mkdir(exist_ok=True)
 (OUT/'ls_06_frontier.json').write_text(json.dumps({'seed':SEED,'counts':counts,'asymptotics':asym,'illustration':{'p':p,'S':S,'r':r,'d0':d0,'a':a,'Bc':Bc,'K_per_pc3':K,'k':k,'kappa':kap,'r1_pc':r1,'l_s_first_step_pc':ls_first,'rows':rows}},indent=2))
 print(f'LS6 FRONTIER seed={SEED}')
 for k_,(p_,f_) in counts.items():print(f'{k_}: PASS={p_} FAIL={f_}')
 for x in asym:print(f"l_s=1e{int(np.log10(x['l_s']))}: R_stop={x['R_stop']:.6g} ratio={x['ratio']:.5f} (limit {k/kap})")
 print(f'home expands only for l_s > {ls_first:.6f} pc')
 for x in rows[::16]:print(f"illustration l_s={x['l_s_pc']:.4g} pc: R_stop={x['R_stop_pc']:.5g} pc realized={x['R_realized_pc']:.5g} pc")
 assert not sum(x[1] for x in counts.values())
if __name__=='__main__':run()
