"""Integration-compatible latency containment and the containment budget B_lat.
Units: budget in b_ref=1e32 J per Julian-year horizon (v8 dictionary); tau in subjective
seconds; delta per subjective year; radiative flux limit sigma*T^4. Parameters are scenarios.
Brute-force clock scans are compared with the endpoint criterion; c is varied to test that it
cancels at the fastest admissible clock.
"""
from pathlib import Path
import json
import numpy as np
from scipy import constants as sc
from scipy.optimize import brentq
SEED=2026101410
OUT=Path(__file__).resolve().parents[1]/'results'
YR=365.25*86400;AU=149597870700.;PC=648000/np.pi*AU;BREF=1e32;LSUN=3.828e26

def Lmin(B,T):return np.sqrt(B*BREF/YR/(np.pi*sc.sigma*T**4))
def rho(B,p,S,r,d0,a):return (1+B)*((r+d0*S)/a+S*np.log1p(B))/p
def lhs(B,T,r1,kap,dtau):return kap*dtau/2*r1/Lmin(B,T)   # dtau = delta*tau in units of (1/yr)*(s)/YR

def Blat(T,r1,kap,dtau,par):
 f=lambda lb:lhs(np.exp(lb),T,r1,kap,dtau)-np.log(rho(np.exp(lb),*par))
 lo,hi=0.,0.
 while f(lo)<=0:lo-=10.
 while f(hi)>=0:hi+=10.
 return np.exp(brentq(f,lo,hi,xtol=1e-13,rtol=1e-14))

def run():
 rng=np.random.default_rng(SEED);counts={k:[0,0] for k in ['criterion_vs_scan','c_cancels','unique_root','monotone','no_latency_floor','endpoint_test','examples']}
 def check(k,c):counts[k][0 if c else 1]+=1
 for _ in range(400):
  par=(10**rng.uniform(-.5,.5),10**rng.uniform(-3.5,-1.5),10**rng.uniform(-2,0),1+10**rng.uniform(-2,0),10**rng.uniform(-.5,2.5))
  T,r1,kap=10**rng.uniform(2,3.3),PC*10**rng.uniform(-.3,.7),rng.uniform(1,3)
  tau,delta=10**rng.uniform(-3,8),10**rng.uniform(-10,0);B=10**rng.uniform(-2,9);dtau=delta*tau/YR
  crit=lhs(B,T,r1,kap,dtau)>=np.log(rho(B,*par))
  for cfac in [1.,.37,4.2]:
   c=sc.c*cfac;wmax=c*tau/(2*Lmin(B,T))
   W=np.append(np.geomspace(wmax*1e-8,wmax,3000),wmax)
   D=delta*W*kap*r1/(c*YR);scan=bool(np.any(np.exp(-D)*rho(B,*par)<=1))
   near=abs(lhs(B,T,r1,kap,dtau)-np.log(rho(B,*par)))<1e-9*(1+abs(np.log(rho(B,*par))))
   check('criterion_vs_scan' if cfac==1. else 'c_cancels',scan==crit or near)
  bl=Blat(T,r1,kap,dtau,par);G=np.geomspace(bl*1e-4,bl*1e4,400)
  fv=lhs(G,T,r1,kap,dtau)-np.log(rho(G,*par))
  check('unique_root',np.all(np.diff(fv)<0) and np.count_nonzero(np.diff(np.sign(fv))!=0)==1)
  for name,args in [('r1',(T,r1*1.5,kap,dtau)),('kappa',(T,r1,kap*1.3,dtau)),('dtau',(T,r1,kap,dtau*3)),('T',(T*1.4,r1,kap,dtau))]:
   check('monotone',Blat(*args,par)>bl)
  b0=np.exp(brentq(lambda lb:np.log(rho(np.exp(lb),*par)),np.log(1e-12),np.log(1e40))) if np.log(rho(1e-12,*par))<0 else 0.
  check('no_latency_floor',bl>=b0*(1-1e-12))
  # Coupled efficiency a=Xi*omega^-alpha: the admissible minimum of g is attained at an endpoint.
  alpha=rng.uniform(0,1);Xi=par[4]*10**rng.uniform(0,3);wlo=10**rng.uniform(-3,0);wmax=sc.c*tau/(2*Lmin(B,T))
  if wmax>wlo:
   p,S,r,d0,_=par;H=(1+B)*np.log1p(B)*S;A=(1+B)*(r+d0*S);beta=delta*kap*r1/(sc.c*YR)
   g=lambda w:np.exp(-beta*w)*(H+A/Xi*w**alpha)
   W=np.geomspace(wlo,wmax,3000);scan=bool(np.any(g(W)<=p));ep=min(g(wlo),g(wmax))<=p
   check('endpoint_test',scan==ep or abs(min(g(wlo),g(wmax))-p)<1e-9*p)
 # Worked examples with the v8 one-channel coefficients.
 par=(1.,.002,.05,1.,50.);r1=1.30*PC;kap=2.;table=[]
 for T in [300.,1000.]:
  for dtau in [1e-8,1e-6,1e-4,1e-2,1.]:
   bl=Blat(T,r1,kap,dtau,par);table.append({'T_K':T,'delta_tau':dtau,'B_lat':float(bl),'B_lat_solar_years':float(bl/(LSUN*YR/BREF)),'Lambda_min_au':float(Lmin(bl,T)/AU)})
  thr=10**brentq(lambda ld:Blat(T,r1,kap,10**ld,par)-LSUN*YR/BREF,-12,2,xtol=1e-12)
  table.append({'T_K':T,'delta_tau_for_one_solar_luminosity':float(thr)})
 B0=np.exp(brentq(lambda lb:np.log(rho(np.exp(lb),*par)),np.log(1e-6),np.log(1e6)))
 check('examples',all(row['B_lat']>=B0 for row in table if 'B_lat' in row))
 check('examples',all(np.diff([row['B_lat'] for row in table if row.get('T_K')==300. and 'B_lat' in row])>0))
 OUT.mkdir(exist_ok=True)
 (OUT/'ls_10_containment.json').write_text(json.dumps({'seed':SEED,'counts':counts,'coefficients':{'p':1.,'S':.002,'r':.05,'d0':1.,'a':50.},'r1_pc':1.30,'kappa':kap,'B_no_latency':B0,'examples':table},indent=2))
 print(f'LS10 CONTAINMENT seed={SEED}')
 for k,(p_,f_) in counts.items():print(f'{k}: PASS={p_} FAIL={f_}')
 print(f'no-latency budget B0={B0:.6f}')
 for row in table:print({k:float(v) for k,v in row.items()})
 assert not sum(x[1] for x in counts.values())
if __name__=='__main__':run()
