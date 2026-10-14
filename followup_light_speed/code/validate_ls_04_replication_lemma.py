"""Corner tests in replication directions certify no replication (Lemma, CITH-LS).
The restricted face is solved exactly through its multiplier; the v8 no-replication
inequalities are then evaluated at that face optimum, and a sample is compared with
full constrained optimization. Instances, value curvatures and production shapes are synthetic.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import brentq,minimize
SEED=2026101404
OUT=Path(__file__).resolve().parents[1]/'results'

def production(kind,th):
 if kind=='log':return (lambda z:np.log1p(z)),(lambda z:1/(1+z))
 if kind=='power':return (lambda z:(1+z)**th-1),(lambda z:th*(1+z)**(th-1))
 return (lambda z:-np.expm1(-th*z)/th),(lambda z:np.exp(-th*z))

def run():
 rng=np.random.default_rng(SEED);counts={k:[0,0] for k in ['face_kkt','multiplier_bound','certified_norep','full_optimizer','not_necessary_witness']}
 def check(k,c):counts[k][0 if c else 1]+=1
 witnesses=0
 for it in range(1200):
  kind=rng.choice(['log','power','sat']);th=rng.uniform(.2,.9) if kind=='power' else 10**rng.uniform(-2,-.5)
  ell,dell=production(kind,th);eta=rng.uniform(0,2.5)
  Fp=lambda C:(1+C)**-eta
  F=(lambda C:np.log1p(C)) if abs(eta-1)<1e-12 else (lambda C:((1+C)**(1-eta)-1)/(1-eta))
  B,C0,a=10**rng.uniform(-.5,2),10**rng.uniform(-1,1),10**rng.uniform(0,3)
  rn,rh=10**rng.uniform(-2,1,2);q=rng.uniform(.5,1);eR,eI=10**rng.uniform(-3,-.5,2)
  AE,AR,AI=10**rng.uniform(-1,1,3);sE,sR,sI=10**rng.uniform(0,1.5,3);pm,pb=10**rng.uniform(-2,0,2)
  dE,dR,dI=AE/sE,AR/sR,AI/sI
  CL=C0+a*ell(B);lamL=Fp(CL)*a*dell(B)
  Sn=(1-q)*dE+eR*dR+eI*dI;Sh=eR*dR+eI*dI
  Pn=(Fp(CL)*rn+Sn)/lamL;Ph=(Fp(CL)*rh+Sh)/lamL
  certified=rng.uniform()<.6
  if certified:pn,ph=Pn*(1+10**rng.uniform(-3,0)),Ph*(1+10**rng.uniform(-3,0))
  else:pn,ph=Pn*10**rng.uniform(-1.5,.5),Ph*10**rng.uniform(-1.5,.5)
  # Exact restricted face: observation and backups are optimized through the multiplier.
  def zof(lam):
   g=lambda z:Fp(C0+a*ell(z))*a*dell(z)-lam
   if g(0)<=0:return 0.
   hi=1.
   while g(hi)>0:hi*=2
   return brentq(g,0,hi,xtol=1e-14,rtol=1e-14)
  wof=lambda lam:sI*max(0.,np.log(eI*dI/(pm*lam)))
  bof=lambda lam:sR*max(0.,np.log(eR*dR/(pb*lam)))
  G=lambda ll:zof(np.exp(ll))+pm*wof(np.exp(ll))+pb*bof(np.exp(ll))-B
  lo,hi=np.log(lamL)-1,np.log(lamL)+1
  while G(lo)<0:lo-=2
  while G(hi)>0:hi+=2
  ll=brentq(G,lo,hi,xtol=1e-14,rtol=1e-15);lam=np.exp(ll);z,w,b=zof(lam),wof(lam),bof(lam)
  C=C0+a*ell(z)
  check('face_kkt',abs(z+pm*w+pb*b-B)<=1e-9*(1+B) and (z<=1e-12 or abs(Fp(C)*a*dell(z)-lam)<=1e-9*lam))
  check('multiplier_bound',lam>=lamL*(1-1e-10))
  gx=Fp(C)*rn+(1-q)*dE+eR*dR*np.exp(-b/sR)+eI*dI*np.exp(-w/sI)-pn*lam
  gy=Fp(C)*rh+eR*dR*np.exp(-b/sR)+eI*dI*np.exp(-w/sI)-ph*lam
  if certified:check('certified_norep',gx<0 and gy<0)
  elif gx<0 and gy<0 and (pn<Pn or ph<Ph):witnesses+=1
  if it%12==0:
   price=np.array([pn,ph,pm,pb])
   def U(s):
    zz=B*max(0.,1-s.sum());x,y,ww,bb=B*s/price
    cc=C0+a*ell(zz)+rn*x+rh*y
    return F(cc)+(1-q)*AE*-np.expm1(-x/sE)+eR*AR*-np.expm1(-(x+y+bb)/sR)+eI*AI*-np.expm1(-(x+y+ww)/sI)
   s0=np.array([0.,0.,pm*w/B,pb*b/B]);u0=U(s0);best=-np.inf
   for start in [s0,np.full(4,.2),np.array([.3,.3,.1,.1]),rng.dirichlet(np.ones(5))[:4]]:
    op=minimize(lambda s:-U(s),start,method='SLSQP',bounds=[(0,1)]*4,constraints={'type':'ineq','fun':lambda s:1-s.sum()},options={'ftol':1e-12,'maxiter':2000})
    best=max(best,-op.fun)
   if gx<0 and gy<0:check('full_optimizer',best-u0<=2e-7*(1+abs(u0)))
   else:check('full_optimizer',best-u0>-1e-9)
 check('not_necessary_witness',witnesses>0)
 OUT.mkdir(exist_ok=True);(OUT/'ls_04_replication_lemma.json').write_text(json.dumps({'seed':SEED,'counts':counts,'uncertified_no_replication_witnesses':witnesses},indent=2))
 print(f'LS4 REPLICATION LEMMA seed={SEED}')
 for k,(p_,f_) in counts.items():print(f'{k}: PASS={p_} FAIL={f_}')
 print(f'no-replication instances that fail a replication corner test: {witnesses}')
 assert not sum(x[1] for x in counts.values())
if __name__=='__main__':run()
