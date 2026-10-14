"""Latency radius: the nearest opportunity binds; frozen-design counts under a common delay.
Reads archive/v7_3/data/allocation/design.npz without changing it. The design rows are
implementation checks, not a prior over civilizations. Delays are expressed through the
dimensionless subjective exponent D; D=kappa*r/l_s for a constant clock.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import brentq,minimize
SEED=2026101405
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parents[1]/'results'

def face_norep(row,D):
 """v8 restricted observation face with every external coefficient discounted by exp(-D)."""
 B,C0,a,rn,rh,AE,AR,AI,se,sr,si,pn,ph,pm=row;d0=1+C0;f=np.exp(-D)
 der=lambda w:.05*AI*f/si*np.exp(-w/si)-pm*a/((1+B-pm*w)*(d0+a*np.log1p(B-pm*w)))
 if der(0)<=0:w=0.
 elif der(B/pm)>=0:w=B/pm
 else:w=brentq(der,0,B/pm,xtol=1e-13,rtol=1e-13)
 z=max(0.,B-pm*w);C=C0+a*np.log1p(z)
 im=.05*AI*f/si*np.exp(-w/si);rm=.01*AR*f/sr;em=.01*AE*f/se
 lam=a/((1+z)*(1+C)) if z>1e-10 else im/pm
 return bool(rn*f/(1+C)+em+rm+im<=pn*lam+1e-11 and rh*f/(1+C)+rm+im<=ph*lam+1e-11)

def run():
 rng=np.random.default_rng(SEED);counts={k:[0,0] for k in ['nearest_binds','radius_formula','design_reproduction','design_nesting','design_monotone']}
 def check(k,c):counts[k][0 if c else 1]+=1
 for _ in range(250):
  B,a,r,S=10**rng.uniform(-1,2),10**rng.uniform(-1,2.5),10**rng.uniform(-2,.5),10**rng.uniform(-3.5,-1);d0=1+10**rng.uniform(-2,1)
  n=rng.integers(2,5);rs=np.sort(10**rng.uniform(-1,1.5,n));ls=10**rng.uniform(-.5,1.5);kap=rng.uniform(1,3);T0=rng.uniform(0,2);dw=10**rng.uniform(-1,0)
  Ds=dw*T0+kap*rs/ls;ps=np.sort(10**rng.uniform(-1,1,n))
  L=np.log1p(B);P=(1+B)*((r+d0*S)/a+S*L)
  def U(x):
   z=max(B-np.dot(ps,x),0.);return np.log(d0+a*np.log1p(z)+np.dot(r*np.exp(-Ds),x))+np.dot(S*np.exp(-Ds),x)
  cons={'type':'ineq','fun':lambda x:B-np.dot(ps,x)-1e-12}
  best=-np.inf
  for start in [np.zeros(n),np.full(n,.2*B/ps.sum()),rng.uniform(0,.5,n)*B/ps/n]:
   op=minimize(lambda x:-U(x),start,method='SLSQP',bounds=[(0,None)]*n,constraints=cons,options={'ftol':1e-13,'maxiter':3000})
   if op.success or np.isfinite(op.fun):best=max(best,-op.fun)
  inward=ps[0]>=np.exp(-Ds[0])*P
  allpass=all(ps[i]>=np.exp(-Ds[i])*P for i in range(n))
  check('nearest_binds',inward==allpass)
  gap=best-U(np.zeros(n))
  check('nearest_binds',gap<=1e-8*(1+abs(U(np.zeros(n)))) if inward else gap>-1e-12)
  rho=P/ps[0];rlat=ls/kap*max(0.,np.log(rho)-dw*T0)
  check('radius_formula',(rs[0]>=rlat)==(Ds[0]>=np.log(rho)) or abs(rs[0]-rlat)<1e-9*(1+rlat))
 P=np.load(ROOT/'archive/v7_3/data/allocation/design.npz')['parameters']
 B,C0,a,rn,rh,AE,AR,AI,se,sr,si,pn,ph,pm=P.T
 CL=C0+a*np.log1p(B);lam=a/((1+B)*(1+CL));em,rm,im=.01*AE/se,.01*AR/sr,.05*AI/si
 rho=np.vstack([(rn/(1+CL)+em+rm+im)/(pn*lam),(rh/(1+CL)+rm+im)/(ph*lam),im/(pm*lam)])
 lmax=np.log(rho.max(0));lrep=np.log(rho[:2].max(0))
 Ds=np.round(np.concatenate([np.arange(0,10.01,.25)]),10)
 allin=[int((lmax<=D).sum()) for D in Ds];cert=[int((lrep<=D).sum()) for D in Ds]
 exact=[];
 for D in Ds:exact.append(int(sum(face_norep(row,D) for row in P)))
 check('design_reproduction',allin[0]==219 and exact[0]==421)
 check('design_reproduction',abs(lmax.max()-8.53937557)<1e-6 and allin[-1]==1024)
 v8=np.load(ROOT/'results/v8_no_replication_classifications.npz')['no_replication']
 check('design_nesting',int(((lrep<=0)&~v8).sum())==0)
 for i,D in enumerate(Ds):
  nr=np.array([face_norep(row,D) for row in P]) if i%8==0 else None
  check('design_nesting',allin[i]<=cert[i]<=exact[i])
  if nr is not None:check('design_nesting',int(((lrep<=D)&~nr).sum())==0 and int(((lmax<=D)&~nr).sum())==0)
 check('design_monotone',all(np.diff(allin)>=0) and all(np.diff(cert)>=0) and all(np.diff(exact)>=0))
 q=lambda frac:float(np.quantile(lmax,frac))
 out={'seed':SEED,'counts':counts,'D_grid':Ds.tolist(),'all_inward':allin,'certified_no_replication':cert,'exact_no_replication_common_delay':exact,
  'ln_rho_max_quantiles':{'min':float(lmax.min()),'q50':q(.5),'q90':q(.9),'q99':q(.99),'max':float(lmax.max())},'ln_rho_rep_max':float(lrep.max())}
 OUT.mkdir(exist_ok=True);(OUT/'ls_05_latency_radius.json').write_text(json.dumps(out,indent=2))
 print(f'LS5 LATENCY RADIUS seed={SEED}')
 for k,(p_,f_) in counts.items():print(f'{k}: PASS={p_} FAIL={f_}')
 for D,x,y,z in zip(Ds,allin,cert,exact):
  if D in (0,1,2,3,4,5,6,7,8,9,10):print(f'D={D:g}: all-inward={x} certified={y} exact no-replication={z}')
 print(f'max ln rho={lmax.max():.8f}; replication-only max={lrep.max():.8f}')
 assert not sum(x[1] for x in counts.values())
if __name__=='__main__':run()
