"""New no-replication classification on unchanged archived reference inputs.
Original arrays and original 219 all-inward classifications are never rewritten.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import brentq,minimize
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results'

def run():
 P=np.load(ROOT/'archive/v7_3/data/allocation/design.npz')['parameters']
 counts={'archived_all_inward':[0,0],'face_optimizer':[0,0],'full_optimizer':[0,0]};nr=[];allin=[];ws=[]
 def check(k,c):counts[k][0 if c else 1]+=1
 for rownum,row in enumerate(P):
  B,C0,a,rn,rh,AE,AR,AI,se,sr,si,pn,ph,pm=row;d0=1+C0
  def obs_der(w):
   z=max(B-pm*w,0);return .05*AI/si*np.exp(-w/si)-pm*a/((1+z)*(d0+a*np.log1p(z)))
  if obs_der(0)<=0:w=0.
  elif obs_der(B/pm)>=0:w=B/pm
  else:w=brentq(obs_der,0,B/pm,xtol=1e-12)
  z=max(0,B-pm*w);C=C0+a*np.log1p(z)
  lm=a/((1+z)*(1+C)) if z>1e-10 else .05*AI/si*np.exp(-w/si)/pm
  im=.05*AI/si*np.exp(-w/si);rm=.01*AR/sr;em=.01*AE/se
  passnr=(rn/(1+C)+em+rm+im<=pn*lm+1e-11 and rh/(1+C)+rm+im<=ph*lm+1e-11)
  nr.append(bool(passnr));ws.append(w)
  CL=C0+a*np.log1p(B);lam=a/((1+B)*(1+CL));ss=np.array([em+rm+.05*AI/si,rm+.05*AI/si,.05*AI/si]);rr=np.array([rn,rh,0.]);price=np.array([pn,ph,pm])
  ai=bool(np.all(rr/(1+CL)+ss<=price*lam));allin.append(ai)
  check('archived_all_inward',not ai or passnr)
  check('face_optimizer',(w>1e-10 or obs_der(0)<=1e-9) and (w<B/pm-1e-10 or obs_der(B/pm)>=-1e-9) and (w<=1e-10 or w>=B/pm-1e-10 or abs(obs_der(w))<1e-8))
  if rownum%16==0:
   def utility(s):
    z=B*max(0,1-s.sum());x,y,ww=B*s/price
    c=C0+a*np.log1p(z)+rn*x+rh*y
    return np.log1p(c)+.01*AE*(-np.expm1(-x/se))+.01*AR*(-np.expm1(-(x+y)/sr))+.05*AI*(-np.expm1(-(x+y+ww)/si))
   s0=np.array([0.,0.,pm*w/B]);u0=utility(s0)
   opt=minimize(lambda s:-utility(s),s0,method='SLSQP',bounds=[(0,1)]*3,constraints={'type':'ineq','fun':lambda s:1-s.sum()},options={'ftol':1e-11,'maxiter':1500})
   gap=-opt.fun-u0
   check('full_optimizer',gap<2e-7 if passnr else gap>0)
 # Count reproductions, not a new prior over actual civilizations.
 check('archived_all_inward',sum(allin)==219)
 out={'counts':counts,'n':len(P),'all_inward':sum(allin),'no_replication':sum(nr),'observation_only_no_replication':int(sum(np.array(nr)&~np.array(allin))),'replication_preferred':len(P)-sum(nr)}
 OUT.mkdir(exist_ok=True);(OUT/'v8_6_design.json').write_text(json.dumps(out,indent=2));np.savez(OUT/'v8_no_replication_classifications.npz',no_replication=nr,all_inward=allin,observation=ws)
 print('6 DESIGN additive classifications on frozen v7.3 reference')
 for k,(p,f) in counts.items():print(f'{k}: PASS={p} FAIL={f}')
 for k,v in out.items():
  if k!='counts':print(f'{k}={v}')
 assert not sum(v[1] for v in counts.values())
if __name__=='__main__':run()
