"""Constructive no-replication persistence despite failure of all-inwardness."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import brentq
SEED=2026101323;OUT=Path(__file__).resolve().parents[1]/'results'
def run():
 rng=np.random.default_rng(SEED);counts={k:[0,0] for k in ['budget','interior_services','strict_replication_test','all_inward_fails','asymptotic_share']}
 def check(k,c):counts[k][0 if c else 1]+=1
 for _ in range(400):
  B=10**rng.uniform(6,9);a=B*B;d0=1+rng.uniform(0,10)
  JI,JR,si,sr,pm,pb=10**rng.uniform(-.5,.5,6);pn=pm+pb+rng.uniform(.1,3);ph=pm+pb+rng.uniform(.1,3);rn,rh=rng.uniform(0,10,2)
  def face(z):
   lam=a/((1+z)*(d0+a*np.log1p(z)))
   w=si*max(0,np.log(JI/(si*pm*lam)));b=sr*max(0,np.log(JR/(sr*pb*lam)))
   return lam,w,b
  z=brentq(lambda z:z+pm*face(z)[1]+pb*face(z)[2]-B,0,B,xtol=1e-8)
  lam,w,b=face(z)
  check('budget',abs(z+pm*w+pb*b-B)<1e-6*(1+B))
  check('interior_services',w>0 and b>0)
  check('strict_replication_test',rn*(1+z)/a<pn-pm-pb and rh*(1+z)/a<ph-pm-pb)
  check('all_inward_fails',(1+B)*np.log1p(B)*(JI/si)>pm)
  check('asymptotic_share',z/B>.999)
 out={'seed':SEED,'counts':counts};OUT.mkdir(exist_ok=True);(OUT/'v8_2c_growing.json').write_text(json.dumps(out,indent=2))
 print(f'2c GROWING NO-REPLICATION seed={SEED}')
 for k,(p,f) in counts.items():print(f'{k}: PASS={p} FAIL={f}')
 assert not sum(v[1] for v in counts.values())
if __name__=='__main__':run()
