"""Light-cone selection: expected visible-but-not-arrived fronts relative to arrived fronts.
Static Euclidean space, homogeneous Poisson origins in space-time, fronts at constant speed v.
Monte Carlo counts are compared with (pi/3)*zeta*v^3*t^4 and (pi/3)*zeta*(c^3-v^3)*t^4.
"""
from pathlib import Path
import json
import numpy as np
SEED=2026101412
OUT=Path(__file__).resolve().parents[1]/'results'

def run():
 rng=np.random.default_rng(SEED);counts={k:[0,0] for k in ['means','independence','ratio','null_conditional']}
 def check(k,c):counts[k][0 if c else 1]+=1
 c,t0=1.,1.;rows=[]
 for beta in [.1,.3,.5,.7,.9,.99]:
  v=beta*c;zeta=60/(np.pi/3*c**3*t0**4);arr=[];vis=[]
  for _ in range(12000):
   n=rng.poisson(zeta*4*np.pi/3*(c*t0)**3*t0);rr=c*t0*rng.uniform(size=n)**(1/3);s=rng.uniform(0,t0,n)
   arr.append(np.sum(rr<=v*(t0-s)));vis.append(np.sum((rr>v*(t0-s))&(rr<=c*(t0-s))))
  arr,vis=np.array(arr),np.array(vis);ma,mv=np.pi/3*zeta*v**3*t0**4,np.pi/3*zeta*(c**3-v**3)*t0**4
  check('means',abs(arr.mean()-ma)<=5*np.sqrt(ma/arr.size)+1e-12 and abs(vis.mean()-mv)<=5*np.sqrt(mv/vis.size))
  cov=np.cov(arr,vis)[0,1];check('independence',abs(cov)<=5*np.sqrt(ma*mv/arr.size)+1e-12)
  R=vis.mean()/arr.mean();seR=R*np.sqrt(1/(arr.size*ma)+1/(vis.size*mv))
  check('ratio',abs(R-((c/v)**3-1))<=5*seR)
  # Independence also means that conditioning on no arrival leaves the visible mean unchanged.
  z=arr==0
  if z.sum()>=200:check('null_conditional',abs(vis[z].mean()-mv)<=5*np.sqrt(mv/z.sum()))
  L=-np.log(.05)
  rows.append({'v_over_c':beta,'visible_to_arrived':(c/v)**3-1,'max_expected_visible_at_P0.05':((c/v)**3-1)*L})
 OUT.mkdir(exist_ok=True);(OUT/'ls_12_visibility.json').write_text(json.dumps({'seed':SEED,'counts':counts,'rows':rows},indent=2))
 print(f'LS12 VISIBILITY seed={SEED}')
 for k,(p_,f_) in counts.items():print(f'{k}: PASS={p_} FAIL={f_}')
 for r in rows:print(f"v/c={r['v_over_c']}: ratio={r['visible_to_arrived']:.6g}, expected visible <= {r['max_expected_visible_at_P0.05']:.6g} at P*=0.05")
 assert not sum(x[1] for x in counts.values())
if __name__=='__main__':run()
