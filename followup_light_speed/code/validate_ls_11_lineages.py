"""Decentralized lineages: first-step dichotomy and percolation bounds for autonomous colonies.
Monte Carlo on regular trees and on the cubic lattice uses fixed seeds; assertions use
binomial or standard-error margins. These are model checks, not estimates of colonization.
"""
from pathlib import Path
import json
from collections import deque
import numpy as np
from scipy.optimize import brentq
SEED=2026101411
OUT=Path(__file__).resolve().parents[1]/'results'
NB=np.array([(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)])

def tree_lineage(rng,Delta,pi,cap):
 """Settled sites from an expanding root on the Delta-regular tree (no overlaps possible)."""
 settled=1+Delta;front=rng.binomial(Delta,pi)
 while front and settled<cap:
  settled+=front*(Delta-1);front=rng.binomial(front*(Delta-1),pi)
 return settled,front>0

def lattice_lineage(rng,pi,radius,cap=400000,root=True):
 settled={(0,0,0)};q=deque([(0,0,0)]) if root else deque();reach=0
 while q and len(settled)<cap:
  s=q.popleft()
  for d in NB:
   t=(s[0]+d[0],s[1]+d[1],s[2]+d[2])
   if t in settled:continue
   settled.add(t);reach=max(reach,max(abs(t[0]),abs(t[1]),abs(t[2])))
   if reach>=radius:return len(settled),True
   if rng.uniform()<pi:q.append(t)
 return len(settled),False

def run():
 rng=np.random.default_rng(SEED);counts={k:[0,0] for k in ['homogeneous_dichotomy','tree_mean','tree_survival','lattice_bound','lattice_regimes']}
 def check(k,c):counts[k][0 if c else 1]+=1
 res={}
 # Identical deterministic colonies: nothing beyond the home, or the whole box.
 for expands in [False,True]:
  n,hit=lattice_lineage(rng,1. if expands else 0.,12,root=expands)
  check('homogeneous_dichotomy',(n==1 and not hit) if not expands else hit)
 for Delta in [3,4,6,12]:
  for frac in [.3,.6,.9]:
   pi=frac/(Delta-1);m=(Delta-1)*pi;exact=1+Delta+Delta*(Delta-1)*pi/(1-m)
   sizes=np.array([tree_lineage(rng,Delta,pi,10**9)[0] for _ in range(20000)])
   se=sizes.std()/np.sqrt(sizes.size);check('tree_mean',abs(sizes.mean()-exact)<=5*se+1e-9)
   res[f'tree_D{Delta}_m{frac}']={'mean':float(sizes.mean()),'exact':exact}
  for frac in [1.3,2.]:
   pi=min(.95,frac/(Delta-1))
   q=brentq(lambda s:(1-pi+pi*s)**(Delta-1)-s,0,1-1e-12) if (Delta-1)*pi>1 else 1.
   theta=1-(1-pi+pi*q)**Delta;surv=np.mean([tree_lineage(rng,Delta,pi,10**6)[1] for _ in range(4000)])
   check('tree_survival',abs(surv-theta)<=5*np.sqrt(theta*(1-theta)/4000)+.01)
   res[f'tree_D{Delta}_supercritical_{frac}']={'pi':pi,'theta':theta,'mc':float(surv)}
 # Cubic lattice, Delta=6: the tree value bounds the mean when 5*pi<1; regimes either side of p_c.
 for pi in [.1,.15,.19]:
  sizes=np.array([lattice_lineage(rng,pi,10**6)[0] for _ in range(4000)]);bound=1+6+30*pi/(1-5*pi)
  check('lattice_bound',sizes.mean()<=bound);res[f'lattice_pi{pi}']={'mean':float(sizes.mean()),'tree_bound':bound}
 hits_low=np.mean([lattice_lineage(rng,.2,30)[1] for _ in range(400)])
 hits_high=np.mean([lattice_lineage(rng,.45,30)[1] for _ in range(400)])
 check('lattice_regimes',hits_low<.02 and hits_high>.3);res['lattice_reach30']={'pi_0.2':float(hits_low),'pi_0.45':float(hits_high)}
 OUT.mkdir(exist_ok=True);(OUT/'ls_11_lineages.json').write_text(json.dumps({'seed':SEED,'counts':counts,'results':res},indent=2))
 print(f'LS11 LINEAGES seed={SEED}')
 for k,(p_,f_) in counts.items():print(f'{k}: PASS={p_} FAIL={f_}')
 for k,v in res.items():print(k,v)
 assert not sum(x[1] for x in counts.values())
if __name__=='__main__':run()
