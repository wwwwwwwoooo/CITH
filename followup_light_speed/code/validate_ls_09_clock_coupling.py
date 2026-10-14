"""Clock speed versus efficiency: discounted required price g(omega)=exp(-beta*omega)(H+A'omega^alpha).
Checks unimodality, the two-interval corner set, the slow-clock tipping limit, the adiabatic
peak formula and agreement with direct maximization of the one-channel objective.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import brentq,minimize_scalar
SEED=2026101409
OUT=Path(__file__).resolve().parents[1]/'results'

def run():
 rng=np.random.default_rng(SEED);counts={k:[0,0] for k in ['unimodal','corner_set','direct_optimizer','slow_limit','adiabatic_peak','no_coupling']}
 def check(k,c):counts[k][0 if c else 1]+=1
 for _ in range(500):
  B,p,r,S,Xi=10**rng.uniform(-.5,2),10**rng.uniform(-.5,.5),10**rng.uniform(-2,0),10**rng.uniform(-3.5,-1),10**rng.uniform(-1,2)
  d0=1+10**rng.uniform(-2,0);alpha=rng.choice([0.,rng.uniform(0,1),1.]);beta=10**rng.uniform(-3,0);wlo=10**rng.uniform(-4,-1)
  L=np.log1p(B);H=(1+B)*L*S;A=(1+B)*(r+d0*S);Ap=A/Xi
  g=lambda w:np.exp(-beta*w)*(H+Ap*w**alpha)
  h=lambda w:alpha*Ap*w**(alpha-1)-beta*(H+Ap*w**alpha)
  W=np.geomspace(wlo,60/beta+1e3*wlo,4000);hv=h(W);gv=g(W)
  check('unimodal',np.all(np.diff(hv)<=1e-12*np.abs(hv[:-1])+1e-300) and np.count_nonzero(np.diff(np.sign(hv))!=0)<=1)
  ok=gv<=p;ch=np.count_nonzero(np.diff(ok.astype(int))!=0)
  # Allowed patterns: all pass; pass-fail-pass; fail-pass; pass-fail-pass with first part empty.
  check('corner_set',ok[-1] and (ch==0 or (ch==1 and not ok[0]) or (ch==2 and ok[0])))
  for w in rng.choice(W,3):
   a=Xi*w**-alpha;f=np.exp(-beta*w)
   obj=lambda x:np.log(d0+a*np.log1p(B-p*x)+r*f*x)+S*f*x
   op=minimize_scalar(lambda x:-obj(x),bounds=(0,B/p),method='bounded',options={'xatol':1e-13*B})
   gain=max(-op.fun,obj(B/p))-obj(0)
   check('direct_optimizer',(gain<=1e-10*(1+abs(obj(0)))) if g(w)<=p*(1-1e-9) else (gain>0 if g(w)>=p*(1+1e-9) else True))
  if alpha>0:check('slow_limit',abs(g(1e-12)-H)<=1e-6*(H+Ap)+Ap*1e-12**alpha*2)
  Btip=brentq(lambda b:(1+b)*np.log1p(b)*S-p,1e-12,1e15)
  check('slow_limit',(p>H)==(B<Btip) or abs(B-Btip)<1e-9*Btip)
  if alpha==1.:
   wh=1/beta-H/Ap
   if wh>0:
    gmax=Ap/beta*np.exp(beta*H/Ap-1);check('adiabatic_peak',np.isclose(g(wh),gmax,rtol=1e-11) and gmax>=gv.max()*(1-1e-9))
   else:check('adiabatic_peak',np.all(np.diff(gv)<0))
  if alpha==0.:check('no_coupling',np.all(np.diff(gv)<0))
 # Illustration used in the figure: v8 one-channel parameters at B=0.9*B_tip.
 p,S,r,d0=1.,.002,.05,1.;Btip=brentq(lambda b:(1+b)*np.log1p(b)*S-p,1,1e4);B=.9*Btip;L=np.log1p(B)
 H=(1+B)*L*S;A=(1+B)*(r+d0*S);beta=1e-3;Xi=1.5e3;rows={}
 for alpha in [0.,.5,1.]:
  W=np.geomspace(1e-1,2e4,600);rows[str(alpha)]={'omega':W.tolist(),'g':(np.exp(-beta*W)*(H+A/Xi*W**alpha)).tolist()}
 g1=lambda w:np.exp(-beta*w)*(H+A/Xi*w)-p;wh=1/beta-H/(A/Xi)
 band=[brentq(g1,1e-9,wh,xtol=1e-12),brentq(g1,wh,1e6,xtol=1e-12)];rows['band_alpha_1']=band
 check('adiabatic_peak',band[0]<wh<band[1] and g1(wh)>0)
 OUT.mkdir(exist_ok=True)
 (OUT/'ls_09_clock_coupling.json').write_text(json.dumps({'seed':SEED,'counts':counts,'illustration':{'B':B,'Btip':Btip,'H':H,'A':A,'beta':beta,'Xi':Xi,'p':p,'curves':rows}},indent=2))
 print(f'LS9 CLOCK COUPLING seed={SEED}')
 for k,(p_,f_) in counts.items():print(f'{k}: PASS={p_} FAIL={f_}')
 print(f'illustration: H={H:.6f} (p=1), B/B_tip=0.9, beta={beta}, Xi={Xi}; alpha=1 fails on ({band[0]:.4f}, {band[1]:.4f})')
 assert not sum(x[1] for x in counts.values())
if __name__=='__main__':run()
