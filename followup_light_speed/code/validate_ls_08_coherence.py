"""Coherence diameter, hull geometry and the integrated heat ceiling.
Constants: SciPy CODATA values; IAU 2015 nominal solar luminosity 3.828e26 W; IAU au and pc;
Earth mean radius 6371.0084 km (JPL Planetary Physical Parameters). Tolerances, clock speeds,
radiator temperatures and flux limits are scenario values, not engineering specifications.
"""
from pathlib import Path
import json
import numpy as np
from scipy import constants as sc
from scipy.spatial import ConvexHull,cKDTree
SEED=2026101408
OUT=Path(__file__).resolve().parents[1]/'results'
YR=365.25*86400;AU=149597870700.;PC=648000/np.pi*AU;LSUN=3.828e26;REARTH=6371.0084e3

def diameter(pts):
 h=ConvexHull(pts);v=pts[h.vertices];d=np.sqrt(((v[:,None,:]-v[None,:,:])**2).sum(-1));return d.max(),h

def run():
 rng=np.random.default_rng(SEED);counts={k:[0,0] for k in ['hull_area_bound','sphere_limit','cauchy_formula','planar_isodiametric','landauer_carnot','scale_identity','table_values']}
 def check(k,c):counts[k][0 if c else 1]+=1
 shapes=['cube','gauss','needle','disk','shell','cluster']
 for it in range(600):
  s=shapes[it%6];n=rng.integers(8,400)
  if s=='cube':P=rng.uniform(-1,1,(n,3))
  elif s=='gauss':P=rng.normal(size=(n,3))*rng.uniform(.1,3,3)
  elif s=='needle':P=rng.normal(size=(n,3))*np.array([5,.05,.05])
  elif s=='disk':P=rng.normal(size=(n,3))*np.array([2,2,.02])
  elif s=='shell':P=rng.normal(size=(n,3));P/=np.linalg.norm(P,axis=1)[:,None]
  else:P=np.concatenate([rng.normal(size=(n//2,3))*.1,rng.normal(size=(n-n//2,3))*.1+rng.uniform(-3,3,3)])
  D,h=diameter(P);check('hull_area_bound',h.area<=np.pi*D*D*(1+1e-12))
  if it%20==0:
   us=rng.normal(size=(3000,3));us/=np.linalg.norm(us,axis=1)[:,None];areas=[]
   for u in us:
    e1=np.cross(u,[1.,0,0] if abs(u[0])<.9 else [0,1.,0]);e1/=np.linalg.norm(e1);e2=np.cross(u,e1)
    Q=np.c_[P@e1,P@e2]
    try:areas.append(ConvexHull(Q).volume)
    except Exception:areas.append(0.)
   m,se=np.mean(areas),np.std(areas)/np.sqrt(len(areas))
   check('cauchy_formula',abs(4*m-h.area)<=5*4*se+1e-9*h.area)
 for _ in range(300):
  n=rng.integers(5,300);Q=rng.normal(size=(n,2))*rng.uniform(.1,3,2)
  hq=ConvexHull(Q);v=Q[hq.vertices];d=np.sqrt(((v[:,None,:]-v[None,:,:])**2).sum(-1)).max()
  check('planar_isodiametric',hq.volume<=np.pi*d*d/4*(1+1e-12))
 for n in [500,5000,40000]:
  # Points on a centred sphere: the farthest partner of p is the point nearest to -p.
  P=rng.normal(size=(n,3));P/=2*np.linalg.norm(P,axis=1)[:,None];h=ConvexHull(P)
  _,j=cKDTree(P).query(-P);D=np.linalg.norm(P-P[j],axis=1).max()
  check('sphere_limit',h.area<=np.pi*D*D and (n<40000 or h.area>.995*np.pi*D*D))
 for _ in range(500):
  Tr=10**rng.uniform(.5,3.5);Tc=Tr*10**rng.uniform(-2,1)
  heat=sc.k*Tc*np.log(2)*(Tr/Tc if Tc<Tr else 1.)
  check('landauer_carnot',heat>=sc.k*Tr*np.log(2)*(1-1e-12))
  dl,tau,om=10**rng.uniform(-9,0),10**rng.uniform(-2,8),10**rng.uniform(0,8)
  ls=sc.c/(dl*om/YR);Lam=sc.c*tau/(2*om);check('scale_identity',np.isclose(ls/Lam,2/(dl*tau/YR),rtol=1e-12))
 table=[]
 for om in [1.,1e3,1e6]:
  for tau,name in [(.1,'0.1 s'),(86400.,'1 day'),(YR,'1 Julian year')]:
   Lam=sc.c*tau/(2*om);P300=np.pi*Lam**2*sc.sigma*300**4;P1000=np.pi*Lam**2*sc.sigma*1000**4
   table.append({'omega':om,'tau':name,'Lambda_m':Lam,'Lambda_au':Lam/AU,'Pmax_300K_W':P300,'Pmax_1000K_W':P1000,'erasure_300K_per_s':np.pi*sc.sigma*300**3*Lam**2/(sc.k*np.log(2))})
 crit=[]
 for T in [300.,1000.]:
  Lmin=np.sqrt(LSUN/(np.pi*sc.sigma*T**4));crit.append({'T_K':T,'Lambda_min_m':Lmin,'Lambda_min_au':Lmin/AU,'omega_crit_tau_1yr':sc.c*YR/(2*Lmin)})
 em=sc.c*.1/(2*1e3);check('table_values',abs(em-14989.6229)<1e-3)
 check('table_values',np.isclose(sc.c/1e6,299.792458) and 3.44<crit[0]['Lambda_min_au']<3.45 and 9180<crit[0]['omega_crit_tau_1yr']<9190)
 subj=[{'target':name,'distance_m':d,'one_way_subjective_years_at_1e6':d/sc.c*1e6/YR} for name,d in [('Earth mean diameter',2*REARTH),('1 au',AU),('1.30 pc',1.30*PC)]]
 check('table_values',15.8<subj[1]['one_way_subjective_years_at_1e6']<15.82 and 4.239e6<subj[2]['one_way_subjective_years_at_1e6']<4.241e6)
 scales=[{'omega':om,'subjective_light_speed_m_s':sc.c/om,'au_one_way_subjective_s':AU/sc.c*om,'nearest_one_way_subjective_yr':1.30*PC/sc.c*om/YR,
  'earth_diameter_one_way_subjective_s':2*REARTH/sc.c*om} for om in [1.,1e3,1e6]]
 check('table_values',np.isclose(scales[0]['au_one_way_subjective_s'],499.00478383615643,rtol=1e-12))
 OUT.mkdir(exist_ok=True)
 (OUT/'ls_08_coherence.json').write_text(json.dumps({'seed':SEED,'counts':counts,'coherence_table':table,'solar_critical_clock':crit,'subjective_light_times':subj,'subjective_scales':scales,'kilo_clock_reaction_diameter_m':em,'subjective_light_speed_at_1e6_m_per_s':sc.c/1e6},indent=2))
 print(f'LS8 COHERENCE seed={SEED}')
 for k,(p_,f_) in counts.items():print(f'{k}: PASS={p_} FAIL={f_}')
 for row in table:print(f"omega={row['omega']:g}, tau={row['tau']}: Lambda={row['Lambda_m']:.4g} m, Pmax(300 K)={row['Pmax_300K_W']:.4g} W")
 for row in crit:print(f"T={row['T_K']:g} K: Lambda_min(L_sun)={row['Lambda_min_au']:.4f} au; omega_crit(1 yr)={row['omega_crit_tau_1yr']:.1f}")
 for row in subj:print(f"{row['target']}: {row['one_way_subjective_years_at_1e6']:.6g} subjective years one way at omega=1e6")
 assert not sum(x[1] for x in counts.values())
if __name__=='__main__':run()
