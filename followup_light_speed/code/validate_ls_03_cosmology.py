"""Finite affectable comoving region in flat Lambda-CDM.
Parameters: Planck 2018 base-LCDM late-universe values H0=67.4 km/s/Mpc, Omega_m=0.315,
flat with Omega_L=1-Omega_m; radiation is neglected except in the stated sensitivity check.
Lengths use the Julian light-year and the IAU 2015 parsec. Quadrature is checked against
closed forms; this is a kinematic bound, not a model of expansion strategy.
"""
from pathlib import Path
import json
import numpy as np
from scipy import constants as sc
from scipy.integrate import quad
SEED=2026101403
OUT=Path(__file__).resolve().parents[1]/'results'
YR=365.25*86400;LY=sc.c*YR;AU=149597870700.;PC=648000/np.pi*AU

def horizons(H0,Om,Or=0.,ae=1.):
 OL=1-Om-Or;DH=sc.c/(H0*1e3/(1e6*PC))
 E=lambda a:np.sqrt(Or*a**-4+Om*a**-3+OL)
 eh=DH*quad(lambda a:1/(a*a*E(a)),ae,np.inf,epsabs=0,epsrel=1e-12,limit=400)[0]
 ph=DH*quad(lambda a:1/(a*a*E(a)),0,1,epsabs=0,epsrel=1e-12,limit=400)[0]
 return DH,eh,ph,E

def run():
 rng=np.random.default_rng(SEED);counts={k:[0,0] for k in ['closed_forms','upper_bound','monotone_emission','probes','planck_numbers','radiation_sensitivity']}
 def check(k,c):counts[k][0 if c else 1]+=1
 # de Sitter: chi_EH(a_e)=c/(H0 a_e). Einstein-de Sitter: particle horizon 2c/H0, no event horizon.
 for ae in [1.,2.,5.]:
  DH=sc.c/(70e3/(1e6*PC));eh=DH*quad(lambda a:1/(a*a),ae,np.inf,epsrel=1e-12)[0];check('closed_forms',np.isclose(eh,DH/ae,rtol=1e-10))
 DH=sc.c/(70e3/(1e6*PC));E=lambda a:a**-1.5
 ph=DH*quad(lambda a:1/(a*a*E(a)),0,1,epsrel=1e-12)[0];check('closed_forms',np.isclose(ph,2*DH,rtol=1e-10))
 grow=[DH*quad(lambda a:1/(a*a*E(a)),1,A,epsrel=1e-12)[0] for A in [10.,1e3,1e5]]
 check('closed_forms',np.allclose(grow,[2*DH*(np.sqrt(A)-1) for A in [10.,1e3,1e5]],rtol=1e-10))
 for _ in range(200):
  H0,Om=rng.uniform(50,90),rng.uniform(.05,.9);DH,eh,ph,E=horizons(H0,Om)
  check('upper_bound',eh<=DH/np.sqrt(1-Om)*(1+1e-12))
  aes=np.sort(rng.uniform(1,20,5));ehs=[horizons(H0,Om,0.,ae)[1] for ae in aes]
  check('monotone_emission',np.all(np.diff(ehs)<0))
  # Probes with conserved comoving momentum u=gamma*beta*a never reach the light horizon.
  for u0 in [.01,.3,3.,30.]:
   chi=DH*quad(lambda a:(u0/a)/np.sqrt(1+(u0/a)**2)/(a*a*E(a)),1,np.inf,epsrel=1e-11,limit=400)[0]
   check('probes',0<chi<eh)
 DH,eh,ph,_=horizons(67.4,.315)
 frac=(eh/ph)**3;check('planck_numbers',eh<ph)
 check('planck_numbers',16.6<eh/LY/1e9<16.8 and 46.9<ph/LY/1e9<47.1 and .044<frac<.046)
 check('planck_numbers',np.isclose(DH/LY/1e9,14.5073,rtol=1e-4))
 sens=[]
 for Or in [0.,5e-5,9.2e-5,1.2e-4]:
  _,e2,p2,_=horizons(67.4,.315,Or);sens.append({'Omega_r':Or,'event_Gly':e2/LY/1e9,'particle_Gly':p2/LY/1e9,'volume_fraction':(e2/p2)**3})
 check('radiation_sensitivity',max(abs(s['event_Gly']-eh/LY/1e9) for s in sens)<2e-3)
 check('radiation_sensitivity',all(.044<s['volume_fraction']<.048 for s in sens))
 out={'seed':SEED,'counts':counts,'H0':67.4,'Omega_m':.315,'Hubble_distance_Gly':DH/LY/1e9,'event_horizon_Gly':eh/LY/1e9,'event_horizon_Gpc':eh/PC/1e9,'particle_horizon_Gly':ph/LY/1e9,'affectable_volume_fraction':frac,'radiation_sensitivity':sens,'dS_bound_Gly':DH/np.sqrt(.685)/LY/1e9}
 OUT.mkdir(exist_ok=True);(OUT/'ls_03_cosmology.json').write_text(json.dumps(out,indent=2))
 print(f'LS3 COSMOLOGY seed={SEED}')
 for k,(p_,f_) in counts.items():print(f'{k}: PASS={p_} FAIL={f_}')
 print(f'event horizon={eh/LY/1e9:.4f} Gly ({eh/PC/1e9:.4f} Gpc); particle horizon={ph/LY/1e9:.4f} Gly; fraction={frac:.5f}; bound={DH/np.sqrt(.685)/LY/1e9:.4f} Gly')
 for s in sens:print(f"Omega_r={s['Omega_r']:g}: event={s['event_Gly']:.4f} particle={s['particle_Gly']:.4f} fraction={s['volume_fraction']:.5f}")
 assert not sum(x[1] for x in counts.values())
if __name__=='__main__':run()
