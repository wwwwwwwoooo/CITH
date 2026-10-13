"""Physical units and thermal illustration, not an engineering feasibility test.
L_sun: IAU 2015 nominal. Masses: JPL Planetary Physical Parameters, 2026-10-13.
Sigma=1 kg/m^2, horizon=Julian year, clock speeds and returns are assumptions.
"""
from pathlib import Path
import json
import numpy as np
from scipy import constants as sc
from scipy.optimize import brentq
SEED=2026101370
OUT=Path(__file__).resolve().parents[1]/'results'

def run():
 rng=np.random.default_rng(SEED);counts={k:[0,0] for k in ['radiator_root','dimensional_scaling','discount_threshold','edges']}
 def check(k,c):counts[k][0 if c else 1]+=1
 for _ in range(500):
  L,M,Sigma=10**rng.uniform(20,28),10**rng.uniform(20,29),10**rng.uniform(-2,2)
  T=(L*Sigma/(sc.sigma*M))**.25
  root=brentq(lambda t:sc.sigma*(M/Sigma)*t**4/L-1,1e-4,1e6,xtol=1e-10)
  check('radiator_root',np.isclose(root,T,rtol=1e-10))
  scale=10**rng.uniform(-3,3)
  check('dimensional_scaling',np.isclose((L*scale*Sigma/(sc.sigma*M*scale))**.25,T))
  ratio=10**rng.uniform(0,50);omega=10**rng.uniform(0,8);lag=10**rng.uniform(-1,3);delta=np.log(ratio)/(omega*lag)
  check('discount_threshold',np.isclose(np.log(ratio)-delta*omega*lag,0,atol=1e-12))
 solar=3.828e26;rad={}
 for name,mass in [('Mercury',.330103e24),('Jupiter',1898.125e24)]:
  T=(solar/(sc.sigma*mass))**.25;rad[name]={'mass_kg':mass,'temperature_K':T,'Wien_peak_micron':sc.Wien/T*1e6}
 check('edges',rad['Mercury']['temperature_K']>100 and rad['Jupiter']['temperature_K']<100)
 check('edges',1/(sc.k*300*np.log(2))>0)
 year=365.25*86400;budget=solar*year/1e32
 table=[]
 for ratio in [1e3,1e10,1e20,1e40]:
  for omega in [1.,1e3,1e6]:table.append({'ratio':ratio,'omega':omega,'lag_year':40,'delta_per_subjective_year':np.log(ratio)/(omega*40)})
 OUT.mkdir(exist_ok=True);(OUT/'v8_7_physical.json').write_text(json.dumps({'seed':SEED,'counts':counts,'radiators':rad,'solar_budget_Eref1e32':budget,'discount_table':table,'sigma_SB':sc.sigma},indent=2))
 print(f'7 PHYSICAL seed={SEED}')
 for k,(p,f) in counts.items():print(f'{k}: PASS={p} FAIL={f}')
 print(f'Nominal solar-year energy / 1e32 J = {budget:.9f}')
 for k,v in rad.items():print(f'{k}: T={v["temperature_K"]:.6f} K; peak={v["Wien_peak_micron"]:.6f} micron')
 print(f'Illustrative delta = {np.log(1e20)/(1e6*40):.12g} per subjective year')
 assert not sum(x[1] for x in counts.values())
if __name__=='__main__':run()
