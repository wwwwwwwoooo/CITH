"""Reconstruct all archived allocation classifications and certificates.

Does not import either of the original generating/verification implementations.
The frozen input array, not an integer seed alone, defines the published design.
"""
from __future__ import annotations
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.stats import qmc

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'/'allocation'
OUT=ROOT/'results'
OUT.mkdir(exist_ok=True)


def sha_array(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a,dtype='<f8').tobytes()).hexdigest()


def threshold(p: np.ndarray, s: dict) -> tuple:
    B,C0,a,rn,rh,AE,AR,AI,e,r,i,*_=p.T
    a=a*float(s['amult'])
    st,sc=float(s['st']),float(s['sc'])
    prices=p[:,11:14]*(float(s['floor'])+1/st+.2/sc)
    remote=np.column_stack([rn,rh,np.zeros_like(B)])*float(s['rmult'])*np.exp(-.1/st)
    slopes=np.column_stack([(1-float(s['q']))*AE/e+float(s['er'])*AR/r+float(s['ei'])*AI/i,
                            float(s['er'])*AR/r+float(s['ei'])*AI/i,float(s['ei'])*AI/i])
    H=(1+B[:,None])*np.log1p(B[:,None])*slopes
    A=(1+B[:,None])*(remote+(1+C0[:,None])*slopes)
    cp=H+A/a[:,None]
    flags=np.all(cp<=prices,axis=1)
    possible=np.all((prices>H)|(A==0),axis=1)
    ac=np.max(np.divide(A,prices-H,out=np.zeros_like(A),where=prices>H),axis=1)
    ac=np.where(possible,ac,np.inf)
    return flags,cp,prices,H,A,ac,remote,slopes,a


def objective_and_gradient(p: np.ndarray, s: dict, v: np.ndarray) -> tuple:
    _,_,prices,_,_,_,remote,_,a=threshold(p,s)
    B,C0,_,_,_,AE,AR,AI,e,r,i,*_=p.T
    z=B*v[:,0]
    x,y,w=(B[:,None]*v[:,1:]/prices).T
    C=C0+a*np.log1p(z)+remote[:,0]*x+remote[:,1]*y
    ee=np.exp(-x/e);er=np.exp(-(x+y)/r);ei=np.exp(-(x+y+w)/i)
    U=np.log1p(C)+(1-s['q'])*AE*(-np.expm1(-x/e)) \
      +s['er']*AR*(-np.expm1(-(x+y)/r))+s['ei']*AI*(-np.expm1(-(x+y+w)/i))
    inv=1/(1+C)
    derivatives=np.column_stack([remote[:,0]*inv+(1-s['q'])*AE/e*ee+s['er']*AR/r*er+s['ei']*AI/i*ei,
                                 remote[:,1]*inv+s['er']*AR/r*er+s['ei']*AI/i*ei,
                                 s['ei']*AI/i*ei])
    grad=np.column_stack([B*inv*a/(1+z),B[:,None]*derivatives/prices])
    return U,grad


def main() -> None:
    design=np.load(DATA/'design.npz',allow_pickle=False)
    p,u,bounds=design['parameters'],design['sobol_vectors'],design['ranges']
    summary=json.loads((DATA/'summary.json').read_text())
    ans={'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
         'input_file_sha256':hashlib.sha256((DATA/'design.npz').read_bytes()).hexdigest(),
         'array_hash_format':'row-major little-endian float64 bytes',
         'sobol_array_sha256':sha_array(u),'parameter_array_sha256':sha_array(p),
         'scenarios':{},'generator_checks':[]}
    total=0;boundary_tests=0;all_gaps=[];all_errors=[]
    for name,entry in summary['scenarios'].items():
        s=entry['settings']; saved=np.load(DATA/(name+'.npz'))
        flags,cp,prices,H,A,ac,*_=threshold(p,s)
        assert np.array_equal(flags,saved['all_local'])
        assert int(flags.sum())==entry['all_local_count']
        v=saved['weights'];U,g=objective_and_gradient(p,s,v)
        assert np.all(v>=0) and np.max(np.abs(v.sum(1)-1))<1e-12
        gaps=np.maximum(0,np.max(g,axis=1)-np.sum(g*v,axis=1))
        ng=gaps/(1+np.abs(U))
        err=float(np.max(np.abs(U-saved['values'])))
        assert np.max(ng)<=2e-6 and err<1e-11
        finite=np.isfinite(ac)&(ac>0)
        aa=ac[finite,None]
        assert np.all(np.any(H[finite]+A[finite]/(aa*(1-1e-6))>prices[finite],axis=1))
        assert np.all(np.all(H[finite]+A[finite]/(aa*(1+1e-6))<prices[finite],axis=1))
        boundary_tests+=2*int(finite.sum())
        rec={'all_inward':int(flags.sum()),'after_price_cut_10pct':int(np.all(cp<=.9*prices,axis=1).sum()),
             'finite_productivity_thresholds':int(np.isfinite(ac).sum()),
             'infeasible_productivity_thresholds':int(np.isinf(ac).sum()),
             'classification_mismatches':int(np.sum(flags!=saved['all_local'])),
             'max_normalized_gap':float(np.max(ng)),'max_value_reconstruction_error':err}
        ans['scenarios'][name]=rec
        total+=len(flags);all_gaps.append(float(np.max(ng)));all_errors.append(err)
    ans['total_archived_solutions']=total
    ans['threshold_sidedness_tests']=boundary_tests
    ans['max_normalized_gap']=max(all_gaps)
    ans['max_value_reconstruction_error']=max(all_errors)
    lo=np.log([v[0] for v in bounds]); hi=np.log([v[1] for v in bounds])
    instant=summary['scenarios']['instantaneous']['settings']
    perfect=summary['scenarios']['instant_perfect_substitution']['settings']
    for label,kwargs in [('seed_integer',{'seed':20260905}),
                         ('rng_integer',{'rng':20260905}),
                         ('rng_generator',{'rng':np.random.default_rng(20260905)})]:
        engine=qmc.Sobol(d=14,scramble=True,**kwargs)
        v=engine.random_base2(m=10)
        draws=np.exp(lo+v*(hi-lo))
        ans['generator_checks'].append({
            'interface':label,'integer_seed':20260905,'d':14,'m':10,
            'scramble':True,'bits':int(engine.bits),'optimization':None,
            'rng_type':type(engine.rng).__name__,'bit_generator':type(engine.rng.bit_generator).__name__,
            'sobol_matches_archive':bool(np.array_equal(v,u)),
            'parameters_match_archive':bool(np.array_equal(draws,p)),
            'max_sobol_abs_difference':float(np.max(np.abs(v-u))),
            'max_parameter_relative_difference':float(np.max(np.abs(draws-p)/p)),
            'sobol_array_sha256':sha_array(v),
            'instantaneous_count':int(threshold(draws,instant)[0].sum()),
            'perfect_substitution_count':int(threshold(draws,perfect)[0].sum())})
    assert ans['generator_checks'][0]['sobol_matches_archive']
    # Floating-point transforms can differ in final bits by layout/platform.
    assert np.allclose(np.exp(lo+u*(hi-lo)),p,rtol=5e-15,atol=0)
    sensitivity=[]
    for measure in ['log_uniform','linear_uniform']:
        for seed in range(20260911,20260919):
            v=qmc.Sobol(d=14,scramble=True,seed=seed).random_base2(m=12)
            draws=np.exp(lo+v*(hi-lo)) if measure=='log_uniform' else bounds[:,0]+v*(bounds[:,1]-bounds[:,0])
            flags=threshold(draws,instant)[0]
            sensitivity.append({'measure':measure,'seed':seed,'n':len(flags),
                                'all_inward_count':int(flags.sum()),'fraction':float(flags.mean())})
    ans['measure_sensitivity']=sensitivity
    old=json.loads((ROOT/'data'/'measure_sensitivity_reference.json').read_text())
    assert sensitivity==old['measure_sensitivity']
    ans['measure_sensitivity_matches_archive']=True
    np.savez_compressed(OUT/'seed_vs_rng_designs.npz',archived_u=u,
                        rng_u=qmc.Sobol(d=14,scramble=True,rng=20260905).random_base2(m=10))
    (OUT/'archive_audit.json').write_text(json.dumps(ans,indent=2,allow_nan=False)+'\n')
    print(json.dumps(ans,indent=2,allow_nan=False))

if __name__=='__main__':main()
