"""Self-contained checks of the single-replicating-channel basin theorem.

Specified benchmark: B_loc(t)=barB, a_t=a0*exp(g_a*t), X_{t+1}=X_t+x_t,
B_t=barB+g_X*X_t. The residual is linear S*x. This is a NEW experiment,
not a reproduction of the unavailable check_trap.py.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq, minimize_scalar

SEED=2026101306
ROOT_TOL=1e-12
BRACKET_TOL=1e-6
OUT=Path(__file__).resolve().parents[1]/'results'

@dataclass(frozen=True)
class P:
    p: float=1.
    r: float=.05
    S: float=.002
    d0: float=1.


def tip(v=P()):return brentq(lambda b:(1+b)*np.log1p(b)-v.p/v.S,0.,v.p/v.S,xtol=ROOT_TOL)

def acrit(B,v=P()):
    den=v.p/(1+B)-v.S*np.log1p(B)
    return np.inf if den<=0 else (v.r+v.d0*v.S)/den

def deriv_sign(x,B,a,v=P()):
    z=max(0.,B-v.p*x)
    # Positive rescaling of the full objective derivative, stable at large a.
    return (v.r+v.S*v.d0+v.S*v.r*x)/a+v.S*np.log1p(z)-v.p/(1+z)

def xstar(B,a,v=P()):
    if deriv_sign(0.,B,a,v)<=0:return 0.
    if deriv_sign(B/v.p,B,a,v)>=0:return B/v.p
    return brentq(lambda x:deriv_sign(x,B,a,v),0.,B/v.p,xtol=ROOT_TOL)

def trajectory(a0,ga,gx,bar,steps=2000,v=P()):
    B=bar;bt=tip(v)
    for t in range(steps):
        loga=np.log(a0)+ga*t
        if loga>650:return dict(regime='UNRESOLVED',t=t,B=B)
        a=np.exp(loga)
        if B>bt+1e-9:return dict(regime='OUT',t=t,B=B,a=float(a))
        x=xstar(B,a,v)
        if x==0 and deriv_sign(0.,B,a,v)<-1e-12:
            return dict(regime='IN',t=t,B=B,a=float(a))
        B+=gx*x
    return dict(regime='UNRESOLVED',t=steps,B=B)

def boundary(ga,gx,bar,v=P()):
    lo=1e-8;hi=acrit(bar,v)*1.01
    assert trajectory(lo,ga,gx,bar,v=v)['regime']=='OUT'
    assert trajectory(hi,ga,gx,bar,v=v)['regime']=='IN'
    while hi-lo>BRACKET_TOL:
        mid=(lo+hi)/2;ans=trajectory(mid,ga,gx,bar,v=v)
        if ans['regime']=='OUT':lo=mid
        elif ans['regime']=='IN':hi=mid
        else:raise RuntimeError('Unresolved trajectory; refusing a forced classification')
    return lo,hi

def run():
    rng=np.random.default_rng(SEED);v=P();bt=tip(v)
    counts={'static_optimizer':[0,0],'comparative_statics':[0,0],
            'runaway_bound':[0,0],'basin_brackets':[0,0],
            'sufficient_certificates':[0,0],'stay_certificate':[0,0],'edge_dynamics':[0,0]}
    def ck(g,ok):counts[g][not bool(ok)]+=1
    for _ in range(200):
        B=rng.uniform(.05,2)*bt;a=10**rng.uniform(-2,3)
        x=xstar(B,a)
        def U(q):return np.log(v.d0+a*np.log1p(B-v.p*q)+v.r*q)+v.S*q
        fit=minimize_scalar(lambda q:-U(q),bounds=(0.,B/v.p),method='bounded',options={'xatol':1e-11})
        best=max(float(U(0)),float(U(B/v.p)),float(-fit.fun))
        ck('static_optimizer',fit.success and abs(float(U(x))-best)<2e-8*(1+abs(best)))
        ck('comparative_statics',xstar(B*1.1,a)>=x-1e-10 and xstar(B,a*1.1)<=x+1e-10)
        B=bt+rng.uniform(.01,2)*bt;a=10**rng.uniform(-2,8)
        x=xstar(B,a)
        ck('runaway_bound',x>(B-bt)/v.p)
    rows=[]
    for ga in (.01,.03,.05,.1,.2):
        lo,hi=boundary(ga,.5,.6*bt)
        low=trajectory(lo,ga,.5,.6*bt);high=trajectory(hi,ga,.5,.6*bt)
        ck('basin_brackets',low['regime']=='OUT' and high['regime']=='IN' and hi-lo<=BRACKET_TOL)
        rows.append(dict(g_a=ga,a_low=lo,a_high=hi,low=low,high=high))
    ck('comparative_statics',all(rows[i+1]['a_high']<rows[i]['a_low'] for i in range(len(rows)-1)))
    # Monotonicity with feedback and with the ceiling ratio.
    gx_bounds=[boundary(.05,gx,.6*bt) for gx in (.1,.5,1.)]
    br_bounds=[boundary(.05,.5,b*bt) for b in (.3,.6,.9)]
    ck('comparative_statics',all(gx_bounds[i][1]<gx_bounds[i+1][0] for i in range(2)))
    ck('comparative_statics',all(br_bounds[i][1]<br_bounds[i+1][0] for i in range(2)))
    # Finite-time sufficient capture certificate B_t <= bar*(1+gX/p)^t.
    for barfrac in (.1,.3,.6):
        for ga in (.1,.3,1.):
            gx=.1;bar=barfrac*bt
            for T in range(4):
                upper=bar*(1+gx/v.p)**T
                if upper>=bt:continue
                a0=1.01*np.exp(-ga*T)*acrit(upper)
                ans=trajectory(a0,ga,gx,bar)
                ck('sufficient_certificates',ans['regime']=='IN' and ans['t']<=T)
    # Exact supremum for a local budget that ramps for N periods, then is constant.
    # After N the budget is constant and exp(-g_a*k) strictly decreases.
    for _ in range(120):
        N=int(rng.integers(1,16));ga=float(rng.uniform(.02,.4));gx=.5
        inherited=float(rng.uniform(0,.1)*bt)
        low=float(rng.uniform(.05,.2)*bt)
        high=float(rng.uniform(.3,.8)*bt)
        frozen=np.linspace(low,high,N+1)+gx*inherited
        needed=np.array([np.exp(-ga*k)*acrit(b) for k,b in enumerate(frozen)])
        threshold=float(needed.max())
        for factor in (.99,1.01):
            a0=factor*threshold;X=inherited;positive=False
            for k in range(N+3):
                bloc=float(np.linspace(low,high,N+1)[min(k,N)])
                B=bloc+gx*X;a=a0*np.exp(ga*k)
                x=xstar(B,a)
                if x>1e-8:positive=True
                X+=x
            ck('stay_certificate',positive == (factor<1))
    scenarios=[]
    for a0,gx,barfrac in ((.5,.5,.6),(50.,.5,.6),(.5,0.,.6),(50.,.5,1.5)):
        ans=trajectory(a0,.05,gx,barfrac*bt)
        scenarios.append(dict(a0=a0,g_X=gx,bar_fraction=barfrac,**ans))
    ck('edge_dynamics',[x['regime'] for x in scenarios]==['OUT','IN','IN','OUT'])
    ck('edge_dynamics',scenarios[2]['t']>0)  # g_X=0 is NOT inward from the start here.
    # g_a=0: below the static threshold crosses, above it remains inward.
    ck('edge_dynamics',trajectory(.8*acrit(.6*bt),0.,.5,.6*bt)['regime']=='OUT')
    ck('edge_dynamics',trajectory(1.2*acrit(.6*bt),0.,.5,.6*bt)['regime']=='IN')
    result=dict(seed=SEED,root_tolerance=ROOT_TOL,bracket_tolerance=BRACKET_TOL,
                counts=counts,Btip=bt,boundaries=rows,scenarios=scenarios,
                feedback_boundaries=gx_bounds,ceiling_boundaries=br_bounds)
    OUT.mkdir(exist_ok=True)
    (OUT/'v8_1e_basin.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'1e BASIN: seed={SEED}; bracket width <= {BRACKET_TOL:g}')
    for g,(passed,failed) in counts.items():print(f'{g}: PASS={passed} FAIL={failed}')
    for row in rows:print(f'g_a={row["g_a"]:.2f}: a_sep in [{row["a_low"]:.9f}, {row["a_high"]:.9f}]')
    for row in scenarios:print(f'a0={row["a0"]:g}, g_X={row["g_X"]:g}, bar/Btip={row["bar_fraction"]:g}: '
                               f'{row["regime"]} at t={row["t"]}, B={row["B"]:.9f}')
    if any(v[1] for v in counts.values()):raise SystemExit(1)

if __name__=='__main__':run()
