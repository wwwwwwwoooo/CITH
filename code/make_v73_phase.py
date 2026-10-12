"""Render Figure 3 only: exact phase geometry and schematic comparative statics.

No population data are generated or plotted. Numeric arrow endpoints are
computed from explicit illustrative starts and transformation factors.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]

def main() -> None:
    out=ROOT/'results'; figdir=ROOT/'manuscript/figures'
    x=np.linspace(0.0,1.0,501)
    boundary=1.0-x
    productivity_start=np.array([0.17,1.04])
    productivity_factor=3.0
    productivity_end=productivity_start/np.array([1.0,productivity_factor])
    price_start=np.array([0.40,0.24])
    price_factor=3.0
    price_end=price_start*price_factor
    assert np.isclose(productivity_start[0],productivity_end[0])
    assert np.isclose(productivity_end[1],productivity_start[1]/productivity_factor)
    assert np.allclose(price_end,price_factor*price_start)
    assert np.allclose(x+boundary,1.0)
    params={
        'basis':'analytic formula; schematic coordinates, not archived civilization samples',
        'boundary':'xi + zeta = 1',
        'productivity':{'factor':productivity_factor,'start':productivity_start.tolist(),
                        'end':productivity_end.tolist()},
        'price_reduction':{'factor':price_factor,'start':price_start.tolist(),'end':price_end.tolist()},
        'unrecoverable_region':'xi >= 1, A > 0, finite positive productivity',
        'archived_point_overlay':False,
        'population_frequency_estimate':False,
        'boundary_points':len(x),
    }
    np.savez_compressed(out/'v73_phase_diagram.npz',xi=x,boundary_zeta=boundary,
                        productivity_start=productivity_start,productivity_end=productivity_end,
                        productivity_factor=productivity_factor,price_start=price_start,
                        price_end=price_end,price_factor=price_factor)
    matplotlib.rcParams['hatch.linewidth']=0.35
    fig,ax=plt.subplots(figsize=(7.2,4.8),layout='constrained')
    ax.fill_between(x,0.0,boundary,alpha=.16)
    ax.axvspan(1.0,1.62,alpha=.09,hatch='/')
    ax.plot(x,boundary,linewidth=1.8)
    ax.axvline(1.0,linestyle=':',linewidth=1.1)
    # Redundant hatch and direct labels allow monochrome reading.
    for start,end in [(productivity_start,productivity_end),(price_start,price_end)]:
        ax.annotate('',xy=end,xytext=start,
                    arrowprops={'arrowstyle':'-|>','lw':1.8,'mutation_scale':13})
    ax.text(.24,1.19,'(a) Higher productivity',fontsize=10)
    ax.text(.24,1.10,r'$a\mapsto\lambda a:\quad(\xi,\zeta)\mapsto(\xi,\zeta/\lambda)$',fontsize=9.5)
    ax.text(.45,.85,'(b) Lower external price',fontsize=10)
    ax.text(.45,.76,r'$p\mapsto p/\kappa:\quad(\xi,\zeta)\mapsto(\kappa\xi,\kappa\zeta)$',fontsize=9.5)
    ax.text(1.31,1.16,'(c) Not recoverable\nby productivity alone',ha='center',va='center',fontsize=9.5)
    ax.text(.40,.11,'Inward corner',ha='center',fontsize=11)
    ax.text(.40,.57,r'$\xi_j+\zeta_j=1$',fontsize=11,rotation=-30)
    ax.text(1.29,.28,r'$\xi_j\geq1$',ha='center',fontsize=11)
    ax.set(xlim=(0.0,1.62),ylim=(0.0,1.35))
    ax.set_xlabel(r'Residual floor relative to price, $\xi_j=H_j/p_j$',fontsize=10.5)
    ax.set_ylabel(r'Finite-productivity term, $\zeta_j=A_j/(p_j a)$',fontsize=10.5)
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(labelsize=9.5)
    fig.savefig(figdir/'fig3.pdf',bbox_inches='tight',pad_inches=.12,
                metadata={'Title':'CITH v7.3 Figure 3: exact phase geometry'})
    fig.savefig(figdir/'fig3.png',dpi=200,bbox_inches='tight',pad_inches=.12)
    plt.close(fig)
    text=json.dumps(params,indent=2)+'\n'
    (out/'v73_phase_checks.json').write_text(text)
    print(text,end='')

if __name__=='__main__':
    main()
