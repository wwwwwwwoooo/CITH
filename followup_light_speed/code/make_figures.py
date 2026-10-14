"""Draw the four CITH-LS figures from validator outputs and closed-form curves.
Run after code/run_all_checks.py. Colours are the first three slots of a validated
categorical palette; reference lines use ink tones, never a series colour.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import constants as sc
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[1];RES=ROOT/'results';FIG=ROOT/'figures'
S1,S2,S3='#2a78d6','#eb6834','#1baf7a'
INK,INK2,MUTED,GRID,AXIS='#0b0b0b','#52514e','#898781','#e1e0d9','#c3c2b7'
YR=365.25*86400;AU=149597870700.;PC=648000/np.pi*AU;LSUN=3.828e26;REARTH=6371.0084e3;BREF=1e32
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.edgecolor':AXIS,'axes.linewidth':.6,
 'axes.labelcolor':INK,'xtick.color':INK2,'ytick.color':INK2,'xtick.major.width':.6,'ytick.major.width':.6,
 'xtick.minor.width':.4,'ytick.minor.width':.4,'axes.grid':True,'grid.color':GRID,'grid.linewidth':.5,'grid.linestyle':'-',
 'legend.frameon':False,'legend.fontsize':8,'lines.linewidth':1.5,'lines.solid_capstyle':'round','pdf.fonttype':42,
 'axes.spines.top':False,'axes.spines.right':False})

def load(name):return json.loads((RES/name).read_text())

def endlabel(ax,x,y,text,dy=0,ha='left'):
 ax.annotate(text,(x,y),xytext=(4,dy),textcoords='offset points',va='center',ha=ha,color=INK2,fontsize=7.5)

def refline(ax,y,text,xpos,ha='right'):
 ax.axhline(y,color=MUTED,lw=.6,zorder=1);ax.annotate(text,(xpos,y),xytext=(0,3),textcoords='offset points',color=MUTED,fontsize=7,ha=ha,va='bottom')

def fig1():
 om=np.geomspace(1,1e8,400);fig,axs=plt.subplots(1,2,figsize=(7.0,2.9),constrained_layout=True)
 ax=axs[0]
 for (tau,lab),col in zip([(.1,r'$\tau=0.1$ s'),(86400.,r'$\tau=1$ day'),(YR,r'$\tau=1$ yr')],[S1,S2,S3]):
  L=sc.c*tau/(2*om);ax.loglog(om,L,color=col,label=lab);endlabel(ax,om[-1],L[-1],{.1:'0.1 s',86400.:'1 day',YR:'1 yr'}[tau])
 for y,t,x in [(2*REARTH,'Earth diameter',3e3),(AU,'1 au',3e6),(1.30*PC,'1.30 pc',3e6)]:refline(ax,y,t,x,ha='center')
 ax.set_xlabel(r'Subjective clock $\omega$ (subjective yr per yr)');ax.set_ylabel(r'Coherence diameter $\Lambda=c\tau/(2\omega)$ (m)')
 ax.set_xlim(1,1e8);ax.set_ylim(1e-1,1e17);ax.legend(loc='lower left');ax.set_title('(a) Integration region',loc='left',fontsize=9,color=INK)
 ax=axs[1]
 for T,col in [(300.,S1),(1000.,S2)]:
  P=np.pi*(sc.c*YR/(2*om))**2*sc.sigma*T**4;ax.loglog(om,P,color=col,label=f'$T_{{\\rm rad}}={T:g}$ K')
  wc=sc.c*YR/(2*np.sqrt(LSUN/(np.pi*sc.sigma*T**4)));ax.plot([wc],[LSUN],'o',ms=5.5,color=col,mec='white',mew=1.2,zorder=5)
  endlabel(ax,wc,LSUN,rf'$\omega_\odot={wc:,.0f}$'.replace(',','{,}'),dy=-9 if T==300. else 9,ha='right')
 refline(ax,LSUN,r'$L_\odot^{\rm N}$',3e7)
 ax.set_xlabel(r'Subjective clock $\omega$ (subjective yr per yr)');ax.set_ylabel(r'Integrated power ceiling $\pi\Lambda^2\sigma T^4$ (W)')
 ax.set_xlim(1,1e8);ax.legend(loc='lower left');ax.set_title(r'(b) Heat ceiling, $\tau=1$ yr',loc='left',fontsize=9,color=INK)
 fig.savefig(FIG/'ls_fig1_scales.pdf');plt.close(fig)

def fig2():
 d=load('ls_05_latency_radius.json');D=np.array(d['D_grid'])
 fig,ax=plt.subplots(figsize=(4.6,3.0),constrained_layout=True)
 for key,lab,col in [('all_inward','All-inward corner',S1),('certified_no_replication','No replication, certified by replication tests',S2),('exact_no_replication_common_delay','No replication, exact face test',S3)]:
  y=np.array(d[key]);ax.step(D,y,where='post',color=col,label=lab);endlabel(ax,0,y[0],f'{y[0]}',ha='right',dy=0)
 ax.set_xlim(-.6,10);ax.set_ylim(0,1080);ax.axhline(1024,color=MUTED,lw=.6)
 ax.annotate('1,024 rows',(10,1024),xytext=(0,3),textcoords='offset points',ha='right',va='bottom',color=MUTED,fontsize=7)
 ax.set_xlabel(r'Common subjective delay exponent $D=\delta\omega T$')
 ax.set_ylabel('Frozen design rows (of 1,024)');ax.legend(loc='lower right')
 fig.savefig(FIG/'ls_fig2_design.pdf');plt.close(fig)

def fig3():
 d=load('ls_06_frontier.json')['illustration'];rows=d['rows']
 ls=np.array([r['l_s_pc'] for r in rows]);Rb=np.array([np.nan if r['R_stop_beyond_r1_pc'] is None else r['R_stop_beyond_r1_pc'] for r in rows]);Rr=np.array([r['R_realized_pc'] for r in rows])
 fig,ax=plt.subplots(figsize=(4.6,3.0),constrained_layout=True)
 ax.loglog(ls,Rb,color=S1,label=r'$R_{\rm stop}$: largest radius $\geq r_1$ where expansion pays')
 m=Rr>0;ax.loglog(ls[m][::2],Rr[m][::2],'o',ms=4.5,color=S2,mec='white',mew=1.,label='Radius reached from the home',zorder=4)
 g=ls>3;ax.loglog(ls[g],1.5*ls[g]*np.log(ls[g]),color=INK2,lw=.9,ls=(0,(4,3)),label=r'$(k/\kappa)\,\ell_s\ln\ell_s$, $k=3$, $\kappa=2$')
 first=d['l_s_first_step_pc'];ax.axvline(first,color=MUTED,lw=.6)
 ax.annotate(f'first step pays only\nfor $\\ell_s>{first:.2f}$ pc',(first,.25),xytext=(-4,0),textcoords='offset points',ha='right',va='bottom',color=MUTED,fontsize=7)
 ax.axhline(d['r1_pc'],color=MUTED,lw=.6);ax.annotate('nearest site, $r_1=1.30$ pc',(ls.max(),d['r1_pc']),xytext=(0,-3),textcoords='offset points',ha='right',va='top',color=MUTED,fontsize=7)
 ax.set_xlim(.3,ls.max());ax.set_ylim(.2,2e4)
 ax.set_xlabel(r'Latency length $\ell_s=c/(\delta\omega)$ (pc)');ax.set_ylabel('Frontier radius (pc)');ax.legend(loc='upper left',fontsize=7)
 fig.savefig(FIG/'ls_fig3_frontier.pdf');plt.close(fig)

def fig4():
 d=load('ls_09_clock_coupling.json')['illustration'];fig,axs=plt.subplots(1,2,figsize=(7.0,2.9),constrained_layout=True)
 ax=axs[0]
 for a,col in [('0.0',S1),('0.5',S2),('1.0',S3)]:
  c=d['curves'][a];ax.semilogx(c['omega'],np.array(c['g'])/d['p'],color=col,label=rf'$\alpha={float(a):g}$')
 lo,hi=d['curves']['band_alpha_1'];ax.axvspan(lo,hi,color=S3,alpha=.08,lw=0)
 ax.annotate('corner fails\n'+r'for $\alpha=1$',(lo*1.15,1.45),ha='left',color=INK2,fontsize=7)
 ax.axhline(1,color=MUTED,lw=.6);ax.annotate('price $p$',(1.5e4,1),xytext=(0,3),textcoords='offset points',ha='right',va='bottom',color=MUTED,fontsize=7)
 ax.set_xlim(1e-1,2e4);ax.set_ylim(0,1.7);ax.set_xlabel(r'Subjective clock $\omega$');ax.set_ylabel(r'Discounted required price $g(\omega)/p$')
 ax.legend(loc='lower left');ax.set_title(r'(a) Efficiency $a=\Xi\omega^{-\alpha}$',loc='left',fontsize=9,color=INK)
 c=load('ls_10_containment.json');par=c['coefficients'];r1=c['r1_pc']*PC;kap=c['kappa'];ax=axs[1]
 def rho(B):return (1+B)*((par['r']+par['d0']*par['S'])/par['a']+par['S']*np.log1p(B))/par['p']
 def Lmin(B,T):return np.sqrt(B*BREF/YR/(np.pi*sc.sigma*T**4))
 def blat(T,dt):
  f=lambda lb:kap*dt/2*r1/Lmin(np.exp(lb),T)-np.log(rho(np.exp(lb)));lo,hi=0.,0.
  while f(lo)<=0:lo-=10.
  while f(hi)>=0:hi+=10.
  return np.exp(brentq(f,lo,hi,xtol=1e-12))
 dts=np.geomspace(1e-9,1,200);sol=LSUN*YR/BREF
 for T,col in [(300.,S1),(1000.,S2)]:
  y=np.array([blat(T,x) for x in dts])/sol;ax.loglog(dts,y,color=col,label=f'$T_{{\\rm rad}}={T:g}$ K')
 ax.axhline(c['B_no_latency']/sol,color=MUTED,lw=.6);ax.annotate(f"no-latency budget, {c['B_no_latency']/sol:.2f}",(1,c['B_no_latency']/sol),xytext=(0,-3),textcoords='offset points',ha='right',va='top',color=MUTED,fontsize=7)
 ax.axhline(1,color=MUTED,lw=.6);ax.annotate(r'one $L_\odot^{\rm N}$ for one year',(1,1),xytext=(0,3),textcoords='offset points',ha='right',va='bottom',color=MUTED,fontsize=7)
 ax.set_ylim(.25,1e10)
 ax.set_xlabel(r'Subjective constant $\delta\tau$');ax.set_ylabel(r'$B_{\rm lat}$ (solar luminosity-years)')
 ax.legend(loc='upper left');ax.set_title(r'(b) Containment budget, $r_1=1.30$ pc',loc='left',fontsize=9,color=INK)
 fig.savefig(FIG/'ls_fig4_clock_containment.pdf');plt.close(fig)

if __name__=='__main__':
 FIG.mkdir(exist_ok=True);fig1();fig2();fig3();fig4();print('figures written to',FIG)
