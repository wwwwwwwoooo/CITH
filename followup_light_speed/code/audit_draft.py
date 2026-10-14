"""Mechanical audit of the CITH-LS draft: references, LaTeX warnings, and quoted numbers.
Every number below is recomputed from the validator JSON outputs and must appear in the
source with the stated rounding. It checks transcription, not mathematical truth.
"""
from pathlib import Path
import json,re,subprocess
import fitz
ROOT=Path(__file__).resolve().parents[1];RES=ROOT/'results'

def load(name):return json.loads((RES/name).read_text())

def main():
 tex=(ROOT/'CITH_LS_draft.tex').read_text();log=(ROOT/'CITH_LS_draft.log').read_text()
 labels=set(re.findall(r'\\label\{([^}]+)\}',tex))
 refs={x.strip() for m in re.finditer(r'\\(?:[cC]ref|ref|eqref)\*?\{([^}]+)\}',tex) for x in m[1].split(',')}
 cites={x.strip() for m in re.finditer(r'\\(?:[cC]ite[pt]?|citealp)\{([^}]+)\}',tex) for x in m[1].split(',')}
 bib=set(re.findall(r'\\bibitem\[[^\]]*\]\{([^}]+)\}',tex))
 warnings=[s for s in log.splitlines() if re.search(r'Warning|Overfull|Underfull|undefined|multiply defined',s)]
 c2=load('ls_02_capped_basin.json');c3=load('ls_03_cosmology.json');c5=load('ls_05_latency_radius.json')
 c6=load('ls_06_frontier.json');c8=load('ls_08_coherence.json');c9=load('ls_09_clock_coupling.json')
 c10=load('ls_10_containment.json');c11=load('ls_11_lineages.json')['results'];c12=load('ls_12_visibility.json')
 crit={r['T_K']:r for r in c8['solar_critical_clock']};ex={(r['T_K'],r.get('delta_tau')):r for r in c10['examples']}
 thr={r['T_K']:r['delta_tau_for_one_solar_luminosity'] for r in c10['examples'] if 'delta_tau_for_one_solar_luminosity' in r}
 cap={r['e3']:r for r in c2['capped_table']}
 vis={r['v_over_c']:r for r in c12['rows']}
 claims={
  'uncapped separator':(f"{c2['a_sep_uncapped']:.6f}",'6.632159'),
  'capped e3=1e-3':(f"{cap[1e-3]['a_sep_capped']:.6f}",'4.950497'),
  'capped e3=1e-5':(f"{cap[1e-5]['a_sep_capped']:.6f}",'4.334291'),
  'capped t_E 1e-5':(f"{cap[1e-5]['t_E']:.2f}",'161.85'),
  'event horizon Gly':(f"{c3['event_horizon_Gly']:.2f}",'16.68'),
  'event horizon Gpc':(f"{c3['event_horizon_Gpc']:.2f}",'5.11'),
  'dS bound':(f"{c3['dS_bound_Gly']:.2f}",'17.53'),
  'particle horizon':(f"{c3['particle_horizon_Gly']:.1f}",'47.0'),
  'particle horizon radiation':(f"{c3['radiation_sensitivity'][2]['particle_Gly']:.1f}",'46.1'),
  'max ln rho':(f"{c5['ln_rho_max_quantiles']['max']:.3f}",'8.539'),
  'rep max ln rho':(f"{c5['ln_rho_rep_max']:.3f}",'7.950'),
  'design D=0':(f"{c5['all_inward'][0]}/{c5['certified_no_replication'][0]}/{c5['exact_no_replication_common_delay'][0]}",'219/333/421'),
  'first step':(f"{c6['illustration']['l_s_first_step_pc']:.3f}",'7.754'),
  'ratio 1e1':(f"{c6['asymptotics'][0]['ratio']:.2f}",'2.12'),
  'ratio 1e14':(f"{c6['asymptotics'][-1]['ratio']:.2f}",'1.69'),
  'omega_sun 300':(f"{crit[300.0]['omega_crit_tau_1yr']:,.0f}".replace(',','{,}'),'9{,}184'),
  'omega_sun 1000':(f"{crit[1000.0]['omega_crit_tau_1yr']:,.0f}".replace(',','{,}'),'102{,}044'),
  'Lambda_min 300':(f"{crit[300.0]['Lambda_min_au']:.2f}",'3.44'),
  'Lambda_min 1000':(f"{crit[1000.0]['Lambda_min_au']:.3f}",'0.310'),
  'clock band low':(f"{c9['illustration']['curves']['band_alpha_1'][0]:.1f}",'51.5'),
  'clock band high':(f"{c9['illustration']['curves']['band_alpha_1'][1]:,.0f}".replace(',','{,}'),'2{,}042'),
  'B0':(f"{c10['B_no_latency']:.2f}",'96.95'),
  'Blat 300 1e-4':(f"{ex[(300.0,1e-4)]['B_lat']:.4g}",'1003'),
  'Blat 1000 1e-2':(f"{ex[(1000.0,1e-2)]['B_lat']:.4g}",'4.43e+07'),
  'threshold 300':(f"{thr[300.0]:.1e}",'3.3e-06'),
  'threshold 1000':(f"{thr[1000.0]:.1e}",'3.0e-07'),
  'lattice means':(f"{c11['lattice_pi0.1']['mean']:.1f}/{c11['lattice_pi0.15']['mean']:.1f}/{c11['lattice_pi0.19']['mean']:.1f}",'11.8/17.1/25.7'),
  'lattice reach':(f"{100*c11['lattice_reach30']['pi_0.45']:.1f}",'94.5'),
  'visible 0.9c':(f"{vis[0.9]['max_expected_visible_at_P0.05']:.1f}",'1.1'),
  'visible 0.99c':(f"{vis[0.99]['max_expected_visible_at_P0.05']:.3f}",'0.092'),
 }
 failures=[k for k,(computed,quoted) in claims.items() if computed!=quoted]
 # Quoted strings that should literally appear in the source (TeX formatting of the same values).
 literal=['6.632159','4.950497','4.334291','16.68','5.11','17.53','47.0','46.1','8.539','7.950','219','333','421','7.754',
  '2.12','1.69','9{,}184','102{,}044','3.44','0.310','51.5','2{,}042','96.95','1.003\\times10^3','4.430\\times10^7',
  '3.3\\times10^{-6}','3.0\\times10^{-7}','11.8','17.1','25.7','94.5','1.1','0.092']
 missing=[x for x in literal if x not in tex]
 pdf=fitz.open(ROOT/'CITH_LS_draft.pdf')
 def count(text):
  tmp=RES/'_count.tex';tmp.write_text(text)
  try:return subprocess.run(['texcount','-brief',str(tmp)],capture_output=True,text=True).stdout.strip().split(':')[0]
  finally:tmp.unlink()
 body=tex.split('\\section{Introduction}',1)[1].split('\\section*{Declarations}',1)[0]
 abstract=re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}',tex,re.S).group(1)
 words={'main_text':count('\\section{Introduction}'+body),'abstract':count(abstract)}
 out={'pages':pdf.page_count,'labels':len(labels),'undefined_refs':sorted(refs-labels),'unused_bibitems':sorted(bib-cites),
  'missing_bibitems':sorted(cites-bib),'latex_warnings':warnings,'quoted_number_mismatches':failures,'missing_literals':missing,
  'claims':{k:{'computed':v[0],'quoted':v[1]} for k,v in claims.items()},'texcount_text+headers+captions':words}
 (RES/'draft_audit.json').write_text(json.dumps(out,indent=2))
 print(json.dumps({k:v for k,v in out.items() if k!='claims'},indent=1))
 assert not (out['undefined_refs'] or out['unused_bibitems'] or out['missing_bibitems'] or warnings or failures or missing)
if __name__=='__main__':main()
