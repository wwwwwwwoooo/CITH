"""Assemble one canonical LaTeX source without altering archived numerical inputs."""
from pathlib import Path
import re,json,hashlib
R=Path(__file__).resolve().parents[1]
O=(R/'archive/v7_3/manuscript/body.tex').read_text()
FULL=(R/'archive/CITH_V7_3_full.tex').read_text()
S=R/'latex'
def take(start,end):return FULL[FULL.index(start):FULL.index(end,FULL.index(start))]
def stripcomments(t):return re.sub(r'^%.*\n','',t,flags=re.M)
# Matched original pieces; descriptive wording changes only where numbering is now structural.
phase=take(r'\subsection{Analytical phase boundary}',r'\subsection{Paired sensitivity design}')
phase=phase.replace(r'\includegraphics[width=.91\textwidth]{fig3.pdf}',r'\includegraphics[width=.91\textwidth]{fig3.pdf}')
# Keep original figure number 3; other inherited illustrations are supplementary.
phase=phase.replace(r'\begin{figure}[!t]',r'\setcounter{figure}{2}'+'\n'+r'\begin{figure}[!t]')
rho=take('A dimensionless robustness diagnostic is',r'\paragraph{Curvature generalization.}')
(S/'phase.tex').write_text(rho+phase)
workload=take('Four controlled examples use',r'\section{Allocation between inward growth and external projects}')
joint=take(r'\subsection{A joint finite-interval witness}',r'\subsection{Budget growth and shared technological progress}')
pure=take(r'\subsection{Budget growth and shared technological progress}',r'\paragraph{Persistence when budgets grow.}')
# Equations with short direct proofs already covered in Stage 1 appendix. Use supplement for full original plot.
dyn=take(r'\section{Future resources can reverse a present inward preference}',r'\section{From an inward civilization to a quiet galaxy}')
dyn=dyn.replace(r'\section{Future resources can reverse a present inward preference}',r'\section{Archived two-period construction}')
dyn=dyn.replace(r'\label{sec:dynamic}',r'\label{supp:dynamic}')
dyn=re.sub(r'\\begin\{proposition\}\[A future-adjusted inward condition\].*?\\end\{proof\}',r'''The archived specialization of \\cref{prop:dynamic} is
\\begin{equation}
 \\frac{r_1-p_1a_1/(1+B_1)}{1+C_{01}+a_1\\ln(1+B_1)}
 +\\beta V_2'(B_{20}+g_LB_1)(g_X-p_1g_L)\\leq0.
 \\label{eq:dynamiccriterion}
\\end{equation}
The derivative proof in the main text applies because both current and
continuation objectives here are concave.
''',dyn,flags=re.S)
# Replace obsolete statements that there is no substantive infinite-horizon theorem.
pos=dyn.find('Thus static inwardness is neither sufficient nor necessary')
if pos>=0:dyn=dyn[:pos]+'''The two reversals are retained as finite-horizon examples. The infinite-horizon
benchmark in \\cref{prop:infinite} now supplies a separate verification theorem,
with its own production and growth assumptions. The nonlinear examples here
are exact special cases of \\cref{prop:dynamic}, not of the linear benchmark.
'''
paired=take(r'\subsection{Paired sensitivity design}',r'\section{Future resources can reverse a present inward preference}')
num=take(r'\section{Numerical specification and paired results}',r'% END expanded body.tex')
num=num.replace(r'\section{Numerical specification and paired results}',r'\subsection{Frozen design and retained paired results}').replace(r'\label{app:numerical}',r'\label{supp:legacy-numerical}')
# Fix symbol conflicts in displayed legacy benefit family, while preserving arrays and parameter order.
paired=paired.replace('e^{-x/e}','\\exp(-x/s_E)').replace('e^{-u/r}','\\exp(-u/s_R)').replace('e^{-v/i}','\\exp(-v/s_I)')
num=num.replace('A_I,e,r,i,p_{n0}', 'A_I,s_E,s_R,s_I,p_{n0}').replace('Values of $e,r,i$ in this appendix are return scales, not energy, radius, or individual identifiers.','The stored columns named e, r, and i are written $s_E,s_R,s_I$ here; their values and ordering are unchanged.')
num=num.replace('$e,r$ &','$s_E,s_R$ &').replace('$i$ &','$s_I$ &')
# Retain archival numerical assertions but clarify that additional historical environments were not rerun now.
num=num.replace('under SciPy~1.17.0 and 1.18.1;', 'in the historical checks under SciPy~1.17.0 and 1.18.1;')
num=num.replace('Supplementary Material section~S2 records the generator specification\nand file checksums.', 'The frozen package and the accompanying integrity manifest record the generator and checksums. The present rerun uses SciPy~1.17.0; the 1.18.1 result is retained archival evidence.')
num=num.replace('in Supplementary Material section~S1','in the archived Supplementary Material and audit outputs')
paired=paired.replace('Table~\\ref{tab:scenarios}',r'\Cref{tab:scenarios}')
supp=(r'''\clearpage
\typeout{CITH-SUPPLEMENT-START}
\label{page:supplement-start}
\part*{Supplementary Material}
\suppressfloats[t]
\setcounter{section}{0}
\renewcommand{\thesection}{S\arabic{section}}
\renewcommand{\theHsection}{supp.\arabic{section}}
\renewcommand{\theHfigure}{supp.\arabic{figure}}
\renewcommand{\theHtable}{supp.\arabic{table}}
\setcounter{figure}{0}\renewcommand{\thefigure}{S\arabic{figure}}
\setcounter{table}{0}\renewcommand{\thetable}{S\arabic{table}}
\section{Retained workload and budget illustrations}
\label{supp:workload}
The numerical examples below are retained from v7.3. Their figures and frozen
arrays are unchanged; only their placement and notation have been updated.
'''+workload+joint+pure+r'\section{Paired sensitivity design and numerical specification}'+'\n'+r'\label{supp:paired}\label{sec:access}'+'\n'+paired+num+dyn)
# Recursion labels retained from archived relations, now supplementary identities.
supp+=r'''
\section{Additional identities and new verification files}
\label{supp:newchecks}
The workload ratios used in the main text are, explicitly,
\begin{equation}
 \frac{D_k(d+1)}{D_k(d)}=b+\frac1{W_d}\longrightarrow b,
 \label{eq:recursion-cost}
\end{equation}
\begin{equation}
 F(C_{d+1})-F(C_d)=\ln\frac{1+c(bW_d+1)}{1+cW_d}\longrightarrow\ln b.
 \label{eq:recursion-value}
\end{equation}
Their proofs divide by $W_d\to\infty$. A nonconcave example retained from
v7.3 is
\begin{equation}
 U(x)=\ln\{1+2\ln(2-x)+0.1x\}+3x^2,\quad 0\leq x\leq1,
 \label{eq:nonconcave}
\end{equation}
with $U'(0)\simeq-0.3772$, $U(0)\simeq0.8697$, and
$U(1)\simeq3.0953$. It violates the concavity assumption, not the corner theorem.

Each new script in \texttt{code/validate\_v8\_*.py} is self-contained and
prints explicit pass and fail counts. The accompanying README maps files
to propositions and distinguishes assertions from independent randomized
cases; each script records its seed and tolerances. The new no-replication classification reads
\texttt{archive/v7\_3/data/allocation/design.npz} without changing it.
It solves observation on the restricted face, then evaluates the two
replication reduced derivatives; 64 of the 1,024 cases are additionally
compared with full constrained optimization. General restricted-face tests
use 320 separate randomized instances with optional backups. Thus the full
optimization check is sampled, whereas every frozen row receives the exact
face criterion.

The original archive scripts have been rerun in a working copy, not in the
immutable archive. Their outputs are included under
\texttt{results/archive\_rerun/}. No figures or archived classification
counts have been revised. Mechanical source and PDF checks are recorded in
\texttt{results/assembly\_audit.json} and
\texttt{results/final\_audit.json}; physical-parameter sources are documented
in the comments of \texttt{code/validate\_v8\_7\_physical.py}.
The finite-efficiency tests are in \texttt{code/validate\_v8\_8\_ceiling.py}
and \texttt{results/v8\_8\_ceiling.json}.
'''
(S/'supplement.tex').write_text(supp)
# Retained macros, theorem environments, and packages.
preamble=FULL[:FULL.index(r'\title{')]
preamble=preamble.replace(r'\graphicspath{{manuscript/figures/}{figures/}}',r'\graphicspath{{figures/}{archive/v7_3/manuscript/figures/}}')
preamble=preamble.replace(r'\usepackage{graphicx,booktabs,tabularx,array}',r'\usepackage{graphicx,booktabs,tabularx,array,longtable}')
header=r'''\title{Why Should They Stop Expanding?\\[3pt]
\large Computational Growth, Tipping Budgets, and the Uniformity Problem}
\author{Hyeonjun Joo\\[3pt]\normalsize Independent Researcher, Republic of Korea\\
\normalsize\href{mailto:jonejoo0805@gmail.com}{jonejoo0805@gmail.com}}
\date{13 October 2026}
\hypersetup{pdftitle={Why Should They Stop Expanding? CITH v8},pdfauthor={Hyeonjun Joo}}
\begin{document}\maketitle
\begin{abstract}
Computational abundance does not by itself remove the incentive to expand.
We develop a resource-allocation theory that separates physical budget
growth, computational efficiency, and the creation of valued virtual worlds.
With logarithmic value, concave local production, bounded external prices,
and persistent residual external motives, unbounded budgets exclude
permanent all-inward allocation regardless of productivity. A tight
production-independent floor yields a unique tipping budget; an explicit
resource-acquisition process has a unique boundary between eventual inward
investment and outward lock-in. Optimizing observation and dormant backups
changes the relevant condition from all-inward investment to no replication.
An expected-catastrophe model shows why protecting a growing internal value
can preserve a seed-backup incentive. A free-world-number lemma connects
specified population ethics to value curvature, and a verified
infinite-horizon benchmark separates current computational returns from
future resource growth. Subjective-delay discounting can weaken remote
returns, but must account for the timing of descendant production and the
risk cost of delayed protection. Finally, a marked-origin model quantifies
the small expansionist tail compatible with a chosen null-arrival likelihood
when potential exposure is large. These results specify the resource,
service-substitution, and timing tests for an inward scenario while identifying
alternative routes and failures. They do not estimate extraterrestrial
preferences, establish indefinite local engineering feasibility, or infer
low detectability from inward allocation.
\end{abstract}
\noindent\textbf{Keywords:} Fermi paradox; computational inward transition;
resource allocation; replication; population ethics; technosignatures.
'''
order=['introduction','model_support','stage1_model','stage3_worlds','stage1_corner','phase','stage1_shared_tech','stage1_persistence','stage4_infinite','stage2_residual','stage5_latency','stage6_population','verification_conclusion']
main='\n'.join(stripcomments((S/(n+'.tex')).read_text()) for n in order)
# World lemma placed after complete basic model, before hub; no forward circularity.
(S/'main_text.tex').write_text(main)
# Load declarations, notation, and curated bibliography prepared separately.
proofs='\n'.join(stripcomments((S/n).read_text()) for n in ['stage1_proofs.tex','proofs_stage23.tex','proofs_stage4.tex','proofs_stage5.tex','proofs_stage67.tex'])
extra=(S/'declarations.tex').read_text() if (S/'declarations.tex').exists() else ''
notation=(S/'numerical_appendix.tex').read_text()+'\n'+(S/'notation.tex').read_text()
bib=FULL[FULL.index(r'\begin{thebibliography}'):FULL.index(r'\end{thebibliography}')+len(r'\end{thebibliography}')]
newbib=(S/'new_bibliography.tex').read_text() if (S/'new_bibliography.tex').exists() else ''
bib=bib.replace(r'\end{thebibliography}',newbib+'\n'+r'\end{thebibliography}')
# Only cited references in the active v8 document remain, with the full original retained in archive.
cited=set()
for m in re.finditer(r'\\cite\w*\*?(?:\[[^\]]*\])?\{([^}]+)\}',main+supp+proofs+notation):cited.update(m.group(1).split(','))
parts=re.split(r'(?=\\bibitem)',bib)
kept=[parts[0]]
for part in parts[1:]:
 m=re.search(r'\\bibitem(?:\[.*?\])?\{([^}]+)\}',part,re.S)
 if m and m.group(1) in cited:kept.append(part.replace(r'\end{thebibliography}',''))
kept=[kept[0]]+sorted(kept[1:],key=lambda s: re.search(r'\\bibitem(?:\[.*?\])?\{([^}]+)\}',s,re.S).group(1).lower())
bib=''.join(kept)+'\n'+r'\end{thebibliography}'
bib=bib.replace(r'2014{\natexlab{b}}', '2014')
(S/'references_used.tex').write_text(bib)
complete=preamble+header+'\n% BEGIN MAIN TEXT\n'+main+'\n% END MAIN TEXT\n'+extra+'\n\\appendix\n'+proofs+'\n'+notation+'\n\\clearpage\n'+bib+'\n'+supp+'\n\\end{document}\n'
(R/'CITH_V8_full.tex').write_text(complete)
labels=re.findall(r'\\label\{([^}]+)\}',complete);refs=[]
for m in re.finditer(r'\\(?:[cC]ref|eqref|ref)\*?\{([^}]+)\}',complete):refs.extend(m.group(1).split(','))
audit={'duplicate_labels':sorted({x for x in labels if labels.count(x)>1}),'unresolved_labels':sorted(set(refs)-set(labels)),'label_count':len(labels),'reference_uses':len(refs),'bibliography_missing':sorted(cited-set(re.findall(r'\\bibitem(?:\[.*?\])?\{([^}]+)\}',bib,re.S)))}
(R/'results/assembly_audit.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit,indent=2))
