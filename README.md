# CITH v8

**Why Should They Stop Expanding? Computational Growth, Tipping Budgets, and the Uniformity Problem**  
Hyeonjun Joo · 13 October 2026

The computational inward transition hypothesis (CITH) is a conditional model of resource allocation between internal computational growth and external investment. Version 8 derives tipping-budget and growth conditions, distinguishes all-inward allocation from no further replication, adds an infinite-horizon benchmark, and examines endogenous risk, population ethics, subjective latency and the expansionist tail compatible with a chosen null-arrival likelihood.

This is a research manuscript, not a peer-reviewed journal publication. It does not estimate extraterrestrial preferences, establish indefinite local engineering feasibility, or infer low detectability from inward allocation. The infinite-horizon result applies to its stated homogeneous benchmark.

## Manuscript and materials

- [Article](CITH_V8_Manuscript.pdf): 44 pages, including proofs, numerical specification and references.
- [Supplement](CITH_V8_Supplement.pdf): 7 pages of retained examples and numerical details.
- [Combined PDF](CITH_V8_full.pdf): 51 pages with cross-document references.
- [Assembled LaTeX](CITH_V8_full.tex): compile from this directory.
- `latex/`: editable modular sources.
- `code/`: 16 v8 numerical validators and reproduction runners.
- `figures/`: the four preserved figure PDFs.
- `results/`: machine-readable numerical outputs and mechanical source/PDF checks.
- `archive/v7_3/`: frozen public v7.3 sources, synthetic inputs and outputs needed for reproduction.

All-inward allocation and no replication are separate classifications. On the frozen 1,024-row synthetic design, v8 retains 219 all-inward cases and finds 421 no-replication cases, including 202 observation-only cases; 603 cases prefer replication. These counts describe the specified design, not civilization prevalence. Every row receives the restricted-face criterion; 64 rows also receive full constrained-optimizer comparisons.

## Follow-up draft: light speed and civilizational growth

[`followup_light_speed/`](followup_light_speed/) contains a draft follow-up, **Fast Minds, Slow Light: Causal Limits on the Growth of Computational Civilizations** ([PDF](followup_light_speed/CITH_LS_draft.pdf), 22 pages). It makes the speed of light explicit in the v8 model, through three limits:

- the light cone bounds resource acquisition;
- remote returns arrive late, and so are discounted on the latency length c/(δω);
- integration within a subjective round trip τ bounds a computation's size, and with it its power, to a ceiling falling as ω⁻².

The draft also shows why none of these limits contains autonomous descendants. It is a draft for author review, not part of the v8 record. It has its own twelve validators, figures, audit and `SHA256SUMS`, and reads the frozen v7.3 design without modifying it. See its README.

## Reproduce the numerical checks

```sh
python -m pip install -r requirements.txt
python code/run_all_checks.py
python code/rerun_archive.py
```

The v8 runner executes all 16 validators and stops on a nonzero exit or failed assertion. The reference execution reports **36,206 passed assertions and zero failures**. Assertions include identities, endpoints, root checks and repeated inequalities; they are not 36,206 independent random draws. Seeds, sample sizes, tolerances and optimizer settings are explicit in each validator. Numerical checks supplement the manuscript proofs.

The archive runner executes four v7.3 commands in a temporary working copy and writes results under `results/archive_rerun/`. It does not modify the archived inputs. The frozen arrays, rather than a seed alone, define the reference design. Its SHA-256 is `ce45ecd0435f2abf834445957db369d1d4a5bc28e9c4e7b686561df50a14194b`.

The reference environment used Python 3.13.5, NumPy 2.3.5 and SciPy 1.17.0. Dependencies are pinned in `requirements.txt`; regenerated floating-point results can vary in their final bits across platforms.

## Validator-to-statement map

The table describes checks actually implemented in each script. A check of a
specialization or necessary inequality is not verification of every hypothesis
or every possible path of the cited result; the analytic proofs remain primary.
All filenames below are under `code/`.

| Validator | Statements and implemented scope |
| --- | --- |
| `validate_v8_1a_hub.py` | `prop:hub`, `prop:corner`, `cor:feasibility`: corner derivatives, logarithmic reductions, constrained optima, zero/tie cases. |
| `validate_v8_1b_floor.py` | `cor:production-growth`: concave-production floor, linear equality witnesses, strict boundary correction. |
| `validate_v8_1c_persistence.py` | `prop:growth-trichotomy`; growing-budget necessary bounds related to `prop:persistence`, subunit-curvature examples, and the exact `eta=0` condition. Not a general arbitrary-path test of `eq:pathwise`. |
| `validate_v8_1d_exponential.py` | `prop:exponential-necessity` production identities; `lem:exponential-endpoints` global comparisons and smooth counterexamples. The necessity proof is analytic. |
| `validate_v8_1e_basin.py` | `prop:basin`, `lem:stay-certificate`: finite basin brackets, policy monotonicity, runaway and capture certificates. |
| `validate_v8_1e_tipping.py` | `prop:tipping`, `prop:curvature`: tipping roots, threshold monotonicity, corner optima, bounded productivity intervals and ties. |
| `validate_v8_2a_replication.py` | `prop:no-replication`, `cor:bundle`: restricted-face KKT, full five-allocation comparisons, shadow-price identities, zero services and ties. |
| `validate_v8_2b_resilience.py` | `prop:endogenous-risk`: exact derivatives, necessary local risk conditions and asymptotics; no global risk-optimality claim. |
| `validate_v8_2c_growing.py` | `cor:growing-no-replication`: the specific unbounded-productivity construction `a=B^2`, including its optimized services. |
| `validate_v8_3_worlds.py` | `lem:world-curvature`: optimized world number, logarithmic welfare, normalization and integer correction. |
| `validate_v8_4_infinite.py` | `prop:infinite`, `cor:infinite-reversals`, `cor:dynamic-residual`, `eq:capped-transition`: Bellman limits, transversality, policy bounds, deviations and capped transitions. |
| `validate_v8_5_latency.py` | `prop:latency`, `eq:discounted-replication`, `prop:depletion`: clock quadrature, corner comparisons, discount thresholds, chain series and depletion edges. |
| `validate_v8_6_design.py` | The frozen 1,024-row classification: 219 all-inward, 421 no-replication, 202 observation-only, 603 replication-preferred; 64 full-optimizer comparisons. |
| `validate_v8_6_population.py` | `prop:uniformity`, `prop:preemption`: marked-Poisson bounds, simulations, strategic endpoint tests and counterexamples. |
| `validate_v8_7_physical.py` | `eq:stellar-budget`, `tab:latency-sensitivity`, `eq:radiator`: dimensional identities and the stated illustrative parameters, not engineering feasibility. |
| `validate_v8_8_ceiling.py` | A1-A4 and shared technology: `cor:ceiling-budget`, `cor:capped-basin`, `app:ceiling-growth-proof`, `app:norep-ceiling-proof`, `cor:shared-technology`. Includes empty/equality cases, finite trajectories, separator brackets, counterexamples, the zero-local-investment boundary, and archived Table S2 checks. |

The ceiling-basin theorem assumes that at least one initial productivity has a
finite crossing. Strict ordering alone is insufficient; the appendix gives an
all-captured family. For a hard cap, a boundary path can instead be captured in
finite time. The validator's finite brackets do not numerically prove an
infinite-time limit. The ceiling no-replication test is necessary; the exact
criterion uses actual productivity at the optimized restricted-face point.

## Build and audit PDFs

Install pdfLaTeX, `texcount`, the packages in the source preamble, and the Python dependencies above. From the repository root:

```sh
python qa/assemble.py
pdflatex -interaction=nonstopmode -halt-on-error CITH_V8_full.tex
pdflatex -interaction=nonstopmode -halt-on-error CITH_V8_full.tex
python qa/final_audit.py
```

The audit checks labels, unresolved references, compilation warnings, word counts, figure identity and PDF text bounds, then splits the combined PDF. The main prose contains 9,592 words excluding the abstract, declarations, appendices, references and supplement; the abstract contains 209 words. Figures 1, 2 and 4 appear as supplementary Figures S1–S3; the main phase diagram retains number 3.

## Integrity, license and citation

Before regeneration, run `sha256sum -c SHA256SUMS`. Generated logs, numerical outputs and PDF metadata may change on rerun; scientific reproduction uses the stated classifications and tolerances.

Original Python and shell code is MIT-licensed. The manuscript, supplementary text, original figures, documentation, original synthetic data and numerical outputs are CC BY 4.0. Third-party components retain their existing licenses, including the LPPL notices in `archive/v7_3/manuscript/cithplainnat.bst`. See `LICENSE` and `LICENSE-CODE`.

For v8, cite this manuscript version and the repository commit you used. `CITATION.cff` contains v8 metadata. A v8 Zenodo DOI has not yet been assigned in this repository update.

The previous **v7.3** release is archived at [10.5281/zenodo.23317930](https://doi.org/10.5281/zenodo.23317930). The [concept DOI](https://doi.org/10.5281/zenodo.23317929) identifies that Zenodo record across versions; it should not be treated as confirmation that v8 is already archived there.
