# CITH v8

**Why Should They Stop Expanding? Computational Growth, Tipping Budgets, and the Uniformity Problem**  
Hyeonjun Joo · 13 October 2026

The computational inward transition hypothesis (CITH) is a conditional model of resource allocation between internal computational growth and external investment. Version 8 derives tipping-budget and growth conditions, distinguishes all-inward allocation from no further replication, adds an infinite-horizon benchmark, and examines endogenous risk, population ethics, subjective latency and the expansionist tail compatible with a chosen null-arrival likelihood.

This is a research manuscript, not a peer-reviewed journal publication. It does not estimate extraterrestrial preferences, establish indefinite local engineering feasibility, or infer low detectability from inward allocation. The infinite-horizon result applies to its stated homogeneous benchmark.

## Manuscript and materials

- [Article](CITH_V8_Manuscript.pdf): 40 pages, including proofs, numerical specification and references.
- [Supplement](CITH_V8_Supplement.pdf): 7 pages of retained examples and numerical details.
- [Combined PDF](CITH_V8_full.pdf): 47 pages with cross-document references.
- [Assembled LaTeX](CITH_V8_full.tex): compile from this directory.
- `latex/`: editable modular sources.
- `code/`: 15 v8 numerical validators and reproduction runners.
- `figures/`: the four preserved figure PDFs.
- `results/`: machine-readable numerical outputs and mechanical source/PDF checks.
- `archive/v7_3/`: frozen public v7.3 sources, synthetic inputs and outputs needed for reproduction.

All-inward allocation and no replication are separate classifications. On the frozen 1,024-row synthetic design, v8 retains 219 all-inward cases and finds 421 no-replication cases, including 202 observation-only cases; 603 cases prefer replication. These counts describe the specified design, not civilization prevalence. Every row receives the restricted-face criterion; 64 rows also receive full constrained-optimizer comparisons.

## Reproduce the numerical checks

```sh
python -m pip install -r requirements.txt
python code/run_all_checks.py
python code/rerun_archive.py
```

The v8 runner executes all 15 validators and stops on a nonzero exit or failed assertion. The reference execution reports **33,017 passed assertions and zero failures**. Assertions include identities, endpoints, root checks and repeated inequalities; they are not 33,017 independent random draws. Seeds, sample sizes, tolerances and optimizer settings are explicit in each validator. Numerical checks supplement the manuscript proofs.

The archive runner executes four v7.3 commands in a temporary working copy and writes results under `results/archive_rerun/`. It does not modify the archived inputs. The frozen arrays, rather than a seed alone, define the reference design. Its SHA-256 is `ce45ecd0435f2abf834445957db369d1d4a5bc28e9c4e7b686561df50a14194b`.

The reference environment used Python 3.13.5, NumPy 2.3.5 and SciPy 1.17.0. Dependencies are pinned in `requirements.txt`; regenerated floating-point results can vary in their final bits across platforms.

## Build and audit PDFs

Install pdfLaTeX, `texcount`, the packages in the source preamble, and the Python dependencies above. From the repository root:

```sh
python qa/assemble.py
pdflatex -interaction=nonstopmode -halt-on-error CITH_V8_full.tex
pdflatex -interaction=nonstopmode -halt-on-error CITH_V8_full.tex
python qa/final_audit.py
```

The audit checks labels, unresolved references, compilation warnings, word counts, figure identity and PDF text bounds, then splits the combined PDF. The main prose contains 9,091 words excluding the abstract, declarations, appendices, references and supplement; the abstract contains 209 words. Figures 1, 2 and 4 appear as supplementary Figures S1–S3; the main phase diagram retains number 3.

## Integrity, license and citation

Before regeneration, run `sha256sum -c SHA256SUMS`. Generated logs, numerical outputs and PDF metadata may change on rerun; scientific reproduction uses the stated classifications and tolerances.

Original Python and shell code is MIT-licensed. The manuscript, supplementary text, original figures, documentation, original synthetic data and numerical outputs are CC BY 4.0. Third-party components retain their existing licenses, including the LPPL notices in `archive/v7_3/manuscript/cithplainnat.bst`. See `LICENSE` and `LICENSE-CODE`.

For v8, cite this manuscript version and the repository commit you used. `CITATION.cff` contains v8 metadata. A v8 Zenodo DOI has not yet been assigned in this repository update.

The previous **v7.3** release is archived at [10.5281/zenodo.23317930](https://doi.org/10.5281/zenodo.23317930). The [concept DOI](https://doi.org/10.5281/zenodo.23317929) identifies that Zenodo record across versions; it should not be treated as confirmation that v8 is already archived there.
