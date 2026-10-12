# CITH v7.3

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23317929.svg)](https://doi.org/10.5281/zenodo.23317929)

**Why Should They Expand? Computational Growth, Resource Sufficiency, and Conditional Interstellar Expansion**  
Hyeonjun Joo · 12 October 2026

The computational inward transition hypothesis (CITH) is a conditional model of how a civilization allocates resources between internal computational growth and external investment. The manuscript derives conditions for an all-inward optimum, constraints on its persistence as budgets grow, and examples in which future resources reverse a present allocation preference.

This is a research manuscript, not a peer-reviewed journal publication. The model does not establish the prevalence of inward civilizations, indefinite physical feasibility, or population-level silence. Allocation and detectability are distinct claims.

## Manuscript and materials

- [Article](CITH_V7_3_revised.pdf)
- [Line-numbered article](CITH_V7_3_line_numbered.pdf)
- [Supplementary material](CITH_V7_3_supplement.pdf)
- `manuscript/`: modular LaTeX sources, bibliography and figures.
- `CITH_V7_3_full.tex`: expanded single-file LaTeX source; compile from this directory.
- `sections/`: section-by-section source copies.
- `code/`: numerical checks, input audits, Figure 3 generation and PDF build tools.
- `data/`: normalized synthetic inputs and archived allocation solutions.
- `results/`: recorded numerical outputs and figure arrays.

## Reproduce the numerical checks

The archived numerical environment used Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0 and Matplotlib 3.10.8. Package versions are pinned in `requirements.txt`.

```sh
python -m pip install -r requirements.txt
sh reproduce.sh
```

The script checks the growth, curvature and production results, the two-period examples, archived allocation classifications, optimality certificates and measure sensitivity. It regenerates Figure 3. Figures 1, 2 and 4 are supplied unchanged; this script does not regenerate them or rerun the population Monte Carlo experiment.

The frozen arrays in `data/allocation/design.npz`, rather than the seed alone, define the reference experiment. The Sobol generator is specified in the article and supplement.

## Build PDFs

A TeX installation providing `pdflatex` and the packages used by the manuscript is required. Supplied `.bbl` files provide the bibliography.

```sh
python code/build_v73.py
```

This builds the article, line-numbered article, supplement and expanded-source PDF, and synchronizes the expanded source and section copies. Build logs are written under `results/build_logs/`.

## Integrity

Before regenerating files, check the distribution with:

```sh
sha256sum -c SHA256SUMS
```

Regeneration can change floating-point final bits and PDF metadata. Scientific reproduction depends on the stated classifications and tolerances, rather than regenerated PDF byte identity.

## License and citation

Original Python and shell code is licensed under MIT. The manuscript, supplementary text, original figures, documentation and original synthetic data are licensed under CC BY 4.0. Third-party components retain their existing licenses; `manuscript/cithplainnat.bst` retains its LaTeX Project Public License notices. See `LICENSE` and `LICENSE-CODE`.

When citing v7.3, use its version DOI:

> Joo, H. (2026). *Why Should They Expand? Computational Growth, Resource Sufficiency, and Conditional Interstellar Expansion* (Version 7.3). Zenodo. https://doi.org/10.5281/zenodo.23317930

The concept DOI [10.5281/zenodo.23317929](https://doi.org/10.5281/zenodo.23317929) refers to all versions and resolves to the latest. Citation metadata is also provided in `CITATION.cff`.
