# CITH follow-up draft: light speed and civilizational growth

**Fast Minds, Slow Light: Causal Limits on the Growth of Computational Civilizations**
Hyeonjun Joo · draft of 14 October 2026 · follow-up to CITH v8

This directory holds a draft follow-up to the CITH v8 manuscript in the parent directory. The draft is for author review. It has not been peer reviewed, and its declarations still contain items marked *[Author confirmation required]*. It was produced with an AI system at the author's request; see the declarations in the PDF.

CITH v8 treats the delay of remote returns as a parameter. This draft makes the speed of light explicit and follows it through three channels of growth. Acquisition is limited by what a light cone contains. Returns from remote projects are delayed, which on a fast subjective clock means heavily discounted. Integration is limited by how large a region can still compute as one unit. A fourth part shows what light speed cannot contain: autonomous descendants.

## Main results

| Channel | Statement in the draft | What it says |
| --- | --- | --- |
| Acquisition | Lemma 3, Corollary 4 | Acquired resources are bounded by the light cone's content, so growth is at most polynomial. Any exponential phase has an explicit maximum duration. |
| | Proposition 5 | A causal cap on v8's acquisition process weakly lowers the basin separator (6.632 → as low as 4.33 in v8's example). Lock-in after a crossing persists. |
| | Lemma 6 | In flat ΛCDM, anything emitted from now on reaches at most 16.68 billion light-years (comoving). The affectable budget is finite. |
| Latency | Lemma 7 | Corner tests in the replication directions alone certify no replication. On v8's frozen design, all 333 certified rows are among v8's 421 no-replication rows. |
| | Proposition 8 | The latency length ℓ_s = c/(δω) turns v8's robustness ratio ρ into a distance. The corner holds exactly when the nearest opportunity lies beyond (ℓ_s/κ)·ln ρ. |
| | Proposition 9 | Centrally evaluated frontier expansion stops at a finite radius, asymptotically (k/κ)·ℓ_s·ln ℓ_s. |
| | Corollary 10 | A backup at distance r lags the home by ωr/c in subjective time, so it insures the past. |
| Integration | Lemma 12, Proposition 13 | A computation that must complete subjective round trips within τ fits in diameter cτ/(2ω). Its power ceiling therefore falls as ω⁻². |
| | Proposition 15 | Resolves v8's open coupling between computational efficiency and clock speed. |
| | Theorem 16 | There is a unique latency-containment budget, and c cancels from it at the fastest admissible clock. |
| Decentralization | Proposition 17 | Lineages of autonomous colonies grow as percolation clusters, so restraint must be inherited colony by colony. |
| | Remark 18 | Near-light expansion fronts are hard to see before they arrive. |

These are conditional results. They do not estimate clocks, discount rates, integration tolerances or heat-flux limits for any civilization; every such value in the examples is labelled as a scenario.

## Files

- [`CITH_LS_draft.pdf`](CITH_LS_draft.pdf): the draft, 22 pages, including proofs and the numerical specification.
- `CITH_LS_draft.tex`: single-file LaTeX source.
- `code/validate_ls_*.py`: twelve self-contained validators, one or more per statement.
- `code/run_all_checks.py`: the runner; it stops on a nonzero exit or a failed assertion.
- `code/make_figures.py`: draws the four figures from the validator outputs.
- `code/audit_draft.py`: a mechanical audit of the draft. It checks references, citations and LaTeX warnings, and checks every quoted number against the JSON outputs.
- `figures/`: the four figure PDFs.
- `results/`: machine-readable outputs, literal stdout logs, `ALL_LS_CHECKS_SUMMARY.json` and `draft_audit.json`.

## Reproduce

From the repository root, with the parent repository's pinned dependencies (`python -m pip install -r requirements.txt`), a pdfLaTeX installation and `texcount`:

```sh
sh followup_light_speed/reproduce.sh
```

The script runs all validators, redraws the figures, compiles the draft three times, and runs the audit. The reference run reports **30,260 passed assertions and zero failures**. As in v8, assertions include identities, endpoints, root checks and repeated inequalities. They are not 30,260 independent random draws. Monte Carlo checks use margins of five standard errors. The full run takes about three minutes, mostly for the cubic-lattice percolation simulations.

`validate_ls_05_latency_radius.py` reads the frozen v7.3 design in `archive/v7_3/data/allocation/design.npz` and v8's `results/v8_no_replication_classifications.npz`. It never writes to either. At zero delay it reproduces v8's 219 all-inward and 421 no-replication counts.

## Validator-to-statement map

| Validator | Statements and implemented scope |
| --- | --- |
| `validate_ls_01_envelope.py` | Lemma 3, Corollary 4: random site sets and causal schedules; explicit duration bound against root finding; truncation of v8's geometric runaway. |
| `validate_ls_02_capped_basin.py` | Proposition 5: policy against direct maximization; domination, monotonicity, capture inclusion, crossing delay, bounded caps, lock-in; v8's separator 6.632159 and Table 2. |
| `validate_ls_03_cosmology.py` | Lemma 6: de Sitter and Einstein–de Sitter closed forms, the 1/√Ω_Λ bound, monotonicity in emission time, probes with decaying peculiar momentum, Planck 2018 numbers, radiation sensitivity. |
| `validate_ls_04_replication_lemma.py` | Lemma 7: exact restricted faces through the multiplier, the bound λ ≥ λ_L, certified cases against v8's exact inequalities, sampled full optimization, witnesses that the converse fails. |
| `validate_ls_05_latency_radius.py` | Proposition 8: multi-distance instances against constrained optimization; frozen-design counts under a common delay, including the exact face test; nesting and monotonicity. |
| `validate_ls_06_frontier.py` | Proposition 9: the general bound on random admissible instances; the benchmark stopping radius against scalar optimization; asymptotic ratios to ℓ_s = 10¹⁴; the no-latency obstruction; Figure 3 data. |
| `validate_ls_07_stale_backup.py` | Corollary 10: complex-step derivatives of the product objective, floor and ceiling equivalences, an explicit staleness model. |
| `validate_ls_08_coherence.py` | Lemmas 2 and 12, Proposition 13: hull areas of random point clouds, Cauchy's formula by Monte Carlo projection, the planar isodiametric inequality, the Landauer–Carnot cost, Table 1. |
| `validate_ls_09_clock_coupling.py` | Proposition 15: unimodality, corner-set structure, agreement with direct optimization, the slow-clock limit, the adiabatic peak, Figure 4a data. |
| `validate_ls_10_containment.py` | Theorem 16: clock scans against the containment criterion, cancellation of c, uniqueness and monotonicity of B_lat, the endpoint test under coupling, Table 3. |
| `validate_ls_11_lineages.py` | Proposition 17: deterministic dichotomy, tree means against the bound, survival probabilities, cubic-lattice regimes. |
| `validate_ls_12_visibility.py` | Remark 18: Poisson space–time simulation of means, independence, ratios, conditioning on no arrival. |

A check of a specialization or a necessary inequality is not a verification of every hypothesis of the cited statement. The proofs in Appendix A remain primary.

## Before any submission

- Confirm or revise the declarations marked *[Author confirmation required]* (funding, competing interests, approval).
- Approve every mathematical claim, interpretation and reference; the generating and checking code is AI-written and is not independent peer review.
- Decide how this draft should relate to v8: a companion paper, a later version, or a section to merge.

## License

As in the parent repository: original code is MIT-licensed; text, figures and numerical outputs are CC BY 4.0.
