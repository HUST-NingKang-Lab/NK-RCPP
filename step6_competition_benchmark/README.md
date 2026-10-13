# Step 6 · Exploratory external NK2R benchmark

**Goal** — place the study's NK2R candidates next to an external set of 10 competition peptides that carry
reported EC50 values, and ask how docking scores, neighbour-score fractions and length-adjusted z-scores
behave. The benchmark provides **context, not independent pharmacological validation**.

## `input/`

| Path | Contents |
|:--|:--|
| `COMP_TOP10_NK2R_peptide.csv` | the 10 competition peptides, docked against NK2R |
| `comp_top10_pos_inputs/`, `comp_top10_neg_inputs/` | the same peptides split per batch for the on-target and neighbour campaigns (`comp01…comp10`, plus `known_peptide` and `known_missing`) |
| `LENGTH_CONTROL_NK2R_peptide.csv`, `length_control_inputs/` | 5 random controls at each of 6, 9, 10 and 15 residues (20 controls) |
| `COMP_TOP10_negative_docking_input.STALE-offby-one.csv` | kept deliberately — see *Watch out* |

The 10-residue random-control pool was built from NK2R-01 and the NKA-derived control, **not** from the full
study-candidate composition; the other lengths use length-specific residue pools.

## `output/`

| File | Rows | Contents |
|:--|--:|:--|
| `competition_positive_scores.csv` | 10 | on-target poses (`rank1…rank10`, `best_score`) for the competition peptides |
| `competition_top10.csv` | 10 | code, name, sequence, workbook EC50, workbook rank |
| `length_control.csv` | 20 | on-target score per length control with its length |
| `COMP_TOP10_summary.csv` | 10 | the full per-peptide row: on-target score, neighbour statistics, `frac_safer`, `delta_mean`, joint-rule flag, `z_vs_length` |
| `NK2R_comparison.csv` | 43 | competition (10) + study NK2R candidates (32) + control, with `source` and `EC50_nM` |

`z_vs_length` = (peptide score − mean of the matching length controls) / (sample s.d. of those controls).

## Results

| Endpoint | Competition peptides (n = 10) | Study NK2R candidates |
|:--|:--:|:--|
| Length-adjusted z < 0 | **9/10** | 11/31 study 10-mers |
| ΔPos < 0 | **5/10** | — |
| Neighbour-score fraction > 0.50 | **8/10** | 14/18 evaluable |
| Joint on-target + fraction rule | **3/10** (3/5 of the competition 9-mers) | 10/18 evaluable |

* Score–potency agreement was **uncertain**: Spearman rho = +0.42, P = 0.23 (n = 10).
* Exploratory two-sided Mann-Whitney: **P = 0.2204** raw scores (10 vs 32 peptides), **P = 0.1699**
  neighbour-score fraction (10 vs 18), **P = 5.09 × 10⁻⁴** length-adjusted z (10 vs 31).
* Raw REF15 values vary with length across the composition-specific control groups — the comparison does not
  isolate a pure length effect.
* NK2R-02 (11 residues) is excluded from the z analysis because no 11-residue baseline exists; it stays in the
  raw-score and fraction comparisons.

> **All of these p values are exploratory, two-sided and uncorrected for multiplicity.** The groups differ in
> length, composition and sampling frame and are not exchangeable validation cohorts. Favourable relative
> scores can support computational ranking; their calibration to potency requires independently documented
> assays. In [`Fig. S6`](../figures/README.md), the overlaid study trace is the lowest-on-target-score
> candidate NK2R-27 — **not** a validated best binder.

## Watch out

* **`COMP_TOP10_negative_docking_input.STALE-offby-one.csv` is intentionally versioned.** The merged
  competition input file carries pre-fix labels shifted by one relative to the per-batch inputs, so looking a
  peptide's sequence up by complex name there silently returns its neighbour's sequence. The per-batch files
  under `comp_top10_neg_inputs/` are authoritative, and `scripts/harvest_docking_scores.py` emits a warning on
  every disagreement. All 243 such warnings in the harvest come from this campaign.
* The `COMP_TOP10_neg` results were stored flat (`<root>/<complex>/`) rather than batched, unlike the receptor
  campaigns; the harvest script handles both layouts.
* Competition sequences and EC50 values were checked against the supplied workbook, but **the originating assay
  protocol and raw experimental records were not available** for independent verification.
