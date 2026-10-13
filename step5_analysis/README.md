# Step 5 · Aggregation and prioritization rules

**Goal** — turn the paired scores from steps 2–4 into the audit's decision layer: neighbour-score fractions,
count proxies, the primary / strict / legacy rules, threshold sensitivity, missingness bounds and the
equal-size shortlist comparison.

All of it is deterministic arithmetic on the archived CSVs. **No p values are computed here** — the
receptor-specific thresholds are prioritization rules, not hypothesis tests.

## `output/`

| File | Rows | Contents | Produced by |
|:--|--:|:--|:--|
| `NK{1,2,3}R_summary.csv` | 33 each | per candidate: on-target score, neighbour score statistics (`neg_mean`, `neg_median`, `neg_min`, `neg_max`, `n_neg`) and the legacy `selectivity` field | `scripts/summarize_results.py` |
| `all_receptors_summary.csv` | 99 | the three summaries stacked, including the 3 controls | `scripts/summarize_results.py` |
| `NK{1,2,3}R_pairwise_selectivity.csv` | 32 / 18 / 32 | the matched aggregates: `n_common`, `n_safer_delta_pos`, `frac_safer`, `delta_mean`, `delta_median`, `delta_min` | `scripts/pairwise_selectivity.py` |
| `NK{1,2,3}R_final_ranking.csv` | 33 each | summary ⊕ aggregates, plus the legacy booleans `better_binding`, `safer`, `both` | `scripts/final_ranking.py` |
| `paired_neighbour_deltas.csv` | 4,518 | **one row per (receptor, candidate, neighbour)** — `delta`, `higher_than_control` | `scripts/paired_neighbour_deltas.py` |
| `paired_deltas_summary.csv` | 82 | recomputed per-candidate counts next to the shipped aggregates, for cross-checking | `scripts/paired_neighbour_deltas.py` |
| `peptide_codes.csv` | 99 | the `NK1R-01…NK3R-32` code system ↔ PepMLM id ↔ sequence (TableS1 source) | `scripts/make_peptide_codes.py` |
| `umap_input.csv` | 1,048 | the comparison set: 949 leaderboard-only + 99 study sequences | `scripts/make_umap_data.py` |
| `umap_coords.csv` | 1,048 | archived two-dimensional coordinates plus `nn_dist_lb` (nearest retained leaderboard entry) | `scripts/make_fig8_umap.py` |

`scripts/embed_esm.py` documents the archived embedding step (esm2_t33_650M_UR50D, layer 33, mean pooling over
residue tokens); the audit **reuses the archived coordinates** rather than regenerating them.

## Column translation

The archive predates the manuscript's vocabulary. Two legacy names are worth flagging:

* **`selectivity`** is `neg_mean − pos_score` — a difference of two model scores, **not** a selectivity ratio
  and not a safety margin. The manuscript does not use the term.
* **`ΔNeg`** is identical to `delta_mean`, the mean candidate-minus-control difference over shared neighbours.

Full mapping: see the terminology table in the [top-level README](../README.md#terminology-manuscript--archive).

## The rules

```
primary rule      ΔPos < 0  AND  frac_safer > 0.50  AND  n_common / (Table S3 rows for that receptor) ≥ 0.90
strict rule       as primary, with frac_safer ≥ 0.75
legacy audit rule ΔPos < 0  AND  frac_safer > 0.50                        (no count proxy)
```

Each rule is applied **independently within each receptor**, and the primary and legacy rules agree on this
archive because every evaluable candidate has a count proxy of exactly **1.00**.

A separate exploratory grid counts candidates at `frac_safer > cutoff` for cutoff = 0.50, 0.55 … 0.75; this grid
is *not* the strict rule.

## Results

| | NK1R | NK2R | NK3R |
|:--|--:|--:|--:|
| Candidates with neighbour summaries | 32 | 18 | 32 |
| Matched comparisons (`n_common` × candidates) | 2,112 | 486 | 1,920 |
| Primary-rule shortlist | 3 | 10 | 3 |
| Strict rule (≥ 0.75) | 0 | 2 | 0 |

* **82 evaluable candidates / 4,518 comparisons** in total; 14 NK2R candidates have no aggregates and are
  classified **missing**, not failed.
* **43** candidates had ΔPos < 0; **33** of them were evaluable; the joint rule retained **16/33**, a
  **51.5%** reduction.
* At matched shortlist sizes (top 3 / 10 / 3 by ascending ΔPos), joint and on-target-only selections shared
  **0/3, 7/10 and 1/3 → 8/16**.
* Raising the cutoff from > 0.50 to **> 0.55 removes all NK1R candidates**; the counts sit close to the
  boundary, so candidates just above and below a cutoff are neighbouring priorities rather than distinct
  biological classes.
* **NK2R missingness bounds:** 14 missing candidates, 10 of them with ΔPos < 0 → **10–20** primary-rule and
  **2–12** strict-rule candidates. These are unconstrained upper bounds; they do not predict the missing
  outcomes.
* All six strict-rule/primary-rule leaders, the three discordant candidates (NK1R-18 29/66, NK2R-27 11/27,
  NK1R-14 41/66 with ΔPos +50.23) and the per-receptor shortlists are listed in the
  [top-level README](../README.md#results-in-the-archive).

## Reproducing the pair-level evidence

```bash
python scripts/paired_neighbour_deltas.py --repo ..
```

```
[ok] 4518 pair-level comparisons over 82 candidates
[check] NK1R: 32 candidates compared, 0 disagree with the shipped aggregates
[check] NK2R: 18 candidates compared, 0 disagree with the shipped aggregates
[check] NK3R: 32 candidates compared, 0 disagree with the shipped aggregates
```

The script joins candidate and control **by sequence**, never by batch label.
