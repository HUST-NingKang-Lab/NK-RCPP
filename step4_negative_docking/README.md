# Step 4 · Matched structural-neighbour docking

**Goal** — score the candidate **and** its receptor-matched control against every member of the step-3 panel,
under one identical protocol, so that each candidate–control comparison is *paired on the same structure*.

This is the step that produces the neighbour-score fraction. It is a **computational stress test of score
behaviour across matched contexts** — it is not an off-target binding panel and not a selectivity ratio.

## `input/`

| Path | Contents |
|:--|:--|
| `NK1R_negative_docking_input.csv` | 32 candidates + `known` (Substance P) × 66 neighbours |
| `NK2R_negative_docking_input.csv` | 18 candidates + `peptide19` (NKA-derived control) × 27 neighbours |
| `NK3R_negative_docking_input.csv` | 32 candidates + `known` (neurokinin B) × 60 neighbours |
| `nk1r_neg_inputs/`, `nk2r_neg_inputs/`, `nk3r_neg_inputs/` | the same rows split one file per batch — **the authoritative copy** |
| *(step 6)* `comp_top10_neg_inputs/` | competition peptides and length controls, same scheme |

Each CSV is the RAPiDock `--protein_peptide_csv` schema: `complex_name, protein_description, peptide_description`.
Every batch carries its own control row, docked through the same code path as the candidates.

## `output/`

| File | Contents |
|:--|:--|
| `negative_scores_long.csv` | one row per (campaign, batch, neighbour): `rank1_score`, `best_score`, `n_poses` — 4,968 rows (NK1R 2,178, NK3R 1,980, NK2R 513, COMP_TOP10 297) |
| `negative_scores_per_batch.csv` | one row per batch with `n_targets` and `mean_score` |

Both are regenerated from the raw run folders by
[`scripts/harvest_docking_scores.py`](../scripts/harvest_docking_scores.py), which handles the two on-disk
layouts (batched `<root>/<batch>/<complex>/` and flat `<root>/<complex>/`).

## Scoring pipeline

`scripts/pairwise_selectivity.py` turns the per-complex scores into the aggregates shipped in
[`../step5_analysis/output/`](../step5_analysis/output/):

```
for each structural neighbour s of the receptor panel
    delta(s)   = candidate_score(s) − control_score(s)
    n_common   = #{ s : both docked }
    n_higher   = #{ s : delta(s) > 0 }
    frac_safer = n_higher / n_common
```

which reproduces the manuscript's **4,518 matched comparisons over 82 evaluable candidates** (66 × 32 = 2,112
for NK1R, 27 × 18 = 486 for NK2R, 60 × 32 = 1,920 for NK3R).

## Watch out

* **The NK2R batch labels are not the canonical candidate codes.** `nk2r_neg_inputs/peptide01…peptide18` map
  to the NK2R-01/03/04/07/09/10/11/13/14/19/22/23/24/26/27/28/30/32 subset, and `peptide19` is the **control**.
  Step 4 NK2R therefore covers **18 of 32** candidates — exactly the 18 with neighbour summaries in the
  manuscript, and the reason 14 NK2R candidates are reported as missing.
* **Never join on the batch label — join on sequence.** `peptideNN` *does* line up with the canonical codes for
  NK1R and NK3R (control batch `known`), which makes the NK2R mismatch easy to miss.
* **Absolute scores in this campaign are not comparable with step 2's.** They come from a different scored system
  (panel protein + peptide rather than receptor + peptide) and can even change sign between receptors. Only the
  within-campaign candidate-minus-control difference, and the fraction of positive deltas, are meaningful.
* **`rank1` is used throughout**; `best_score` is shipped alongside it so the pose choice can be varied.
