# Step 2 · On-target docking (RAPiDock)

**Goal** — score every candidate against its own receptor, and against the receptor-matched control ligand,
under one identical protocol.

## `input/`

| Path | Contents |
|:--|:--|
| `pdbs/NK1R.pdb`, `pdbs/NK2R.pdb`, `pdbs/NK3R.pdb` | The receptor structures handed to RAPiDock. NK1R = PDB 6E59 chain A; NK2R = AlphaFold2 model of UniProt P21452; NK3R = predicted model (residues 76–369). |
| `NK1R_docking_input.csv` | 32 candidates + `Substance_P` control |
| `NK2R_docking_input.csv` | 32 candidates + `pep_known` control |
| `NK3R_docking_input.csv` | 32 candidates + `Neurokinin_B` control |

Each CSV is the RAPiDock `--protein_peptide_csv` schema: `complex_name, protein_description, peptide_description`.
The control row is what makes the comparison paired — it is docked through the same code path as the
candidates, never carried over from a published structure.

## How it was run

```bash
# conda env: RAPiDock   (see scripts/run_docking_gpu.sh for the driver script)
python inference.py --config default_inference_args.yaml \
                    --protein_peptide_csv NK2R_peptide.csv \
                    --output_dir results/NK2R_peptide_len10
```

Key settings from `RAPiDock/default_inference_args.yaml`:

| Setting | Value |
|:--|:--|
| `N` (poses per complex) | 10 |
| `inference_steps` / `actual_steps` | 16 / 16 |
| `scoring_function` | `ref2015` |
| `ckpt` | `rapidock_local.pt` (local docking model) |
| `conformation_partial` | `1:1:1` (polyproline / helical / extended peptide starts) |
| `batch_size` | 4 |

## `output/`

| File | Contents |
|:--|:--|
| `positive_scores_long.csv` | one row per (campaign, complex, pose) — 1290 rows = 129 complexes × 10 poses |
| `positive_scores_wide.csv` | one row per complex, `rank1…rank10` plus `best_score` |

`campaign` ∈ `{NK1R, NK2R, NK3R, COMP_TOP10, LENGTH_CONTROL}`; `is_control` flags the receptor-matched
ligand. Both files are regenerated from the raw run folders by
[`scripts/harvest_docking_scores.py`](../scripts/harvest_docking_scores.py).

## Result

| Receptor | Control (REF15) | Best candidate | Best Delta-on | Candidates better than control |
|:--|--:|:--|--:|:--:|
| NK1R | 574.39 | NK1R-18 (517.29) | −57.1 | 14 / 32 |
| NK2R | 771.13 | NK2R-27 (636.35) | −134.8 | 23 / 32 |
| NK3R | 2240.54 | NK3R-21 (1965.25) | −275.3 | 6 / 32 |

## Watch out

* **REF15 is a Rosetta all-atom energy, not a binding free energy**, and its absolute value depends on the
  system (protein size, cofactors, protonation) — hence the control ligand in every comparison. Nothing in
  this repository compares raw scores *across* receptors.
* `rank1` is the pose RAPiDock ranked first; the table also carries `best_score` (minimum over the 10
  poses) so the two choices can be compared. All downstream analysis uses `rank1`.
* Diffusional docking is sensitive to the peptide starting conformation; the three `conformation_partial`
  starts are what give the 10 poses their spread.
