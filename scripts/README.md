# Scripts

Grouped by pipeline step. Many still carry their original descriptive filenames (`make_fig8.py` produces
manuscript Fig. 5) — the [figures index](../figures/README.md) gives the authoritative file ↔ figure mapping.

**Paths.** Every script resolves its inputs relative to the checkout (`Path(__file__).resolve().parent.parent`), so a
clone runs anywhere. Raw RAPiDock run directories, the competition workbook and the Foldseek work directories are
**not** shipped — big binaries — so point the scripts at your own copies with environment variables:

| Variable | Meaning | Default |
|:--|:--|:--|
| `NK_RCPP_BASE` | repository root | the script's grandparent directory |
| `NK_RCPP_RAW` | raw docking runs, `peptide_docking/` | `<repo>/raw` |
| `NK_RCPP_COMP_XLSX` | competition workbook | `<repo>/step6_competition_benchmark/input/…xlsx` |

The `*.sh` drivers follow the same idea: `RAPIDOCK_DIR` (default `$HOME/RAPiDock`), `CONDA_ENV`, `PDB_DIR`,
`DOCKING_DIR`, `PEPTIDES_DIR`. The audit scripts under *Step 5* also take `--repo` / `--source` arguments.

## Step 1 · Candidate generation

| Script | Purpose |
|:--|:--|
| `make_candidate_tables.py` | builds `step1_peptide_generation/output/*_pepmlm_candidates.csv` from the archived Colab exports, asserting that the campaign numbering matches the docking inputs |
| `make_peptide_codes.py` | the `NK1R-01…NK3R-32` code ↔ PepMLM id ↔ sequence map (TableS1 source) |
| `extract_chains.py` | extracts the receptor chains from the PDBs into `receptor_chains.fasta` |

## Step 2 · On-target docking

| Script | Purpose |
|:--|:--|
| `run_docking_gpu.sh`, `run_full_pipeline_master.sh`, `launch_workers.sh`, `launch_workerD.sh` | RAPiDock runs on a GPU machine (per receptor, then batched) |
| `harvest_docking_scores.py` | harvests raw RAPiDock output folders into the step 2 / 4 / 6 score tables; handles both on-disk layouts and warns when per-batch sequences disagree with a merged input file |
| `summarize_results.py` | per-receptor summaries (`*_summary.csv`, `all_receptors_summary.csv`) |

## Step 3 · Structural-neighbour retrieval

| Script | Purpose |
|:--|:--|
| `foldseek_search.py` | Foldseek web-API query, combined 3Di + amino-acid mode |
| `parse_and_download.py` | parses the PDB100 result table, redundancy-reduces, downloads structures |
| `cluster_sequences.py` | pure-Python fallback clustering (greedy 3-mer Jaccard 0.50) used when CD-HIT is unavailable |
| `download_pdbs_parallel.sh` | parallel structure download (8 threads) |
| `harvest_foldseek_hits.py` | condenses the alignment table to `pdb100_best_hits.csv` + companions |

## Step 4 · Matched structural-neighbour docking

| Script | Purpose |
|:--|:--|
| `create_neg_input.py` | builds the per-batch neighbour-docking inputs |
| `run_fix_batches.sh`, `comp_top10_worker.sh`, `run_comp_top10_gpu.sh`, `run_length_control.sh` | the neighbour-docking campaigns on a GPU machine |
| `pairwise_selectivity.py` | the candidate–control aggregates (`n_common`, `n_safer_delta_pos`, `frac_safer`, `delta_mean/median/min`) |
| `paired_neighbour_deltas.py` | **rebuilds the pair-level comparisons** from the per-complex scores — 4,518 over 82 candidates, zero disagreement with the shipped aggregates |
| `fix_shifted_labels.py` | repairs the COMP-campaign label shift; the stale merged input is kept in step 6 on purpose |

## Step 5 · Prioritization and figures

| Script | Purpose |
|:--|:--|
| `final_ranking.py` | merges on-target and neighbour evidence into `*_final_ranking.csv` with the legacy booleans |
| `make_supp_tables.py` | TableS1–S3 |
| `make_ms_fig2.py`, `make_ms_fig3.py`, `make_ms_fig4.py` | manuscript Fig. 2, Fig. 3, Fig. 4 |
| `make_fig8.py`, `make_fig8_umap.py`, `make_umap_data.py`, `embed_esm.py` | manuscript Fig. 5 (pocket views + archived ESM-2/UMAP projection) |
| `make_fig2.py`, `make_fig3.py`, `make_fig4.py` | Fig. S8, Fig. S5, Fig. S9 (earlier, wider-scale versions) |
| `make_supp_figs.py` | Fig. S1–S3 |
| `make_binding_modes.py`, `pdb_render.py` | Fig. S4 pose renderings (pure Bio.PDB + matplotlib, no PyMOL) |
| `fig_style.py` | shared house style and `save_pub()` (SVG + PDF + 600-dpi TIFF) |

## Step 6 · External NK2R benchmark

| Script | Purpose |
|:--|:--|
| `make_comp_top10_inputs.py` | extracts the top 10 competition peptides from the workbook and builds docking inputs |
| `make_length_control.py` | generates the composition-matched random controls (5 each at 6/9/10/15 residues) |
| `summarize_comp_top10.py` | summarises competition + study NK2R results (`NK2R_comparison.csv`) |
| `make_comp_table.py` | TableS4 |
| `make_fig5_comparison.py`, `make_fig6_headtohead.py` | Fig. S6, Fig. S7 |
