# Step 3 · Structural-neighbour retrieval (Foldseek)

**Goal** — for each receptor, retrieve a redundancy-reduced panel of structurally similar proteins whose
sequences can be docked against candidate *and* control, so that every comparison in step 4 is matched.

> These panels are **structural neighbours, not validated pharmacological off-targets.** They are heterogeneous
> by construction — receptor structures, unrelated enzymes and, for NK1R, NK1R structures themselves. A high
> neighbour-score fraction is a computational scoring statement, not a safety or selectivity claim.

## How it was run

The bundled search script (`scripts/foldseek_search.py`) queries the **Foldseek web service** in combined
**3Di + amino-acid** mode against a list of nine database identifiers:

```
afdb50, BFVD, afdb-swissprot, afdb-proteome, bfmd, mgnify_esm30, pdb100, cath50, gmgcl_id
```

Its command-line default selects four of them — `pdb100,afdb50,afdb-swissprot,afdb-proteome` — and the
downstream tables in this directory are the **PDB100** result table only. Redundancy reduction used a CD-HIT
sequence-identity threshold of 0.50 when CD-HIT was available, and otherwise a greedy 3-mer Jaccard-similarity
threshold of 0.50: **these are different measures and are not interchangeable identity criteria.** Chain
selection took the Foldseek-matched chain, falling back to the first chain when the match could not be
resolved.

## `input/`

Empty. The search ran through the web service against the receptor chains in
[`../step1_peptide_generation/input/receptor_chains.fasta`](../step1_peptide_generation/input/receptor_chains.fasta);
no local query files were archived. The receptor structures themselves live in
[`../step2_positive_docking/input/pdbs/`](../step2_positive_docking/input/pdbs/).

## `output/<RECEPTOR>/`

| File | Contents |
|:--|:--|
| `pdb100_best_hits.csv` | One row per retrieved structural neighbour — best alignment per neighbour by bitscore, parsed from the 21-column Foldseek alignment table. Columns: `query, neighbour, pdb, chain, description, seq_identity_pct, aln_len, evalue, bits, taxname, in_negative_panel`. |
| `alis_pdb100_report.m8` | The taxonomy report shipped alongside the alignment table (percent / count / rank / taxid / taxon-name rows). It is named `.m8` but is a taxon report, not an alignment table. |
| `pdb_ids.txt` | The retrieved PDB identifiers, used to fetch the structures. **Absent for NK2R** — the NK2R search predates the per-receptor naming convention. |
| `alis_pdb100_cdhit50.fasta` | The redundancy-reduced representative sequences (CD-HIT-50, or the greedy 3-mer fallback). |

Condensed by [`scripts/harvest_foldseek_hits.py`](../scripts/harvest_foldseek_hits.py). A Foldseek target
string looks like `6e59-assembly1.cif.gz_A`; a repeated chain is suffixed (`…_A-2`) and the suffix is stripped
so the key matches the panel folder name.

## Result

| Receptor | Retrieved neighbours | In the step-4 panel | Panel FASTA | `pdb_ids.txt` |
|:--|--:|--:|--:|:--|
| NK1R | 892 | **66** | 255 seqs | 254 ids |
| NK2R | 728 | **27** | 142 seqs | — |
| NK3R | 869 | **60** | 215 seqs | 214 ids |

The `in_negative_panel` flag marks the neighbours that survived filtering into the docking panel. Those
counts — 66 / 27 / 60 — are the denominators of the reported comparison-count proxy defined in
[step 5](../step5_analysis/).

## Watch out

* **The raw Foldseek result archives and filtering logs were not preserved,** so the exact database selection,
  the intermediate exclusions and the effect of panel composition cannot be independently assessed. The
  scripts document the intended processing, not a verified execution history.
* **`pdb_ids.txt` counts differ from the FASTA counts** (254 vs 255 for NK1R) because one CD-HIT-50
  representative has no single downloadable PDB id. Neither number equals the 66/27/60 panel sizes — the panel
  was filtered further.
* **Annotation composition (Fig. S3) describes the supplied structure annotations** and does not establish
  pharmacological relevance.
