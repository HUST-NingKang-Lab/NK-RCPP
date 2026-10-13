# Step 1 · Candidate peptide generation (PepMLM)

**Goal** — produce a receptor-specific candidate library: 32 ten-residue peptides per receptor.

## How it was run

PepMLM was used through **its authors' public Colab notebook** (the online version), not a local
install — the receptors were pasted in as the conditioning input, peptide length was fixed at **10**,
and the notebook returned a ranked list of designs with a pseudo-perplexity score per sequence.

The Colab environment is not pinned here because nothing in this repository depends on it: the deliverable
is the sequence list, and everything downstream is a deterministic function of that list.

## `input/`

| File | Contents |
|:--|:--|
| `receptor_chains.fasta` | The target chain of each receptor, extracted from the structure files used everywhere downstream. NK1R: PDB 6E59 chain A (construct numbering 28–1196). NK2R: AlphaFold2 model of UniProt **P21452**, chain A, residues 1–398. NK3R: predicted model, chain B, residues 76–369. |

> These sequences were re-extracted from the receptor files for the record. The exact string that was
> pasted into the Colab notebook was not archived.

## `output/`

| File | Rows | Contents |
|:--|--:|:--|
| `NK1R_pepmlm_candidates.csv` | 33 | 32 candidates + control, with `pseudo_perplexity` |
| `NK2R_pepmlm_candidates.csv` | 33 | 32 candidates + control, **`pseudo_perplexity` empty** — the NK2R Colab export was not archived, so the sequences were recovered from the docking input file instead |
| `NK3R_pepmlm_candidates.csv` | 33 | 32 candidates + control, with `pseudo_perplexity` |

Columns: `receptor, peptide_id, sequence, control_name, pseudo_perplexity, role`.

Built by [`scripts/make_candidate_tables.py`](../scripts/make_candidate_tables.py), which also asserts
that the campaign numbering here (`pep1…pep32`) is identical to the numbering used by the docking input
files — for NK1R and NK3R it checks the Colab export against `NK1R_peptide.csv` / `NK3R_peptide.csv`
sequence by sequence.

## Receptor-matched controls

Every receptor carries its natural ligand through the pipeline as the reference point:

| Receptor | Control | Sequence | Source |
|:--|:--|:--|:--|
| NK1R | Substance P | `RPKPQQFFGLM` | endogenous tachykinin |
| NK2R | Neurokinin A derived | `HKTDSFVGLM` | NKA (10-mer) |
| NK3R | Neurokinin B | `DMHDFFVGLM` | endogenous tachykinin |

## Watch out

* **This numbering is load-bearing.** `pep1…pep32` here is the `NK2R-01…NK2R-32` code system used in the
  figures, in TableS1 and in the positive-docking results. The negative-docking stage uses a *different*
  NK2R labelling — see [step 4](../step4_negative_docking/README.md).
* Pseudo-perplexity is a language-model confidence, not an affinity estimate. It was not used as a filter.
