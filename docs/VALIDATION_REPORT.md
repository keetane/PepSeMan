# MDM2 linear baseline validation

Run date: 2026-10-28. Runtime: Python 3.12.9, RDKit 2026.09.1, macOS Apple Silicon.

## Source and roundtrip

| Measure | Result |
|---|---:|
| CSV measurement rows | 492 |
| Columns | 49 |
| Unique Molecule ChEMBL IDs | 404 |
| Target | CHEMBL5023 |
| RDKit-readable source SMILES | 492 / 492 |
| CSV→2D SDF records written/read | 492 / 492 |
| Canonical isomeric SMILES or non-empty property mismatches | 0 |

## Formal conservative decomposition

| Measure | Result |
|---|---:|
| Unique structures | 404 |
| Accepted unique 5mer+ structures | 49 |
| Accepted measurement rows | 50 |
| Strict structure MMP pairs | 29 |
| Pairs with at least one exact same-assay/same-endpoint pChEMBL comparison | 9 |
| First-pass elapsed time | 0.854 s |

Measurement-row statuses were: 50 accepted, 138 short peptide, 232 non-peptide, 1 non-linear, 64 ambiguous/unsupported backbone, and 7 multiple-fragment. The earlier 68-row heuristic is intentionally not treated as ground truth: the formal parser accepted 50 rows and refused ambiguous candidates rather than selecting a path arbitrarily.

Accepted unique length distribution: 5mer 1, 6mer 1, 8mer 2, 9mer 1, 10mer 1, 12mer 25, 13mer 1, 14mer 1, 15mer 1, and 17mer 15.

The first pass computed 404 unique decomposition cache entries and reused 88 repeated measurement-row structures during the same import. It created 712 indexed site contexts. Activity-only fixture updates recomputed zero structures and zero context buckets; identical file reimports queued no structural work. Full regression commands and assertions are in `tests/`.

## Unresolved scope

The baseline above predates the cyclic parser. It remains a historical linear result, not a cyclic classification benchmark. Ambiguous branched paths, automatic salt parent selection, residue insertion/deletion MMPs, and arbitrary noncanonical D/L reference classes remain unsupported.

## Cyclic implementation checks (2026-10-29)

Synthetic 5mer+ tests cover head-to-tail, thioether, disulfide, ring C=C, repeated and mixed bridges, head-to-tail plus side-chain bridges, cap–side-chain thioether, E/Z preservation, cyclic rotation, strict one-residue MMPs, and ambiguity/false-positive rejection. The 2026-10-29 check passed 40 non-CSV Python tests and 7 browser E2E tests. The three CSV baseline tests could not run at that time because the root-level CSV files had been moved.

A **read-only** comparison of the 2,535 previously stored `backbone-v1` structures against the new parser found all 712 previously accepted structures retained the same status, length, and residue variant hashes. It newly accepted 51 formerly unsupported/non-linear structures with inferred cyclic topologies: 28 disulfide-only, 7 ring C=C, 3 thioether, 3 head-to-tail only, and 10 head-to-tail plus disulfide. These are structural classifications, **not human-validated ground truth or accuracy measurements**. The remaining status changes were 4 structures newly classified as short peptides and 3 as unsupported. No live project database was migrated during that earlier check.

The “olefin metathesis” UI category is inferred from a ring C=C bridge; it does not establish the synthesis route. Hydrogenated bridges are not assigned that category. Human-reviewed labels for the user's real cyclic dataset are still needed before any classification-accuracy claim.

## Local source import and analysis (2026-10-30)

The eight user-supplied CSVs in `data/source/` were named `20261030_filename.csv`. `20261030_oxitosin_SS.csv` is byte-identical to `20261030_Oxytocin_S.csv`, so the two were retained locally but counted only once in the `ring-Peptide` import. The five distinct literature files were combined into one 119-row snapshot with each originating filename retained in `Source File`. A consistent SQLite backup was saved before changing imports, and the preexisting MDM2 and GLP-1 import history was retained. All three revisions reported in this section used `backbone-v2` and `strict-topology-v2` at the time of measurement.

| Project | Source rows | Unique structures | Accepted unique structures | MMP pairs |
| --- | ---: | ---: | ---: | ---: |
| `20261028_MDM2` | 492 | 404 | 56 | 32 |
| `20261028_GLP-1` | 2,131 | 2,131 | 707 | 1,466 |
| `ring-Peptide` | 119 | 89 | 88 | 31 |

The `ring-Peptide` rows comprise 13 Grb2-SH2, 16 MEGA-A4, 35 Oxytocin, 7 RGD, and 48 UbE2 measurements. R decomposition accepted 118 rows; one RGD entry was excluded as a peptide shorter than five residues. The accepted rows contained 6 head-to-tail, 28 thioether, 35 disulfide, and 43 ring C=C classifications. The 31 distinct structure MMPs were 10 head-to-tail, 16 thioether, and 5 disulfide. The 43 accepted ring C=C structures had distinct topology signatures, so no same-topology one-residue pair was eligible for an MMP in this sample. This is not evidence of a broken MMP engine or of the absence of experimental analogues.

The 43 Python tests passed, including the three local CSV checks. In-process API checks returned HTTP 200 for decomposition and pair lists across all three projects, all four cyclic-type filters, a cyclic pair detail, and MMP sorting. SQLite integrity was `ok`; all three reanalysis/import jobs completed. Real-data classifications were **not** compared with human-reviewed reference labels, so classification accuracy is not established.

## Local RAP1 triazole check (2026-10-30)

The additional local `data/source/RAP1.csv` contains 25 ChEMBL activity rows from document `CHEMBL3751807`. On the `backbone-v3` / `strict-topology-v3` implementation, a fresh isolated SQLite import accepted all 25 structures. Twelve had an N-substituted 1,2,3-triazole bridge, and thirteen had no detected cyclization link. The strict MMP engine produced 22 pairs: 16 within the triazole topology and 6 among the linear structures. It made no cross-topology pair. The source file was not copied into the Git-tracked tree. These are parser outputs against unlabeled real structures, **not** a human-validated classification-accuracy result. The already-running local app retains its previous revision until its queued reanalysis is processed by the worker.

Synthetic regression tests also check a triazole-bridged 5mer, its one-residue MMP partner, an otherwise similar pendant triazole that must not be classified as cyclic, and invariance to SMILES atom order. The structural category does not prove a click-chemistry synthesis route, and the current SMARTS is limited to an N-substituted aromatic 1,2,3-triazole motif.

