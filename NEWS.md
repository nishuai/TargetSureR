# TargetSureR 1.1.0

## New function
- `build_annodb()` builds the transcript sequence index that
  `annotate_sites()` needs in sequence mode, for one feature class
  (`"3UTR"`, `"5UTR"` or `"CDS"`). This step previously existed only as
  code to be copied out of the vignette, so a fresh installation could
  not run Mode A without it. The four transcript filters that define the
  scanned set — protein-coding biotype, sequence names restricted to
  `1:22`/`X`/`Y`/`MT`, longest sequence per gene within the class, and a
  minimum length of 8 nt — are documented and adjustable.
  With Ensembl 113 the result is 19,022 3'UTR, 18,973 5'UTR and 19,355
  coding sequences, byte-identical to the databases used for the
  published analyses.

## Documentation
- The vignette gains a **"What you need before you start"** section: the
  four external resources with their sizes and sources, a suggested
  directory layout, and a table of what remains usable when one of them
  is absent. Retrieval is now documented step by step, including the
  manual-download route for an `EnsDb` whose AnnotationHub download
  fails partway, GTEx GCT conversion, and a runnable BioMart ortholog
  download.
- `gtex_example` and `ortholog_table_example` are now described as what
  they are: 20 genes of **simulated** values for documentation only.
  Scoring real off-targets against them returns an unquantified table
  rather than an error, so the help pages for `enrich_annotations()`,
  `species_match()` and `load_gtex_example()` now say so at the point of
  use. The README previously called `gtex_example` a GTEx subset, which
  it is not.
- `reference_set`'s documented format no longer lists an `expr_summary`
  element, which the object does not contain.
- Guidance added on substituting your own reference cohort via
  `rank_sirna(reference = "custom")`, and on `SeedMatchR` versions: CRAN
  carries 1.1.1, while 2.0.0 (used for the published analyses) comes
  from the authors' GitHub repository.

## Development scripts
- Scripts under `data-raw/` no longer carry absolute paths or a
  hard-coded R library version. Locations resolve through
  `data-raw/_paths.R` from the environment variables
  `TARGETSURER_HOME`, `TARGETSURER_DATA` and `TARGETSURER_INPUT`, each
  with a documented default.

# TargetSureR 1.0.0

Initial public release.

## Core pipeline
- `annotate_sites()` dispatcher supporting two input modes:
  - Mode A: siRNA sequence (via `SeedMatchR` + `Biostrings`)
  - Mode B: `(transcript, position)` (via `ensembldb`)
- `enrich_annotations()` appends three curated gene-list flags and a
  GTEx tissue TPM block to the core table.
- `rank_sirna()` produces count endpoints, cohort percentiles and an
  optional four-dimension composite score with tiered
  Low/Medium/High/Critical grading, normalised either within the supplied
  set (`reference = "self"`) or against a reference cohort
  (`reference = "builtin"` / `"custom"`).
- `tissue_safety()` computes a per-gene `bio_attention` score.
- `species_match()` recommends animal models by Layer-1 ortholog
  conservation across mouse, rat, cyno, rhesus, rabbit and dog.
- `normalize_guide()` converts a reported siRNA strand to the canonical
  19 nt guide: strips a written-out `dTdT`/`UU` 3' overhang,
  reverse-complements a sense strand when asked, and caps the length from
  the 3' end so the seed region (guide positions 2-8) is preserved. This
  is the same normalisation applied to the built-in reference cohort, so
  user siRNAs can be made like-for-like comparable before scanning.
- `generate_report()` bundles CSV outputs.

## Strand contract
`annotate_sites()` searches the transcriptome for
`reverseComplement(sequence)` and therefore requires the **guide
(antisense)** strand. A sense strand raises no error but silently yields
seed-only hits with no on-target match; the documentation now states the
contract explicitly and gives a self-check (a true guide produces a
`full_complementarity` hit with `n_mismatch == 0` on its own target gene).

## Configuration
- `ensdb_path` argument lets users pin any future Ensembl release without
  updating the package.
- REST API fallback for transcripts newer than the local EnsDb (`httr2`
  Suggest). When `rest_fallback = TRUE` (default), unresolved transcript
  IDs are re-queried against `https://rest.ensembl.org/lookup/id` with
  session-level caching.
- User-tunable composite weights via `default_weights()` (four equal
  dimensions by default).
- `with_lookup_source` in `annotate_from_position()` for provenance
  tracking (ensdb / rest / unresolved).

## Data
- `reference_set`: 94 experimentally reported human siRNAs curated from
  the MIT/ICBP siRNA Database (144 human-reactive entries filtered to the
  94 synthetic duplexes, excluding 50 long-hairpin shRNA constructs), each
  scanned across 3'UTR / 5'UTR / CDS against Ensembl v113. All guides are
  19 nt. `composite_score` is self-normalised within the cohort, i.e. on
  the same scale that `rank_sirna(reference = "builtin")` places a new
  siRNA on.
- Built-in gene lists: `ae_genes`, `cancer_genes_v2`,
  `immune_genes_immport`.
- Example resources: `gtex_example`, `ortholog_table_example`.

## Validation
- 0 errors / 0 warnings on `R CMD check --as-cran` (R 4.5.3).
- Regression tests pin the reference cohort's internal consistency: the
  shipped `composite_score` must equal the value recomputed from
  `tier_thresholds`, and those thresholds must be the min/max of the
  stored metrics.
- End-to-end real-data tests (see `data-raw/realtest_*.R`).
