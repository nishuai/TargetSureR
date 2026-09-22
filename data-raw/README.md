# data-raw/

How the datasets bundled with TargetSureR were produced. These scripts are
excluded from the installed package (`.Rbuildignore`); they exist so that
each `data/*.rda` can be traced to its source rather than taken on trust.

They are provenance records, not a turnkey build. Two of them read curated
source tables that are not redistributed with the package, and all of them
need reference data that must be downloaded first — see
`vignette("01_quickstart")`, section *"What you need before you start"*.

## Paths

No script contains an absolute path. Locations resolve through `_paths.R`
from three environment variables, each with a documented fallback:

| Variable | Holds | Default |
|---|---|---|
| `TARGETSURER_HOME` | package source root | parent of `data-raw/` |
| `TARGETSURER_DATA` | prepared reference data (EnsDb, annodb, expression matrix, orthologs) | `~/TargetSureR_data` |
| `TARGETSURER_INPUT` | curated source tables (siRNA cohort, gene-list sources) | `../siRNA_prediction_model` |

```r
Sys.setenv(TARGETSURER_DATA = "~/TargetSureR_data")
source("data-raw/_paths.R")   # defines pkg_path(), data_path(), input_path()
```

## Scripts

| Script | Produces | Needs |
|---|---|---|
| `_paths.R` | — (sourced by the others) | — |
| `make_internal_data.R` | the bundled `.rda` files | — |
| `rebuild_gene_lists.R` | `ae_genes`, `cancer_genes_v2`, `immune_genes_immport` | source gene-list tables in `TARGETSURER_INPUT` |
| `build_annodb_all.R` | three cached `annodb` objects | BSgenome + EnsDb |
| `build_reference_set.R` | `reference_set` (the 94-siRNA cohort) | annodb, expression matrix, ortholog table, siRNA source table |
| `rebuild_reference_current.R` | refreshed `reference_set` from cached sites | cached build intermediates, ortholog table |
| `reflag_reference.R` | re-applied gene-list flags on cached sites | cached build intermediates |

### Bundled data and its status

`ae_genes` (35 symbols), `cancer_genes_v2` (979), `immune_genes_immport`
(1,794) and `reference_set` (94 siRNAs) are complete and used for analysis.

`gtex_example` and `ortholog_table_example` are 20 genes of **simulated**
values, present so that examples and tests run without a multi-gigabyte
download. They are not measured data and must be replaced for any real
analysis; the help pages for `enrich_annotations()`, `species_match()` and
`load_gtex_example()` say so at the point of use.

### Reference cohort

`build_reference_set.R` scans 94 synthetic siRNA duplexes curated from the
MIT/ICBP siRNA Database, normalised to 19 nt guides by `normalize_guide()`.
End-to-end runtime is roughly an hour on one workstation, dominated by
`SeedMatchR::SeedMatchR()` on the larger feature classes. Its input sequence
table is not redistributed here, so treat the script as a record of how the
bundled cohort was built rather than a runnable recipe. To normalise against
a different cohort, scan your own guides and pass the result via
`rank_sirna(reference = "custom")`; the vignette section *"Using a different
reference cohort"* gives the pattern.

### Not tracked

Cached build artefacts stay local (`.gitignore`): the three `annodb` `.rds`
files (~24 MB) and `reference_set_build/` intermediates (~137 MB). Rebuild
them with `build_annodb_all.R` and `build_reference_set.R`.
