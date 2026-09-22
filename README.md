# TargetSureR

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![R >= 4.1.0](https://img.shields.io/badge/R-%3E%3D%204.1.0-1f425f.svg)](https://www.r-project.org/)
<!-- badges: end -->

**TargetSureR** is a modular R toolkit for evaluating small interfering
RNA (siRNA) off-target effects end-to-end: from raw sequence or
predicted sites, through curated biological annotation, to risk
ranking, tissue-level safety profiling, and cross-species animal-model
recommendation — all in a single consistent pipeline.

## Why TargetSureR

| | Existing tools | TargetSureR |
|---|---|---|
| Accepts siRNA sequence | Some (SeedMatchR, siDirect) | ✅ via SeedMatchR backend |
| Accepts pre-computed `(transcript, position)` | Rarely | ✅ direct input |
| Transcript-relative region (5'UTR / CDS / 3'UTR) | Usually not | ✅ via Ensembl |
| Curated gene-list enrichment (oncogenes, AE genes, immune) | Manual | ✅ built-in |
| GTEx tissue-level safety scoring | No | ✅ 68 tissues |
| Cross-species animal-model recommendation | No | ✅ Layer-1 ortholog |
| CRAN-ready, pure R | Mixed | ✅ no shell tools required |

## Installation

```r
# install.packages("remotes")
remotes::install_github("nishuai/TargetSureR")
```

Heavy dependencies are *optional*. Install only what you need:

```r
# For Mode B (transcript + position input):
BiocManager::install(c("ensembldb", "AnnotationFilter",
                        "IRanges", "S4Vectors", "BiocGenerics"))

# For Mode A (siRNA sequence input):
BiocManager::install(c("SeedMatchR", "Biostrings",
                        "AnnotationHub", "GenomicRanges"))

# For REST API fallback (resolves newer transcript IDs not in local EnsDb):
install.packages("httr2")
```

## Pipeline at a glance

```
┌─────────────────────────────────────────────────────────────────┐
│  INPUT                                                           │
│                                                                  │
│  Mode A: siRNA sequence           Mode B: (transcript, position) │
│     list(siRNA_name,                 data.frame(siRNA_name,      │
│          sequence)                                transcript,    │
│                                                    position)     │
└──────────┬──────────────────────────────┬───────────────────────┘
           │                              │
           │      annotate_sites()        │
           ▼                              ▼
    ╔══════════════════════════════════════════╗
    ║  10-column core off-target table          ║
    ║  (siRNA_name, strand, transcript, ...)   ║
    ╚═════════════════╤════════════════════════╝
                      │
            enrich_annotations()
                      │
    ╔═════════════════▼════════════════════════╗
    ║  81-column enriched table                 ║
    ║  + 3 gene-list flags                      ║
    ║  + 68 GTEx tissue TPM columns             ║
    ╚═════════════════╤════════════════════════╝
                      │
    ┌─────────────────┼─────────────────┐
    ▼                 ▼                 ▼
rank_sirna()   tissue_safety()    species_match()
per-siRNA      per-gene           per-siRNA ×
composite      bio-attention      per-species
risk score     score              conservation
    │                 │                 │
    └────────┬────────┴─────────┬───────┘
             ▼                  ▼
       generate_report()   CSV / plots
```

## Two input modes — same output shape

**Mode A: from siRNA sequence**

```r
library(TargetSureR)

sites <- annotate_sites(
  input   = list(siRNA_name = "hsiR_PMP22",
                 sequence   = "UCCUCUCCAUUUGAGAUAGUU"),
  mode    = "sequence",
  species = "human"     # maps to hg38 for SeedMatchR
)
```

**Mode B: from `(transcript, position)` pairs**

```r
sites <- annotate_sites(
  input = data.frame(
    siRNA_name = "hsiR22",
    transcript = "ENST00000426385.4",
    position   = 379
  ),
  mode      = "position",
  ensdb_path = "~/ensdb/ensembl_113_hsapiens.sqlite"
)
```

**Both modes return the same 10-column schema:**

```
siRNA_name   siRNA   strand   transcript    gene_name
biotype      position region   match_type    risk_score
```

## Choose your Ensembl release

You are never locked to an old Ensembl version. Any `.sqlite` file that
`ensembldb::EnsDb()` can read — from AnnotationHub, Bioconductor data
packages, or future releases — can be dropped in via `ensdb_path`:

```r
# Pin to a specific release you downloaded (recommended for reproducibility)
annotate_sites(...,
  ensdb_path = "~/ensdb/ensembl_113_hsapiens.sqlite")

# Let AnnotationHub fetch the latest (needs network)
annotate_sites(..., species = "human")
```

For transcripts newer than the local EnsDb (e.g., Ensembl 114+ IDs),
TargetSureR automatically queries the Ensembl REST API as a fallback
(requires `httr2`). Typical resolution rate after fallback: >95%.

## End-to-end example

```r
library(TargetSureR)

# 1. Annotate (Mode B with REST fallback for newer transcripts)
sites <- annotate_sites(
  data.frame(
    siRNA_name = c("hsiR22", "hsiR22", "hsiR21"),
    transcript = c("ENST00000426385.4",
                   "ENST00000335508.11",
                   "ENST00000380152.7"),
    position   = c(379L, 3136L, 500L)),
  mode       = "position",
  ensdb_path = "~/ensdb/ensembl_113_hsapiens.sqlite")

# 2. Enrich with gene-list flags and tissue expression
enriched <- enrich_annotations(sites)

# 3. Rank siRNAs
rank_sirna(enriched, target_sirnas = "hsiR22")

# 4. Tissue safety
tissue_safety(enriched, primary_tissue = "Nerve_Tibial")

# 5. Cross-species recommendation (6 model species)
species_match(sites,
              candidates = c("mouse", "rat", "cyno",
                             "rhesus", "rabbit", "dog"))

# 6. One-call bundle
generate_report(enriched, output_dir = "results/",
                target_sirnas = "hsiR22")
```

## Modules and public API

| Function                    | Role                                             |
|-----------------------------|--------------------------------------------------|
| `build_annodb()`            | Build the transcript sequence index for Mode A   |
| `normalize_guide()`         | Reported strand -> canonical 19 nt guide (antisense) |
| `annotate_sites()`          | Dispatcher: picks Mode A or Mode B              |
| `annotate_from_sequence()`  | Mode A — scan transcriptome via SeedMatchR      |
| `annotate_from_position()`  | Mode B — resolve metadata via ensembldb         |
| `enrich_annotations()`      | + gene-list flags + GTEx tissue TPM             |
| `compute_risk_score()`      | Rule-based `risk_score` (region × match_type)   |
| `rank_sirna()`              | count endpoints, cohort percentiles, optional composite |
| `tissue_safety()`           | Per-gene `bio_attention`                        |
| `species_match()`           | Per-species conservation score                  |
| `generate_report()`         | Bundle CSVs                                     |
| `default_weights()`         | four equal composite weights (0.25 each)        |
| `load_gene_lists()`         | Return built-in gene lists                      |
| `load_gtex_example()`       | Simulated example expression matrix (not real data) |
| `core_columns()` / `enriched_columns()` | Column-name specifications             |

## Built-in reference data

| Data        | Contents                                                     |
|-------------|--------------------------------------------------------------|
| `ae_genes`              | Adverse-event gene symbols                        |
| `cancer_genes_v2`       | Curated cancer gene panel                         |
| `immune_genes_immport`  | ImmPort-derived immune-related panel              |
| `reference_set`          | 94-siRNA MIT/ICBP normalisation cohort (19 nt guides) |
| `gtex_example`          | **Simulated** 20-gene × 10-tissue matrix — examples only |
| `ortholog_table_example`| **Simulated** 20-gene table — examples only        |

The three gene lists and `reference_set` are complete and ready for
analysis. The two `*_example` datasets are **not real data**: the numbers
are invented, so that examples run without a download. Scoring a real
off-target set against them returns a table in which almost every gene is
unquantified rather than an error. Replace both before acting on any
`tissue_safety()` or `species_match()` output.

All bundled data total ~18 KB (CRAN-friendly).

## External data you must supply

TargetSureR ships no genome, transcriptome or expression data. Four
external resources are needed for full functionality, assembled once:

| Resource | Size | Needed for | Source |
|---|---|---|---|
| `BSgenome.Hsapiens.UCSC.hg38` | ~800 MB | Mode A scanning | Bioconductor |
| `EnsDb` SQLite (Ensembl 113 = `AH119325`) | ~510 MB | Mode A and Mode B | AnnotationHub |
| Tissue median-expression matrix | ~20 MB | `tissue_safety()` | GTEx portal or in-house |
| Ortholog table | ~200 KB | `species_match()` | Ensembl BioMart |

The transcript sequence index (`annodb`) is then built locally with
`build_annodb()`, about 25 s per feature class, and cached.

Partial setups work: cohort percentiles from `rank_sirna()` need only the
genome, EnsDb and `annodb`, because `reference_set` is bundled. See the
**"What you need before you start"** section of
`vignette("01_quickstart")` for step-by-step retrieval, a suggested
directory layout, and what each missing resource costs you.


## Related tools and citations

- **SeedMatchR** — Cazares T. et al. (2024). *SeedMatchR: identify
  off-target effects mediated by siRNA seed regions in RNA-seq
  experiments*. Bioinformatics Advances, 4(1), vbae002.
- **ensembldb** — Rainer J. et al. (2019). *ensembldb: an R package to
  create and use Ensembl-based annotation resources*. Bioinformatics,
  35(17), 3151-3153.
- **GTEx** — GTEx Consortium (2020). *The GTEx Consortium atlas of
  genetic regulatory effects across human tissues*. Science, 369(6509),
  1318-1330.

## License

MIT © 2026 Shuai Ni

## Contact

Shuai Ni &lt;nishuai@wakerbio.com&gt; · Waker Bio

Bug reports and feature requests:
[github.com/nishuai/TargetSureR/issues](https://github.com/nishuai/TargetSureR/issues)
