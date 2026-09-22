# Build internal datasets shipped with TargetSureR.
#
# This script is a template. It builds four .rda files into data/:
#   ae_genes.rda             - adverse-event gene symbols
#   cancer_genes_v2.rda      - curated cancer gene panel
#   immune_genes_immport.rda - ImmPort-curated immune panel
#   gtex_example.rda         - small GTEx v10 example subset
#   ortholog_table_example.rda - tiny ortholog example for species_match()
#
# Run from the package root with: Rscript data-raw/make_internal_data.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

out_dir <- "data"
dir.create(out_dir, showWarnings = FALSE)

# -- gene lists -----------------------------------------------------------
#
# Read from the curated source tables rather than an inline panel, so the
# shipped lists match the files documented in the manuscript. Override the
# location with TARGETSURE_INPUTDATA if the sources sit elsewhere.

input_dir <- Sys.getenv("TARGETSURE_INPUTDATA", unset = "")
if (!nzchar(input_dir))
  input_dir <- file.path("..", "siRNA_prediction_model", "inputdata")
if (!dir.exists(input_dir))
  stop("gene-list source directory not found: ", input_dir, call. = FALSE)

read_symbols <- function(file, column) {
  path <- file.path(input_dir, file)
  if (!file.exists(path)) stop("missing gene-list source: ", path, call. = FALSE)
  x <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  if (!column %in% names(x))
    stop(sprintf("column '%s' not found in %s", column, file), call. = FALSE)
  s <- toupper(trimws(as.character(x[[column]])))
  s <- s[!is.na(s) & nzchar(s) & s != "NA"]
  sort(unique(s))
}

# Adverse-event genes: reported insertional-mutagenesis / AE associations.
ae_genes <- read_symbols("AE_GeneList.txt", "Gene")

# Cancer gene panel (v2), keyed on the queried input symbol.
cancer_genes_v2 <- read_symbols("Cancer_gene_IDv2.txt", "Input")

# ImmPort immune gene list (Gene Summary export; 18 functional categories).
immune_genes_immport <- read_symbols("Immune_GeneList.txt", "Symbol")

message(sprintf(
  "gene lists: AE %d, cancer %d, immune %d (union %d unique)",
  length(ae_genes), length(cancer_genes_v2), length(immune_genes_immport),
  length(unique(c(ae_genes, cancer_genes_v2, immune_genes_immport)))))

save(ae_genes,             file = file.path(out_dir, "ae_genes.rda"),             compress = "xz")
save(cancer_genes_v2,      file = file.path(out_dir, "cancer_genes_v2.rda"),      compress = "xz")
save(immune_genes_immport, file = file.path(out_dir, "immune_genes_immport.rda"), compress = "xz")

# -- GTEx example (small subset) ------------------------------------------

set.seed(42)
example_genes <- toupper(c("PMP22", "SF3B1", "CTNNB1", "PALB2", "TP53",
                            "KRAS", "EGFR", "MYC", "IL6", "TNF",
                            "MSH6", "NRAS", "RELA", "JAK2", "STAT3",
                            "BRCA1", "PTEN", "APC", "BRAF", "PIK3CA"))
example_tissues <- c("Nerve_Tibial", "Brain_Cortex", "Liver",
                     "Heart_Left_Ventricle", "Muscle_Skeletal",
                     "Kidney_Cortex", "Lung", "Whole_Blood",
                     "Spleen", "Adipose_Subcutaneous")

gtex_example <- tibble::tibble(gene_symbol = example_genes)
for (t in example_tissues) {
  gtex_example[[t]] <- round(stats::rlnorm(length(example_genes),
                                            meanlog = 1, sdlog = 1), 3)
}
save(gtex_example, file = file.path(out_dir, "gtex_example.rda"), compress = "xz")

# -- ortholog example -----------------------------------------------------

n_g <- length(example_genes)
ortholog_table_example <- tibble::tibble(
  human_symbol         = example_genes,
  mouse_ortholog_type  = rep("one2one", n_g),
  rat_ortholog_type    = rep("one2one", n_g),
  cyno_ortholog_type   = rep("one2one", n_g),
  rhesus_ortholog_type = rep("one2one", n_g),
  rabbit_ortholog_type = rep("one2one", n_g),
  dog_ortholog_type    = rep("one2one", n_g),
  mouse_identity_pct   = round(runif(n_g, 70, 99), 1),
  rat_identity_pct     = round(runif(n_g, 70, 99), 1),
  cyno_identity_pct    = round(runif(n_g, 85, 99), 1),
  rhesus_identity_pct  = round(runif(n_g, 85, 99), 1),
  rabbit_identity_pct  = round(runif(n_g, 75, 95), 1),
  dog_identity_pct     = round(runif(n_g, 75, 95), 1)
)
save(ortholog_table_example,
     file = file.path(out_dir, "ortholog_table_example.rda"),
     compress = "xz")

message("Wrote: ", paste(list.files(out_dir), collapse = ", "))
