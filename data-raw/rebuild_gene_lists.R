# Rebuild only the three curated gene-list datasets from their source tables.
# The example GTEx and ortholog datasets are left untouched.
suppressPackageStartupMessages(library(readr))

input_dir <- Sys.getenv("TARGETSURE_INPUTDATA", unset = "")
if (!nzchar(input_dir))
  input_dir <- file.path("..", "siRNA_prediction_model", "inputdata")
stopifnot(dir.exists(input_dir))

read_symbols <- function(file, column) {
  x <- read_tsv(file.path(input_dir, file), show_col_types = FALSE,
                progress = FALSE)
  stopifnot(column %in% names(x))
  s <- toupper(trimws(as.character(x[[column]])))
  sort(unique(s[!is.na(s) & nzchar(s) & s != "NA"]))
}

ae_genes             <- read_symbols("AE_GeneList.txt", "Gene")
cancer_genes_v2      <- read_symbols("Cancer_gene_IDv2.txt", "Input")
immune_genes_immport <- read_symbols("Immune_GeneList.txt", "Symbol")

save(ae_genes,             file = "data/ae_genes.rda",             compress = "xz")
save(cancer_genes_v2,      file = "data/cancer_genes_v2.rda",      compress = "xz")
save(immune_genes_immport, file = "data/immune_genes_immport.rda", compress = "xz")

cat(sprintf("AE %d | cancer %d | immune %d | union %d\n",
            length(ae_genes), length(cancer_genes_v2),
            length(immune_genes_immport),
            length(unique(c(ae_genes, cancer_genes_v2, immune_genes_immport)))))
