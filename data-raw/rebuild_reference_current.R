# Rebuild reference_set from cached Mode-A sites after downstream-method updates.
# No transcriptome rescan is required; raw site calls remain unchanged.
# Paths come from _paths.R; see that file for the environment variables.
.a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
source(file.path(if (length(.a)) dirname(sub("^--file=", "", .a[1])) else "data-raw",
                 "_paths.R"))
load_pkg_source()
setwd(PKG_HOME)
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr) })

cache <- readRDS("data-raw/reference_set_build/build_intermediates.rds")
sites_all <- cache$sites_all
# Gene flags are sequence-independent and retained from the cached enrichment.
enriched <- cache$enriched
ortholog <- read_csv(data_path("ortholog_biomart.csv"), show_col_types = FALSE)

ranking <- rank_sirna(enriched, reference = "self")
species <- species_match(enriched,
  candidates = c("mouse","rat","cyno","rhesus","rabbit","dog"),
  ortholog_table = ortholog)

dims <- names(default_weights())
tier_thresholds <- lapply(dims, function(d) {
  v <- ranking[[d]]
  c(min = min(v), max = max(v),
    q25 = quantile(v, 0.25, names = FALSE),
    q50 = quantile(v, 0.50, names = FALSE),
    q75 = quantile(v, 0.75, names = FALSE))
})
names(tier_thresholds) <- dims
tier_thresholds <- as.data.frame(do.call(rbind, tier_thresholds))
tier_thresholds$dimension <- rownames(tier_thresholds)
tier_thresholds <- tier_thresholds[, c("dimension","min","max","q25","q50","q75")]

gene_summary <- enriched |>
  group_by(.data$siRNA_name) |>
  summarise(n_critical_genes = n_distinct(.data$gene_name[
    .data$cancer_gene_v2 | .data$AE_gene | .data$immune_gene_immport]),
    .groups = "drop")

# Preserve metadata and provenance from the shipped object; replace only
# quantities affected by the downstream algorithm and input-table changes.
reference_set$metrics <- ranking
reference_set$tier_thresholds <- tier_thresholds
reference_set$species_score <- species$species_score
reference_set$gene_summary <- gene_summary
reference_set$expr_summary <- NULL
reference_set$risk_distribution <- enriched |>
  count(.data$siRNA_name, .data$risk_score, name = "n") |>
  complete(.data$siRNA_name,
           risk_score = sort(unique(enriched$risk_score)), fill = list(n = 0L))
reference_set$built_at <- as.character(Sys.time())
reference_set$pkg_version <- read.dcf(pkg_path("DESCRIPTION"), "Version")[1, 1]
reference_set$n_siRNAs <- nrow(reference_set$metadata)
attr(reference_set, "species_method") <- paste(
  "unique critical genes equally weighted; Ensembl Genes 116 BioMart",
  "downloaded 2026-09-14; human annotation Ensembl 113")

stopifnot(nrow(reference_set$metrics) == 94L,
          nrow(reference_set$tier_thresholds) == 4L,
          identical(reference_set$tier_thresholds$dimension, dims),
          nrow(reference_set$species_score) == 94L,
          is.null(reference_set$expr_summary))
save(reference_set, file = "data/reference_set.rda", compress = "xz")
saveRDS(list(sites_all = sites_all, enriched = enriched,
             ranking = ranking, species = species),
        "data-raw/reference_set_build/build_intermediates.rds", compress = "xz")
cat("reference_set rebuilt\n")
cat("thresholds:", paste(reference_set$tier_thresholds$dimension, collapse=", "), "\n")
cat("species score medians:\n")
print(vapply(reference_set$species_score[-1], median, numeric(1)))
