# Re-apply gene-list flags to the cached reference-cohort sites after the
# curated lists changed. Site calls themselves are unaffected, so no rescan.
suppressPackageStartupMessages({ library(readr); library(dplyr) })
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)
for (f in list.files("data", pattern = "[.]rda$", full.names = TRUE))
  load(f, envir = .GlobalEnv)

cache <- readRDS("data-raw/reference_set_build/build_intermediates.rds")
sites <- cache$sites_all
gtex <- read_tsv(file.path("..", "siRNA_prediction_model", "paper", "realrun",
                           "rerun19", "gtex_v10_median_tpm.tsv"),
                 show_col_types = FALSE)

enriched <- enrich_annotations(sites, expression_matrix = gtex)
crit <- enriched$cancer_gene_v2 | enriched$AE_gene | enriched$immune_gene_immport
cat(sprintf("reference sites: %d | critical sites: %d | unique critical genes: %d\n",
            nrow(enriched), sum(crit),
            length(unique(toupper(enriched$gene_name[crit])))))

saveRDS(list(sites_all = sites, enriched = enriched),
        "data-raw/reference_set_build/build_intermediates.rds", compress = "xz")
cat("cached enrichment refreshed\n")
