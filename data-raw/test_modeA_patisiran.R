# Mode A test: 2 siRNAs
# 1. Patisiran (ALN-18328 / ONPATTRO) — guide AUGGAAUACUCUUGGUUACTT
#    Source: drugs.com Patisiran monograph; FDA-approved 2018 for hATTR
# 2. siRNA #2 (user-provided)        — guide UCCUGUUGCUGAGUAUCAUTT

# (user library resolved by R itself; no .libPaths() override)

setwd(Sys.getenv("TARGETSURER_HOME", unset = "."))
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
for (f in list.files("data", pattern = "\\.rda$", full.names = TRUE)) {
  load(f, envir = .GlobalEnv)
}

suppressPackageStartupMessages({
  library(dplyr); library(tibble)
})

out_dir <- "demo_output/modeA_test"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sirna_list <- list(
  list(name = "Patisiran",  seq = "AUGGAAUACUCUUGGUUACTT"),
  list(name = "siRNA_user", seq = "UCCUGUUGCUGAGUAUCAUTT")
)

# Pre-load annotation once (SeedMatchR: 3UTR features for hg38)
cat("=== Loading hg38 3UTR annotation (one-time) ===\n")
t0 <- Sys.time()
annodb <- SeedMatchR::load_annotations(
  reference.name    = "hg38",
  feature.type      = "3UTR",
  canonical         = FALSE,
  min.feature.width = 8,
  longest.utr       = TRUE,
  return_gene_name  = FALSE
)
cat(sprintf("annodb loaded in %.1f s; %d transcripts\n",
            as.numeric(Sys.time()-t0, units="secs"),
            length(annodb$seqs)))

all_sites <- list()
for (s in sirna_list) {
  cat(sprintf("\n=== %s ===\n", s$name))
  cat(sprintf("Guide (5'->3'): %s (%d nt)\n", s$seq, nchar(s$seq)))

  t0 <- Sys.time()
  sites <- annotate_sites(
    input  = list(siRNA_name = s$name, sequence = s$seq),
    mode   = "sequence",
    species = "human",
    feature_type = "3UTR",
    seed_name = "mer7m8",
    max_mismatch_full = 3L,
    annodb = annodb
  )
  cat(sprintf("  hits: %d sites, %d transcripts, %d genes (%.1f s)\n",
              nrow(sites), n_distinct(sites$transcript),
              n_distinct(sites$gene_name),
              as.numeric(Sys.time()-t0, units="secs")))
  if (nrow(sites) > 0) {
    cat("  match_type:\n"); print(table(sites$match_type))
  }
  all_sites[[s$name]] <- sites
}

sites_all <- bind_rows(all_sites)
write.csv(sites_all, file.path(out_dir, "01_sites_core.csv"),
          row.names = FALSE)

if (nrow(sites_all) == 0) {
  cat("\nNo hits across both siRNAs. Stopping.\n"); quit()
}

# Step 2: enrich
cat("\n=== Step 2: enrich_annotations() ===\n")
enriched <- enrich_annotations(sites_all)
n_crit <- sum(enriched$cancer_gene_v2 | enriched$AE_gene |
              enriched$immune_gene_immport)
cat(sprintf("enriched: %d rows; critical: %d (%.1f%%)\n",
            nrow(enriched), n_crit, 100*n_crit/nrow(enriched)))
write.csv(enriched, file.path(out_dir, "02_sites_enriched.csv"),
          row.names = FALSE)

# Step 3: rank (joint normalization across both siRNAs)
cat("\n=== Step 3: rank_sirna() ===\n")
ranking <- rank_sirna(enriched,
                      target_sirnas = sapply(sirna_list, `[[`, "name"))
print(ranking[, c("siRNA_name","composite_score","rank","risk_tier",
                  "high_risk_critical","critical_gene_count",
                  "full_comp_count","cds_count","total_offtargets")])
write.csv(ranking, file.path(out_dir, "03_ranking.csv"), row.names = FALSE)

# Step 4: tissue safety
cat("\n=== Step 4: tissue_safety() ===\n")
tissue <- tissue_safety(enriched, primary_tissue = "Nerve_Tibial")
cat("Top 10 by bio_attention:\n")
print(head(tissue[, c("gene_name","n_hits","primary_tpm","max_risk_score",
                      "bio_attention","is_key_gene")], 10))
write.csv(tissue, file.path(out_dir, "04_tissue_safety.csv"),
          row.names = FALSE)

# Step 5: species match
cat("\n=== Step 5: species_match() ===\n")
sp <- species_match(sites_all,
                    candidates = c("mouse","rat","cyno",
                                   "rhesus","rabbit","dog"))
cat("species_score:\n"); print(sp$species_score)
cat("\nrecommendation:\n")
cat(paste0("  ", sp$recommendation, collapse="\n"), "\n")
write.csv(sp$species_score, file.path(out_dir, "05_species_score.csv"),
          row.names = FALSE)
write.csv(sp$site_detail, file.path(out_dir, "05_species_site_detail.csv"),
          row.names = FALSE)

# Per-siRNA summary file
for (s in sirna_list) {
  out_file <- file.path(out_dir, sprintf("summary_%s.csv", s$name))
  r <- ranking[ranking$siRNA_name == s$name, ]
  sp_row <- sp$species_score[sp$species_score$siRNA_name == s$name, ]
  this_sites <- all_sites[[s$name]]
  this_enriched <- enriched[enriched$siRNA_name == s$name, ]
  this_crit <- sum(this_enriched$cancer_gene_v2 | this_enriched$AE_gene |
                   this_enriched$immune_gene_immport)
  this_genes <- this_enriched$gene_name
  this_tissue <- tissue[tissue$gene_name %in% unique(this_genes), ]
  this_tissue <- head(this_tissue[order(-this_tissue$bio_attention), ], 10)

  sink(out_file)
  cat(sprintf("# Mode A Result: %s\n", s$name))
  cat(sprintf("# Guide: 5'-%s-3' (%d nt)\n", s$seq, nchar(s$seq)))
  cat(sprintf("# Sites: %d | Critical: %d | Tier: %s | Composite: %.4f\n\n",
              nrow(this_sites), this_crit,
              as.character(r$risk_tier), r$composite_score))
  cat("## Ranking Metrics\n")
  sink()
  suppressWarnings(write.table(r, out_file, append=TRUE, sep=",",
                               row.names=FALSE, quote=TRUE))
  sink(out_file, append=TRUE); cat("\n## Top Genes by Bio-Attention\n"); sink()
  if (nrow(this_tissue) > 0) {
    suppressWarnings(write.table(this_tissue, out_file, append=TRUE,
                                 sep=",", row.names=FALSE, quote=TRUE))
  }
  sink(out_file, append=TRUE); cat("\n## Species Conservation\n"); sink()
  if (nrow(sp_row) > 0) {
    suppressWarnings(write.table(sp_row, out_file, append=TRUE,
                                 sep=",", row.names=FALSE, quote=TRUE))
  }
  sink(out_file, append=TRUE)
  rec <- grep(s$name, sp$recommendation, value=TRUE)
  if (length(rec) > 0) cat(sprintf("\n# %s\n", rec))
  sink()
}

cat("\n=== ALL DONE ===\n")
cat(sprintf("Output: %s\n", normalizePath(out_dir)))
for (f in list.files(out_dir)) {
  cat(sprintf("  %s\n", f))
}
