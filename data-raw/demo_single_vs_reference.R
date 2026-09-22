# Demo: single siRNA with built-in reference cohort
# Test siRNA: Patisiran (TTR-targeting, FDA-approved)

# (user library resolved by R itself; no .libPaths() override)

setwd(Sys.getenv("TARGETSURER_HOME", unset = "."))
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
for (f in list.files("data", pattern = "\\.rda$", full.names = TRUE)) {
  load(f, envir = .GlobalEnv)
}

suppressPackageStartupMessages({
  library(dplyr); library(tibble)
})

annodb <- readRDS("data-raw/hg38_3UTR_annodb.rds")

out_dir <- "demo_output/single_siRNA_vs_reference"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Run Mode A on a single siRNA ---------------------------------
sites <- annotate_sites(
  input  = list(siRNA_name = "Patisiran",
                sequence   = "AUGGAAUACUCUUGGUUAC"),
  mode   = "sequence", annodb = annodb,
  save_raw = file.path(out_dir, "raw"))
cat(sprintf("Mode A: %d off-target sites\n", nrow(sites)))

# ---- 2. Enrich + downstream ------------------------------------------
enriched <- enrich_annotations(sites)
cat(sprintf("Enriched: %d sites x %d cols\n", nrow(enriched), ncol(enriched)))

# ---- 3. Rank using BUILT-IN REFERENCE (default now) ------------------
ranking <- rank_sirna(enriched, target_sirnas = "Patisiran")
cat("\nRanking with reference='builtin' (default):\n")
print(ranking[, c("siRNA_name","composite_score","rank","risk_tier",
                  "high_risk_critical","critical_gene_count",
                  "full_comp_count","total_offtargets")])

# Compare to legacy self mode
ranking_self <- rank_sirna(enriched, target_sirnas = "Patisiran",
                            reference = "self")
cat("\nLegacy reference='self' (single siRNA -> all 0):\n")
print(ranking_self[, c("siRNA_name","composite_score","rank","risk_tier")])

# Show where Patisiran sits in the reference cohort
ref <- reference_set
cat(sprintf("\nReference cohort (n=%d):\n", ref$n_siRNAs))
cat(sprintf("  composite_score range: [%.3f, %.3f]\n",
            min(ref$metrics$composite_score),
            max(ref$metrics$composite_score)))
pos <- mean(ref$metrics$composite_score < ranking$composite_score)
cat(sprintf("  Patisiran (score=%.3f) percentile: %.0f%% (lower = safer)\n",
            ranking$composite_score, 100 * pos))

# ---- 4. Tissue safety + species --------------------------------------
tissue <- tissue_safety(enriched, primary_tissue = "Nerve_Tibial")
species <- species_match(sites,
              candidates = c("mouse","rat","cyno","rhesus","rabbit","dog"))

# ---- 5. Report bundle (CSVs + figures) -------------------------------
paths <- generate_report(enriched,
                          output_dir     = out_dir,
                          target_sirnas  = "Patisiran",
                          primary_tissue = "Nerve_Tibial",
                          species_result = species,
                          with_figures   = TRUE)

# ---- 6. Generate a ranking_bar with reference overlay ----------------
p_ref <- plot_ranking_bar(ranking, reference = "builtin")
ggplot2::ggsave(file.path(out_dir, "figures", "01b_ranking_vs_reference.pdf"),
                 p_ref, width = 9, height = 5)
ggplot2::ggsave(file.path(out_dir, "figures", "01b_ranking_vs_reference.png"),
                 p_ref, width = 9, height = 5, dpi = 150)

cat("\n=== DONE ===\n")
cat(sprintf("Output: %s\n", normalizePath(out_dir)))
for (f in list.files(out_dir, recursive = TRUE)) {
  cat(sprintf("  %s\n", f))
}
