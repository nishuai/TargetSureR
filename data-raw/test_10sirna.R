# Test TargetSureR with 10 real therapeutic siRNA names
# Each siRNA gets a separate output CSV file

# (user library resolved by R itself; no .libPaths() override)

setwd(Sys.getenv("TARGETSURER_HOME", unset = "."))
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
for (f in list.files("data", pattern = "\\.rda$", full.names = TRUE)) {
  load(f, envir = .GlobalEnv)
}

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
})

out_dir <- "demo_output/test_10sirna"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

v2_path <- file.path(Sys.getenv("TARGETSURER_INPUT", unset = "../siRNA_prediction_model"),
                     "siRNA_offtarget_annotated_final_v2.tsv")
v2 <- read.delim(v2_path, sep = "\t", stringsAsFactors = FALSE)
all_tx <- unique(v2$transcript)
all_mt <- v2$match_type

# 10 real therapeutic siRNAs (FDA-approved or late-stage clinical)
sirna_info <- tibble::tribble(
  ~siRNA_name,    ~target_gene, ~drug_name,
  "Patisiran",    "TTR",        "ONPATTRO",
  "Givosiran",    "ALAS1",      "GIVLAARI",
  "Lumasiran",    "HAO1",       "OXLUMO",
  "Inclisiran",   "PCSK9",      "LEQVIO",
  "Vutrisiran",   "TTR",        "AMVUTTRA",
  "Nedosiran",    "LDHA",       "RIVFLOZA",
  "Fitusiran",    "SERPINC1",   "ALHEMO",
  "Teprasiran",   "TP53",       "QPI-1002",
  "Cosdosiran",   "CASP2",      "QPI-1007",
  "Tivanisiran",  "TRPV1",      "SYL1001"
)

v113 <- file.path(Sys.getenv("LOCALAPPDATA"),
                   "R/AnnotationHub/AH119325.sqlite")

set.seed(42)
cat("=== TargetSureR Pipeline Test: 10 Therapeutic siRNAs ===\n\n")

# Generate all sites together
all_sites_raw <- lapply(seq_len(nrow(sirna_info)), function(i) {
  n_sites <- sample(80:180, 1)
  tibble(
    siRNA_name = sirna_info$siRNA_name[i],
    transcript = sample(all_tx, n_sites, replace = TRUE),
    position   = sample(100:5000, n_sites, replace = TRUE)
  )
}) |> bind_rows()
cat(sprintf("Total synthetic sites: %d\n\n", nrow(all_sites_raw)))

# Step 1: Annotate all at once
cat("Step 1: Annotating...\n")
sites <- annotate_sites(all_sites_raw, mode = "position", ensdb_path = v113,
                        rest_fallback = TRUE, with_lookup_source = TRUE)
sites$match_type <- sample(all_mt, nrow(sites), replace = TRUE)
sites$risk_score <- compute_risk_score(sites$region, sites$match_type)
cat(sprintf("  resolved: %d/%d (%.0f%%)\n",
            sum(sites$lookup_source != "unresolved"), nrow(sites),
            100*mean(sites$lookup_source != "unresolved")))

# Step 2: Enrich all
cat("Step 2: Enriching...\n")
sites_core <- sites[, c(core_columns(), "lookup_source")]
enriched <- enrich_annotations(sites_core[, core_columns()])

# Step 3: Rank all together (proper normalization across 10 siRNAs)
cat("Step 3: Ranking...\n")
ranking <- rank_sirna(enriched, target_sirnas = sirna_info$siRNA_name)
print(ranking[, c("siRNA_name","composite_score","rank","risk_tier",
                  "high_risk_critical","critical_gene_count","total_offtargets")])

# Step 4: Tissue safety (global)
cat("\nStep 4: Tissue safety...\n")
tissue <- tissue_safety(enriched, primary_tissue = "Nerve_Tibial")

# Step 5: Species match (global)
cat("Step 5: Species match...\n")
sp <- species_match(sites_core[, core_columns()],
                    candidates = c("mouse","rat","cyno",
                                   "rhesus","rabbit","dog"))
cat("Recommendations:\n")
cat(paste0("  ", sp$recommendation, collapse = "\n"), "\n")

# Output: one file per siRNA
cat("\n== Writing per-siRNA result files ==\n")
for (i in seq_len(nrow(sirna_info))) {
  sn <- sirna_info$siRNA_name[i]
  r <- ranking[ranking$siRNA_name == sn, ]
  t_genes <- tissue[tissue$gene_name %in%
    unique(enriched$gene_name[enriched$siRNA_name == sn]), ]
  t_genes <- head(t_genes[order(-t_genes$bio_attention), ], 10)
  sp_row <- sp$species_score[sp$species_score$siRNA_name == sn, ]

  out_file <- file.path(out_dir, sprintf("%02d_%s.csv", i, sn))
  sink(out_file)
  cat(sprintf("# TargetSureR Result: %s\n", sn))
  cat(sprintf("# Target: %s | Drug: %s\n",
              sirna_info$target_gene[i], sirna_info$drug_name[i]))
  cat(sprintf("# Total sites: %d | Risk tier: %s | Composite: %.4f | Rank: %d/10\n\n",
              r$total_offtargets, as.character(r$risk_tier),
              r$composite_score, r$rank))
  cat("## Ranking Metrics\n")
  sink()
  suppressWarnings(write.table(r, out_file, append=TRUE, sep=",",
                               row.names=FALSE, quote=TRUE))
  sink(out_file, append=TRUE)
  cat("\n## Top Genes by Bio-Attention (this siRNA's off-targets)\n")
  sink()
  if (nrow(t_genes) > 0) {
    suppressWarnings(write.table(
      t_genes[, c("gene_name","n_hits","primary_tpm","max_risk_score",
                  "bio_attention","is_key_gene")],
      out_file, append=TRUE, sep=",", row.names=FALSE, quote=TRUE))
  }
  sink(out_file, append=TRUE)
  cat("\n## Species Conservation Score\n")
  sink()
  if (nrow(sp_row) > 0) {
    suppressWarnings(write.table(sp_row, out_file, append=TRUE, sep=",",
                                 row.names=FALSE, quote=TRUE))
  }
  sink(out_file, append=TRUE)
  rec <- grep(sn, sp$recommendation, value=TRUE)
  if (length(rec) > 0) cat(sprintf("\n# %s\n", rec))
  sink()
  cat(sprintf("  [%d] %s -> %s\n", i, sn, basename(out_file)))
}

cat("\n=== ALL DONE ===\n")
cat(sprintf("Output: %s\n", normalizePath(out_dir)))
