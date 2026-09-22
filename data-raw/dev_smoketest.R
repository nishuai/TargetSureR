options(warn = 1)

cat("Loading dependencies...\n")
for (p in c("dplyr","tibble","tidyr","readr","purrr","rlang")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("  installing:", p, "\n")
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

cat("\nSourcing package R files...\n")
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  source(f, local = .GlobalEnv)
}

cat("\nLoading internal data into globalenv...\n")
for (f in list.files("data", pattern = "\\.rda$", full.names = TRUE)) {
  load(f, envir = .GlobalEnv)
}

cat("\n--- Test 1: compute_risk_score ---\n")
print(compute_risk_score(
  region     = c("CDS","CDS","3'UTR","3'UTR","5'UTR","ncRNA"),
  match_type = c("full_complementarity","partial_match",
                 "full_complementarity","partial_match",
                 "unknown","unknown")
))

cat("\n--- Test 2: column specs ---\n")
cat("core_columns:    ", length(core_columns()), "\n")
cat("enriched_columns:", length(enriched_columns()), "\n")

cat("\n--- Test 3: weights ---\n")
w <- default_weights()
print(w)
cat("sum:", sum(w), "\n")

cat("\n--- Test 4: full Mode-B pipeline (with pre-filled fields) ---\n")
sites <- tibble::tibble(
  siRNA_name = c("hsiR22","hsiR22","hsiR21","hsiR21","hsiR30"),
  siRNA      = NA_character_,
  strand     = "guide",
  transcript = c("ENST1","ENST2","ENST3","ENST4","ENST5"),
  gene_name  = c("CTNNB1","TP53","MYC","IL6","EGFR"),
  biotype    = "protein_coding",
  position   = c(379L,4687L,100L,250L,900L),
  region     = c("CDS","3'UTR","CDS","3'UTR","5'UTR"),
  match_type = c("full_complementarity","partial_match",
                 "full_complementarity","partial_match","unknown"),
  risk_score = NA_integer_
)
sites$risk_score <- compute_risk_score(sites$region, sites$match_type)

cat("\nsites:\n"); print(sites)

cat("\n--- enrich_annotations ---\n")
enriched <- enrich_annotations(sites)
cat("dim:", dim(enriched), "\n")
cat("AE/cancer/immune flags:\n")
print(enriched[, c("gene_name","AE_gene","cancer_gene_v2","immune_gene_immport")])

cat("\n--- rank_sirna ---\n")
ranking <- rank_sirna(enriched, target_sirnas = "hsiR22")
print(ranking[, c("siRNA_name","composite_score","rank","risk_tier","is_target")])

cat("\n--- tissue_safety ---\n")
ts <- tissue_safety(enriched, primary_tissue = "Nerve_Tibial")
print(ts[, c("gene_name","primary_tpm","max_risk_score","bio_attention")])

cat("\n--- species_match ---\n")
res <- species_match(sites, candidates = c("mouse","rat","cyno"))
print(res$species_score)
cat("\nrecommendation:\n")
cat(paste0("  ", res$recommendation, collapse = "\n"), "\n")

cat("\n--- generate_report ---\n")
od <- tempfile("tsr_report_")
paths <- generate_report(enriched, output_dir = od,
                          target_sirnas = "hsiR22",
                          primary_tissue = "Nerve_Tibial",
                          species_result = res)
cat("written files:\n")
for (p in unlist(paths)) cat("  ", p, "\n")

cat("\nALL OK\n")
