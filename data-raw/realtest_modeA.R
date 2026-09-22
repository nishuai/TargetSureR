# Mode A real test: scan a real published siRNA against the human
# transcriptome using SeedMatchR + Biostrings.
#
# Sequence: PMP22 siRNA guide from the literature (Charcot-Marie-Tooth 1A
# preclinical candidate). 21 nt, 5'->3'. We expect at least one hit on
# PMP22 transcripts (full complementarity), plus some seed-mediated
# partial matches elsewhere.

# (user library resolved by R itself; no .libPaths() override)

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
for (f in list.files("data", pattern = "\\.rda$", full.names = TRUE)) {
  load(f, envir = .GlobalEnv)
}

suppressPackageStartupMessages({
  library(SeedMatchR)
  library(Biostrings)
})

cat("R:", R.version.string, "\n")
cat("SeedMatchR:", as.character(packageVersion("SeedMatchR")), "\n\n")

# ---- 1. load human 3'UTR annotation (SeedMatchR default) --------------
cat("Loading SeedMatchR annotations for hg38 (3'UTR)...\n")
t0 <- Sys.time()
annodb <- SeedMatchR::load_annotations(
  reference.name    = "hg38",
  feature.type      = "3UTR",
  canonical         = FALSE,
  min.feature.width = 8,
  longest.utr       = TRUE,
  return_gene_name  = FALSE     # critical: keep tx_id as seqnames
)
cat(sprintf("Loaded in %.1f s\n", as.numeric(Sys.time() - t0, units = "secs")))
cat("Total 3'UTR features:", length(annodb$seqs), "\n\n")

# ---- 2. real siRNA ---------------------------------------------------
sirna <- list(
  siRNA_name = "hsiR_PMP22_lit",
  sequence   = "UCCUCUCCAUUUGAGAUAGUU"   # 21 nt
)
cat("Test siRNA:\n")
cat("  name    :", sirna$siRNA_name, "\n")
cat("  guide   :", sirna$sequence, "\n")
cat("  length  :", nchar(sirna$sequence), "\n\n")

# ---- 3. run Mode A ----------------------------------------------------
cat("Running annotate_from_sequence()...\n")
t0 <- Sys.time()
result <- annotate_from_sequence(sirna,
                                  species = "hg38",
                                  annodb  = annodb,
                                  max_mismatch_full = 3L)
cat(sprintf("Done in %.1f s\n\n", as.numeric(Sys.time() - t0, units = "secs")))

cat("Result dim:", dim(result), "\n")
cat("Columns:", paste(names(result), collapse = ", "), "\n\n")

cat("Match-type breakdown:\n")
print(table(result$match_type, useNA = "ifany"))

cat("\nRisk-score breakdown:\n")
print(table(result$risk_score, useNA = "ifany"))

cat("\nTop 15 by risk_score:\n")
ord <- order(result$risk_score, decreasing = TRUE)
print(head(result[ord, ], 15))

cat("\nAny hits on PMP22?\n")
pmp22_hits <- subset(result, !is.na(gene_name) & gene_name == "PMP22")
cat("PMP22 hit rows:", nrow(pmp22_hits), "\n")
print(pmp22_hits)

# ---- 4. save ----------------------------------------------------------
out_dir <- "data-raw/modeA_realtest_out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(result, file.path(out_dir, "annotate_from_sequence_result.csv"),
          row.names = FALSE)
cat("\nSaved to:", normalizePath(out_dir), "\n")
