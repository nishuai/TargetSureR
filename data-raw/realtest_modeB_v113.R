# Mode B real test with the LATEST EnsDb (v113, Oct 2024) instead of v86.
# Downloads ~80 MB on first run, cached under R/AnnotationHub afterward.

# (user library resolved by R itself; no .libPaths() override)

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
for (f in list.files("data", pattern = "\\.rda$", full.names = TRUE)) {
  load(f, envir = .GlobalEnv)
}

suppressPackageStartupMessages({
  library(AnnotationHub)
  library(ensembldb)
})

cat("R:", R.version.string, "\n")

# ---- pull latest EnsDb ------------------------------------------------
# v113 sqlite was manually downloaded to this path (AnnotationHub download
# kept failing due to partial file issues). Load directly via EnsDb().
v113_sqlite <- file.path(Sys.getenv("LOCALAPPDATA"),
                          "R/AnnotationHub/AH119325.sqlite")
stopifnot(file.exists(v113_sqlite))
cat("Loading v113 from:", v113_sqlite, "\n")
t0 <- Sys.time()
edb <- ensembldb::EnsDb(v113_sqlite)
cat(sprintf("Loaded EnsDb in %.1f s\n", as.numeric(Sys.time()-t0, units="secs")))
cat("EnsDb organism :", ensembldb::organism(edb), "\n")
cat("EnsDb version  :", ensembldb::ensemblVersion(edb), "\n\n")

# ---- input ------------------------------------------------------------
v2_path <- file.path(Sys.getenv("TARGETSURER_INPUT", unset = "../siRNA_prediction_model"),
                     "siRNA_offtarget_annotated_final_v2.tsv")
v2 <- read.delim(v2_path, sep = "\t", stringsAsFactors = FALSE)
v2$position <- as.integer(v2$position)

minimal <- v2[, c("siRNA_name", "transcript", "position")]
cat("Input rows:", nrow(minimal), " | unique tx:", length(unique(minimal$transcript)),"\n\n")

# ---- run Mode B -------------------------------------------------------
cat("Running annotate_from_position with EnsDb v113...\n")
t0 <- Sys.time()
result <- annotate_from_position(minimal, ensdb = edb)
cat(sprintf("Done in %.1f s\n\n", as.numeric(Sys.time()-t0, units="secs")))

# ---- reconciliation ---------------------------------------------------
v2$tx_key     <- sub("\\..*$", "", v2$transcript)
result$tx_key <- sub("\\..*$", "", result$transcript)

cmp <- data.frame(
  siRNA_name  = v2$siRNA_name,
  tx_key      = v2$tx_key,
  position    = v2$position,
  gene_v2     = v2$gene_name,
  gene_ours   = result$gene_name,
  region_v2   = v2$region,
  region_ours = result$region,
  stringsAsFactors = FALSE
)

cat("=== Gene-name agreement ===\n")
gene_match <- !is.na(cmp$gene_ours) &
              toupper(cmp$gene_v2) == toupper(cmp$gene_ours)
cat(sprintf("  exact match : %d / %d (%.1f%%)\n",
            sum(gene_match), nrow(cmp), 100*mean(gene_match)))
cat(sprintf("  NA in ours  : %d (tx not found in EnsDb v113)\n",
            sum(is.na(cmp$gene_ours))))
cat(sprintf("  disagree    : %d (likely HGNC alias drift)\n",
            sum(!is.na(cmp$gene_ours) & !gene_match)))

cat("\n=== Region agreement (overall) ===\n")
region_match <- !is.na(cmp$region_ours) & cmp$region_v2 == cmp$region_ours
cat(sprintf("  exact match : %d / %d (%.1f%%)\n",
            sum(region_match), nrow(cmp), 100*mean(region_match)))

cat("\n=== Region cross-tab (v2 row x ours col) ===\n")
print(table(cmp$region_v2, cmp$region_ours, useNA="ifany"))

cat("\n=== Resolved subset (tx found in EnsDb v113) ===\n")
resolved <- cmp[!is.na(cmp$gene_ours), ]
cat(sprintf("  size              : %d / %d\n", nrow(resolved), nrow(cmp)))
cat(sprintf("  gene exact match  : %d / %d (%.1f%%)\n",
            sum(toupper(resolved$gene_v2) == toupper(resolved$gene_ours), na.rm=TRUE),
            nrow(resolved),
            100*mean(toupper(resolved$gene_v2)==toupper(resolved$gene_ours), na.rm=TRUE)))
cat(sprintf("  region exact match: %d / %d (%.1f%%)\n",
            sum(resolved$region_v2 == resolved$region_ours, na.rm=TRUE),
            nrow(resolved),
            100*mean(resolved$region_v2==resolved$region_ours, na.rm=TRUE)))

cat("\n=== Sample gene disagreements (alias drift) ===\n")
gd <- resolved[toupper(resolved$gene_v2) != toupper(resolved$gene_ours), ]
print(head(unique(gd[, c("gene_v2","gene_ours")]), 20))

# ---- save -------------------------------------------------------------
out_dir <- "data-raw/modeB_realtest_v113_out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(result, file.path(out_dir, "annotate_from_position_result.csv"),
          row.names = FALSE)
write.csv(cmp,    file.path(out_dir, "reconciliation.csv"),
          row.names = FALSE)
cat("\nSaved to:", normalizePath(out_dir), "\n")
