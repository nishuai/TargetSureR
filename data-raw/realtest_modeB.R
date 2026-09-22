# Real test for Mode B.
# Reads the real off-target table (v2.tsv) from the original project,
# strips down to (siRNA_name, transcript, position), feeds it into
# annotate_from_position() with EnsDb.Hsapiens.v86, then reconciles the
# returned gene_name / region against the original file.

options(repos = c(CRAN = "https://cloud.r-project.org"))
# (user library resolved by R itself; no .libPaths() override)

cat("R:", R.version.string, "\n\n")

# ---- 1. Load TargetSureR (dev mode: source R/) -------------------------
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
for (f in list.files("data", pattern = "\\.rda$", full.names = TRUE)) {
  load(f, envir = .GlobalEnv)
}
cat("Loaded TargetSureR R/ files and data/ rda\n\n")

# ---- 2. Read real v2.tsv (15 columns, 956 hits, 7 siRNAs) -------------
v2_path <- file.path(Sys.getenv("TARGETSURER_INPUT", unset = "../siRNA_prediction_model"),
                     "siRNA_offtarget_annotated_final_v2.tsv")
stopifnot(file.exists(v2_path))

v2 <- readr::read_tsv(v2_path,
                       col_types = readr::cols(.default = "c",
                                                position = "i"))
cat("v2.tsv: ", nrow(v2), "rows,",
    dplyr::n_distinct(v2$siRNA_name), "siRNAs,",
    dplyr::n_distinct(v2$transcript), "transcripts\n\n")

# ---- 3. Build minimal Mode-B input ------------------------------------
minimal <- v2[, c("siRNA_name", "transcript", "position")]
cat("Minimal input (first 6 rows):\n")
print(head(minimal))
cat("\n")

# ---- 4. Load the EnsDb (offline, no AnnotationHub needed) -------------
suppressPackageStartupMessages({
  library(EnsDb.Hsapiens.v86)
})
edb <- EnsDb.Hsapiens.v86
cat("EnsDb metadata:\n")
print(ensembldb::organism(edb))
print(ensembldb::ensemblVersion(edb))
cat("\n")

# ---- 5. Run annotate_from_position ------------------------------------
cat("Running annotate_from_position()...\n")
t0 <- Sys.time()
result <- annotate_from_position(minimal, ensdb = edb)
cat(sprintf("Done in %.1f s\n\n", as.numeric(Sys.time() - t0, units = "secs")))

cat("Result dim:", dim(result), "\n")
cat("Columns:", paste(names(result), collapse = ", "), "\n\n")
cat("First 6 rows:\n")
print(head(result))

# ---- 6. Reconcile with v2.tsv -----------------------------------------
cat("\n=== Reconciliation with v2.tsv ===\n")

v2$tx_key      <- sub("\\..*$", "", v2$transcript)
result$tx_key  <- sub("\\..*$", "", result$transcript)

cmp <- dplyr::tibble(
  siRNA_name      = v2$siRNA_name,
  tx_key          = v2$tx_key,
  position        = v2$position,
  gene_v2         = v2$gene_name,
  gene_ours       = result$gene_name,
  region_v2       = v2$region,
  region_ours     = result$region,
  match_v2        = v2$match_type,
  match_ours      = result$match_type
)

cat("\nGene-name agreement:\n")
gene_match <- !is.na(cmp$gene_ours) &
              toupper(cmp$gene_v2) == toupper(cmp$gene_ours)
cat("  exact match :", sum(gene_match), "/", nrow(cmp),
    sprintf("(%.1f%%)\n", 100 * sum(gene_match) / nrow(cmp)))
cat("  NA in ours  :", sum(is.na(cmp$gene_ours)), "\n")
cat("  disagree    :",
    sum(!is.na(cmp$gene_ours) & !gene_match), "\n")

cat("\nRegion agreement:\n")
region_match <- !is.na(cmp$region_ours) & cmp$region_v2 == cmp$region_ours
cat("  exact match :", sum(region_match), "/", nrow(cmp),
    sprintf("(%.1f%%)\n", 100 * sum(region_match) / nrow(cmp)))
cat("\n  region cross-tab (v2 row x ours col):\n")
print(table(cmp$region_v2, cmp$region_ours, useNA = "ifany"))

cat("\nSample disagreements (region):\n")
disagree <- cmp[!region_match & !is.na(cmp$region_ours), ]
print(head(disagree[, c("tx_key", "position", "gene_v2",
                          "region_v2", "region_ours")], 10))

cat("\nSample disagreements (gene_name):\n")
gd <- cmp[!gene_match & !is.na(cmp$gene_ours), ]
print(head(gd[, c("tx_key", "position", "gene_v2", "gene_ours")], 10))

cat("\nWhich tx_ids had NO gene returned (likely Ensembl version drift):\n")
no_gene <- cmp[is.na(cmp$gene_ours), ]
print(head(no_gene[, c("tx_key", "gene_v2")], 10))
cat(sprintf("  total tx_ids with no gene_ours: %d (%d unique tx_keys)\n",
            nrow(no_gene), dplyr::n_distinct(no_gene$tx_key)))

# ---- 7. Save artifacts for inspection ---------------------------------
out_dir <- "data-raw/modeB_realtest_out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
readr::write_csv(result, file.path(out_dir, "annotate_from_position_result.csv"))
readr::write_csv(cmp,    file.path(out_dir, "reconciliation.csv"))

cat("\n=== DONE ===\n")
cat("Results in:", normalizePath(out_dir), "\n")
