# Mode B real-data test batches (3 scenarios).
#
# Batch 1: full v2.tsv (956 rows, 7 siRNAs)            - coverage + scale
# Batch 2: hsiR20 on PMP22/SNAP25/SF3B1 only          - precise per-site audit
# Batch 3: synthetic edge cases (5'UTR, CDS boundaries, ncRNA, NA pos)

# (user library resolved by R itself; no .libPaths() override)

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
for (f in list.files("data", pattern = "\\.rda$", full.names = TRUE)) {
  load(f, envir = .GlobalEnv)
}

suppressPackageStartupMessages({
  library(EnsDb.Hsapiens.v86)
})
edb <- EnsDb.Hsapiens.v86

out_root <- "data-raw/modeB_batches_out"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

hr <- function(msg) {
  cat("\n\n", strrep("=", 70), "\n", msg, "\n", strrep("=", 70), "\n\n", sep = "")
}

summarise_result <- function(batch, input_df, result, v2_df = NULL) {
  cat(sprintf("input : %d rows\n", nrow(input_df)))
  cat(sprintf("output: %d rows\n", nrow(result)))
  cat("region distribution:\n"); print(table(result$region, useNA = "ifany"))
  cat("match_type distribution:\n"); print(table(result$match_type, useNA = "ifany"))
  cat("risk_score distribution:\n"); print(table(result$risk_score, useNA = "ifany"))
  resolved <- sum(!is.na(result$gene_name))
  cat(sprintf("gene resolution: %d / %d (%.1f%%)\n",
              resolved, nrow(result), 100*resolved/nrow(result)))
  if (!is.null(v2_df) && all(c("gene_name","region") %in% names(v2_df))) {
    gene_match <- !is.na(result$gene_name) &
                  toupper(v2_df$gene_name) == toupper(result$gene_name)
    region_match <- !is.na(result$region) &
                    v2_df$region == result$region
    cat(sprintf("vs v2 gene   : %d / %d (%.1f%%)\n",
                sum(gene_match, na.rm=TRUE), nrow(result),
                100*mean(gene_match, na.rm=TRUE)))
    cat(sprintf("vs v2 region : %d / %d (%.1f%%)\n",
                sum(region_match, na.rm=TRUE), nrow(result),
                100*mean(region_match, na.rm=TRUE)))
  }
  out_csv <- file.path(out_root, sprintf("batch_%s_result.csv", batch))
  write.csv(result, out_csv, row.names = FALSE)
  cat("saved:", normalizePath(out_csv), "\n")
}

# ========================================================================
# Batch 1: Full v2.tsv (size + coverage)
# ========================================================================
hr("Batch 1: Full v2.tsv (956 rows, 7 siRNAs, 910 transcripts)")

v2_path <- file.path(Sys.getenv("TARGETSURER_INPUT", unset = "../siRNA_prediction_model"),
                     "siRNA_offtarget_annotated_final_v2.tsv")
v2 <- read.delim(v2_path, sep = "\t", stringsAsFactors = FALSE)
v2$position <- as.integer(v2$position)

b1_in <- v2[, c("siRNA_name", "transcript", "position")]
t0 <- Sys.time()
b1_out <- annotate_from_position(b1_in, ensdb = edb)
cat(sprintf("elapsed: %.2f s\n\n", as.numeric(Sys.time() - t0, units = "secs")))
summarise_result("1_full_v2", b1_in, b1_out, v2_df = v2)

# ========================================================================
# Batch 2: hsiR20 only (PMP22 target, 47 rows - precision audit)
# ========================================================================
hr("Batch 2: hsiR20 (47 rows, precise per-site audit)")

b2_in <- v2[v2$siRNA_name == "hsiR20", c("siRNA_name", "transcript", "position")]
b2_v2 <- v2[v2$siRNA_name == "hsiR20", ]
t0 <- Sys.time()
b2_out <- annotate_from_position(b2_in, ensdb = edb)
cat(sprintf("elapsed: %.3f s\n\n", as.numeric(Sys.time() - t0, units = "secs")))
summarise_result("2_hsiR20", b2_in, b2_out, v2_df = b2_v2)

cat("\nPer-gene breakdown:\n")
b2_cmp <- data.frame(
  gene_v2      = b2_v2$gene_name,
  gene_ours    = b2_out$gene_name,
  region_v2    = b2_v2$region,
  region_ours  = b2_out$region,
  stringsAsFactors = FALSE
)
for (g in unique(b2_v2$gene_name)) {
  sub <- b2_cmp[b2_cmp$gene_v2 == g, ]
  rm <- mean(sub$region_v2 == sub$region_ours, na.rm = TRUE)
  cat(sprintf("  %-10s  n=%3d   region match=%5.1f%%   gene resolved=%d/%d\n",
              g, nrow(sub), 100*rm,
              sum(!is.na(sub$gene_ours)), nrow(sub)))
}

# ========================================================================
# Batch 3: synthetic edge cases
# ========================================================================
hr("Batch 3: synthetic edge cases")

# Use known PMP22 transcript to craft boundary positions.
flt <- AnnotationFilter::TxIdFilter("ENST00000426385")
txl <- getFromNamespace(".transcriptLengths", "ensembldb")(
  edb, with.utr5_len = TRUE, with.cds_len = TRUE,
  filter = AnnotationFilter::TxIdFilter(
    c("ENST00000426385",  # PMP22 protein-coding
      "ENST00000380152",  # BRCA2 coding
      "ENST00000335508",  # SF3B1 coding
      "ENST00000415707",  # a random pseudogene-like
      "ENST00000361624")))# MT-CO1 mitochondrial protein
cat("Lengths reference table:\n"); print(txl[, c("tx_id","utr5_len","cds_len","tx_len")])

# Construct edge cases
b3_in <- data.frame(
  siRNA_name = c(
    "edge_pmp22_utr5_first", "edge_pmp22_utr5_last", "edge_pmp22_cds_first",
    "edge_pmp22_cds_last",   "edge_pmp22_utr3_first", "edge_pmp22_utr3_last",
    "edge_brca2_mid",        "edge_sf3b1_mid",        "edge_mtco1_mid",
    "edge_nonexistent",      "edge_past_end"
  ),
  transcript = c(
    "ENST00000426385.4", "ENST00000426385.4", "ENST00000426385.4",
    "ENST00000426385.4", "ENST00000426385.4", "ENST00000426385.4",
    "ENST00000380152.7", "ENST00000335508.11", "ENST00000361624.2",
    "ENST99999999999.1", "ENST00000426385.4"
  ),
  position = c(
    1,                                                       # first nt = 5'UTR
    max(1L, txl$utr5_len[txl$tx_id == "ENST00000426385"]),   # last nt of 5'UTR
    txl$utr5_len[txl$tx_id == "ENST00000426385"] + 1L,       # first nt of CDS
    txl$utr5_len[txl$tx_id == "ENST00000426385"] +
      txl$cds_len[txl$tx_id == "ENST00000426385"],           # last nt of CDS
    txl$utr5_len[txl$tx_id == "ENST00000426385"] +
      txl$cds_len[txl$tx_id == "ENST00000426385"] + 1L,      # first nt of 3'UTR
    txl$tx_len[txl$tx_id == "ENST00000426385"],              # last nt of tx
    500L,  # middle-ish of BRCA2
    100L,  # middle-ish of SF3B1
    500L,  # mitochondrial
    100L,  # non-existent transcript -> ncRNA fallback expected
    999999L  # way past tx end
  ),
  stringsAsFactors = FALSE
)
print(b3_in)

t0 <- Sys.time()
b3_out <- annotate_from_position(b3_in, ensdb = edb)
cat(sprintf("\nelapsed: %.3f s\n\n", as.numeric(Sys.time() - t0, units = "secs")))
cat("Detailed result:\n")
print(b3_out, n = Inf)

summarise_result("3_edge_cases", b3_in, b3_out)

cat("\n\nExpected vs observed for edge cases:\n")
expected <- c("5'UTR","5'UTR","CDS","CDS","3'UTR","3'UTR",
              "CDS","CDS","CDS","ncRNA","3'UTR")
observed <- b3_out$region
cat(sprintf("  %-30s  %-8s  %-8s  %s\n", "case", "exp", "obs", "ok?"))
for (i in seq_along(expected)) {
  ok <- identical(expected[i], observed[i])
  cat(sprintf("  %-30s  %-8s  %-8s  %s\n",
              b3_in$siRNA_name[i], expected[i], observed[i],
              if (ok) "OK" else "MISMATCH"))
}

cat("\n=== Batches 1-3 done ===\n")
