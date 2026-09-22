# Mode A real test WITHOUT AnnotationHub/2bit download.
#
# Strategy: build annodb$seqs locally from EnsDb.Hsapiens.v86 +
# BSgenome.Hsapiens.UCSC.hg38 (or from any transcript FASTA if present).
# This bypasses SeedMatchR's load_annotations() network calls entirely.
#
# Required packages (install if missing):
#   - EnsDb.Hsapiens.v86             (already installed)
#   - BSgenome.Hsapiens.UCSC.hg38    (needed only if extracting UTR seqs)
#
# Fallback: if BSgenome is not installed, read a user-supplied FASTA of
# 3'UTR sequences (one entry per transcript, seqnames = tx_id).

# (user library resolved by R itself; no .libPaths() override)

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
for (f in list.files("data", pattern = "\\.rda$", full.names = TRUE)) {
  load(f, envir = .GlobalEnv)
}

suppressPackageStartupMessages({
  library(SeedMatchR)
  library(Biostrings)
  library(ensembldb)
  library(EnsDb.Hsapiens.v86)
  library(GenomicRanges)
  library(BiocGenerics)
})

cat("R:", R.version.string, "\n")
cat("SeedMatchR:", as.character(packageVersion("SeedMatchR")), "\n\n")

edb <- EnsDb.Hsapiens.v86

# ---- 1. try local BSgenome, else skip -------------------------------
use_bsg <- requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
if (!use_bsg) {
  cat("BSgenome.Hsapiens.UCSC.hg38 not installed.\n")
  cat("For now we will test with a MINIMAL custom DNAStringSet built\n")
  cat("from just the PMP22-family 3'UTR sequences (synthetic, for structural test).\n\n")
}

# ---- 2. build minimal annodb by hand for structural test ------------
# We will just create a tiny DNAStringSet that guarantees at least ONE
# full-complementarity hit to our guide, so we verify merge/region/scoring
# path without needing real UTR sequences.
guide  <- "UCCUCUCCAUUUGAGAUAGUU"
target <- as.character(Biostrings::reverseComplement(
            Biostrings::DNAString(gsub("U", "T", guide))))
cat("Guide (U->T) revcomp target:", target, "\n\n")

# Synthesize 3 fake 3'UTRs:
#   - ENST_PMP22_fake: contains exact target at position 80
#   - ENST_NOISE1:     contains seed-only match somewhere
#   - ENST_NOISE2:     no match
noise1_backbone <- paste(sample(c("A","C","G","T"), 300, replace=TRUE), collapse="")
seed_target <- as.character(Biostrings::reverseComplement(
                 Biostrings::DNAString(gsub("U","T", substr(guide, 2, 8))))) # mer7m8
# Insert seed at pos 150 of noise1
substring(noise1_backbone, 150, 150 + nchar(seed_target) - 1) <- seed_target

noise2 <- paste(sample(c("A","C","G","T"), 300, replace=TRUE), collapse="")
pmp22  <- paste(sample(c("A","C","G","T"), 200, replace=TRUE), collapse="")
substring(pmp22, 80, 80 + nchar(target) - 1) <- target  # plant exact hit at 80

seqs <- Biostrings::DNAStringSet(c(
  ENST_PMP22_fake = pmp22,
  ENST_NOISE1     = noise1_backbone,
  ENST_NOISE2     = noise2
))

# build gtf-like GRanges with mcols for tx_meta
gtf <- GenomicRanges::GRanges(
  seqnames = names(seqs),
  ranges   = IRanges::IRanges(start = 1, width = BiocGenerics::width(seqs)),
  strand   = "+"
)
S4Vectors::mcols(gtf)$tx_id      <- names(seqs)
S4Vectors::mcols(gtf)$gene_name  <- c("PMP22", "NOISE1", "NOISE2")
S4Vectors::mcols(gtf)$tx_biotype <- "protein_coding"

annodb <- list(
  seqs = seqs,
  gtf  = gtf,
  txdb = edb
)

cat("Constructed synthetic annodb:\n")
cat("  seqs : ", length(annodb$seqs), "transcripts\n")
cat("  widths:", BiocGenerics::width(annodb$seqs), "\n\n")

# ---- 3. run Mode A ---------------------------------------------------
sirna <- list(siRNA_name = "hsiR_PMP22_lit", sequence = guide)

cat("Running annotate_from_sequence()...\n")
t0 <- Sys.time()
result <- annotate_from_sequence(sirna, species = "human",
                                  annodb = annodb,
                                  max_mismatch_full = 3L)
cat(sprintf("Done in %.3f s\n\n", as.numeric(Sys.time() - t0, units = "secs")))

cat("Result:\n")
print(result)

cat("\n--- Expectations ---\n")
cat("Exact target planted at ENST_PMP22_fake position 80 (full_complementarity expected)\n")
cat("Seed planted at ENST_NOISE1 position 150 (partial_match expected)\n")
cat("Nothing planted on ENST_NOISE2\n")

# ---- 4. sanity checks ------------------------------------------------
pm <- subset(result, transcript == "ENST_PMP22_fake")
stopifnot("ENST_PMP22_fake in result" = nrow(pm) >= 1)
stopifnot("position on PMP22 = 80" = any(pm$position == 80))
stopifnot("match_type on PMP22 = full_complementarity" =
            any(pm$match_type == "full_complementarity"))

n1 <- subset(result, transcript == "ENST_NOISE1")
if (nrow(n1) > 0) {
  stopifnot("NOISE1 hits are partial" = all(n1$match_type == "partial_match"))
}

cat("\nALL STRUCTURAL CHECKS PASSED\n")

# ---- 5. save ---------------------------------------------------------
out_dir <- "data-raw/modeA_synthetic_out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(result, file.path(out_dir, "annotate_from_sequence_result.csv"),
          row.names = FALSE)
cat("Saved:", normalizePath(out_dir), "\n")
