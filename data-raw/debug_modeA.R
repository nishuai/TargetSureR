# (user library resolved by R itself; no .libPaths() override)
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
for (f in list.files("data", pattern = "\\.rda$", full.names = TRUE)) load(f, envir = .GlobalEnv)
suppressPackageStartupMessages({
  library(SeedMatchR); library(Biostrings); library(ensembldb);
  library(EnsDb.Hsapiens.v86); library(GenomicRanges); library(BiocGenerics)
})

guide  <- "UCCUCUCCAUUUGAGAUAGUU"
target <- as.character(Biostrings::reverseComplement(
            Biostrings::DNAString(gsub("U","T", guide))))
pmp22 <- paste(sample(c("A","C","G","T"), 200, replace=TRUE), collapse="")
set.seed(1)
substring(pmp22, 80, 80 + nchar(target) - 1) <- target

seqs <- DNAStringSet(c(ENST_PMP22_fake = pmp22))

cat("Target to search:", target, "\n")
cat("Inserted at positions", 80, "-", 80+nchar(target)-1, "\n")
cat("Actual substring at 80..100:", substring(pmp22, 80, 100), "\n\n")

# Test vmatchPattern directly
cat("=== vmatchPattern (max.mismatch=0) ===\n")
mi0 <- Biostrings::vmatchPattern(DNAString(target), seqs, max.mismatch=0, fixed=TRUE)
print(mi0)

cat("\n=== vmatchPattern (max.mismatch=3) ===\n")
mi3 <- Biostrings::vmatchPattern(DNAString(target), seqs, max.mismatch=3, fixed=TRUE)
print(mi3)
cat("starts:", BiocGenerics::start(mi3)[[1]], "\n\n")

# Now call internal .scan_full
cat("=== .scan_full output ===\n")
annodb <- list(seqs = seqs,
               gtf = {g <- GRanges(); mcols(g) <- DataFrame(tx_id=character(),gene_name=character(),tx_biotype=character()); g})
hits <- .scan_full(guide, annodb, 3)
print(hits)

cat("\n=== debug .scan_full internal variables ===\n")
sequence <- guide
guide_dna <- Biostrings::DNAString(gsub("U","T", as.character(sequence)))
target2    <- Biostrings::reverseComplement(guide_dna)
mi <- Biostrings::vmatchPattern(target2, seqs, max.mismatch = 3, with.indels = FALSE, fixed = TRUE)
cat("mi lens:\n"); print(lengths(mi))
lens <- lengths(mi)
tx_ids <- rep(names(mi), lens)
starts <- unlist(BiocGenerics::start(mi), use.names = FALSE)
cat("tx_ids:", tx_ids, "\n")
cat("starts:", starts, "\n")

exact_mi <- Biostrings::vmatchPattern(target2, seqs[unique(tx_ids)],
                                        max.mismatch = 0L, fixed = TRUE)
cat("\nexact_mi names:", names(exact_mi), "\n")
exact_lens <- lengths(exact_mi)
cat("exact_lens:", exact_lens, "\n")
cat("names(exact_lens):", names(exact_lens), "\n")
cat("exact_set:", names(exact_lens)[exact_lens > 0], "\n")
cat("tx_ids %in% exact_set:", tx_ids %in% names(exact_lens)[exact_lens > 0], "\n")
