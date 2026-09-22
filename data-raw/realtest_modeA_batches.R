# Mode A real-data test batches (3 scenarios).
#
# Batch 1: one literature siRNA vs medium synthetic transcriptome  (accuracy)
# Batch 2: 4 siRNAs in sequence vs same transcriptome               (batch use)
# Batch 3: edge cases: no match, seed near boundary, short tx, IUPAC N
#
# Note: all batches use synthetic DNAStringSet (no hg38 download needed).
# The synthetic sequences PLANT known hits at known positions so we can
# assert exact ground truth.

# (user library resolved by R itself; no .libPaths() override)

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
for (f in list.files("data", pattern = "\\.rda$", full.names = TRUE)) {
  load(f, envir = .GlobalEnv)
}

suppressPackageStartupMessages({
  library(SeedMatchR); library(Biostrings)
  library(GenomicRanges); library(S4Vectors); library(IRanges)
  library(BiocGenerics)
})

cat("R:", R.version.string, "\n")
cat("SeedMatchR:", as.character(packageVersion("SeedMatchR")), "\n\n")

out_root <- "data-raw/modeA_batches_out"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

hr <- function(msg) cat("\n\n", strrep("=", 70), "\n", msg, "\n",
                         strrep("=", 70), "\n\n", sep = "")

# Helper: build an annodb from a named character vector of sequences
build_annodb <- function(named_seqs, gene_names, biotypes = "protein_coding") {
  seqs <- DNAStringSet(named_seqs)
  gtf  <- GRanges(seqnames = names(seqs),
                  ranges   = IRanges(start = 1, width = width(seqs)),
                  strand   = "+")
  mcols(gtf)$tx_id      <- names(seqs)
  mcols(gtf)$gene_name  <- gene_names
  mcols(gtf)$tx_biotype <- biotypes
  list(seqs = seqs, gtf = gtf)
}

# Helper: insert substring at given position (1-based start)
insert_at <- function(backbone, pos, piece) {
  pre  <- substr(backbone, 1, pos - 1)
  post <- substr(backbone, pos + nchar(piece), nchar(backbone))
  paste0(pre, piece, post)
}

set.seed(42)
random_seq <- function(len) paste(sample(c("A","C","G","T"), len, replace=TRUE), collapse="")

revc <- function(rna) as.character(
  reverseComplement(DNAString(gsub("U","T", rna))))

seed_of <- function(rna, name = "mer7m8") {
  if (name == "mer7m8") revc(substr(rna, 2, 8))
  else stop("only mer7m8 supported here")
}

summarise_result <- function(batch, result) {
  cat(sprintf("  rows   : %d\n", nrow(result)))
  cat("  match_type:\n"); print(table(result$match_type, useNA = "ifany"))
  cat("  risk_score:\n"); print(table(result$risk_score, useNA = "ifany"))
  out_csv <- file.path(out_root, sprintf("batch_%s_result.csv", batch))
  write.csv(result, out_csv, row.names = FALSE)
  cat("  saved:", normalizePath(out_csv), "\n")
}

# =======================================================================
# Batch 1: literature PMP22 siRNA vs 20-transcript background (accuracy)
# =======================================================================
hr("Batch 1: PMP22 lit siRNA vs 20-tx synthetic transcriptome")

guide <- "UCCUCUCCAUUUGAGAUAGUU"
target <- revc(guide)                          # exact DNA target
seed   <- seed_of(guide, "mer7m8")             # seed DNA target (for mer7m8)
cat("Guide:", guide, "\n")
cat("Full target (exact):", target, "\n")
cat("Seed target (mer7m8):", seed, "\n\n")

# Plant a full hit on PMP22_1 at pos 80; a partial (1 mm) on PMP22_2 at 250;
# a seed-only on SNAP25 at 150; rest are clean random.
backbones <- setNames(replicate(20, random_seq(600)), paste0("ENST_T", sprintf("%02d", 1:20)))

# Full-comp on ENST_T01 at 80 -> PMP22_1
backbones["ENST_T01"] <- insert_at(backbones["ENST_T01"], 80, target)
# 1-mismatch version on ENST_T02 at 250
target_1mm <- target; substr(target_1mm, 10, 10) <- "A"
backbones["ENST_T02"] <- insert_at(backbones["ENST_T02"], 250, target_1mm)
# 2-mismatch on ENST_T03 at 400
target_2mm <- target; substr(target_2mm, 5, 5) <- "G"; substr(target_2mm, 15, 15) <- "C"
backbones["ENST_T03"] <- insert_at(backbones["ENST_T03"], 400, target_2mm)
# Seed only on ENST_T04 at 150
backbones["ENST_T04"] <- insert_at(backbones["ENST_T04"], 150, seed)

genes <- paste0("GENE_", sprintf("%02d", 1:20))
genes[1:4] <- c("PMP22", "PMP22", "PMP22", "SNAP25")
annodb <- build_annodb(backbones, genes)

cat("Annodb:", length(annodb$seqs), "tx, widths range",
    min(width(annodb$seqs)), "-", max(width(annodb$seqs)), "\n\n")

t0 <- Sys.time()
r1 <- annotate_from_sequence(
  list(siRNA_name = "hsiR_lit_PMP22", sequence = guide),
  annodb = annodb, max_mismatch_full = 3L)
cat(sprintf("elapsed: %.3f s\n", as.numeric(Sys.time()-t0, units="secs")))
summarise_result("1_pmp22_lit", r1)

cat("\nPlanted ground truth:\n")
cat("  ENST_T01 pos 80  : full_complementarity (exact)\n")
cat("  ENST_T02 pos 250 : partial_match (1 mm)\n")
cat("  ENST_T03 pos 400 : partial_match (2 mm)\n")
cat("  ENST_T04 pos 150 : partial_match (seed only)\n")
cat("\nActual hits:\n")
print(r1[, c("transcript","position","match_type","risk_score")], n = Inf)

expected <- list(
  c("ENST_T01", 80,  "full_complementarity"),
  c("ENST_T02", 250, "partial_match"),
  c("ENST_T03", 400, "partial_match"),
  c("ENST_T04", 150, "partial_match")
)
hits_key <- paste(r1$transcript, r1$position, r1$match_type)
ok <- sum(vapply(expected, function(e) paste(e[1],e[2],e[3]) %in% hits_key, logical(1)))
cat(sprintf("\nRecall: %d / %d ground-truth hits recovered\n", ok, length(expected)))

# =======================================================================
# Batch 2: 4 siRNAs × same transcriptome (batch processing)
# =======================================================================
hr("Batch 2: 4 siRNAs sequentially against same annodb")

sirnas <- list(
  list(siRNA_name="hsiR_A", sequence="AACCGGUUAACCGGUUAACCGG"),  # 22 nt
  list(siRNA_name="hsiR_B", sequence="UCCUCUCCAUUUGAGAUAGUU"),   # 21 nt (PMP22 lit)
  list(siRNA_name="hsiR_C", sequence="UGGAUGACCCUGUCGAUAGUU"),   # 21 nt random
  list(siRNA_name="hsiR_D", sequence="GCGCGCGCGCGCGCGCGCGCG")    # low complexity
)

# Plant one exact hit for hsiR_A and hsiR_B; no planted for C/D
bb2 <- setNames(replicate(10, random_seq(500)), paste0("ENST_M", sprintf("%02d", 1:10)))
bb2["ENST_M01"] <- insert_at(bb2["ENST_M01"], 100, revc(sirnas[[1]]$sequence))
bb2["ENST_M02"] <- insert_at(bb2["ENST_M02"], 200, revc(sirnas[[2]]$sequence))
genes2 <- paste0("GENE_M", sprintf("%02d", 1:10))
annodb2 <- build_annodb(bb2, genes2)

all_res <- list()
for (s in sirnas) {
  cat(sprintf("\n--- %s (%d nt) ---\n", s$siRNA_name, nchar(s$sequence)))
  t0 <- Sys.time()
  r <- annotate_from_sequence(s, annodb = annodb2, max_mismatch_full = 3L)
  cat(sprintf("  hits: %d  elapsed: %.3f s\n",
              nrow(r), as.numeric(Sys.time()-t0, units="secs")))
  if (nrow(r) > 0) print(r[, c("transcript","position","match_type","risk_score")], n = 10)
  all_res[[s$siRNA_name]] <- r
}
combined <- do.call(rbind, all_res)
summarise_result("2_multi_sirna", combined)

# =======================================================================
# Batch 3: edge cases
# =======================================================================
hr("Batch 3: edge cases")

guide3 <- "UCCUCUCCAUUUGAGAUAGUU"  # reuse
target3 <- revc(guide3)
seed3   <- seed_of(guide3)

# A short tx (exactly 21 nt = only target)
# A tx full of N
# A tx where seed is right at position 1
# A tx with target spanning the tail (position = width - 20)
short_exact <- target3                             # width 21
all_n       <- strrep("N", 100)
seed_at_1   <- paste0(seed3, random_seq(200 - nchar(seed3)))
tail_target <- paste0(random_seq(200 - nchar(target3)), target3)
no_hit      <- paste(sample(c("A","C","G","T"), 200, replace=TRUE), collapse="")

bb3 <- c(ENST_SHORT = short_exact,
         ENST_ALLN  = all_n,
         ENST_SEED0 = seed_at_1,
         ENST_TAIL  = tail_target,
         ENST_CLEAN = no_hit)
annodb3 <- build_annodb(bb3,
                        gene_names = c("GENE_SHORT","GENE_ALLN","GENE_SEED0",
                                        "GENE_TAIL","GENE_CLEAN"))

cat("Widths:\n"); print(width(annodb3$seqs))

t0 <- Sys.time()
r3 <- annotate_from_sequence(
  list(siRNA_name = "hsiR_edge", sequence = guide3),
  annodb = annodb3, max_mismatch_full = 3L)
cat(sprintf("elapsed: %.3f s\n", as.numeric(Sys.time()-t0, units="secs")))
summarise_result("3_edge_cases", r3)

cat("\nDetailed:\n")
print(r3, n = Inf)

cat("\nExpectations:\n")
cat("  ENST_SHORT at pos 1 : full_complementarity (the only hit)\n")
cat("  ENST_ALLN           : no hit (N != ACGT under fixed=TRUE)\n")
cat("  ENST_SEED0 at pos 1 : partial_match (seed only)\n")
cat("  ENST_TAIL at pos 180: full_complementarity\n")
cat("  ENST_CLEAN          : no hit\n")

cat("\n=== Mode A batches 1-3 done ===\n")
