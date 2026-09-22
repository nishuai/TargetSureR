# Mode A test (offline, no AnnotationHub):
# Build SeedMatchR-compatible annodb from local BSgenome + EnsDb v113.
# Test 2 siRNAs:
#   1. Patisiran (TTR, FDA-approved siRNA drug)
#   2. user-provided  UCCUGUUGCUGAGUAUCAUTT

# (user library resolved by R itself; no .libPaths() override)

setwd(Sys.getenv("TARGETSURER_HOME", unset = "."))
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
for (f in list.files("data", pattern = "\\.rda$", full.names = TRUE)) {
  load(f, envir = .GlobalEnv)
}

suppressPackageStartupMessages({
  library(dplyr); library(tibble)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(ensembldb)
  library(GenomicRanges); library(GenomicFeatures); library(Biostrings)
  library(BiocGenerics); library(S4Vectors); library(IRanges)
  library(SeedMatchR)
})

bsg  <- BSgenome.Hsapiens.UCSC.hg38
v113 <- file.path(Sys.getenv("LOCALAPPDATA"),
                   "R/AnnotationHub/AH119325.sqlite")
edb  <- ensembldb::EnsDb(v113)

out_dir <- "demo_output/modeA_offline"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. extract 3'UTR sequences from EnsDb + BSgenome ----------------
annodb_cache <- "data-raw/hg38_3UTR_annodb.rds"
if (file.exists(annodb_cache)) {
  cat(sprintf("=== Step 0: load cached annodb (%s) ===\n", annodb_cache))
  annodb <- readRDS(annodb_cache)
  cat(sprintf("annodb loaded: %d transcripts\n\n", length(annodb$seqs)))
} else {
  cat("=== Step 0: build offline 3UTR annodb ===\n")
  t0 <- Sys.time()

three_utr <- threeUTRsByTranscript(edb)
cat(sprintf("3UTR ranges loaded: %d transcripts (%.1f s)\n",
            length(three_utr),
            as.numeric(Sys.time()-t0, units="secs")))

# Filter to standard chromosomes and protein-coding genes
tx_meta <- transcripts(edb,
  columns = c("tx_id","tx_biotype","gene_name","seq_name"),
  return.type = "data.frame")
tx_meta <- subset(tx_meta,
  tx_biotype == "protein_coding" &
  seq_name %in% c(as.character(1:22),"X","Y","MT"))
cat(sprintf("Protein-coding tx: %d\n", nrow(tx_meta)))

three_utr <- three_utr[names(three_utr) %in% tx_meta$tx_id]
cat(sprintf("3UTR after filter: %d transcripts\n", length(three_utr)))

# Need UCSC-style chr names for BSgenome
seqlevelsStyle(three_utr) <- "UCSC"

# Extract 3UTR sequences (concat exons)
cat("Extracting 3'UTR sequences from BSgenome...\n")
t1 <- Sys.time()
utr_seqs <- extractTranscriptSeqs(bsg, three_utr)
cat(sprintf("Extracted %d 3UTR seqs in %.1f s\n",
            length(utr_seqs),
            as.numeric(Sys.time()-t1, units="secs")))

# Take longest 3UTR per gene to mirror SeedMatchR's longest.utr=TRUE
tx2gene <- setNames(tx_meta$gene_name, tx_meta$tx_id)
utr_widths <- BiocGenerics::width(utr_seqs)
df <- data.frame(tx = names(utr_seqs),
                  gene = tx2gene[names(utr_seqs)],
                  w = utr_widths,
                  stringsAsFactors = FALSE)
df <- df[order(df$gene, -df$w), ]
df_keep <- df[!duplicated(df$gene) & !is.na(df$gene), ]
cat(sprintf("After longest-per-gene filter: %d (1 transcript per gene)\n",
            nrow(df_keep)))

utr_seqs <- utr_seqs[df_keep$tx]
utr_seqs <- utr_seqs[BiocGenerics::width(utr_seqs) >= 8]
cat(sprintf("Final transcripts: %d\n", length(utr_seqs)))

# Build a gtf-like GRanges with mcols
gtf_meta <- subset(tx_meta, tx_id %in% names(utr_seqs))
gtf <- GRanges(
  seqnames = paste0("chr", gtf_meta$seq_name),
  ranges   = IRanges(start = 1, width = nchar(utr_seqs[gtf_meta$tx_id])),
  strand   = "+"
)
S4Vectors::mcols(gtf)$tx_id      <- gtf_meta$tx_id
S4Vectors::mcols(gtf)$gene_name  <- gtf_meta$gene_name
S4Vectors::mcols(gtf)$tx_biotype <- gtf_meta$tx_biotype

annodb <- list(seqs = utr_seqs, gtf = gtf, txdb = edb)
cat(sprintf("annodb ready: %d transcripts\n\n", length(annodb$seqs)))

# Save annodb for reuse
saveRDS(annodb, "data-raw/hg38_3UTR_annodb.rds", compress = "xz")
cat("Saved annodb to data-raw/hg38_3UTR_annodb.rds\n\n")
}  # end else-build

# ---- 2. run Mode A on 2 siRNAs --------------------------------------
# Note: dTdT 3' overhang is removed; SeedMatchR scans the RNA core only.
sirna_list <- list(
  list(name = "Patisiran",  seq = "AUGGAAUACUCUUGGUUAC"),
  list(name = "siRNA_user", seq = "UCCUGUUGCUGAGUAUCAU")
)

all_sites <- list()
for (s in sirna_list) {
  cat(sprintf("=== %s ===\n", s$name))
  cat(sprintf("Guide (5'->3'): %s\n", s$seq))
  t0 <- Sys.time()
  sites <- annotate_sites(
    input  = list(siRNA_name = s$name, sequence = s$seq),
    mode   = "sequence", species = "human",
    feature_type = "3UTR", seed_name = "mer7m8",
    max_mismatch_full = 3L, annodb = annodb,
    save_raw = file.path(out_dir, "raw")
  )
  cat(sprintf("  hits: %d sites, %d transcripts, %d genes (%.1f s)\n",
              nrow(sites), n_distinct(sites$transcript),
              n_distinct(sites$gene_name),
              as.numeric(Sys.time()-t0, units="secs")))
  if (nrow(sites) > 0) {
    cat("  match_type:\n"); print(table(sites$match_type))
    cat("  n_mismatch distribution:\n"); print(table(sites$n_mismatch, useNA="ifany"))
    cat("  Top 3 hits (with target_site):\n")
    show <- head(sites[order(-sites$risk_score, sites$n_mismatch), ], 3)
    for (k in seq_len(nrow(show))) {
      cat(sprintf("    [%s @ %s:%d] %s mm=%d\n",
                  show$gene_name[k], show$transcript[k],
                  show$position[k], show$match_type[k],
                  show$n_mismatch[k]))
      cat("    alignment:\n")
      cat("    ", gsub("\n", "\n    ", show$alignment[k]), "\n", sep="")
    }
  }
  all_sites[[s$name]] <- sites
  cat("\n")
}

sites_all <- bind_rows(all_sites)
write.csv(sites_all, file.path(out_dir, "01_sites_core.csv"), row.names=FALSE)

if (nrow(sites_all) == 0) {
  cat("No hits, stopping.\n"); quit()
}

# ---- 3. enrich + rank + tissue + species ----------------------------
cat("=== enrich_annotations() ===\n")
enriched <- enrich_annotations(sites_all)
n_crit <- sum(enriched$cancer_gene_v2 | enriched$AE_gene |
              enriched$immune_gene_immport)
cat(sprintf("enriched: %d rows; critical: %d (%.1f%%)\n",
            nrow(enriched), n_crit, 100*n_crit/nrow(enriched)))
write.csv(enriched, file.path(out_dir, "02_sites_enriched.csv"),
          row.names=FALSE)

cat("\n=== rank_sirna() ===\n")
ranking <- rank_sirna(enriched,
                      target_sirnas = sapply(sirna_list, `[[`, "name"))
print(ranking[, c("siRNA_name","composite_score","rank","risk_tier",
                  "high_risk_critical","critical_gene_count",
                  "full_comp_count","cds_count","total_offtargets")])
write.csv(ranking, file.path(out_dir, "03_ranking.csv"), row.names=FALSE)

cat("\n=== tissue_safety() ===\n")
tissue <- tissue_safety(enriched, primary_tissue = "Nerve_Tibial")
cat("Top 10:\n")
print(head(tissue[, c("gene_name","n_hits","primary_tpm","max_risk_score",
                      "bio_attention","is_key_gene")], 10))
write.csv(tissue, file.path(out_dir, "04_tissue_safety.csv"),
          row.names=FALSE)

cat("\n=== species_match() ===\n")
sp <- species_match(sites_all,
                    candidates = c("mouse","rat","cyno","rhesus","rabbit","dog"))
cat("species_score:\n"); print(sp$species_score)
cat("\nrecommendation:\n")
cat(paste0("  ", sp$recommendation, collapse="\n"), "\n")
write.csv(sp$species_score, file.path(out_dir, "05_species_score.csv"),
          row.names=FALSE)

cat("\n=== generate_figures() ===\n")
fig_dir <- file.path(out_dir, "figures")
fig_paths <- generate_figures(enriched, ranking, tissue,
                               species_result = sp,
                               output_dir = fig_dir)
cat(sprintf("Wrote %d figure files to %s\n",
            length(fig_paths), normalizePath(fig_dir)))

# Per-siRNA summary
for (s in sirna_list) {
  out_file <- file.path(out_dir, sprintf("summary_%s.csv", s$name))
  r <- ranking[ranking$siRNA_name == s$name, ]
  this_e <- enriched[enriched$siRNA_name == s$name, ]
  this_crit <- sum(this_e$cancer_gene_v2 | this_e$AE_gene |
                   this_e$immune_gene_immport)
  this_tissue <- tissue[tissue$gene_name %in% unique(this_e$gene_name), ]
  this_tissue <- head(this_tissue[order(-this_tissue$bio_attention), ], 10)
  sp_row <- sp$species_score[sp$species_score$siRNA_name == s$name, ]

  sink(out_file)
  cat(sprintf("# Mode A: %s\n", s$name))
  cat(sprintf("# Guide: 5'-%s-3' (%d nt)\n", s$seq, nchar(s$seq)))
  cat(sprintf("# Sites: %d | Critical: %d | Tier: %s | Composite: %.4f\n\n",
              nrow(this_e), this_crit,
              as.character(r$risk_tier), r$composite_score))
  cat("## Ranking\n"); sink()
  suppressWarnings(write.table(r, out_file, append=TRUE, sep=",",
                                row.names=FALSE, quote=TRUE))
  sink(out_file, append=TRUE); cat("\n## Top genes\n"); sink()
  if (nrow(this_tissue) > 0) {
    suppressWarnings(write.table(this_tissue, out_file, append=TRUE,
                                 sep=",", row.names=FALSE, quote=TRUE))
  }
  sink(out_file, append=TRUE); cat("\n## Species\n"); sink()
  if (nrow(sp_row) > 0) {
    suppressWarnings(write.table(sp_row, out_file, append=TRUE,
                                 sep=",", row.names=FALSE, quote=TRUE))
  }
  sink(out_file, append=TRUE)
  rec <- grep(s$name, sp$recommendation, value=TRUE)
  if (length(rec) > 0) cat(sprintf("\n# %s\n", rec))
  sink()
}

cat("\n=== DONE ===\n")
cat(sprintf("Output: %s\n", normalizePath(out_dir)))
for (f in list.files(out_dir)) {
  fp <- file.path(out_dir, f)
  cat(sprintf("  %-30s  %.1f KB\n", f, file.info(fp)$size/1024))
}
