# Dry-run validator for Mode A using a "shimmed" copy of annotate_from_sequence
# that uses our mock symbols directly (replacing the pkg::name lookups).

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
for (f in list.files("data", pattern = "\\.rda$", full.names = TRUE)) load(f, envir = .GlobalEnv)

# Mock objects (same as before)
mock_load_annotations <- function(reference.name, feature.type,
                                   canonical, min.feature.width,
                                   longest.utr, return_gene_name) {
  cat("  [mock] load_annotations: ref=", reference.name,
      " feature=", feature.type, " return_gene_name=", return_gene_name, "\n", sep="")
  stopifnot(reference.name %in% c("hg38","hg38-old","mm39","mm10","rnor6","rnor7"))
  stopifnot(return_gene_name == FALSE)
  list(
    seqs = list(ENST_A = "AAA", ENST_B = "BBB", ENST_C = "CCC"),
    gtf  = structure(list(),
                     mcols = data.frame(
                       tx_id = c("ENST_A","ENST_B","ENST_C"),
                       gene_name = c("CTNNB1","TP53","MYC"),
                       tx_biotype = "protein_coding",
                       stringsAsFactors = FALSE))
  )
}

mock_SeedMatchR <- function(seqs, sequence, seed.name, res.format) {
  stopifnot(res.format == "granges")
  cat("  [mock] SeedMatchR: seed=", seed.name, "\n", sep="")
  structure(list(seqnames = c("ENST_A","ENST_B"), start = c(120L, 250L)),
            class = "GRanges")
}
mock_seqnames <- function(g) g$seqnames

mock_DNAString <- function(s) structure(s, class = "DNAString")
mock_revComp   <- function(s) structure(paste(rev(strsplit(chartr("ACGT","TGCA", s),"")[[1]]),
                                                collapse=""), class = "DNAString")
mock_vmatchPattern <- function(pattern, subject, max.mismatch,
                                with.indels = FALSE, fixed = TRUE) {
  cat("  [mock] vmatchPattern: max.mismatch=", max.mismatch, "\n", sep="")
  if (max.mismatch == 0L) {
    list(ENST_A = list(start = c(50L)),
         ENST_B = list(start = integer()),
         ENST_C = list(start = integer()))
  } else {
    list(ENST_A = list(start = c(50L)),
         ENST_B = list(start = c(300L)),
         ENST_C = list(start = integer()))
  }
}

# Shimmed annotate_from_sequence: same body, but ::-lookups replaced with
# mock function references in this script's environment.
annotate_from_sequence_SHIM <- function(input,
                                         species = "human",
                                         feature_type = "3UTR",
                                         seed_name = "mer7m8",
                                         max_mismatch_full = 3L,
                                         strand_default = "guide",
                                         annodb = NULL) {
  if (is.null(annodb)) {
    ref <- .map_species_to_seedmatchr(species)
    annodb <- mock_load_annotations(
      reference.name = ref, feature.type = feature_type,
      canonical = FALSE, min.feature.width = 8,
      longest.utr = TRUE, return_gene_name = FALSE)
  }
  region_label <- .feature_to_region(feature_type)
  meta <- tibble::tibble(
    tx_id = attr(annodb$gtf,"mcols")$tx_id,
    gene_name = attr(annodb$gtf,"mcols")$gene_name,
    biotype = attr(annodb$gtf,"mcols")$tx_biotype
  )

  seed_hits <- {
    gr <- mock_SeedMatchR(annodb$seqs, input$sequence, seed_name, "granges")
    if (length(gr$seqnames) == 0) tibble::tibble(transcript=character(),
                                                  position=integer(),
                                                  match_type=character())
    else tibble::tibble(transcript = mock_seqnames(gr),
                        position   = as.integer(gr$start),
                        match_type = "partial_match")
  }
  full_hits <- {
    guide_dna <- mock_DNAString(gsub("U","T",input$sequence))
    target    <- mock_revComp(guide_dna)
    mi <- mock_vmatchPattern(target, annodb$seqs,
                              max.mismatch = max_mismatch_full)
    starts_list <- lapply(mi, function(z) z$start)
    lens <- vapply(starts_list, length, integer(1))
    if (sum(lens) == 0) tibble::tibble(transcript=character(),
                                         position=integer(),
                                         match_type=character())
    else {
      tx_ids <- rep(names(mi), lens)
      starts <- unlist(starts_list, use.names=FALSE)
      exact_mi <- mock_vmatchPattern(target, annodb$seqs[unique(tx_ids)],
                                      max.mismatch = 0L)
      exact_starts <- lapply(exact_mi, function(z) z$start)
      exact_set <- names(exact_starts)[vapply(exact_starts, length, integer(1)) > 0]
      mt <- ifelse(tx_ids %in% exact_set, "full_complementarity", "partial_match")
      tibble::tibble(transcript = tx_ids,
                     position   = as.integer(starts),
                     match_type = mt)
    }
  }
  merged <- .merge_hits(seed_hits, full_hits, meta = meta,
                         region_label = region_label)
  if (nrow(merged) == 0) return(tibble::tibble())
  merged$siRNA_name <- input$siRNA_name
  merged$siRNA      <- as.character(input$sequence)
  merged$strand     <- strand_default
  merged$risk_score <- compute_risk_score(merged$region, merged$match_type)
  merged[, core_columns()]
}

cat("=== Mode A dry-run with shim ===\n\n")
out_a <- annotate_from_sequence_SHIM(
  list(siRNA_name = "hsiR1", sequence = "UUAUAGAGCAAGAACACUGUUUU"),
  species = "human"
)
print(out_a)
stopifnot(all(core_columns() %in% names(out_a)))
stopifnot(nrow(out_a) >= 1)
stopifnot(all(out_a$strand == "guide"))
stopifnot(all(out_a$siRNA_name == "hsiR1"))
stopifnot(all(!is.na(out_a$gene_name)))
cat("\nMode A shim end-to-end: PASS\n\n")

cat("--- downstream enrich + rank ---\n")
enr <- enrich_annotations(out_a)
print(enr[, c("siRNA_name","gene_name","cancer_gene_v2","AE_gene","risk_score")])
rnk <- rank_sirna(enr, target_sirnas = "hsiR1")
print(rnk)

cat("\n=== ALL MODE A CHECKS PASSED ===\n")
