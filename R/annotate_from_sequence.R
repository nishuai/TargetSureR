#' Annotate off-target sites from siRNA sequence (Mode A)
#'
#' Scans a reference transcriptome for seed-mediated (via \pkg{SeedMatchR})
#' and near-full-complementarity off-target sites of a supplied siRNA guide
#' sequence. Returns the 10-column standard off-target table.
#'
#' @details
#' Two passes are performed and merged:
#' \itemize{
#'   \item \pkg{SeedMatchR} canonical seed search
#'     (default \code{"mer7m8"}). Hits are labelled
#'     \code{match_type = "partial_match"}.
#'   \item \code{Biostrings::vmatchPattern} on the reverse-complement of the
#'     full guide with \code{max.mismatch = max_mismatch_full}. Hits with
#'     0 mismatches are \code{"full_complementarity"}; hits with 1-3
#'     mismatches are \code{"partial_match"}.
#' }
#' \pkg{SeedMatchR}'s \code{load_annotations()} accepts only the keys
#' \code{c("hg38", "hg38-old", "mm39", "mm10", "rnor6", "rnor7")}. Common
#' species names (\code{"human"}, \code{"mouse"}, \code{"rat"}) are mapped
#' to a sensible default.
#'
#' By default \code{load_annotations()} extracts only one feature class at a
#' time (3' UTR by default). The \code{region} column reflects whichever
#' feature class was loaded into \code{annodb}.
#'
#' @param input A list with \code{siRNA_name} and \code{sequence}. The
#'   sequence MUST be the \strong{guide (antisense) strand}, oriented
#'   5'->3', as RNA or DNA character. The scan searches the transcriptome
#'   for \code{reverseComplement(sequence)}, so supplying a sense strand
#'   raises no error and instead returns a plausible-looking table of
#'   seed-only hits with no on-target match. To self-check: a true guide
#'   yields a \code{full_complementarity} hit with \code{n_mismatch == 0}
#'   on its own target gene. Use \code{\link{normalize_guide}} to convert a
#'   sense strand, strip a \code{dTdT}/\code{UU} overhang, or match the
#'   19 nt convention of the built-in \code{\link{reference_set}}.
#' @param species Organism key. One of \code{"human"}, \code{"mouse"},
#'   \code{"rat"}, or a key understood by
#'   \code{SeedMatchR::load_annotations()}.
#' @param feature_type Feature class to load
#'   (\code{"3UTR"}, \code{"5UTR"}, \code{"exons"}, \code{"cds"}).
#'   Default \code{"3UTR"}, matching SeedMatchR's primary use case.
#' @param seed_name Seed definition passed to \pkg{SeedMatchR}
#'   (default \code{"mer7m8"}).
#' @param max_mismatch_full Maximum mismatches allowed by the
#'   full-complementarity scan. Default \code{3}.
#' @param strand_default Default strand to tag. Usually \code{"guide"}.
#' @param annodb Optional pre-loaded annotation database (the list returned
#'   by \code{SeedMatchR::load_annotations}). Providing it speeds up
#'   repeated calls and avoids re-downloads.
#' @param save_raw Optional. Path to a directory; if non-NULL, the raw
#'   per-hit output from \pkg{SeedMatchR} (seed scan) and
#'   \pkg{Biostrings::vmatchPattern} (full scan) is written to a
#'   tab-separated file \code{raw_<siRNA_name>.tsv} in that directory,
#'   preserving every raw-engine column for audit.
#' @param ... Reserved for future arguments.
#'
#' @return A tibble with \code{\link{core_columns}} plus three columns
#'   describing the matched site: \code{target_site} (mRNA subsequence
#'   in RNA notation, 5'->3', same width as the guide), \code{n_mismatch}
#'   (integer count of mismatches between guide-revcomp and target_site),
#'   and \code{alignment} (3-line text alignment of guide / target /
#'   match line).
#'
#' @examples
#' \dontrun{
#' annotate_from_sequence(
#'   list(siRNA_name = "hsiR22",
#'        sequence   = "UUAUAGAGCAAGAACACUGUUUU"),
#'   species = "human")
#' }
#' @export
annotate_from_sequence <- function(input,
                                   species = "human",
                                   feature_type = c("3UTR", "5UTR",
                                                     "exons", "cds"),
                                   seed_name = "mer7m8",
                                   max_mismatch_full = 3L,
                                   strand_default = "guide",
                                   annodb = NULL,
                                   save_raw = NULL,
                                   ...) {
  require_pkg("SeedMatchR",    "sequence-mode annotation")
  require_pkg("Biostrings",    "sequence-mode annotation")
  require_pkg("GenomicRanges", "sequence-mode annotation")
  require_pkg("BiocGenerics",  "sequence-mode annotation")

  stopifnot(is.list(input),
            !is.null(input$siRNA_name),
            !is.null(input$sequence))
  feature_type <- match.arg(feature_type)

  if (is.null(annodb)) {
    ref <- .map_species_to_seedmatchr(species)
    annodb <- SeedMatchR::load_annotations(
      reference.name    = ref,
      feature.type      = feature_type,
      canonical         = FALSE,
      min.feature.width = 8,
      longest.utr       = TRUE,
      return_gene_name  = FALSE   # critical: keep tx_id as seqnames
    )
  }

  region_label <- .feature_to_region(feature_type)
  meta <- .build_tx_meta(annodb)

  guide_len <- nchar(as.character(input$sequence))
  guide_rna <- chartr("T", "U", toupper(as.character(input$sequence)))

  seed_hits <- .scan_seed(input$sequence, annodb, seed_name = seed_name,
                          guide_rna = guide_rna, guide_len = guide_len)
  full_hits <- .scan_full(input$sequence, annodb, max_mismatch_full,
                          guide_rna = guide_rna, guide_len = guide_len)

  if (!is.null(save_raw)) {
    .save_raw_scan(seed_hits, full_hits, save_raw, input$siRNA_name,
                   input$sequence, seed_name, max_mismatch_full)
  }

  merged <- .merge_hits(seed_hits, full_hits, meta = meta,
                        region_label = region_label)

  if (nrow(merged) == 0) {
    out <- tibble::tibble(
      siRNA_name = character(), siRNA = character(),
      strand = character(), transcript = character(),
      gene_name = character(), biotype = character(),
      position = integer(), region = character(),
      match_type = character(), risk_score = integer(),
      target_site = character(), n_mismatch = integer(),
      alignment = character()
    )
    return(out)
  }

  merged$siRNA_name <- input$siRNA_name
  merged$siRNA      <- guide_rna
  merged$strand     <- strand_default
  merged$risk_score <- compute_risk_score(merged$region, merged$match_type)
  merged$alignment  <- .format_alignment(rep(guide_rna, nrow(merged)),
                                          merged$target_site)

  merged[, c(core_columns(), "target_site", "n_mismatch", "alignment"),
         drop = FALSE]
}

.map_species_to_seedmatchr <- function(species) {
  key <- tolower(species)
  m <- c(
    human  = "hg38",
    hg38   = "hg38",
    mouse  = "mm39",
    mm39   = "mm39",
    mm10   = "mm10",
    rat    = "rnor7",
    rnor7  = "rnor7",
    rnor6  = "rnor6"
  )
  if (key %in% names(m)) return(unname(m[key]))
  species
}

.feature_to_region <- function(feature_type) {
  switch(feature_type,
         `3UTR`  = "3'UTR",
         `5UTR`  = "5'UTR",
         exons   = "exon",
         cds     = "CDS",
         "3'UTR")
}

#' Build a transcript-to-gene mapping from annodb$gtf
#' @keywords internal
#' @noRd
.build_tx_meta <- function(annodb) {
  gtf <- annodb$gtf
  if (is.null(gtf)) {
    return(tibble::tibble(tx_id = character(),
                          gene_name = character(),
                          biotype = character()))
  }
  m <- S4Vectors::mcols(gtf)
  tx_id <- if ("tx_id"      %in% names(m)) as.character(m$tx_id)      else character(length(gtf))
  gene  <- if ("gene_name"  %in% names(m)) as.character(m$gene_name)
           else if ("gene_id" %in% names(m)) as.character(m$gene_id)
           else rep(NA_character_, length(gtf))
  bt    <- if ("tx_biotype" %in% names(m)) as.character(m$tx_biotype)
           else if ("biotype" %in% names(m)) as.character(m$biotype)
           else rep(NA_character_, length(gtf))
  tibble::tibble(tx_id = tx_id, gene_name = gene, biotype = bt)
}

.scan_seed <- function(sequence, annodb, seed_name,
                       guide_rna, guide_len) {
  gr <- SeedMatchR::SeedMatchR(
    seqs = annodb$seqs,
    sequence = sequence,
    seed.name = seed_name,
    res.format = "granges"
  )
  if (length(gr) == 0) {
    return(tibble::tibble(transcript = character(),
                          position   = integer(),
                          match_type = character(),
                          target_site = character(),
                          n_mismatch  = integer()))
  }
  tx <- as.character(GenomicRanges::seqnames(gr))
  pos <- as.integer(BiocGenerics::start(gr))
  ts <- .extract_target_site(annodb$seqs, tx, pos, guide_len)
  guide_rc <- as.character(Biostrings::reverseComplement(
                Biostrings::RNAString(guide_rna)))
  nmm <- .count_mismatches(rep(guide_rc, length(tx)), ts)
  tibble::tibble(
    transcript  = tx,
    position    = pos,
    match_type  = "partial_match",
    target_site = ts,
    n_mismatch  = nmm
  )
}

.scan_full <- function(sequence, annodb, max_mismatch,
                       guide_rna, guide_len) {
  guide_dna <- Biostrings::DNAString(gsub("U", "T", as.character(sequence)))
  target    <- Biostrings::reverseComplement(guide_dna)

  mi <- Biostrings::vmatchPattern(target, annodb$seqs,
                                   max.mismatch = max_mismatch,
                                   with.indels = FALSE,
                                   fixed = TRUE)
  lens <- lengths(mi)
  if (sum(lens) == 0) {
    return(tibble::tibble(transcript = character(),
                          position   = integer(),
                          match_type = character(),
                          target_site = character(),
                          n_mismatch  = integer()))
  }
  tx_ids <- rep(names(mi), lens)
  starts <- unlist(BiocGenerics::start(mi), use.names = FALSE)

  hit_seqs <- Biostrings::subseq(annodb$seqs[tx_ids],
                                  start = starts,
                                  width = length(target))
  nmm <- as.integer(Biostrings::neditAt(target, hit_seqs, fixed = TRUE))
  mt  <- ifelse(nmm == 0L, "full_complementarity", "partial_match")

  ts <- chartr("T", "U", as.character(hit_seqs))

  tibble::tibble(
    transcript  = tx_ids,
    position    = as.integer(starts),
    match_type  = mt,
    target_site = ts,
    n_mismatch  = as.integer(nmm)
  )
}

.merge_hits <- function(seed_hits, full_hits, meta, region_label) {
  m <- dplyr::bind_rows(seed_hits, full_hits)
  if (nrow(m) == 0) return(m)

  priority <- c(full_complementarity = 1L, partial_match = 2L, unknown = 3L)
  m$.pri <- priority[m$match_type]
  m <- m[order(m$transcript, m$position, m$.pri), , drop = FALSE]
  m <- m[!duplicated(m[, c("transcript", "position")]), , drop = FALSE]
  m$.pri <- NULL

  m$gene_name <- meta$gene_name[match(m$transcript, meta$tx_id)]
  m$biotype   <- meta$biotype[match(m$transcript, meta$tx_id)]
  m$region    <- region_label
  m
}

# Save raw, low-level scan output (one TSV per siRNA) before merging/dedup.
# Includes BOTH SeedMatchR seed hits and Biostrings full-comp hits, with
# scan_engine to distinguish them.
.save_raw_scan <- function(seed_hits, full_hits, save_dir, siRNA_name,
                           guide_seq, seed_name, max_mismatch_full) {
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

  rows <- list()
  if (nrow(seed_hits) > 0) {
    rows[[1]] <- tibble::tibble(
      siRNA_name  = siRNA_name,
      scan_engine = sprintf("SeedMatchR (%s)", seed_name),
      transcript  = seed_hits$transcript,
      position    = seed_hits$position,
      match_type  = seed_hits$match_type,
      target_site = seed_hits$target_site,
      n_mismatch  = seed_hits$n_mismatch
    )
  }
  if (nrow(full_hits) > 0) {
    rows[[2]] <- tibble::tibble(
      siRNA_name  = siRNA_name,
      scan_engine = sprintf("Biostrings::vmatchPattern (mm<=%d)",
                            max_mismatch_full),
      transcript  = full_hits$transcript,
      position    = full_hits$position,
      match_type  = full_hits$match_type,
      target_site = full_hits$target_site,
      n_mismatch  = full_hits$n_mismatch
    )
  }
  raw <- if (length(rows) > 0) dplyr::bind_rows(rows) else
    tibble::tibble(siRNA_name = character(), scan_engine = character(),
                   transcript = character(), position = integer(),
                   match_type = character(), target_site = character(),
                   n_mismatch = integer())

  fn <- file.path(save_dir, sprintf("raw_%s.tsv",
                                     gsub("[^A-Za-z0-9_.-]", "_",
                                          siRNA_name)))
  header <- sprintf(
    "# TargetSureR raw scan output\n# siRNA_name: %s\n# guide: %s\n# seed_name: %s\n# max_mismatch_full: %d\n# n_seed_hits: %d\n# n_full_hits: %d\n",
    siRNA_name, guide_seq, seed_name, max_mismatch_full,
    nrow(seed_hits), nrow(full_hits))
  writeLines(header, fn)
  suppressWarnings(
    write.table(raw, fn, append = TRUE, sep = "\t",
                row.names = FALSE, quote = FALSE)
  )
  invisible(fn)
}
