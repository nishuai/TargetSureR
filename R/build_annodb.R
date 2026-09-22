#' Build a transcript sequence index for sequence-mode scanning
#'
#' Constructs the \code{annodb} object that \code{\link{annotate_sites}}
#' requires in sequence mode (Mode A): a set of transcript sequences for one
#' feature class, paired with the coordinate metadata needed to report hits.
#' This is the slowest preparatory step in the workflow, so the result should
#' be cached with \code{\link[base]{saveRDS}} and reloaded in later sessions.
#'
#' @section Prerequisites:
#' Two external resources must be installed before this function can run.
#' Neither ships with TargetSureR; both are one-time downloads.
#'
#' \describe{
#'   \item{\code{BSgenome.Hsapiens.UCSC.hg38}}{The GRCh38 genome sequence,
#'     installed as a Bioconductor annotation package (about 800 MB on disk).
#'     Install with
#'     \code{BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")}.}
#'   \item{An \code{EnsDb} SQLite file}{An Ensembl annotation snapshot, about
#'     510 MB. The analyses published with this package used Ensembl release
#'     113 (AnnotationHub record \code{AH119325}). See the
#'     \dQuote{Introduction to TargetSureR} vignette for retrieval, including
#'     the manual-download route for networks where
#'     \code{AnnotationHub::cache()} fails partway through.}
#' }
#'
#' @section Transcript selection:
#' Four filters are applied, in this order. They define the transcript set
#' that every downstream count is relative to, so a published site count is
#' only reproducible when all four match.
#'
#' \enumerate{
#'   \item \code{tx_biotype == "protein_coding"}.
#'   \item Sequence name in \code{1:22}, \code{X}, \code{Y}, \code{MT} —
#'     scaffolds and patch contigs are dropped.
#'   \item One transcript per gene: the longest sequence \emph{within this
#'     feature class}. A gene may therefore contribute a different
#'     representative transcript to each of the three classes, and a per-gene
#'     site count summed across classes is an upper bound on the number of
#'     distinct positions.
#'   \item Sequences shorter than \code{min_width} (default 8 nt) are dropped,
#'     since a 7-nt seed cannot be located in them.
#' }
#'
#' With Ensembl 113 these filters yield 19,022 3'UTR, 18,973 5'UTR and 19,355
#' coding sequences.
#'
#' @param feature Feature class to build: \code{"3UTR"}, \code{"5UTR"} or
#'   \code{"CDS"}. Case-insensitive.
#' @param ensdb_path Path to the \code{EnsDb} SQLite file. Ignored when
#'   \code{edb} is supplied.
#' @param edb An \code{EnsDb} object, as an alternative to
#'   \code{ensdb_path}. Supply this to build several feature classes from one
#'   open connection.
#' @param bsgenome A \code{BSgenome} object supplying the genome sequence.
#'   Defaults to \code{BSgenome.Hsapiens.UCSC.hg38}, loaded on demand.
#' @param seq_names Sequence names to retain. Defaults to the 22 autosomes
#'   plus \code{X}, \code{Y} and \code{MT}. Pass \code{NULL} to disable this
#'   filter.
#' @param min_width Minimum sequence length in nucleotides. Default 8.
#' @param verbose Print progress and the resulting transcript count.
#'   Default \code{TRUE}.
#'
#' @return A list of three elements, the structure
#'   \code{\link{annotate_sites}} expects as its \code{annodb} argument:
#'   \describe{
#'     \item{\code{seqs}}{A \code{DNAStringSet} of transcript sequences, named
#'       by Ensembl transcript ID.}
#'     \item{\code{gtf}}{A \code{GRanges} of transcript-relative coordinates
#'       carrying \code{tx_id}, \code{gene_name} and \code{tx_biotype}.}
#'     \item{\code{txdb}}{The \code{EnsDb} used, retained so that a result
#'       records its own annotation source.}
#'   }
#'
#' @seealso \code{\link{annotate_sites}} for scanning against the result.
#'
#' @examples
#' \dontrun{
#' # Point these at your own copies; see the vignette for how to obtain them.
#' ensdb <- "~/TargetSureR_data/EnsDb.Hsapiens.v113.sqlite"
#'
#' annodb_3utr <- build_annodb("3UTR", ensdb_path = ensdb)
#' annodb_5utr <- build_annodb("5UTR", ensdb_path = ensdb)
#' annodb_cds  <- build_annodb("CDS",  ensdb_path = ensdb)
#'
#' # Cache once, reload thereafter
#' saveRDS(annodb_3utr, "hg38_3UTR_annodb.rds", compress = "xz")
#' annodb_3utr <- readRDS("hg38_3UTR_annodb.rds")
#' }
#' @export
build_annodb <- function(feature,
                         ensdb_path = NULL,
                         edb = NULL,
                         bsgenome = NULL,
                         seq_names = c(as.character(1:22), "X", "Y", "MT"),
                         min_width = 8L,
                         verbose = TRUE) {
  require_pkg("ensembldb",       "annodb construction")
  require_pkg("GenomicFeatures", "annodb construction")
  require_pkg("GenomicRanges",   "annodb construction")
  require_pkg("GenomeInfoDb",    "annodb construction")
  require_pkg("Biostrings",      "annodb construction")
  require_pkg("IRanges",         "annodb construction")
  require_pkg("S4Vectors",       "annodb construction")
  require_pkg("BiocGenerics",    "annodb construction")

  feature <- .match_feature(feature)
  edb <- .resolve_ensdb(edb, ensdb_path)
  bsgenome <- .resolve_bsgenome(bsgenome)

  if (verbose) {
    message(sprintf("Building %s annodb from Ensembl %s ...",
                    feature, ensembldb::ensemblVersion(edb)))
  }
  t0 <- Sys.time()

  feat <- switch(feature,
    `3UTR` = ensembldb::threeUTRsByTranscript(edb),
    `5UTR` = ensembldb::fiveUTRsByTranscript(edb),
    CDS    = ensembldb::cdsBy(edb, by = "tx"))

  tx_meta <- ensembldb::transcripts(edb,
    columns     = c("tx_id", "tx_biotype", "gene_name", "seq_name"),
    return.type = "data.frame")
  tx_meta <- tx_meta[tx_meta$tx_biotype == "protein_coding", , drop = FALSE]
  if (!is.null(seq_names)) {
    tx_meta <- tx_meta[tx_meta$seq_name %in% seq_names, , drop = FALSE]
  }
  if (nrow(tx_meta) == 0) {
    stop("No protein-coding transcripts remain after filtering; ",
         "check 'seq_names' against the EnsDb's sequence naming.",
         call. = FALSE)
  }

  feat <- feat[names(feat) %in% tx_meta$tx_id]
  if (length(feat) == 0) {
    stop(sprintf("No %s regions remain for the selected transcripts.", feature),
         call. = FALSE)
  }
  GenomeInfoDb::seqlevelsStyle(feat) <- "UCSC"
  seqs <- GenomicFeatures::extractTranscriptSeqs(bsgenome, feat)

  # One transcript per gene: longest within this feature class.
  tx2gene <- tx_meta$gene_name
  names(tx2gene) <- tx_meta$tx_id
  keep <- data.frame(
    tx    = names(seqs),
    gene  = tx2gene[names(seqs)],
    width = BiocGenerics::width(seqs),
    stringsAsFactors = FALSE)
  keep <- keep[order(keep$gene, -keep$width), , drop = FALSE]
  keep <- keep[!duplicated(keep$gene) & !is.na(keep$gene), , drop = FALSE]
  seqs <- seqs[keep$tx]
  seqs <- seqs[BiocGenerics::width(seqs) >= min_width]

  gtf_meta <- tx_meta[tx_meta$tx_id %in% names(seqs), , drop = FALSE]
  gtf <- GenomicRanges::GRanges(
    seqnames = paste0("chr", gtf_meta$seq_name),
    ranges   = IRanges::IRanges(
      start = 1L,
      width = BiocGenerics::width(seqs[gtf_meta$tx_id])),
    strand   = "+")
  S4Vectors::mcols(gtf)$tx_id      <- gtf_meta$tx_id
  S4Vectors::mcols(gtf)$gene_name  <- gtf_meta$gene_name
  S4Vectors::mcols(gtf)$tx_biotype <- gtf_meta$tx_biotype

  if (verbose) {
    message(sprintf("  %d transcripts retained (%.1f s)",
                    length(seqs),
                    as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  list(seqs = seqs, gtf = gtf, txdb = edb)
}

.match_feature <- function(feature) {
  if (length(feature) != 1L || !is.character(feature)) {
    stop("'feature' must be a single string: \"3UTR\", \"5UTR\" or \"CDS\".",
         call. = FALSE)
  }
  key <- toupper(gsub("[^A-Za-z0-9]", "", feature))
  switch(key,
    "3UTR" = "3UTR", "THREEUTR" = "3UTR",
    "5UTR" = "5UTR", "FIVEUTR"  = "5UTR",
    "CDS"  = "CDS",  "CODING"   = "CDS",
    stop(sprintf("Unknown feature '%s'; expected \"3UTR\", \"5UTR\" or \"CDS\".",
                 feature), call. = FALSE))
}

.resolve_ensdb <- function(edb, ensdb_path) {
  if (!is.null(edb)) return(edb)
  if (is.null(ensdb_path)) {
    stop("Supply either 'ensdb_path' (path to an EnsDb SQLite file) or 'edb' ",
         "(an EnsDb object). See the vignette for how to obtain one.",
         call. = FALSE)
  }
  if (!file.exists(ensdb_path)) {
    stop(sprintf("EnsDb file not found: '%s'", ensdb_path), call. = FALSE)
  }
  ensembldb::EnsDb(ensdb_path)
}

.resolve_bsgenome <- function(bsgenome) {
  if (!is.null(bsgenome)) return(bsgenome)
  pkg <- "BSgenome.Hsapiens.UCSC.hg38"
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf(paste0("Package '%s' is required to extract transcript ",
                        "sequences. Install with: BiocManager::install('%s') ",
                        "(about 800 MB), or pass a different genome via ",
                        "'bsgenome'."), pkg, pkg), call. = FALSE)
  }
  getExportedValue(pkg, pkg)
}
