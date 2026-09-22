#' Annotate off-target sites from (transcript, position) input (Mode B)
#'
#' Given a minimal data frame of off-target positions, resolves gene symbol,
#' transcript biotype, and transcript region using an \pkg{ensembldb}
#' \code{EnsDb} object. Transcripts not found in the local EnsDb are
#' re-queried against the Ensembl REST API
#' (\url{https://rest.ensembl.org/lookup/id}) when \code{rest_fallback = TRUE}
#' (the default), so newer transcript IDs released after the local annotation
#' snapshot can still be resolved. The result conforms to
#' \code{\link{core_columns}}.
#'
#' @param input A data frame with required columns \code{siRNA_name},
#'   \code{transcript}, \code{position}. Optional columns \code{siRNA} and
#'   \code{strand} are preserved if present.
#' @param species Organism key. Used only when both \code{ensdb} and
#'   \code{ensdb_path} are \code{NULL}, in which case
#'   \code{\link[AnnotationHub]{AnnotationHub}} is queried.
#' @param ensdb An \pkg{ensembldb} \code{EnsDb} object. If \code{NULL},
#'   the function falls back to \code{ensdb_path}, then to
#'   \pkg{AnnotationHub} based on \code{species}.
#' @param ensdb_path Path to a locally stored EnsDb SQLite file. Loaded
#'   via \code{\link[ensembldb]{EnsDb}}. Useful when
#'   \pkg{AnnotationHub} downloads are unreliable or when you want to pin
#'   the annotation to a fixed release.
#' @param strand_default Default strand assignment when the input has no
#'   \code{strand} column. Defaults to \code{"guide"}.
#' @param rest_fallback Logical. If \code{TRUE} (default), transcripts not
#'   found in the local EnsDb are re-queried against the Ensembl REST API.
#'   Requires the \pkg{httr2} package and network access. When the network
#'   is unavailable a warning is emitted and unresolved transcripts retain
#'   \code{NA} for \code{gene_name} / \code{biotype}.
#' @param with_lookup_source Logical. If \code{TRUE}, append an extra
#'   column \code{lookup_source} with values \code{"ensdb"}, \code{"rest"},
#'   or \code{"unresolved"} indicating where each annotation came from.
#'   Useful for diagnostics; default \code{FALSE} preserves the
#'   10-column \code{\link{core_columns}} contract.
#'
#' @return A tibble with \code{\link{core_columns}}, optionally followed by
#'   a \code{lookup_source} column when \code{with_lookup_source = TRUE}.
#'
#' @examples
#' \dontrun{
#' # Pin to a specific Ensembl release via a local SQLite file (recommended)
#' annotate_from_position(
#'   data.frame(siRNA_name = "hsiR22",
#'              transcript = "ENST00000426385",
#'              position   = 379),
#'   ensdb_path = "~/ensdb/ensembl_113_hsapiens.sqlite")
#'
#' # Or let AnnotationHub fetch the latest EnsDb (network)
#' annotate_from_position(
#'   data.frame(siRNA_name = "hsiR22",
#'              transcript = "ENST00000426385",
#'              position   = 379),
#'   species = "human")
#' }
#' @export
annotate_from_position <- function(input,
                                   species = "human",
                                   ensdb = NULL,
                                   ensdb_path = NULL,
                                   strand_default = "guide",
                                   rest_fallback = TRUE,
                                   with_lookup_source = FALSE) {
  stopifnot(is.data.frame(input))
  abort_missing_cols(input, c("siRNA_name", "transcript", "position"),
                     "annotate_from_position")

  x <- tibble::as_tibble(input)
  x <- coalesce_col(x, "siRNA",  NA_character_)
  x <- coalesce_col(x, "strand", strand_default)

  x$tx_key <- sub("\\..*$", "", x$transcript)

  if (is.null(ensdb)) {
    if (!is.null(ensdb_path)) {
      require_pkg("ensembldb", "EnsDb loading from file")
      if (!file.exists(ensdb_path)) {
        stop("ensdb_path does not exist: ", ensdb_path, call. = FALSE)
      }
      ensdb <- ensembldb::EnsDb(ensdb_path)
    } else {
      ensdb <- resolve_ensdb(species)
    }
  }

  meta <- .lookup_transcript_meta(x$tx_key, ensdb)
  ensdb_resolved <- meta$tx_id

  rest_meta <- NULL
  if (rest_fallback) {
    unresolved <- setdiff(unique(stats::na.omit(x$tx_key)), ensdb_resolved)
    if (length(unresolved) > 0 &&
        requireNamespace("httr2", quietly = TRUE)) {
      rest_meta <- .lookup_via_rest(unresolved)
      if (nrow(rest_meta) > 0) {
        meta <- dplyr::bind_rows(
          meta,
          tibble::tibble(tx_id      = rest_meta$tx_id,
                         gene_name  = rest_meta$gene_name,
                         tx_biotype = rest_meta$tx_biotype)
        )
      }
    }
  }

  x$gene_name <- meta$gene_name[match(x$tx_key, meta$tx_id)]
  x$biotype   <- meta$tx_biotype[match(x$tx_key, meta$tx_id)]

  x$region <- .classify_region(x$tx_key, x$position, ensdb,
                               rest_lens = rest_meta)

  x$match_type <- .infer_match_type(x)
  x$risk_score <- compute_risk_score(x$region, x$match_type)

  if (with_lookup_source) {
    src <- rep("unresolved", nrow(x))
    src[x$tx_key %in% ensdb_resolved] <- "ensdb"
    if (!is.null(rest_meta) && nrow(rest_meta) > 0) {
      src[x$tx_key %in% rest_meta$tx_id] <- "rest"
    }
    x$lookup_source <- src
    x <- x[, c(core_columns(), "lookup_source"), drop = FALSE]
  } else {
    x <- x[, core_columns(), drop = FALSE]
  }
  x
}

#' Resolve an EnsDb from a species key
#'
#' @param species Character. Common name or Latin name understood here.
#' @return An \code{EnsDb} object.
#' @keywords internal
#' @noRd
resolve_ensdb <- function(species) {
  require_pkg("AnnotationHub", "transcript annotation in Mode B")
  require_pkg("ensembldb",     "transcript annotation in Mode B")

  key <- tolower(species)
  species_map <- c(
    human  = "Homo sapiens",
    mouse  = "Mus musculus",
    rat    = "Rattus norvegicus",
    cyno   = "Macaca fascicularis",
    monkey = "Macaca fascicularis",
    rhesus = "Macaca mulatta",
    rabbit = "Oryctolagus cuniculus",
    dog    = "Canis lupus familiaris",
    beagle = "Canis lupus familiaris"
  )
  latin <- species_map[key]
  if (is.na(latin)) latin <- species

  ah <- AnnotationHub::AnnotationHub(ask = FALSE)
  q <- AnnotationHub::query(ah, c("EnsDb", latin))
  if (length(q) == 0) {
    stop(sprintf("No EnsDb found for '%s' in AnnotationHub", species),
         call. = FALSE)
  }
  q[[length(q)]]  # latest available
}

#' Look up transcript metadata via ensembldb
#' @keywords internal
#' @noRd
.lookup_transcript_meta <- function(tx_ids, ensdb) {
  require_pkg("ensembldb",        "transcript annotation")
  require_pkg("AnnotationFilter", "transcript annotation")
  tx_ids <- unique(stats::na.omit(tx_ids))
  if (length(tx_ids) == 0) {
    return(tibble::tibble(tx_id = character(),
                          gene_name = character(),
                          tx_biotype = character()))
  }
  flt <- AnnotationFilter::TxIdFilter(tx_ids)
  df <- ensembldb::transcripts(
    ensdb, filter = flt,
    columns = c("tx_id", "tx_biotype", "gene_name"),
    return.type = "data.frame"
  )
  tibble::as_tibble(df[, c("tx_id", "gene_name", "tx_biotype"), drop = FALSE])
}

#' Classify each (transcript, position) into a region label
#'
#' Uses \code{ensembldb}'s internal transcript-length table to fetch
#' \code{utr5_len}, \code{cds_len}, \code{tx_len} for each transcript in a
#' single query, then classifies each (transcript, position) by direct
#' numerical comparison against transcript-relative coordinates.
#'
#' Logic per transcript (position is 1-based, relative to 5' end of the
#' mature mRNA):
#' \itemize{
#'   \item tx_id not in EnsDb:                      \code{NA} -> \code{"ncRNA"}.
#'   \item No CDS recorded (non-coding transcript): \code{"ncRNA"}.
#'   \item position <= utr5_len:                    \code{"5'UTR"}.
#'   \item position <= utr5_len + cds_len:          \code{"CDS"}.
#'   \item otherwise:                               \code{"3'UTR"}.
#' }
#'
#' Note: \code{tx_cds_seq_start/end} returned by
#' \code{ensembldb::transcripts(columns = ...)} are GENOMIC coordinates, not
#' transcript-relative; this function deliberately avoids them.
#'
#' @keywords internal
#' @noRd
.classify_region <- function(tx_ids, positions, ensdb, rest_lens = NULL) {
  require_pkg("ensembldb", "region classification")

  n <- length(tx_ids)
  if (n == 0) return(character())

  unique_tx <- unique(stats::na.omit(tx_ids))
  if (length(unique_tx) == 0) return(rep("ncRNA", n))

  # ensembldb:::.transcriptLengths returns utr5_len, cds_len, tx_len
  tx_lens <- tryCatch(
    {
      fn <- getFromNamespace(".transcriptLengths", "ensembldb")
      fn(ensdb, with.utr5_len = TRUE, with.cds_len = TRUE,
         filter = AnnotationFilter::TxIdFilter(unique_tx))
    },
    error = function(e) NULL
  )
  if (is.null(tx_lens)) {
    tx_lens <- data.frame(tx_id = character(),
                          utr5_len = integer(),
                          cds_len = integer(),
                          tx_len = integer(),
                          stringsAsFactors = FALSE)
  }

  if (!is.null(rest_lens) && nrow(rest_lens) > 0) {
    extra <- rest_lens[!rest_lens$tx_id %in% tx_lens$tx_id,
                       c("tx_id", "utr5_len", "cds_len", "tx_len"),
                       drop = FALSE]
    if (nrow(extra) > 0) {
      tx_lens <- dplyr::bind_rows(
        tx_lens[, c("tx_id", "utr5_len", "cds_len", "tx_len"), drop = FALSE],
        extra
      )
    }
  }

  if (nrow(tx_lens) == 0) return(rep("ncRNA", n))

  idx <- match(tx_ids, tx_lens$tx_id)
  utr5 <- tx_lens$utr5_len[idx]
  cdsl <- tx_lens$cds_len[idx]

  out <- rep(NA_character_, n)
  not_in_db <- is.na(idx)
  no_cds    <- !not_in_db & (is.na(cdsl) | cdsl == 0)

  out[not_in_db] <- "ncRNA"
  out[no_cds]    <- "ncRNA"

  coding <- !not_in_db & !no_cds
  u <- ifelse(is.na(utr5[coding]), 0L, utr5[coding])
  c <- cdsl[coding]
  p <- positions[coding]
  out[coding] <- ifelse(p <= u,         "5'UTR",
                 ifelse(p <= u + c,     "CDS",
                                         "3'UTR"))
  out
}

.infer_match_type <- function(x) {
  rep("unknown", nrow(x))
}
