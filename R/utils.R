#' Internal helpers
#' @keywords internal
#' @importFrom rlang .data
#' @importFrom utils getFromNamespace write.table
#' @noRd
NULL

#' Pipe operator re-exported from dplyr
#' @importFrom dplyr %>%
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
NULL

require_pkg <- function(pkg, reason = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- sprintf("Package '%s' is required", pkg)
    if (!is.null(reason)) msg <- sprintf("%s for %s", msg, reason)
    msg <- sprintf("%s. Install with: BiocManager::install('%s')", msg, pkg)
    stop(msg, call. = FALSE)
  }
  invisible(TRUE)
}

abort_missing_cols <- function(df, required, context) {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(sprintf("%s: missing required column(s): %s",
                 context, paste(missing, collapse = ", ")),
         call. = FALSE)
  }
  invisible(TRUE)
}

coalesce_col <- function(df, col, default = NA) {
  if (!col %in% names(df)) df[[col]] <- default
  df
}

safe_toupper <- function(x) {
  if (is.null(x)) return(character(0))
  toupper(as.character(x))
}

normalize_01 <- function(x) {
  x <- as.numeric(x)
  r <- range(x, na.rm = TRUE)
  if (!is.finite(diff(r)) || diff(r) == 0) return(rep(0, length(x)))
  (x - r[1]) / diff(r)
}

# Extract a transcript subsequence at (transcript, start, width), return as
# RNA character (T -> U), 5'->3'. Returns NA when out of bounds.
.extract_target_site <- function(seqs, transcript, start, width) {
  out <- rep(NA_character_, length(transcript))
  for (i in seq_along(transcript)) {
    s <- seqs[[transcript[i]]]
    if (is.null(s)) next
    en <- start[i] + width - 1L
    if (start[i] < 1L || en > length(s)) next
    sub <- as.character(Biostrings::subseq(s, start[i], en))
    out[i] <- chartr("T", "U", sub)
  }
  out
}

.count_mismatches <- function(a, b) {
  n <- length(a)
  out <- integer(n)
  for (i in seq_len(n)) {
    if (is.na(a[i]) || is.na(b[i])) { out[i] <- NA_integer_; next }
    sa <- strsplit(a[i], "")[[1]]
    sb <- strsplit(b[i], "")[[1]]
    if (length(sa) != length(sb)) { out[i] <- NA_integer_; next }
    out[i] <- sum(sa != sb)
  }
  out
}

# Build a 3-line alignment string. guide is the siRNA RNA sequence (5'->3'),
# target is the mRNA subsequence at the hit (5'->3', RNA notation, same length).
# The match line uses '|' where guide-revcomp matches target, ' ' otherwise.
.format_alignment <- function(guide, target) {
  out <- character(length(target))
  for (i in seq_along(target)) {
    if (is.na(guide[i]) || is.na(target[i])) { out[i] <- NA_character_; next }
    g <- guide[i]
    t <- target[i]
    if (nchar(g) != nchar(t)) {
      out[i] <- sprintf("guide  5'-%s-3'\ntarget 5'-%s-3'\n# length mismatch (%d vs %d)",
                        g, t, nchar(g), nchar(t))
      next
    }
    g_rc <- as.character(Biostrings::reverseComplement(
              Biostrings::RNAString(g)))
    sg <- strsplit(g_rc, "")[[1]]
    st <- strsplit(t,    "")[[1]]
    match_line <- ifelse(sg == st, "|", " ")
    out[i] <- sprintf("guide  5'-%s-3'\ntarget 5'-%s-3'\n          %s",
                      g, t, paste(match_line, collapse = ""))
  }
  out
}

# Session-scoped cache for Ensembl REST lookups.
# Keyed by transcript id (without version suffix). Values are 1-row tibbles
# matching the schema of `.lookup_via_rest()`.
.tsr_rest_cache <- new.env(parent = emptyenv())

.tsr_rest_cache_clear <- function() {
  rm(list = ls(.tsr_rest_cache, all.names = TRUE), envir = .tsr_rest_cache)
  invisible(NULL)
}

# Low-level POST helper. Isolated so tests can mock it.
.tsr_rest_post <- function(path, body,
                           base_url = "https://rest.ensembl.org") {
  require_pkg("httr2", "REST API fallback")
  req <- httr2::request(base_url)
  req <- httr2::req_url_path_append(req, path)
  req <- httr2::req_method(req, "POST")
  req <- httr2::req_headers(req,
                            `Content-Type` = "application/json",
                            Accept         = "application/json")
  req <- httr2::req_body_json(req, body)
  req <- httr2::req_throttle(req, rate = 15 / 1)
  req <- httr2::req_retry(req, max_tries = 3)
  req <- httr2::req_timeout(req, 30)
  resp <- httr2::req_perform(req)
  httr2::resp_body_json(resp, simplifyVector = FALSE)
}

# Walk exons to derive transcript-relative utr5_len / cds_len / tx_len from a
# single REST `/lookup/id?expand=1` payload. Returns NULL if the payload is
# missing required fields. ncRNA (no Translation) -> cds_len = NA.
.tx_lens_from_rest <- function(entry) {
  exons <- entry$Exon
  if (is.null(exons) || length(exons) == 0) return(NULL)

  starts <- vapply(exons, function(e) as.integer(e$start), integer(1))
  ends   <- vapply(exons, function(e) as.integer(e$end),   integer(1))
  widths <- ends - starts + 1L
  tx_len <- sum(widths)
  strand <- as.integer(entry$strand %||% 1L)

  ord <- if (strand >= 0) order(starts) else order(-starts)
  starts_o <- starts[ord]; ends_o <- ends[ord]; widths_o <- widths[ord]

  tr <- entry$Translation
  if (is.null(tr)) {
    return(list(utr5_len = NA_integer_, cds_len = NA_integer_, tx_len = tx_len))
  }
  tr_start_g <- as.integer(tr$start)
  tr_end_g   <- as.integer(tr$end)

  cds_g_lo <- if (strand >= 0) tr_start_g else tr_end_g
  cds_g_hi <- if (strand >= 0) tr_end_g   else tr_start_g

  pos_in_tx <- function(g_pos) {
    cum <- 0L
    for (i in seq_along(starts_o)) {
      s <- starts_o[i]; e <- ends_o[i]
      if (strand >= 0) {
        if (g_pos >= s && g_pos <= e) return(cum + (g_pos - s + 1L))
      } else {
        if (g_pos >= s && g_pos <= e) return(cum + (e - g_pos + 1L))
      }
      cum <- cum + widths_o[i]
    }
    NA_integer_
  }

  cds_start_tx <- pos_in_tx(cds_g_lo)
  cds_end_tx   <- pos_in_tx(cds_g_hi)
  if (is.na(cds_start_tx) || is.na(cds_end_tx)) {
    return(list(utr5_len = NA_integer_, cds_len = NA_integer_, tx_len = tx_len))
  }
  utr5_len <- cds_start_tx - 1L
  cds_len  <- cds_end_tx - cds_start_tx + 1L
  list(utr5_len = utr5_len, cds_len = cds_len, tx_len = tx_len)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Resolve a vector of transcript ids via Ensembl REST. Returns a tibble:
#   tx_id, gene_name, tx_biotype, utr5_len, cds_len, tx_len
# Empty tibble on any failure (network, parse, throttle exhaustion).
.lookup_via_rest <- function(tx_ids, chunk_size = 500L) {
  empty <- tibble::tibble(
    tx_id      = character(),
    gene_name  = character(),
    tx_biotype = character(),
    utr5_len   = integer(),
    cds_len    = integer(),
    tx_len     = integer()
  )
  tx_ids <- unique(stats::na.omit(as.character(tx_ids)))
  if (length(tx_ids) == 0) return(empty)

  cached <- list()
  todo <- character()
  for (id in tx_ids) {
    if (exists(id, envir = .tsr_rest_cache, inherits = FALSE)) {
      cached[[id]] <- get(id, envir = .tsr_rest_cache, inherits = FALSE)
    } else {
      todo <- c(todo, id)
    }
  }

  if (length(todo) > 0) {
    fetched <- tryCatch(
      .tsr_rest_lookup_chunked(todo, chunk_size = chunk_size),
      error = function(e) {
        warning("REST fallback unavailable: ", conditionMessage(e),
                call. = FALSE)
        NULL
      }
    )
    if (!is.null(fetched)) {
      for (i in seq_len(nrow(fetched))) {
        id <- fetched$tx_id[i]
        assign(id, fetched[i, , drop = FALSE], envir = .tsr_rest_cache)
        cached[[id]] <- fetched[i, , drop = FALSE]
      }
    }
  }
  if (length(cached) == 0) return(empty)
  out <- dplyr::bind_rows(cached)
  tibble::as_tibble(out)
}

.tsr_rest_lookup_chunked <- function(tx_ids, chunk_size) {
  chunks <- split(tx_ids, ceiling(seq_along(tx_ids) / chunk_size))
  rows <- list()
  gene_id_to_name <- list()
  for (chk in chunks) {
    body <- list(ids = as.list(chk), expand = 1L)
    resp <- .tsr_rest_post("/lookup/id", body)
    for (id in chk) {
      entry <- resp[[id]]
      if (is.null(entry)) next
      lens <- .tx_lens_from_rest(entry)
      parent_gid <- entry$Parent
      rows[[id]] <- tibble::tibble(
        tx_id      = id,
        gene_name  = NA_character_,
        tx_biotype = entry$biotype %||% NA_character_,
        utr5_len   = if (is.null(lens)) NA_integer_ else as.integer(lens$utr5_len),
        cds_len    = if (is.null(lens)) NA_integer_ else as.integer(lens$cds_len),
        tx_len     = if (is.null(lens)) NA_integer_ else as.integer(lens$tx_len),
        gene_id    = parent_gid %||% NA_character_
      )
    }
  }
  if (length(rows) == 0) {
    return(tibble::tibble(tx_id = character(), gene_name = character(),
                          tx_biotype = character(), utr5_len = integer(),
                          cds_len = integer(), tx_len = integer()))
  }
  tx_tbl <- dplyr::bind_rows(rows)

  gene_ids <- unique(stats::na.omit(tx_tbl$gene_id))
  if (length(gene_ids) > 0) {
    gchunks <- split(gene_ids, ceiling(seq_along(gene_ids) / chunk_size))
    for (gchk in gchunks) {
      gresp <- .tsr_rest_post("/lookup/id", list(ids = as.list(gchk)))
      for (gid in gchk) {
        ge <- gresp[[gid]]
        if (is.null(ge)) next
        gene_id_to_name[[gid]] <- ge$display_name %||% ge$external_name %||% NA_character_
      }
    }
    gn_vec <- unlist(gene_id_to_name)
    tx_tbl$gene_name <- unname(gn_vec[tx_tbl$gene_id])
  }
  tx_tbl$gene_id <- NULL
  tx_tbl
}

# Robust loader: works whether the package is installed (utils::data())
# or sourced during development (object already in globalenv).
.load_pkg_data <- function(name) {
  if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
    return(get(name, envir = .GlobalEnv))
  }
  e <- new.env()
  ok <- tryCatch({
    utils::data(list = name, package = "TargetSureR", envir = e)
    TRUE
  }, error = function(err) FALSE,
     warning = function(err) FALSE)
  if (ok && exists(name, envir = e, inherits = FALSE)) return(get(name, envir = e))
  rda <- file.path("data", paste0(name, ".rda"))
  if (file.exists(rda)) {
    load(rda, envir = e)
    if (exists(name, envir = e, inherits = FALSE)) return(get(name, envir = e))
  }
  stop(sprintf("Internal dataset '%s' not found", name), call. = FALSE)
}
