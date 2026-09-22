#' Enrich the core off-target table with gene-list flags and tissue expression
#'
#' Left-joins three curated gene-list boolean flags
#' (\code{AE_gene}, \code{cancer_gene_v2}, \code{immune_gene_immport}) and a
#' tissue expression matrix onto a core off-target table produced by
#' \code{\link{annotate_sites}}.
#'
#' @param sites A data frame with at least \code{\link{core_columns}}.
#' @param gene_lists A named list with character vectors for
#'   \code{AE_gene}, \code{cancer_gene_v2}, \code{immune_gene_immport}.
#'   If \code{NULL}, the package defaults from \code{\link{load_gene_lists}}
#'   are used.
#' @param expression_matrix A data frame whose first column is
#'   \code{gene_symbol} (matched case-insensitively against
#'   \code{gene_name}), with one numeric column per tissue. If \code{NULL},
#'   the \strong{simulated} 20-gene example
#'   (\code{\link{load_gtex_example}}) is used, which leaves nearly every
#'   real off-target gene unquantified; supply a measured matrix for
#'   analysis. Check the \code{"genes_missing_expression"} attribute of the
#'   result to see how many genes the matrix did not cover.
#' @param fill_missing_tpm Value used for genes absent from
#'   \code{expression_matrix}. Default \code{NA_real_}, which keeps
#'   \dQuote{expression unknown} distinct from \dQuote{measured as zero};
#'   pass \code{0} to restore the previous fill behaviour.
#'
#' @return A tibble whose columns are \code{\link{enriched_columns}} followed
#'   by one numeric column per tissue. The number of off-target genes absent
#'   from \code{expression_matrix} is reported as the
#'   \code{"genes_missing_expression"} attribute.
#' @examples
#' \dontrun{
#' enriched <- enrich_annotations(sites)
#' }
#' @export
enrich_annotations <- function(sites,
                               gene_lists = NULL,
                               expression_matrix = NULL,
                               fill_missing_tpm = NA_real_) {
  abort_missing_cols(sites, core_columns(), "enrich_annotations")

  if (is.null(gene_lists))        gene_lists        <- load_gene_lists()
  if (is.null(expression_matrix)) expression_matrix <- load_gtex_example()

  x <- tibble::as_tibble(sites)
  gene_upper <- safe_toupper(x$gene_name)

  x$AE_gene             <- gene_upper %in% safe_toupper(gene_lists$AE_gene)
  x$cancer_gene_v2      <- gene_upper %in% safe_toupper(gene_lists$cancer_gene_v2)
  x$immune_gene_immport <- gene_upper %in% safe_toupper(gene_lists$immune_gene_immport)

  tpm <- .prepare_tpm(expression_matrix)
  tissues <- setdiff(names(tpm), "gene_symbol")
  idx <- match(gene_upper, tpm$gene_symbol)
  tpm_block <- tpm[idx, tissues, drop = FALSE]
  missing <- is.na(idx)
  if (any(missing)) tpm_block[missing, ] <- fill_missing_tpm

  out <- dplyr::bind_cols(x[, enriched_columns(), drop = FALSE], tpm_block)
  out <- tibble::as_tibble(out)
  attr(out, "genes_missing_expression") <-
    length(unique(gene_upper[missing]))
  out
}

.prepare_tpm <- function(em) {
  if (!"gene_symbol" %in% names(em)) {
    stop("expression_matrix must have a 'gene_symbol' column",
         call. = FALSE)
  }
  em$gene_symbol <- toupper(as.character(em$gene_symbol))
  em[!duplicated(em$gene_symbol), , drop = FALSE]
}
