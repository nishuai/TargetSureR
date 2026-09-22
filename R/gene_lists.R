#' Load built-in curated gene lists
#'
#' Returns the three gene lists used by \code{\link{enrich_annotations}}:
#' adverse-event genes (\code{AE_gene}), a curated cancer gene panel
#' (\code{cancer_gene_v2}), and a curated immune-related gene panel
#' (\code{immune_gene_immport}).
#'
#' Lists are matched case-insensitively against \code{gene_name}.
#'
#' @return A named list of three character vectors (uppercase symbols):
#'   \code{AE_gene}, \code{cancer_gene_v2}, \code{immune_gene_immport}.
#' @examples
#' gl <- load_gene_lists()
#' lapply(gl, length)
#' @export
load_gene_lists <- function() {
  list(
    AE_gene             = .load_pkg_data("ae_genes"),
    cancer_gene_v2      = .load_pkg_data("cancer_genes_v2"),
    immune_gene_immport = .load_pkg_data("immune_genes_immport")
  )
}

#' Load the example GTEx median-TPM matrix shipped with the package
#'
#' Returns \code{\link{gtex_example}}, a 20-gene by 10-tissue matrix of
#' \strong{simulated} median TPM values shaped like a GTEx matrix. It exists so
#' that examples and tests run without a large download.
#'
#' \strong{The values are invented, not measured expression.} For analysis,
#' download a median-expression matrix from the GTEx portal (or use in-house
#' RNA-seq) and pass it to \code{\link{enrich_annotations}} via the
#' \code{expression_matrix} argument. See
#' \code{vignette("01_quickstart")}, step 4.
#'
#' @return A data frame with a \code{gene_symbol} column and one numeric column
#'   per tissue.
#' @examples
#' head(load_gtex_example())
#' @export
load_gtex_example <- function() {
  .load_pkg_data("gtex_example")
}
