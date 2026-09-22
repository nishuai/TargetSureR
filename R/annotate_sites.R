#' Annotate off-target sites: dispatcher
#'
#' Produces the 10-column standard off-target table. Supports two input modes:
#'
#' \describe{
#'   \item{\code{mode = "position"}}{Input is a data frame with at least
#'     \code{siRNA_name}, \code{transcript}, and \code{position}. Optional
#'     columns: \code{siRNA}, \code{strand}. Uses \pkg{ensembldb} to resolve
#'     \code{gene_name}, \code{biotype}, and \code{region}.}
#'   \item{\code{mode = "sequence"}}{Input is a list with \code{siRNA_name}
#'     and \code{sequence}. Delegates to \pkg{SeedMatchR} to search the
#'     transcriptome.}
#' }
#'
#' Both paths emit the same column set; see \code{\link{core_columns}}.
#'
#' @param input Mode-dependent; see Details.
#' @param mode One of \code{"position"} or \code{"sequence"}.
#' @param species Organism key understood by \pkg{SeedMatchR}
#'   (e.g. \code{"human"}, \code{"mouse"}, \code{"rat"}).
#' @param ... Passed to the mode-specific backend
#'   (\code{\link{annotate_from_position}} or \code{\link{annotate_from_sequence}}).
#'
#' @return A tibble with \code{\link{core_columns}}.
#' @seealso \code{\link{annotate_from_position}}, \code{\link{annotate_from_sequence}}
#' @examples
#' \dontrun{
#' # Mode B - default (uses AnnotationHub or a supplied EnsDb)
#' annotate_sites(
#'   data.frame(siRNA_name = "hsiR22",
#'              transcript = "ENST00000426385.4",
#'              position   = 379),
#'   mode = "position", species = "human")
#'
#' # Mode B - pinned to a local SQLite release
#' annotate_sites(
#'   data.frame(siRNA_name = "hsiR22",
#'              transcript = "ENST00000426385.4",
#'              position   = 379),
#'   mode = "position",
#'   ensdb_path = "~/ensdb/ensembl_113_hsapiens.sqlite")
#'
#' # REST fallback resolves newer transcripts automatically (default ON)
#' annotate_sites(
#'   data.frame(siRNA_name = "hsiR22",
#'              transcript = "ENST00000856451.1",
#'              position   = 200),
#'   mode = "position",
#'   ensdb_path = "~/ensdb/ensembl_113_hsapiens.sqlite",
#'   rest_fallback = TRUE)
#' }
#' @export
annotate_sites <- function(input,
                           mode = c("position", "sequence"),
                           species = "human",
                           ...) {
  mode <- match.arg(mode)
  switch(
    mode,
    position = annotate_from_position(input, species = species, ...),
    sequence = annotate_from_sequence(input, species = species, ...)
  )
}
