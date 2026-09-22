#' Bundle a report directory from an enriched off-target table
#'
#' Writes the enriched table, the ranking table, the tissue-safety table,
#' and (optionally) the species-match result to CSV files under
#' \code{output_dir}. With \code{with_figures = TRUE}, also writes the
#' six standard figures via \code{\link{generate_figures}}.
#'
#' @param enriched Output of \code{\link{enrich_annotations}}.
#' @param output_dir Directory to write files into. Created if absent.
#' @param target_sirnas Passed to \code{\link{rank_sirna}}.
#' @param primary_tissue Passed to \code{\link{tissue_safety}}.
#' @param species_result Optional output of \code{\link{species_match}}.
#' @param with_figures Logical. If \code{TRUE}, also writes the standard
#'   6-figure bundle (PDF + PNG) into \code{output_dir/figures/}.
#'   Requires the optional plotting Suggests (\pkg{ggplot2}, \pkg{tidyr},
#'   \pkg{ggrepel}, \pkg{RColorBrewer}). Default \code{FALSE} for
#'   backward compatibility.
#'
#' @return Invisibly, a list of written file paths.
#' @examples
#' \dontrun{
#' generate_report(enriched, output_dir = tempfile("report_"),
#'                 with_figures = TRUE)
#' }
#' @export
generate_report <- function(enriched,
                            output_dir,
                            target_sirnas = character(),
                            primary_tissue = NULL,
                            species_result = NULL,
                            with_figures = FALSE) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  paths <- list()

  paths$enriched <- file.path(output_dir, "enriched_offtargets.csv")
  readr::write_csv(enriched, paths$enriched)

  ranking <- rank_sirna(enriched, target_sirnas = target_sirnas)
  paths$ranking <- file.path(output_dir, "sirna_ranking.csv")
  readr::write_csv(ranking, paths$ranking)

  tissue <- tryCatch(
    tissue_safety(enriched, primary_tissue = primary_tissue),
    error = function(e) NULL
  )
  if (!is.null(tissue)) {
    paths$tissue <- file.path(output_dir, "gene_tissue_safety.csv")
    readr::write_csv(tissue, paths$tissue)
  }

  if (!is.null(species_result)) {
    paths$species_score <- file.path(output_dir, "species_score_matrix.csv")
    paths$site_detail   <- file.path(output_dir, "site_conservation_detail.csv")
    paths$gene_detail   <- file.path(output_dir, "gene_conservation_detail.csv")
    readr::write_csv(species_result$species_score, paths$species_score)
    readr::write_csv(species_result$site_detail,   paths$site_detail)
    if (!is.null(species_result$gene_detail))
      readr::write_csv(species_result$gene_detail, paths$gene_detail)
  }

  if (with_figures && !is.null(tissue)) {
    fig_dir <- file.path(output_dir, "figures")
    fig_paths <- generate_figures(enriched, ranking, tissue,
                                   species_result = species_result,
                                   output_dir = fig_dir)
    paths$figures <- fig_paths
  }

  invisible(paths)
}
