#' TargetSureR: siRNA off-target annotation, risk ranking, tissue and species analysis
#'
#' @description
#' A modular toolkit for evaluating small interfering RNA (siRNA) off-target
#' effects. The package has four layers:
#'
#' \describe{
#'   \item{Annotation}{\code{\link{annotate_sites}} produces a 10-column
#'     standard off-target table from either siRNA sequences (Mode A, via
#'     \pkg{SeedMatchR}) or (transcript, position) tuples (Mode B).}
#'   \item{Enrichment}{\code{\link{enrich_annotations}} appends curated
#'     gene-list flags and GTEx tissue expression to the core table.}
#'   \item{Analysis}{\code{\link{rank_sirna}}, \code{\link{tissue_safety}},
#'     and \code{\link{species_match}}.}
#'   \item{Reporting}{\code{\link{generate_report}} bundles CSVs and plots.}
#' }
#'
#' @section Core table schema:
#' The 10-column standard off-target table:
#' \code{siRNA_name}, \code{siRNA}, \code{strand}, \code{transcript},
#' \code{gene_name}, \code{biotype}, \code{position}, \code{region},
#' \code{match_type}, \code{risk_score}. See \code{\link{core_columns}}.
#'
#' @keywords internal
"_PACKAGE"
