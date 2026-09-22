#' Column specifications for the TargetSureR pipeline
#'
#' Named character vectors describing the expected columns at each pipeline
#' stage. Use these in validation and downstream code rather than hard-coding
#' column names.
#'
#' @return A character vector of column names.
#' @examples
#' core_columns()
#' enriched_columns()
#' @name column_specs
NULL

#' @rdname column_specs
#' @export
core_columns <- function() {
  c("siRNA_name", "siRNA", "strand", "transcript", "gene_name",
    "biotype", "position", "region", "match_type", "risk_score")
}

#' @rdname column_specs
#' @export
enriched_columns <- function() {
  c(core_columns(),
    "AE_gene", "cancer_gene_v2", "immune_gene_immport")
}
