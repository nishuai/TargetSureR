#' Default composite-score weights for \code{rank_sirna}
#'
#' Four equally weighted count dimensions used by \code{\link{rank_sirna}};
#' the weights sum to 1. Two describe where sites fall, and therefore which
#' silencing mechanism is in play: \code{cds_count} (coding-sequence sites,
#' cleavage) and \code{seed_region_count} (3'UTR seed sites, miRNA-like
#' translational repression). These are the least correlated pair among the
#' available count dimensions (Spearman rho 0.58 across the built-in
#' \code{\link{reference_set}}). Two describe engagement of curated
#' cancer / adverse-event / immune genes: \code{critical_gene_count} (hits)
#' and \code{high_risk_critical} (those in high-risk regions).
#'
#' Equal weights are deliberate: any other split would encode an unvalidated
#' claim about relative importance. Pass a named vector to
#' \code{rank_sirna(weights = ...)} to override.
#'
#' Three dimensions available in \code{reference_set$metrics} are excluded.
#' \code{high_risk_count} is collinear with \code{cds_count} (rho 1.000).
#' \code{full_comp_count} and \code{full_comp_cds} take only the values 0-2,
#' so a single site rescales to 0.5 and, at their former combined weight of
#' 0.27, dominated the score; that site is normally the intended on-target
#' match rather than an off-target liability.
#'
#' @return A named numeric vector of length 4.
#' @examples
#' default_weights()
#' sum(default_weights())  # 1
#' @export
default_weights <- function() {
  c(
    cds_count           = 0.25,
    seed_region_count   = 0.25,
    critical_gene_count = 0.25,
    high_risk_critical  = 0.25
  )
}
