#' Rank siRNAs by a composite off-target risk score
#'
#' Computes seven per-siRNA metrics from an enriched off-target table, scales
#' each metric to \code{[0, 1]}, applies \code{\link{default_weights}} (or
#' user-supplied weights) to produce a \code{composite_score}, then ranks and
#' bins siRNAs into \code{Low / Medium / High / Critical} tiers.
#'
#' Lower composite score is better (fewer / less dangerous off-targets).
#'
#' @section Normalization references:
#' Three modes via \code{reference}:
#' \describe{
#'   \item{\code{"builtin"}}{(default) Normalize against the package's
#'     built-in \code{\link{reference_set}} of 94 known siRNAs. Single
#'     siRNA gets a meaningful absolute score (its position relative to
#'     the reference cohort distribution).}
#'   \item{\code{"self"}}{Normalize each metric across the siRNAs in
#'     \code{enriched}. Backward-compatible legacy behavior; works well
#'     when ranking a batch of >= 5 siRNAs against each other. Single-
#'     siRNA input gets all-zero scores.}
#'   \item{\code{"custom"}}{Use \code{reference_set}, a user-supplied
#'     reference object with the same shape as \code{\link{reference_set}}.}
#' }
#'
#' @param enriched A data frame from \code{\link{enrich_annotations}}.
#' @param target_sirnas Character vector of siRNA names to highlight as
#'   \code{is_target = TRUE}. Optional.
#' @param weights Named numeric vector matching
#'   \code{names(default_weights())}. Re-normalised to sum 1 internally.
#' @param tier_breaks Numeric breakpoints for tier assignment on the composite
#'   score. Default \code{c(0.25, 0.50, 0.75)}.
#' @param reference One of \code{"builtin"} (default), \code{"self"},
#'   \code{"custom"}. See Details.
#' @param reference_set When \code{reference = "custom"}, the user-supplied
#'   reference object (a list with at least \code{tier_thresholds}).
#'   Defaults to the built-in \code{\link{reference_set}}.
#'
#' @return A tibble with per-siRNA metrics, \code{composite_score},
#'   \code{rank}, \code{risk_tier}, and \code{is_target}.
#' @examples
#' \dontrun{
#' # Default: compare to built-in 94-siRNA reference cohort
#' rank_sirna(enriched)
#'
#' # Self-normalize within the input batch (legacy behavior)
#' rank_sirna(enriched, reference = "self")
#' }
#' @export
rank_sirna <- function(enriched,
                       target_sirnas = character(),
                       weights = default_weights(),
                       tier_breaks = c(0.25, 0.50, 0.75),
                       reference = c("builtin", "self", "custom"),
                       reference_set = NULL) {
  reference <- match.arg(reference)
  abort_missing_cols(enriched,
                     c("siRNA_name", "region", "match_type", "risk_score",
                       "cancer_gene_v2", "AE_gene", "immune_gene_immport"),
                     "rank_sirna")

  x <- tibble::as_tibble(enriched)
  x$is_critical <- x$cancer_gene_v2 | x$AE_gene | x$immune_gene_immport
  x$is_high_risk <- x$risk_score >= 5
  x$is_full_comp <- x$match_type == "full_complementarity"
  x$is_cds       <- x$region == "CDS"
  x$is_seed_utr  <- x$region == "3'UTR" & x$match_type == "partial_match"

  metrics <- dplyr::summarise(
    dplyr::group_by(x, .data$siRNA_name),
    high_risk_critical  = sum(.data$is_high_risk & .data$is_critical, na.rm = TRUE),
    critical_gene_count = sum(.data$is_critical, na.rm = TRUE),
    full_comp_cds       = sum(.data$is_full_comp & .data$is_cds, na.rm = TRUE),
    full_comp_count     = sum(.data$is_full_comp, na.rm = TRUE),
    cds_count           = sum(.data$is_cds, na.rm = TRUE),
    high_risk_count     = sum(.data$is_high_risk, na.rm = TRUE),
    seed_region_count   = sum(.data$is_seed_utr, na.rm = TRUE),
    total_offtargets    = dplyr::n(),
    .groups = "drop"
  )

  dims <- names(default_weights())
  if (!all(dims %in% names(weights))) {
    stop("weights must contain: ", paste(dims, collapse = ", "), call. = FALSE)
  }
  w <- weights[dims] / sum(weights[dims])

  normed <- metrics
  if (reference == "self") {
    for (d in dims) normed[[d]] <- normalize_01(metrics[[d]])
  } else {
    ref <- if (reference == "builtin") .load_pkg_data("reference_set") else
           reference_set
    if (is.null(ref) || !"tier_thresholds" %in% names(ref)) {
      stop("reference_set must be a list containing 'tier_thresholds'",
           call. = FALSE)
    }
    th <- ref$tier_thresholds
    for (d in dims) {
      mn <- th$min[th$dimension == d]
      mx <- th$max[th$dimension == d]
      if (length(mn) == 0 || mx == mn) {
        normed[[d]] <- 0
      } else {
        v <- (metrics[[d]] - mn) / (mx - mn)
        v[v < 0] <- 0; v[v > 1] <- 1
        normed[[d]] <- v
      }
    }
  }

  composite <- as.matrix(normed[, dims, drop = FALSE]) %*% as.numeric(w)
  metrics$composite_score <- as.numeric(composite)

  metrics <- metrics[order(metrics$composite_score), , drop = FALSE]
  metrics$rank <- seq_len(nrow(metrics))
  metrics$risk_tier <- cut(metrics$composite_score,
                           breaks = c(-Inf, tier_breaks, Inf),
                           labels = c("Low", "Medium", "High", "Critical"),
                           right = TRUE)
  metrics$is_target <- metrics$siRNA_name %in% target_sirnas
  attr(metrics, "reference") <- reference
  tibble::as_tibble(metrics)
}
