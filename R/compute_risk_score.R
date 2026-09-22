#' Compute the rule-based risk score for off-target sites
#'
#' Assigns an integer risk score in \code{[1, 10]} based on the region of the
#' off-target hit and the type of match. Designed to work in the absence of
#' thermodynamic (\code{mfe}) information.
#'
#' @section Formula:
#' \preformatted{
#'   base  <- 5 if region == "CDS"
#'            3 if region == "3'UTR"
#'            1 otherwise
#'   bonus <- 2 if match_type == "full_complementarity"
#'            0 otherwise
#'   risk_score <- pmin(base + bonus, 10)
#' }
#'
#' @param region Character vector of region labels
#'   (\code{"CDS"}, \code{"3'UTR"}, \code{"5'UTR"}, \code{"ncRNA"}, ...).
#' @param match_type Character vector of match types
#'   (\code{"full_complementarity"}, \code{"partial_match"}, \code{"unknown"}).
#'
#' @return Integer vector, same length as \code{region}.
#' @examples
#' compute_risk_score(
#'   region     = c("CDS", "3'UTR", "5'UTR"),
#'   match_type = c("full_complementarity", "partial_match", "unknown")
#' )
#' @export
compute_risk_score <- function(region, match_type) {
  stopifnot(length(region) == length(match_type))

  region <- as.character(region)
  match_type <- as.character(match_type)

  base <- ifelse(region == "CDS",   5L,
          ifelse(region == "3'UTR", 3L, 1L))

  bonus <- ifelse(match_type == "full_complementarity", 2L, 0L)

  pmin(base + bonus, 10L)
}
