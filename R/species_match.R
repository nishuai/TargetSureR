#' Rank candidate animal models by cross-species conservation
#'
#' For each siRNA in a core off-target table, estimates how well each
#' candidate species preserves the human off-target profile. Each unique
#' off-target gene contributes once, irrespective of the number of predicted
#' sites or feature-class transcripts assigned to that gene. Scoring is
#' restricted to \dQuote{critical} genes by default.
#'
#' Requires a pre-computed ortholog table; see \code{data-raw/} for the
#' builder template. Layer 1 only (no sequence-level alignment).
#'
#' @param sites A data frame with \code{\link{core_columns}} or the output of
#'   \code{\link{enrich_annotations}}. If the three gene flags are missing,
#'   they will be computed via \code{\link{load_gene_lists}}.
#' @param candidates Character vector of species keys to evaluate. Defaults
#'   to the six commonly used preclinical model species:
#'   \code{"mouse"}, \code{"rat"}, \code{"cyno"} (\emph{Macaca fascicularis}),
#'   \code{"rhesus"} (\emph{Macaca mulatta}),
#'   \code{"rabbit"} (\emph{Oryctolagus cuniculus}),
#'   \code{"dog"} (\emph{Canis lupus familiaris}).
#' @param ortholog_table A data frame keyed by \code{human_symbol}, with
#'   per-species \code{\{species\}_identity_pct} and
#'   \code{\{species\}_ortholog_type} columns. If \code{NULL}, the
#'   \strong{simulated} 20-gene \code{\link{ortholog_table_example}} is used:
#'   genes outside those 20 resolve to no ortholog and score zero, so a real
#'   result is dragged towards zero regardless of true orthology. Supply a
#'   BioMart download covering your own off-target genes for analysis; see
#'   \code{vignette("01_quickstart")}, step 5.
#' @param critical_only Logical. If \code{TRUE} (default), only critical
#'   genes contribute to the species score.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{species_score}}{siRNA x species score matrix (tibble).}
#'     \item{\code{site_detail}}{Long tibble with per-site, per-species
#'       conservation.}
#'     \item{\code{gene_detail}}{Long tibble with one row per unique
#'       siRNA, gene and species; this is the table used for scoring.}
#'     \item{\code{recommendation}}{Character vector with textual
#'       recommendations.}
#'   }
#' @examples
#' \dontrun{
#' species_match(sites,
#'               candidates = c("mouse", "rat", "cyno",
#'                              "rhesus", "rabbit", "dog"))
#' }
#' @export
species_match <- function(sites,
                          candidates = c("mouse", "rat", "cyno",
                                         "rhesus", "rabbit", "dog"),
                          ortholog_table = NULL,
                          critical_only = TRUE) {
  abort_missing_cols(sites, c("siRNA_name", "gene_name"),
                     "species_match")

  x <- tibble::as_tibble(sites)

  if (!all(c("AE_gene", "cancer_gene_v2", "immune_gene_immport") %in% names(x))) {
    gl <- load_gene_lists()
    gu <- safe_toupper(x$gene_name)
    x$AE_gene             <- gu %in% safe_toupper(gl$AE_gene)
    x$cancer_gene_v2      <- gu %in% safe_toupper(gl$cancer_gene_v2)
    x$immune_gene_immport <- gu %in% safe_toupper(gl$immune_gene_immport)
  }
  x$is_critical <- x$cancer_gene_v2 | x$AE_gene | x$immune_gene_immport

  if (critical_only) x <- x[x$is_critical, , drop = FALSE]
  if (nrow(x) == 0) {
    return(list(species_score = tibble::tibble(siRNA_name = character()),
                site_detail    = tibble::tibble(),
                gene_detail    = tibble::tibble(),
                recommendation = character()))
  }

  if (is.null(ortholog_table)) {
    ortholog_table <- .load_pkg_data("ortholog_table_example")
  }

  detail <- .conservation_per_site(x, ortholog_table, candidates)

  # Weight each unique off-target gene once. Counting rows here would
  # implicitly give more weight to genes represented by several sites or
  # feature-class transcripts.
  gene_detail <- detail |>
    dplyr::distinct(.data$siRNA_name, .data$gene_name, .data$species,
                    .keep_all = TRUE)
  agg <- dplyr::summarise(
    dplyr::group_by(gene_detail, .data$siRNA_name, .data$species),
    species_score = mean(.data$conservation),
    n_genes       = dplyr::n(),
    .groups = "drop"
  )
  mat <- tidyr::pivot_wider(agg, id_cols = "siRNA_name",
                            names_from = "species",
                            values_from = "species_score")
  rec <- .recommend(agg)
  list(species_score  = tibble::as_tibble(mat),
       site_detail    = tibble::as_tibble(detail),
       gene_detail    = tibble::as_tibble(gene_detail),
       recommendation = rec)
}

.conservation_per_site <- function(x, ot, species) {
  ot$human_symbol <- safe_toupper(ot$human_symbol)
  gu <- safe_toupper(x$gene_name)
  idx <- match(gu, ot$human_symbol)

  out <- list()
  for (sp in species) {
    id_col    <- sprintf("%s_identity_pct", sp)
    type_col  <- sprintf("%s_ortholog_type", sp)
    if (!all(c(id_col, type_col) %in% names(ot))) {
      warning(sprintf("ortholog_table missing columns for species '%s'; skipped", sp),
              call. = FALSE)
      next
    }
    ident <- ot[[id_col]][idx]
    otype <- ot[[type_col]][idx]
    cons <- .identity_to_conservation(ident, otype)
    out[[sp]] <- tibble::tibble(
      siRNA_name    = x$siRNA_name,
      gene_name     = x$gene_name,
      species       = sp,
      identity_pct  = ident,
      ortholog_type = otype,
      conservation  = cons
    )
  }
  dplyr::bind_rows(out)
}

.identity_to_conservation <- function(identity_pct, ortholog_type) {
  ortholog_type <- as.character(ortholog_type)
  cons <- rep(0, length(identity_pct))
  no   <- is.na(ortholog_type) | ortholog_type == "none"
  cons[!no & identity_pct >= 90]                                        <- 1.0
  cons[!no & identity_pct >= 80 & identity_pct < 90]                    <- 0.6
  cons[!no & identity_pct <  80 & !is.na(identity_pct)]                 <- 0.3
  cons[no] <- 0
  cons
}

.recommend <- function(agg) {
  best <- dplyr::slice_max(dplyr::group_by(agg, .data$siRNA_name),
                            order_by = .data$species_score,
                            n = 1, with_ties = FALSE)
  vapply(seq_len(nrow(best)), function(i) {
    sprintf("%s: highest conservation in %s (score %.2f, %d critical genes)",
            best$siRNA_name[i], best$species[i],
            best$species_score[i], best$n_genes[i])
  }, character(1))
}
