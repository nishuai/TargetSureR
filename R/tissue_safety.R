#' Summarise tissue-level expression risk of off-target genes
#'
#' For every off-target \code{gene_name} in an enriched table, computes a
#' per-tissue mean TPM and a bio-attention score that combines tissue
#' expression, maximum observed \code{risk_score}, CDS full-complementarity
#' hit rate, and critical-gene membership. The result is a gene-level tibble
#' ordered by decreasing \code{bio_attention}.
#'
#' @param enriched A data frame from \code{\link{enrich_annotations}} whose
#'   right-hand columns are tissue TPMs.
#' @param primary_tissue Column name of the tissue treated as the intended
#'   target tissue (e.g. \code{"Nerve_Tibial"}). Used to emphasise on-target
#'   expression vs. systemic expression in the bio-attention score.
#' @param weights Named numeric vector; components
#'   \code{expr, risk, fc_cds, key}. Default
#'   \code{c(0.35, 0.25, 0.25, 0.15)}.
#' @param log_expr Logical. Compress \code{primary_tpm} with \code{log1p}
#'   before rescaling \code{score_expr}. Default \code{TRUE}. TPM spans four
#'   to five orders of magnitude, so a plain min-max is dominated by the
#'   single highest-expressed gene and pushes every other gene to near zero.
#' @param exclude_genes Regular expression matched against
#'   \code{gene_name}; matching genes are dropped before scoring. Defaults to
#'   mitochondrial genes (\code{"^(MT-|MT[ARN]?[0-9])"}), which top the TPM
#'   range in every tissue but are not plausible siRNA liabilities. Pass
#'   \code{NULL} to keep every gene.
#' @param on_target_genes Character vector of intended target gene symbols
#'   (matched case-insensitively), dropped before scoring. Default
#'   \code{NULL}, which keeps them. A target gene enters this table via its
#'   own perfect CDS match - on-target success, not an off-target liability -
#'   and then holds the maximum of both \code{max_risk_score} and
#'   \code{prop_fc_cds}. Since those are rescaled by min-max, one such gene
#'   flattens every genuine off-target to a single value across half the
#'   score's weight. Supply the query's target gene(s) to score off-targets
#'   only; leave \code{NULL} to see the on-target hit as an internal control.
#'
#' @return A tibble with one row per \code{gene_name} and columns for the
#'   component scores and the final \code{bio_attention}. The logical column
#'   \code{expression_missing} marks genes absent from the expression matrix,
#'   whose \code{primary_tpm} is \code{NA} and which are scored as zero
#'   expression for ranking purposes only.
#' @examples
#' \dontrun{
#' tissue_safety(enriched, primary_tissue = "Nerve_Tibial")
#' }
#' @export
tissue_safety <- function(enriched,
                          primary_tissue = NULL,
                          weights = c(expr = 0.35, risk = 0.25,
                                      fc_cds = 0.25, key = 0.15),
                          log_expr = TRUE,
                          exclude_genes = "^(MT-|MT[ARN]?[0-9])",
                          on_target_genes = NULL) {
  abort_missing_cols(enriched,
                     c("gene_name", "region", "match_type", "risk_score",
                       "AE_gene", "cancer_gene_v2", "immune_gene_immport"),
                     "tissue_safety")

  tissues <- setdiff(names(enriched), enriched_columns())
  if (length(tissues) == 0) {
    stop("No tissue TPM columns found. Did you call enrich_annotations()?",
         call. = FALSE)
  }
  if (is.null(primary_tissue)) {
    primary_tissue <- tissues[1]
  }
  if (!primary_tissue %in% tissues) {
    stop(sprintf("primary_tissue '%s' not in expression matrix",
                 primary_tissue), call. = FALSE)
  }

  x <- tibble::as_tibble(enriched)

  # Drop genes matched by exclude_genes before scoring. Mitochondrial
  # transcripts sit at the top of every tissue's TPM range yet are not
  # plausible siRNA liabilities, and because score_expr is rescaled across
  # whatever genes remain, leaving one in flattens every other gene.
  if (!is.null(exclude_genes) && nzchar(exclude_genes)) {
    drop <- grepl(exclude_genes, x$gene_name, ignore.case = TRUE)
    if (any(drop)) x <- x[!drop, , drop = FALSE]
    if (nrow(x) == 0)
      stop("exclude_genes removed every row", call. = FALSE)
  }

  # An intended target gene reaches this table through its own perfect CDS
  # match, which is on-target success rather than an off-target liability. It
  # then owns the maximum of both max_risk_score and prop_fc_cds, and because
  # those are rescaled by min-max, every genuine off-target collapses to a
  # single value on half the score's weight. Removing it restores that spread.
  if (!is.null(on_target_genes) && length(on_target_genes) > 0) {
    keep <- !(safe_toupper(x$gene_name) %in% safe_toupper(on_target_genes))
    if (!all(keep)) x <- x[keep, , drop = FALSE]
    if (nrow(x) == 0)
      stop("on_target_genes removed every row", call. = FALSE)
  }

  gene_mat <- dplyr::summarise(
    dplyr::group_by(x, .data$gene_name),
    n_hits         = dplyr::n(),
    primary_tpm    = mean(.data[[primary_tissue]], na.rm = TRUE),
    max_risk_score = max(.data$risk_score, na.rm = TRUE),
    prop_fc_cds    = mean(.data$match_type == "full_complementarity"
                          & .data$region == "CDS", na.rm = TRUE),
    is_cancer      = any(.data$cancer_gene_v2, na.rm = TRUE),
    is_ae          = any(.data$AE_gene, na.rm = TRUE),
    is_immune      = any(.data$immune_gene_immport, na.rm = TRUE),
    .groups = "drop"
  )
  gene_mat$is_key_gene <- gene_mat$is_cancer | gene_mat$is_ae | gene_mat$is_immune

  w <- weights / sum(weights)
  # TPM spans 4-5 orders of magnitude (mitochondrial and ribosomal genes reach
  # 10^4 while a typical coding gene sits at 10^1), so a plain min-max would be
  # dominated by the single largest gene and compress everything else to ~0.
  # Compress on log1p first, then rescale.
  # An off-target gene absent from the expression matrix carries NA rather
  # than a measured zero. It is flagged, then treated as zero for ranking
  # only, so that "expression unknown" stays distinguishable in the output.
  gene_mat$expression_missing <- is.na(gene_mat$primary_tpm)
  tpm <- tidyr::replace_na(gene_mat$primary_tpm, 0)
  gene_mat$score_expr   <- if (isTRUE(log_expr)) normalize_01(log1p(tpm))
                           else normalize_01(tpm)
  gene_mat$score_risk   <- normalize_01(gene_mat$max_risk_score)
  gene_mat$score_fc_cds <- normalize_01(gene_mat$prop_fc_cds)
  gene_mat$score_key    <- as.numeric(gene_mat$is_key_gene)
  gene_mat$bio_attention <- w["expr"]   * gene_mat$score_expr +
                            w["risk"]   * gene_mat$score_risk +
                            w["fc_cds"] * gene_mat$score_fc_cds +
                            w["key"]    * gene_mat$score_key
  gene_mat <- gene_mat[order(-gene_mat$bio_attention), , drop = FALSE]
  tibble::as_tibble(gene_mat)
}
