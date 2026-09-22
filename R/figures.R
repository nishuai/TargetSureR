#' Plotting helpers for TargetSureR
#'
#' Six figure helpers that consume the standard pipeline outputs and emit
#' \pkg{ggplot2} objects. All return ggplots; they do not write to disk.
#' Use \code{\link{generate_figures}} for a one-call bundle that writes
#' both PDF and PNG to a directory.
#'
#' Heavy plotting dependencies (\pkg{ggplot2}, \pkg{tidyr}, \pkg{ggrepel},
#' \pkg{RColorBrewer}, \pkg{scales}) are declared under Suggests; calling
#' any of these helpers without them throws an actionable error.
#'
#' @name figures
NULL

.theme_tsr <- function() {
  ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "#555555"),
      panel.grid.major.y = ggplot2::element_line(color = "#eeeeee"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA))
}

.require_plot_pkgs <- function() {
  for (p in c("ggplot2","tidyr","RColorBrewer","scales"))
    require_pkg(p, "plotting")
}

# Qualitative palette for user siRNAs, returned NAMED so that
# scale_*_manual() maps colours by siRNA name rather than by position.
.user_palette <- function(names_vec) {
  nm <- sort(unique(as.character(names_vec)))
  pal <- RColorBrewer::brewer.pal(max(3L, min(9L, length(nm))), "Set1")
  pal <- rep(pal, length.out = length(nm))
  stats::setNames(pal, nm)
}

#' @describeIn figures Bar chart of composite_score per siRNA, colored by tier.
#'   When \code{reference} is supplied, switches to a violin layout: the
#'   reference-cohort distribution is shown as a violin with all reference
#'   siRNAs jittered as gray points, the user siRNA(s) is highlighted as
#'   colored diamonds, and the 25/50/75 percentile of the reference cohort
#'   is annotated with dashed lines.
#' @param ranking Output of \code{\link{rank_sirna}}.
#' @param reference Optional. \code{"builtin"} for the package built-in
#'   reference, or a list with the same shape as \code{\link{reference_set}}.
#'   Default \code{NULL} (no reference overlay; classic bar chart).
#' @export
plot_ranking_bar <- function(ranking, reference = NULL) {
  .require_plot_pkgs()
  # Sequential light-to-dark, so that increasing burden reads as increasing
  # intensity. A red-green diverging scheme was previously used here, which both
  # inverted the ordering (Low took the red end) and is unreadable under
  # red-green colour vision deficiency.
  tier_pal <- RColorBrewer::brewer.pal(4, "YlOrRd")
  names(tier_pal) <- c("Low", "Medium", "High", "Critical")

  # Coerce to a factor carrying all four tiers, so that the legend shows the
  # full scale even when the supplied siRNAs occupy only some of it. A
  # character column would leave `drop = FALSE` nothing to preserve.
  ranking$risk_tier <- factor(as.character(ranking$risk_tier),
                              levels = names(tier_pal))

  if (is.null(reference)) {
    ranking$siRNA_name <- factor(ranking$siRNA_name,
      levels = ranking$siRNA_name[order(ranking$composite_score)])
    return(ggplot2::ggplot(ranking,
      ggplot2::aes(x = .data$siRNA_name, y = .data$composite_score,
                   fill = .data$risk_tier)) +
      ggplot2::geom_col(width = 0.7, show.legend = TRUE) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$composite_score)),
                          vjust = -0.4, size = 3.5) +
      # Keep a key for every tier, so the legend shows the full scale even
      # when the supplied siRNAs occupy only part of it.
      ggplot2::scale_fill_manual(
        values = tier_pal, drop = FALSE,
        guide = ggplot2::guide_legend(
          override.aes = list(colour = "black", linewidth = 0.3))) +
      ggplot2::labs(title = "siRNA off-target risk ranking",
         x = NULL, y = "Composite risk score", fill = "Risk tier") +
      .theme_tsr())
  }

  # ---- reference overlay: violin + jitter + highlight + quantile lines
  ref <- if (identical(reference, "builtin"))
    .load_pkg_data("reference_set") else reference
  if (!"metrics" %in% names(ref)) {
    stop("reference must contain 'metrics' (with composite_score)",
         call. = FALSE)
  }
  ref_scores <- ref$metrics$composite_score
  q <- stats::quantile(ref_scores, c(0.25, 0.50, 0.75), na.rm = TRUE)

  ref_df  <- data.frame(group = "Reference cohort", score = ref_scores)
  user_df <- data.frame(
    group = "Reference cohort",                 # plot user dots on ref column
    score = ranking$composite_score,
    siRNA_name = ranking$siRNA_name,
    risk_tier = ranking$risk_tier)

  ggplot2::ggplot() +
    ggplot2::geom_violin(
      data = ref_df,
      ggplot2::aes(x = .data$group, y = .data$score),
      fill = "#dddddd", color = "#888888",
      width = 0.55, alpha = 0.75, trim = FALSE) +
    ggplot2::geom_jitter(
      data = ref_df,
      ggplot2::aes(x = .data$group, y = .data$score),
      width = 0.11, alpha = 0.45, size = 0.9, color = "#666666") +
    ggplot2::geom_hline(yintercept = q, linetype = "dashed",
                         color = "#444444", linewidth = 0.4) +
    ggplot2::annotate("text",
      x = 1.45, y = q,
      label = sprintf("Q%d = %.3f", c(25, 50, 75), q),
      size = 3.0, hjust = 0, color = "#444444") +
    ggplot2::geom_point(
      data = user_df,
      ggplot2::aes(x = .data$group, y = .data$score, fill = .data$risk_tier),
      shape = 23, size = 5, color = "black", stroke = 0.6,
      show.legend = TRUE) +
    ggrepel::geom_text_repel(
      data = user_df,
      ggplot2::aes(x = .data$group, y = .data$score,
                   label = sprintf("%s (%.3f)", .data$siRNA_name, .data$score)),
      size = 3.4, fontface = "bold",
      nudge_x = -0.35, direction = "y", segment.color = "#888888",
      min.segment.length = 0) +
    ggplot2::scale_fill_manual(values = tier_pal, drop = FALSE,
                                name = "Risk tier",
                                # draw a key for every tier, so the legend
                                # shows the full scale even where no siRNA
                                # falls in a tier
                                guide = ggplot2::guide_legend(
                                  override.aes = list(shape = 23, size = 4,
                                                      colour = "black",
                                                      alpha = 1))) +
    # keep the single violin column at a normal width instead of letting it
    # stretch across the whole panel
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(1.1, 1.1))) +
    ggplot2::labs(title = "siRNA off-target risk ranking vs. reference cohort",
       x = NULL, y = "Composite risk score (lower = safer)") +
    .theme_tsr()
}

#' @describeIn figures Stacked bar of off-target site counts by risk_score.
#'   With \code{reference}, draws all 94 reference siRNAs side-by-side
#'   together with the user's siRNA(s); the user's bar(s) get a thick
#'   black outline and a red label so they're easy to locate among the
#'   sorted cohort.
#' @param enriched Output of \code{\link{enrich_annotations}}.
#' @export
plot_risk_distribution <- function(enriched, reference = NULL) {
  .require_plot_pkgs()
  d <- dplyr::count(enriched, .data$siRNA_name, .data$risk_score)

  if (is.null(reference)) {
    return(ggplot2::ggplot(d,
      ggplot2::aes(x = .data$siRNA_name, y = .data$n,
                   fill = factor(.data$risk_score))) +
      ggplot2::geom_col(position = ggplot2::position_stack(reverse = TRUE)) +
      ggplot2::scale_fill_brewer(palette = "OrRd", direction = 1,
        name = "risk_score\n(CDS=5, +full=+2)") +
      ggplot2::labs(title = "Off-target site count by risk_score",
        x = NULL, y = "Off-target sites") +
      .theme_tsr())
  }

  ref <- if (identical(reference, "builtin"))
    .load_pkg_data("reference_set") else reference
  if (!"risk_distribution" %in% names(ref)) {
    stop("reference must contain 'risk_distribution'", call. = FALSE)
  }

  ref_d  <- ref$risk_distribution
  ref_d$is_user  <- FALSE
  user_d <- d
  user_d$is_user <- TRUE

  combined <- dplyr::bind_rows(
    ref_d[,  c("siRNA_name","risk_score","n","is_user")],
    user_d[, c("siRNA_name","risk_score","n","is_user")])

  totals <- dplyr::summarise(
    dplyr::group_by(combined, .data$siRNA_name),
    total = sum(.data$n), is_user = any(.data$is_user),
    .groups = "drop")
  totals <- totals[order(totals$total), ]
  combined$siRNA_name <- factor(combined$siRNA_name,
                                 levels = totals$siRNA_name)
  combined$is_user <- factor(combined$is_user, levels = c(FALSE, TRUE))

  user_levels <- as.character(totals$siRNA_name[totals$is_user])

  ggplot2::ggplot(combined,
    ggplot2::aes(x = .data$siRNA_name, y = .data$n,
                 fill = factor(.data$risk_score))) +
    ggplot2::geom_col(
      ggplot2::aes(color = .data$is_user, linewidth = .data$is_user),
      position = ggplot2::position_stack(reverse = TRUE)) +
    ggplot2::scale_fill_brewer(palette = "OrRd", direction = 1,
       name = "risk_score\n(CDS=5, +full=+2)") +
    ggplot2::scale_color_manual(values = c(`FALSE` = NA, `TRUE` = "black"),
                                 guide = "none") +
    ggplot2::scale_linewidth_manual(values = c(`FALSE` = 0, `TRUE` = 0.7),
                                     guide = "none") +
    ggplot2::labs(title = "Off-target site count by risk_score: cohort view",
       x = "siRNA (sorted by total off-targets)",
       y = "Off-target sites") +
    .theme_tsr() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90, hjust = 1, vjust = 0.5, size = 4,
        color = ifelse(levels(combined$siRNA_name) %in% user_levels,
                       "#d62728", "#999999"),
        face  = ifelse(levels(combined$siRNA_name) %in% user_levels,
                       "bold", "plain")))
}

#' @describeIn figures Scatter of top genes by bio_attention components.
#' @param tissue Output of \code{\link{tissue_safety}}.
#' @param top_n Number of top genes to label (default 30).
#' @param y Which variable to put on the y axis. \code{"risk"} uses
#'   \code{score_risk}, which is derived from \code{risk_score} and therefore
#'   takes at most four distinct values - when the shortlist shares one
#'   region/match combination it collapses to a single row of points.
#'   \code{"hits"} (default) uses \code{log1p(n_hits)}, which is continuous
#'   enough to separate genes in that situation.
#' @param x Which variable to put on the x axis. \code{"tpm"} (default) uses
#'   the gene's median TPM in the primary tissue on a log scale, in reported
#'   units; \code{"score"} uses the rescaled expression component instead.
#' @export
plot_gene_attention <- function(tissue, top_n = 30, y = c("hits", "risk"),
                                x = c("tpm", "score")) {
  y <- match.arg(y); x <- match.arg(x)
  .require_plot_pkgs()
  require_pkg("ggrepel", "plotting")
  scatter_pal <- RColorBrewer::brewer.pal(3, "Set1")[c(2, 1)]
  d <- utils::head(tissue, top_n)
  d$.yval <- if (y == "risk") d$score_risk else d$n_hits
  y_lab <- if (y == "risk") "Normalized max risk_score"
           else "Off-target sites on this gene"
  d$.xval <- if (x == "score") d$score_expr else d$primary_tpm
  x_lab <- if (x == "score") "Normalized expression (primary tissue)"
           else "Median expression in primary tissue (TPM)"
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$.xval, y = .data$.yval,
                                   size = .data$n_hits, color = .data$is_key_gene)) +
    ggplot2::geom_point(alpha = 0.75) +
    ggrepel::geom_text_repel(ggplot2::aes(label = .data$gene_name),
                              size = 3.2, max.overlaps = top_n) +
    ggplot2::scale_color_manual(
      values = c(`FALSE` = scatter_pal[1], `TRUE` = scatter_pal[2]),
      labels = c("Non-critical","Critical"), name = "Gene category") +
    ggplot2::scale_size_continuous(range = c(2, 10), name = "# hits") +
    ggplot2::labs(title = "Off-target gene attention",
       x = x_lab, y = y_lab)
  if (x == "tpm")
    p <- p + ggplot2::scale_x_log10(labels = scales::label_number())
  if (y == "hits")
    p <- p + ggplot2::scale_y_log10(breaks = c(1, 2, 3, 5, 10, 20))
  p + .theme_tsr()
}

#' @describeIn figures Histogram of the number of mismatches between a guide
#'   and each site it binds, overlaid per siRNA so guides can be compared. Reads the CORE table (the output of
#'   \code{\link{annotate_sites}}), not the enriched one, because
#'   \code{\link{enrich_annotations}} drops \code{n_mismatch}. Sites with
#'   \code{n_mismatch == 0} are perfect complements (typically the intended
#'   on-target) and are marked with a dashed rule.
#' @param sites Output of \code{\link{annotate_sites}} (the core table),
#'   which must retain the \code{n_mismatch} column.
#' @export
plot_mismatch_distribution <- function(sites) {
  .require_plot_pkgs()
  abort_missing_cols(sites, c("siRNA_name", "n_mismatch"),
                     "plot_mismatch_distribution")
  d <- tibble::as_tibble(sites)
  d <- d[!is.na(d$n_mismatch), ]
  if (nrow(d) == 0) stop("no sites with a non-missing n_mismatch",
                          call. = FALSE)
  d$n_mismatch <- as.integer(d$n_mismatch)
  pal <- .user_palette(unique(d$siRNA_name))

  ggplot2::ggplot(d, ggplot2::aes(x = .data$n_mismatch,
                                   fill = .data$siRNA_name)) +
    ggplot2::geom_histogram(binwidth = 1, alpha = 0.55,
                             position = "identity", color = "white",
                             linewidth = 0.2) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                         color = "#444444", linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = pal, name = "siRNA") +
    ggplot2::scale_x_continuous(breaks = scales::breaks_width(2)) +
    ggplot2::labs(
      title = "Mismatch distribution: guide vs. bound sites",
      x = "Mismatches vs. reverse-complement of the guide",
      y = "Off-target sites") +
    .theme_tsr()
}

#' @describeIn figures Cumulative view of the same mismatch counts as
#'   \code{\link{plot_mismatch_distribution}}. The histogram is dominated by the
#'   high-mismatch mode, which hides the low-mismatch tail that actually
#'   drives off-target risk; the cumulative curve reads that tail directly
#'   ("what fraction of sites sit at <= k mismatches"). A lower, flatter
#'   left-hand region means a more specific guide.
#'   The y axis is zoomed to the data actually inside the x window: the full
#'   cumulative curve reaches 100% only at high mismatch counts, so a fixed
#'   0-100% axis would squeeze the informative band into ~1% of the plot
#'   height.
#' @param max_mismatch Integer. Right-hand limit of the x axis, so the
#'   informative left tail is not compressed. Default \code{8}.
#' @param y Either \code{"share"} (default, cumulative percentage of that
#'   guide's own sites - comparable across guides with different totals) or
#'   \code{"count"} (cumulative number of sites - shows absolute burden).
#' @export
plot_mismatch_cumulative <- function(sites, max_mismatch = 8L,
                                  y = c("share", "count")) {
  y <- match.arg(y)
  .require_plot_pkgs()
  abort_missing_cols(sites, c("siRNA_name", "n_mismatch"),
                     "plot_mismatch_cumulative")
  d <- tibble::as_tibble(sites)
  d <- d[!is.na(d$n_mismatch), ]
  if (nrow(d) == 0) stop("no sites with a non-missing n_mismatch",
                          call. = FALSE)
  d$n_mismatch <- as.integer(d$n_mismatch)

  # cumulative count and share per siRNA, evaluated on a common grid
  grid <- seq(0L, max(d$n_mismatch), by = 1L)
  cum <- do.call(rbind, lapply(split(d, d$siRNA_name), function(g) {
    data.frame(siRNA_name = g$siRNA_name[1], n_mismatch = grid,
               cum_n = vapply(grid, function(k) sum(g$n_mismatch <= k),
                              numeric(1)),
               share = vapply(grid, function(k) mean(g$n_mismatch <= k),
                              numeric(1)),
               stringsAsFactors = FALSE)
  }))
  lab <- do.call(rbind, lapply(split(cum, cum$siRNA_name), function(g) {
    g[g$n_mismatch == min(max_mismatch, max(g$n_mismatch)), ]
  }))
  pal <- .user_palette(unique(d$siRNA_name))

  cum$yval <- if (y == "share") cum$share else cum$cum_n
  lab$yval  <- if (y == "share") lab$share  else lab$cum_n
  # Zoom y to what is inside the x window. The full curve only reaches 100%
  # at high k, so a 0-100% axis would flatten the informative band.
  y_hi  <- max(cum$yval[cum$n_mismatch <= max_mismatch], na.rm = TRUE)
  y_pad <- if (y == "share") 1.35 else 1.30

  p <- ggplot2::ggplot(cum, ggplot2::aes(x = .data$n_mismatch, y = .data$yval,
                                         color = .data$siRNA_name)) +
    ggplot2::geom_step(linewidth = 0.9, direction = "hv") +
    ggplot2::geom_point(size = 1.8) +
    ggrepel::geom_text_repel(
      data = lab,
      # The siRNA name is already in the legend, so label only the value and
      # keep it inside the panel: anchored left of the end point, right
      # aligned, repelled vertically only.
      ggplot2::aes(label = if (y == "share")
                     sprintf("%d (%.2f%%)", .data$cum_n, 100 * .data$share)
                   else
                     sprintf("%d", .data$cum_n)),
      # sit just ABOVE the end point so the guide's own step line cannot
      # run through the text
      size = 3.1, fontface = "bold", direction = "y", hjust = 0.5,
      nudge_y = y_hi * 0.07, box.padding = 0.4, point.padding = 0.3,
      segment.color = "#aaaaaa", segment.size = 0.3,
      min.segment.length = 0, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = pal, name = "siRNA") +
    # right margin on x so the last step and its label are not flush against
    # the panel edge; y headroom so the topmost label clears the frame
    ggplot2::coord_cartesian(xlim = c(0, max_mismatch + 0.45),
                             ylim = c(0, y_hi * y_pad)) +
    ggplot2::labs(
      title = "Mismatch distribution, cumulative view",
      x = "Mismatches vs. reverse-complement of the guide (k)",
      y = if (y == "share") "Cumulative share of that guide's sites"
          else "Cumulative off-target sites")
  if (y == "share")
    p <- p + ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 0.1))
  p + .theme_tsr()
}

#' @describeIn figures Stacked bar of off-target sites broken down by
#'   transcript region crossed with match type. \code{risk_score} is derived
#'   from exactly these two variables (CDS = 5, 3'UTR = 3, 5'UTR = 1, plus 2
#'   for full complementarity), so this view shows where a guide's risk comes
#'   from rather than the aggregate.
#' @param position Either \code{"stack"} (default, absolute site counts) or
#'   \code{"fill"} (composition as a percentage of each guide's own sites).
#' @export
plot_region_matchtype <- function(enriched, position = c("stack", "fill")) {
  position <- match.arg(position)
  .require_plot_pkgs()
  abort_missing_cols(enriched, c("siRNA_name", "region", "match_type"),
                     "plot_region_matchtype")
  d <- tibble::as_tibble(enriched)
  d$strata <- paste(d$region,
                    ifelse(d$match_type == "full_complementarity",
                           "full", "partial"), sep = " / ")
  lv <- c("CDS / full", "CDS / partial", "3'UTR / full", "3'UTR / partial",
          "5'UTR / full", "5'UTR / partial")
  lv <- c(intersect(lv, unique(d$strata)), setdiff(unique(d$strata), lv))
  d$strata <- factor(d$strata, levels = lv)
  pal <- stats::setNames(
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "OrRd")[3:9])(
      length(lv)), lv)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$siRNA_name,
                                       fill = .data$strata)) +
    ggplot2::geom_bar(position = if (position == "fill") "fill" else "stack",
                      width = 0.7) +
    ggplot2::scale_fill_manual(values = pal, name = "Region / match") +
    ggplot2::labs(
      title = if (position == "fill")
        "Off-target composition by region and match type"
      else "Off-target sites by region and match type",
      x = NULL,
      y = if (position == "fill") "Share of the guide's sites"
          else "Off-target sites") +
    .theme_tsr()
  if (position == "fill")
    p <- p + ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1))
  p
}

#' @describeIn figures Distribution of per-gene ortholog identity per
#'   species, from the gene-level table used for species scoring. Table-style summaries
#'   report one mean per species; this shows the spread behind that mean.
#'   Only genes with a resolved one-to-one ortholog carry an
#'   \code{identity_pct} and are plotted.
#' @export
plot_conservation_spread <- function(species_result) {
  .require_plot_pkgs()
  d <- if (!is.null(species_result$gene_detail))
    species_result$gene_detail else species_result$site_detail
  abort_missing_cols(d, c("siRNA_name", "gene_name", "species", "identity_pct"),
                     "plot_conservation_spread")
  d <- tibble::as_tibble(d) |>
    dplyr::distinct(.data$siRNA_name, .data$gene_name, .data$species,
                    .keep_all = TRUE)
  d <- d[!is.na(d$identity_pct), ]
  if (nrow(d) == 0)
    stop("no rows with a non-missing identity_pct", call. = FALSE)
  lv <- intersect(c("mouse","rat","cyno","rhesus","rabbit","dog"),
                  unique(d$species))
  d$species <- factor(d$species, levels = lv)
  med <- stats::aggregate(identity_pct ~ species, d, stats::median)

  ggplot2::ggplot(d, ggplot2::aes(x = .data$species,
                                  y = .data$identity_pct)) +
    ggplot2::geom_violin(fill = "#dddddd", color = "#888888",
                         width = 0.85, alpha = 0.75, trim = FALSE) +
    ggplot2::geom_jitter(width = 0.14, alpha = 0.35, size = 0.9,
                         color = "#666666") +
    ggplot2::stat_summary(fun = stats::median, geom = "crossbar",
                          width = 0.5, linewidth = 0.35,
                          color = "#c0392b") +
    ggplot2::geom_text(data = med,
      ggplot2::aes(label = sprintf("%.1f", .data$identity_pct)),
      vjust = -0.9, size = 3.0, fontface = "bold", color = "#c0392b") +
    ggplot2::labs(
      title = "Ortholog identity spread across candidate species",
      x = NULL, y = "Ortholog identity (%)") +
    .theme_tsr()
}

#' @describeIn figures Ortholog coverage per species: the share of critical
#'   off-target genes for which a one-to-one ortholog was resolved. This is
#'   the complement of \code{\link{plot_conservation_spread}}, which can only
#'   plot the genes that were resolved; a species may show high identity on a
#'   thin slice of genes, and that trade-off is invisible in an
#'   identity-only view.
#' @export
plot_ortholog_coverage <- function(species_result) {
  .require_plot_pkgs()
  d <- if (!is.null(species_result$gene_detail))
    species_result$gene_detail else species_result$site_detail
  abort_missing_cols(d, c("siRNA_name", "species", "gene_name", "ortholog_type"),
                     "plot_ortholog_coverage")
  d <- tibble::as_tibble(d) |>
    dplyr::distinct(.data$siRNA_name, .data$gene_name, .data$species,
                    .keep_all = TRUE)
  d$resolved <- !is.na(d$ortholog_type) & d$ortholog_type == "one2one"
  agg <- stats::aggregate(resolved ~ species, d,
                          function(z) c(n = length(z), k = sum(z)))
  cov <- data.frame(species = agg$species,
                    n = agg$resolved[, "n"], k = agg$resolved[, "k"])
  cov$share <- cov$k / cov$n
  lv <- intersect(c("mouse","rat","cyno","rhesus","rabbit","dog"), cov$species)
  cov$species <- factor(cov$species, levels = lv)
  cov <- cov[order(cov$species), ]
  hi <- max(cov$share)

  ggplot2::ggplot(cov, ggplot2::aes(x = .data$species, y = .data$share,
                                    fill = .data$share)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.0f%%\n(%d/%d)",
                                   100 * .data$share, .data$k, .data$n)),
      vjust = -0.25, size = 3.0, lineheight = 0.9) +
    ggplot2::scale_fill_gradientn(
      colours = grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(9, "YlGnBu"))(100),
      limits = c(0, 1), guide = "none") +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1)) +
    ggplot2::coord_cartesian(ylim = c(0, hi * 1.25)) +
    ggplot2::labs(
      title = "One-to-one ortholog coverage of critical off-target genes",
      x = NULL, y = "Share of gene-species pairs resolved") +
    .theme_tsr()
}

#' @describeIn figures Faceted comparison of a query siRNA against the
#'   reference cohort on raw count dimensions, one panel per dimension. Unlike
#'   \code{\link{plot_ranking_bar}}, nothing here is normalised, weighted or
#'   combined: each panel shows the cohort's own distribution of an
#'   observable count (sites, CDS sites, seed sites, critical-gene hits) with
#'   the query drawn on top, so the position can be read in the original
#'   units. Counts span orders of magnitude, so panels default to a log10 x
#'   axis.
#' @param dims Character vector of \code{reference_set$metrics} columns to
#'   panel. Default drops \code{high_risk_count} (collinear with
#'   \code{cds_count}) and the full-complementarity counts, which only take
#'   the values 0-2.
#' @param log_x Logical, log10 the count axis. Default \code{TRUE}.
#' @export
plot_dimension_panels <- function(ranking, reference = "builtin",
                                  dims = c("total_offtargets", "cds_count",
                                           "seed_region_count",
                                           "high_risk_critical",
                                           "critical_gene_count"),
                                  log_x = TRUE) {
  .require_plot_pkgs()
  require_pkg("ggrepel", "plotting")
  ref <- if (identical(reference, "builtin"))
    .load_pkg_data("reference_set") else reference
  if (!"metrics" %in% names(ref))
    stop("reference must contain 'metrics'", call. = FALSE)
  rm_ <- ref$metrics
  dims <- intersect(dims, intersect(names(rm_), names(ranking)))
  if (length(dims) == 0)
    stop("none of `dims` are present in both ranking and reference metrics",
         call. = FALSE)

  nice <- c(total_offtargets = "Total off-target sites",
            cds_count = "Sites in CDS",
            seed_region_count = "Seed sites in 3'UTR",
            high_risk_critical = "High-risk sites in critical genes",
            critical_gene_count = "Hits in critical genes",
            high_risk_count = "High-risk sites",
            full_comp_count = "Full-complementarity sites",
            full_comp_cds = "Full-complementarity sites in CDS")
  lab_for <- function(k) ifelse(k %in% names(nice), nice[k], k)

  ref_long <- do.call(rbind, lapply(dims, function(d)
    data.frame(dim = lab_for(d), value = as.numeric(rm_[[d]]),
               stringsAsFactors = FALSE)))
  q_long <- do.call(rbind, lapply(dims, function(d)
    data.frame(dim = lab_for(d), value = as.numeric(ranking[[d]]),
               siRNA_name = ranking$siRNA_name, stringsAsFactors = FALSE)))
  lv <- lab_for(dims)
  ref_long$dim <- factor(ref_long$dim, levels = lv)
  q_long$dim   <- factor(q_long$dim,   levels = lv)

  # a log axis cannot show zero; shift both sets by the same constant
  if (log_x) {
    ref_long$value <- ref_long$value + 1
    q_long$value   <- q_long$value + 1
  }
  pal <- .user_palette(unique(q_long$siRNA_name))

  p <- ggplot2::ggplot() +
    ggplot2::geom_violin(
      data = ref_long,
      ggplot2::aes(x = .data$value, y = 1),
      fill = "#dddddd", color = "#888888", alpha = 0.75, width = 0.9) +
    ggplot2::geom_jitter(
      data = ref_long, ggplot2::aes(x = .data$value, y = 1),
      height = 0.22, alpha = 0.3, size = 0.7, color = "#666666") +
    ggplot2::geom_point(
      data = q_long,
      ggplot2::aes(x = .data$value, y = 1, fill = .data$siRNA_name),
      shape = 23, size = 4, color = "black") +
    ggrepel::geom_text_repel(
      data = q_long,
      ggplot2::aes(x = .data$value, y = 1,
                   label = format(.data$value - as.integer(log_x),
                                  big.mark = ",", trim = TRUE)),
      size = 2.9, fontface = "bold", direction = "y", nudge_y = 0.3,
      segment.color = "#aaaaaa", segment.size = 0.3,
      min.segment.length = 0, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = pal, name = "Query siRNA") +
    ggplot2::facet_wrap(~ dim, ncol = 1, scales = "free_x") +
    ggplot2::labs(
      title = "Query siRNA vs. reference cohort, raw counts",
      x = if (log_x) "Count (log10, zero shifted by 1)" else "Count",
      y = NULL) +
    .theme_tsr() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold"))
  if (log_x) p <- p + ggplot2::scale_x_log10(labels = scales::label_number())
  p
}

#' @describeIn figures Where a single siRNA falls inside the reference cohort
#'   on each raw count dimension, as a percentile with the observed count
#'   printed alongside. One row per dimension; the 50th percentile is marked.
#'   Lower is safer, so a bar reaching the right-hand side means the query
#'   carries more of that liability than most of the cohort.
#' @param siRNA Name of the single siRNA in \code{ranking} to profile. Default
#'   \code{NULL} uses the first row.
#' @export
plot_dimension_percentile <- function(ranking, reference = "builtin",
                                      siRNA = NULL,
                                      dims = c("total_offtargets", "cds_count",
                                               "seed_region_count",
                                               "high_risk_critical",
                                               "critical_gene_count")) {
  .require_plot_pkgs()
  ref <- if (identical(reference, "builtin"))
    .load_pkg_data("reference_set") else reference
  if (!"metrics" %in% names(ref))
    stop("reference must contain 'metrics'", call. = FALSE)
  rm_ <- ref$metrics
  if (is.null(siRNA)) siRNA <- ranking$siRNA_name[1]
  row <- ranking[ranking$siRNA_name == siRNA, , drop = FALSE]
  if (nrow(row) == 0)
    stop(sprintf("siRNA '%s' not found in ranking", siRNA), call. = FALSE)
  dims <- intersect(dims, intersect(names(rm_), names(row)))
  if (length(dims) == 0)
    stop("no usable dims", call. = FALSE)

  nice <- c(total_offtargets = "Total off-target sites",
            cds_count = "Sites in CDS",
            seed_region_count = "Seed sites in 3'UTR",
            high_risk_critical = "High-risk sites in critical genes",
            critical_gene_count = "Hits in critical genes")
  d <- do.call(rbind, lapply(dims, function(k) {
    v <- as.numeric(row[[k]][1]); pool <- as.numeric(rm_[[k]])
    data.frame(dim = ifelse(k %in% names(nice), nice[k], k),
               value = v, pct = 100 * mean(pool < v),
               med = stats::median(pool), stringsAsFactors = FALSE)
  }))
  d$dim <- factor(d$dim, levels = rev(d$dim))

  # Build the annotation and decide, per row, whether it fits beyond the bar.
  # The axis is fixed to 0-100, so a bar ending near 100 has no room outside.
  d$lab <- sprintf("%sth pct - %s observed (cohort median %s)",
                   format(round(d$pct), trim = TRUE),
                   format(d$value, big.mark = ",", trim = TRUE),
                   format(round(d$med), big.mark = ",", trim = TRUE))
  # ~0.62 axis units per character at size 3 on a 0-100 axis.
  d$lab_inside <- d$pct + 0.62 * nchar(d$lab) > 100
  d$lab_hjust <- ifelse(d$lab_inside, 1.02, -0.02)

  ggplot2::ggplot(d, ggplot2::aes(x = .data$pct, y = .data$dim,
                                  fill = .data$pct)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_vline(xintercept = 50, linetype = "dashed",
                        color = "#888888", linewidth = 0.4) +
    # A long bar leaves no room beyond its end, so the annotation sits inside
    # the bar in that case and beyond it otherwise.
    ggplot2::geom_text(
      ggplot2::aes(label = .data$lab, hjust = .data$lab_hjust,
                   colour = .data$lab_inside),
      size = 3.0, show.legend = FALSE) +
    ggplot2::scale_colour_manual(
      values = c(`TRUE` = "white", `FALSE` = "#222222"), guide = "none") +
    ggplot2::scale_fill_gradientn(
      colours = grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(9, "YlOrRd"))(100),
      limits = c(0, 100), guide = "none") +
    ggplot2::coord_cartesian(xlim = c(0, 100)) +
    ggplot2::labs(
      title = sprintf("%s within the reference cohort", siRNA),
      x = "Percentile of the reference cohort", y = NULL) +
    .theme_tsr() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

#' @describeIn figures Heatmap of siRNA x species conservation score.
#' @param species_result Output of \code{\link{species_match}}.
#' @export
plot_species_heatmap <- function(species_result) {
  .require_plot_pkgs()
  sp_pal <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(9, "YlGnBu"))(100)
  sp_long <- tidyr::pivot_longer(species_result$species_score,
    -"siRNA_name", names_to = "species", values_to = "score")
  ggplot2::ggplot(sp_long, ggplot2::aes(x = .data$species, y = .data$siRNA_name,
                                         fill = .data$score)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$score)),
                        size = 3.5) +
    ggplot2::scale_fill_gradientn(colours = sp_pal, limits = c(0, 1),
                                   name = "Conservation\nscore") +
    ggplot2::labs(title = "Animal model recommendation",
       x = NULL, y = NULL) +
    .theme_tsr() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

#' @describeIn figures Boxplot of per-species conservation score in the
#'   reference cohort, with the user's siRNA(s) overlaid as colored
#'   diamonds. Shows how each candidate animal model for the user's
#'   siRNA compares to the distribution observed across known siRNAs.
#' @export
plot_species_vs_reference <- function(species_result, reference = "builtin") {
  .require_plot_pkgs()
  ref <- if (identical(reference, "builtin"))
    .load_pkg_data("reference_set") else reference
  if (!"species_score" %in% names(ref)) {
    stop("reference must contain 'species_score'", call. = FALSE)
  }

  ref_long <- tidyr::pivot_longer(ref$species_score,
    -"siRNA_name", names_to = "species", values_to = "score")
  ref_long <- ref_long[!is.na(ref_long$score), ]
  user_long <- tidyr::pivot_longer(species_result$species_score,
    -"siRNA_name", names_to = "species", values_to = "score")
  user_long <- user_long[!is.na(user_long$score), ]

  species_levels <- intersect(
    c("mouse","rat","cyno","rhesus","rabbit","dog"),
    union(unique(ref_long$species), unique(user_long$species)))
  ref_long$species  <- factor(ref_long$species,  levels = species_levels)
  user_long$species <- factor(user_long$species, levels = species_levels)

  n_user <- max(3, length(unique(user_long$siRNA_name)))
  user_pal <- RColorBrewer::brewer.pal(n_user, "Set1")[
    seq_len(length(unique(user_long$siRNA_name)))]

  ggplot2::ggplot() +
    ggplot2::geom_violin(
      data = ref_long,
      ggplot2::aes(x = .data$species, y = .data$score),
      fill = "#dddddd", color = "#888888", alpha = 0.7,
      width = 0.85, trim = FALSE) +
    ggplot2::geom_jitter(
      data = ref_long,
      ggplot2::aes(x = .data$species, y = .data$score),
      width = 0.15, alpha = 0.25, size = 0.6, color = "#666666") +
    ggplot2::geom_point(
      data = user_long,
      ggplot2::aes(x = .data$species, y = .data$score,
                   color = .data$siRNA_name),
      size = 4.5, shape = 18) +
    ggplot2::scale_color_manual(values = user_pal, name = "Your siRNA") +
    # coord_cartesian (not scale_y_continuous(limits=)) so the violin's
    # untrimmed density tail and the jittered reference points are clipped
    # visually rather than dropped before the stat is computed.
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(title = "Animal model fit: your siRNA vs. reference cohort",
       x = NULL, y = "Conservation score (0-1; higher = better model)") +
    .theme_tsr()
}

#' @describeIn figures Percentile rank of the user siRNA's species
#'   conservation score within the reference cohort, one bar per species.
#'   Higher percentile means \dQuote{this animal model has higher
#'   conservation than X\% of reference siRNAs} — a stronger model fit.
#' @export
plot_species_percentile <- function(species_result, reference = "builtin") {
  .require_plot_pkgs()
  ref <- if (identical(reference, "builtin"))
    .load_pkg_data("reference_set") else reference
  if (!"species_score" %in% names(ref)) {
    stop("reference must contain 'species_score'", call. = FALSE)
  }

  user_long <- tidyr::pivot_longer(species_result$species_score,
    -"siRNA_name", names_to = "species", values_to = "score")
  user_long <- user_long[!is.na(user_long$score), ]

  ref_long <- tidyr::pivot_longer(ref$species_score,
    -"siRNA_name", names_to = "species", values_to = "score")
  ref_long <- ref_long[!is.na(ref_long$score), ]

  user_long$percentile <- vapply(seq_len(nrow(user_long)), function(i) {
    sp <- as.character(user_long$species[i])
    sc <- user_long$score[i]
    pool <- ref_long$score[ref_long$species == sp]
    if (length(pool) == 0) return(NA_real_)
    100 * mean(pool <= sc, na.rm = TRUE)
  }, numeric(1))

  species_levels <- intersect(
    c("mouse","rat","cyno","rhesus","rabbit","dog"),
    unique(user_long$species))
  user_long$species <- factor(user_long$species, levels = species_levels)

  # Sequential rather than red-green diverging: higher percentile means better
  # conservation, and the scheme stays legible under colour vision deficiency.
  bar_pal <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(9, "YlGnBu"))(100)

  p <- ggplot2::ggplot(user_long,
    ggplot2::aes(x = .data$species, y = .data$percentile,
                 fill = .data$percentile)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.0f%%\n(score=%.2f)",
                                    .data$percentile, .data$score)),
      vjust = -0.15, size = 3.0, lineheight = 0.85) +
    ggplot2::scale_fill_gradientn(colours = bar_pal, limits = c(0, 100),
                                   guide = "none") +
    ggplot2::scale_y_continuous(limits = c(0, 118),
                                 breaks = c(0, 25, 50, 75, 100)) +
    ggplot2::geom_hline(yintercept = 50, linetype = "dashed",
                         color = "#888888", linewidth = 0.4) +
    ggplot2::labs(title = "Animal model fit: percentile vs. reference cohort",
       x = NULL, y = "Percentile (0-100)") +
    .theme_tsr()

  if (length(unique(user_long$siRNA_name)) > 1) {
    p <- p + ggplot2::facet_wrap(~ .data$siRNA_name)
  }
  p
}

#' @describeIn figures Heatmap of top off-target genes x tissues (log1p TPM).
#' @param tissues_of_interest Character vector of tissue column names.
#'   Defaults to a 10-tissue panel covering nerve, brain, liver, heart,
#'   muscle, kidney, lung, blood, spleen, adipose.
#' @param top_n Number of top genes (default 20).
#' @export
plot_tissue_heatmap <- function(enriched, tissue,
                                tissues_of_interest = NULL,
                                top_n = 20) {
  .require_plot_pkgs()
  if (is.null(tissues_of_interest)) {
    # Lung and Liver lead: the two tissues most often reported for systemic
    # and hepatically delivered siRNA.
    tissues_of_interest <- c("Lung", "Liver", "Nerve_Tibial", "Brain_Cortex",
      "Heart_Left_Ventricle", "Muscle_Skeletal", "Kidney_Cortex",
      "Whole_Blood", "Spleen", "Adipose_Subcutaneous")
  }
  top_gene_names <- utils::head(tissue$gene_name, top_n)
  use_tissues <- intersect(tissues_of_interest, names(enriched))
  if (length(use_tissues) == 0)
    stop("none of `tissues_of_interest` are columns of `enriched`",
         call. = FALSE)
  d <- dplyr::filter(enriched, .data$gene_name %in% top_gene_names)
  d <- dplyr::group_by(d, .data$gene_name)
  d <- dplyr::summarise(d,
    dplyr::across(dplyr::all_of(use_tissues),
                   ~ mean(.x, na.rm = TRUE)),
    .groups = "drop")
  # A gene with no measurable expression in any displayed tissue contributes
  # an all-blank row, so drop it and keep however many genes remain rather
  # than padding the panel out to top_n.
  keep <- rowSums(as.matrix(d[, use_tissues, drop = FALSE]), na.rm = TRUE) > 0
  d <- d[keep, , drop = FALSE]
  if (nrow(d) == 0)
    stop("every candidate gene has zero expression in the selected tissues",
         call. = FALSE)
  top_gene_names <- top_gene_names[top_gene_names %in% d$gene_name]
  d <- tidyr::pivot_longer(d, -"gene_name",
    names_to = "tissue", values_to = "tpm")
  d$tissue <- factor(d$tissue, levels = use_tissues)
  d$gene_name <- factor(d$gene_name, levels = rev(top_gene_names))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$tissue, y = .data$gene_name,
                                   fill = log1p(.data$tpm))) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_distiller(palette = "YlGnBu", direction = 1,
       name = "log1p(TPM)") +
    ggplot2::labs(title = "Tissue expression of top off-target genes",
       x = NULL, y = NULL) +
    .theme_tsr() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                    panel.grid = ggplot2::element_blank())
}

#' @describeIn figures Bubble chart of siRNA x critical-gene off-targets.
#' @param top_n Number of genes to display, taken in order of
#'   \code{bio_attention}. Default \code{10}; plotting every critical gene
#'   makes the axis unreadable.
#' @export
plot_critical_bubble <- function(enriched, tissue, top_n = 10) {
  .require_plot_pkgs()
  d <- dplyr::filter(enriched,
    .data$cancer_gene_v2 | .data$AE_gene | .data$immune_gene_immport)
  if (nrow(d) == 0) {
    return(ggplot2::ggplot() +
           ggplot2::annotate("text", x = 0.5, y = 0.5,
             label = "No critical-gene off-targets") +
           .theme_tsr())
  }
  # Every critical gene would otherwise be drawn, which makes the axis
  # illegible; keep the highest-ranked genes only.
  keep <- utils::head(tissue$gene_name[order(-tissue$bio_attention)], top_n)
  d <- dplyr::filter(d, .data$gene_name %in% keep)
  if (nrow(d) == 0)
    stop("none of the top-ranked genes carry critical-gene sites", call. = FALSE)
  d <- dplyr::count(d, .data$siRNA_name, .data$gene_name, name = "n_hits")
  d <- dplyr::left_join(d, tissue[, c("gene_name", "bio_attention")],
                          by = "gene_name")
  d$gene_name <- factor(d$gene_name,
    levels = unique(d$gene_name[order(d$bio_attention)]))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$siRNA_name, y = .data$gene_name,
                                   size = .data$n_hits, color = .data$bio_attention)) +
    ggplot2::geom_point() +
    ggplot2::scale_size_continuous(range = c(2, 10), name = "# hits") +
    # a discrete axis expands by only 0.6 units by default, which clips the
    # largest bubbles on the first and last row
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.8)) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::scale_color_distiller(palette = "YlOrRd", direction = 1,
                                    name = "bio_attention") +
    ggplot2::labs(title = "siRNA x critical-gene off-target map",
       x = NULL, y = NULL) +
    .theme_tsr()
}

#' Generate the standard 6-figure bundle and write to disk
#'
#' Convenience wrapper that runs all six figure helpers and writes both
#' PDF and PNG versions to \code{output_dir}.
#'
#' @param enriched Output of \code{\link{enrich_annotations}}.
#' @param ranking Output of \code{\link{rank_sirna}}.
#' @param tissue Output of \code{\link{tissue_safety}}.
#' @param species_result Output of \code{\link{species_match}}.
#' @param output_dir Directory to write to (created if absent).
#' @param formats Character vector of formats. Default \code{c("pdf","png")}.
#'
#' @return Invisibly, a named list of file paths written.
#' @examples
#' \dontrun{
#' generate_figures(enriched, ranking, tissue, species_result,
#'                  output_dir = "results/figures")
#' }
#' @param with_reference Logical. If \code{TRUE} (default), additional
#'   comparison figures using the built-in \code{\link{reference_set}}
#'   are written: \code{01b_ranking_vs_reference} (violin),
#'   \code{04b_species_vs_reference} (boxplot per species),
#'   \code{04c_species_percentile} (percentile bars).
#' @param sites Optional. The CORE table from \code{\link{annotate_sites}},
#'   which still carries \code{n_mismatch}. Supplying it adds the
#'   mismatch-distribution figures; \code{enriched} alone cannot produce them
#'   because \code{\link{enrich_annotations}} drops that column.
#' @param variants Either \code{"primary"} (default; one figure per data
#'   view) or \code{"both"}, which additionally writes the alternative
#'   rendering of each: cumulative mismatch curve (\code{07b}), percentage
#'   composition (\code{08b}) and ortholog coverage (\code{09b}).
#' @export
generate_figures <- function(enriched, ranking, tissue,
                             species_result = NULL,
                             output_dir,
                             formats = c("pdf","png"),
                             with_reference = TRUE,
                             sites = NULL,
                             variants = c("primary", "both")) {
  variants <- match.arg(variants)
  .require_plot_pkgs()
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  plots <- list(
    `01_sirna_ranking_bar`     = list(p = plot_ranking_bar(ranking),    w=6, h=5),
    `02_risk_distribution`     = list(p = plot_risk_distribution(enriched), w=6, h=5),
    `03_gene_attention_scatter`= list(p = plot_gene_attention(tissue),  w=9, h=6.5),
    `05_tissue_heatmap`        = list(p = plot_tissue_heatmap(enriched, tissue), w=10, h=6),
    `06_sirna_gene_bubble`     = list(p = plot_critical_bubble(enriched, tissue), w=9, h=7)
  )
  if (!is.null(species_result)) {
    plots$`04_species_heatmap` <- list(
      p = plot_species_heatmap(species_result), w=9, h=5)
  }
  # --- extra dimensions: mismatch decay, region x match type, ortholog ---
  if (!is.null(sites)) {
    plots$`07a_mismatch_distribution` <- list(
      p = plot_mismatch_distribution(sites), w = 9, h = 5)
    if (variants == "both")
      plots$`07b_mismatch_cumulative` <- list(
        p = plot_mismatch_cumulative(sites), w = 9, h = 5)
  }
  plots$`08a_region_matchtype` <- list(
    p = plot_region_matchtype(enriched), w = 6, h = 5)
  if (variants == "both")
    plots$`08b_region_matchtype_share` <- list(
      p = plot_region_matchtype(enriched, position = "fill"), w = 6, h = 5)
  if (!is.null(species_result) &&
      "site_detail" %in% names(species_result)) {
    plots$`09a_conservation_spread` <- list(
      p = plot_conservation_spread(species_result), w = 9, h = 5)
    if (variants == "both")
      plots$`09b_ortholog_coverage` <- list(
        p = plot_ortholog_coverage(species_result), w = 9, h = 5)
  }

  if (with_reference) {
    plots$`01b_ranking_vs_reference` <- list(
      p = plot_ranking_bar(ranking, reference = "builtin"), w=6, h=6)
    plots$`02b_risk_distribution_vs_reference` <- list(
      p = plot_risk_distribution(enriched, reference = "builtin"),
      w = 16, h = 6)
    if (!is.null(species_result)) {
      plots$`04b_species_vs_reference` <- list(
        p = plot_species_vs_reference(species_result, reference = "builtin"),
        w=10, h=6)
      plots$`04c_species_percentile` <- list(
        p = plot_species_percentile(species_result, reference = "builtin"),
        w=10, h=6)
    }
  }

  paths <- list()
  for (nm in names(plots)) {
    for (fmt in formats) {
      fn <- file.path(output_dir, paste0(nm, ".", fmt))
      ggplot2::ggsave(fn, plots[[nm]]$p,
                       width = plots[[nm]]$w, height = plots[[nm]]$h,
                       dpi = 150)
      paths[[paste0(nm, ".", fmt)]] <- fn
    }
  }
  invisible(paths)
}
