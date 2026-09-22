#' Adverse-event gene panel
#'
#' A character vector of uppercase gene symbols flagged as adverse-event
#' associated. Used as the \code{AE_gene} source in
#' \code{\link{enrich_annotations}}.
#'
#' @format A character vector.
#' @source Curated manually (see \code{data-raw/make_internal_data.R}).
"ae_genes"

#' Cancer gene panel (version 2)
#'
#' Uppercase symbols for a curated cancer gene panel.
#' Replace with the full panel before production use.
#'
#' @format A character vector.
#' @source Curated manually (see \code{data-raw/make_internal_data.R}).
"cancer_genes_v2"

#' Immune-related gene panel (ImmPort-curated)
#'
#' Uppercase symbols for a curated immune-related gene panel derived from
#' ImmPort.
#'
#' @format A character vector.
#' @source ImmPort, curated subset (see \code{data-raw/make_internal_data.R}).
"immune_genes_immport"

#' Simulated tissue-expression matrix for examples
#'
#' A 20-gene by 10-tissue tibble of \strong{simulated} median TPM values,
#' shaped like a GTEx median-TPM matrix. It exists so that examples, tests
#' and the vignette run without a large download.
#'
#' \strong{The numbers are not measured expression data and must not be used
#' for analysis.} The tissue names follow GTEx conventions, but the values are
#' invented. Any gene absent from these 20 receives \code{NA}, so a real
#' off-target set scored against this matrix leaves most genes unquantified
#' and the resulting \code{tissue_safety()} ranking is meaningless.
#'
#' For analysis, download the GTEx median-TPM matrix and pass it via the
#' \code{expression_matrix} argument of \code{\link{enrich_annotations}}; the
#' vignette section \dQuote{Tissue expression} covers retrieval and
#' preparation. The published analyses used GTEx v10 (57,853 gene symbols
#' across 68 tissues).
#'
#' @format A tibble with \code{gene_symbol} plus 10 tissue columns.
#' @source Simulated (see \code{data-raw/make_internal_data.R}).
"gtex_example"

#' Simulated ortholog table for examples
#'
#' A 20-gene table of \strong{simulated} ortholog assignments across six
#' candidate nonclinical species (\code{mouse}, \code{rat}, \code{cyno},
#' \code{rhesus}, \code{rabbit}, \code{dog}). Used by
#' \code{\link{species_match}} when no \code{ortholog_table} is supplied, so
#' that examples run without a BioMart download.
#'
#' \strong{The identity percentages are invented and must not be used for
#' species selection.} Genes outside these 20 resolve to no ortholog and
#' score zero, which drags any real conservation score towards zero
#' regardless of the true orthology.
#'
#' For analysis, download a table covering your own off-target genes from
#' Ensembl BioMart and pass it as \code{ortholog_table}; the vignette section
#' \dQuote{Cross-species orthologs} gives a ready-to-run download script. The
#' published analyses used Ensembl Genes 116 over 2,143 genes.
#'
#' @format A tibble with columns \code{human_symbol},
#'   \code{\{species\}_ortholog_type}, \code{\{species\}_identity_pct}
#'   for each of the six species.
#' @source Simulated (see \code{data-raw/make_internal_data.R}).
"ortholog_table_example"

#' Built-in reference set of 94 known human siRNAs
#'
#' Pre-computed Mode A results for 94 experimentally reported human siRNAs
#' curated from the MIT/ICBP siRNA Database, aggregated to a compact summary
#' suitable for use as a normalization reference in \code{\link{rank_sirna}}
#' and as a comparison background in \code{\link{plot_ranking_bar}}.
#'
#' Used by passing \code{reference = "builtin"} to \code{\link{rank_sirna}}.
#'
#' All 94 guides are 19 nt, normalized with
#' \code{\link{normalize_guide}(target_len = 19)}. Passing a user siRNA
#' through the same function before scanning keeps the comparison
#' like-for-like; in practice the effect of a 21 nt versus 19 nt guide on
#' the site counts is small (under 1% in the authors' tests), so this is a
#' reproducibility convention rather than a large numerical correction.
#'
#' \strong{Scale of \code{composite_score}.} The stored scores are
#' self-normalized \emph{within} this cohort (each dimension rescaled to
#' its own observed min-max across the 94 siRNAs). This is the same scale
#' that \code{rank_sirna(reference = "builtin")} places a new siRNA on,
#' because that path rescales using \code{tier_thresholds}, whose
#' \code{min}/\code{max} are taken from this cohort. User scores and the
#' stored reference scores are therefore directly comparable, as are the
#' \code{risk_tier} labels derived from them.
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{\code{metadata}}{tibble (94 x 6): siRNA_name, guide,
#'     guide_len, target_gene, refseq, source.}
#'   \item{\code{metrics}}{tibble (94 x 13): the raw count endpoints +
#'     composite_score + rank + risk_tier for each siRNA, after
#'     self-normalization on the 94-siRNA cohort (see above).}
#'   \item{\code{tier_thresholds}}{data.frame (7 x 6): per-dimension
#'     min/max/q25/q50/q75 thresholds. These are what \code{rank_sirna}
#'     uses to normalize a NEW siRNA against the reference cohort.}
#'   \item{\code{species_score}}{tibble (94 x 7): siRNA x species
#'     conservation matrix.}
#'   \item{\code{gene_summary}}{tibble (94 x 2): per-siRNA count of
#'     critical (cancer / AE / immune) gene hits.}
#'   \item{\code{risk_distribution}}{tibble (376 x 3): per-siRNA count of
#'     off-target sites at each integer \code{risk_score} value, used by
#'     \code{\link{plot_risk_distribution}} for reference overlay.}
#'   \item{\code{built_at}, \code{pkg_version}, \code{n_siRNAs}}{
#'     audit metadata.}
#' }
#' @source Curated from the MIT/ICBP siRNA Database (Sharp Laboratory, Koch
#'   Institute / MIT Center for Cancer Research,
#'   \url{https://web.mit.edu/sirna/}, accessed 2026-05-15). The 144
#'   human-reactive entries were filtered to the 94 synthetic siRNA duplexes
#'   (19-30 nt), excluding 50 long-hairpin shRNA constructs. Built by
#'   \code{data-raw/build_reference_set.R}; off-target sites scanned against
#'   Ensembl v113, guides normalized by \code{\link{normalize_guide}}.
"reference_set"
