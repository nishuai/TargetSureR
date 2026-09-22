make_fake_sites <- function() {
  tibble::tibble(
    siRNA_name = c("hsiR22", "hsiR22", "hsiR21", "hsiR21", "hsiR30"),
    siRNA      = NA_character_,
    strand     = "guide",
    transcript = c("ENST1", "ENST2", "ENST3", "ENST4", "ENST5"),
    gene_name  = c("CTNNB1", "TP53", "MYC", "IL6", "EGFR"),
    biotype    = "protein_coding",
    position   = c(379L, 4687L, 100L, 250L, 900L),
    region     = c("CDS", "3'UTR", "CDS", "3'UTR", "5'UTR"),
    match_type = c("full_complementarity", "partial_match",
                   "full_complementarity", "partial_match", "unknown"),
    risk_score = NA_integer_
  ) |>
    dplyr::mutate(risk_score = compute_risk_score(.data$region, .data$match_type))
}

test_that("enrich_annotations produces enriched_columns + tissues", {
  sites <- make_fake_sites()
  enriched <- enrich_annotations(sites)

  expect_true(all(enriched_columns() %in% names(enriched)))
  expect_equal(nrow(enriched), nrow(sites))

  tissues <- setdiff(names(enriched), enriched_columns())
  expect_true(length(tissues) >= 1)
  for (t in tissues) expect_type(enriched[[t]], "double")
})

test_that("enrich_annotations flags critical genes case-insensitively", {
  sites <- make_fake_sites()
  enriched <- enrich_annotations(sites)
  expect_true(enriched$cancer_gene_v2[enriched$gene_name == "TP53"])
})

test_that("rank_sirna returns one row per siRNA with tier", {
  sites <- make_fake_sites()
  enriched <- enrich_annotations(sites)
  ranking <- rank_sirna(enriched, target_sirnas = "hsiR22",
                        reference = "self")

  expect_equal(nrow(ranking), dplyr::n_distinct(sites$siRNA_name))
  expect_true(all(c("composite_score", "rank", "risk_tier", "is_target")
                  %in% names(ranking)))
  expect_true(ranking$is_target[ranking$siRNA_name == "hsiR22"])
})

test_that("rank_sirna with reference='builtin' gives non-zero scores for single siRNA", {
  sites <- make_fake_sites()
  enriched <- enrich_annotations(sites)
  ranking <- rank_sirna(enriched, reference = "builtin")
  expect_equal(nrow(ranking), dplyr::n_distinct(sites$siRNA_name))
  expect_true("composite_score" %in% names(ranking))
  expect_equal(attr(ranking, "reference"), "builtin")
})

test_that("tissue_safety returns per-gene attention scores", {
  sites <- make_fake_sites()
  enriched <- enrich_annotations(sites)
  ts <- tissue_safety(enriched)
  expect_true(all(c("gene_name", "bio_attention") %in% names(ts)))
  expect_equal(nrow(ts), dplyr::n_distinct(sites$gene_name))
})

test_that("species_match emits score matrix and recommendation", {
  sites <- make_fake_sites()
  res <- species_match(sites,
                       candidates = c("mouse", "rat", "cyno",
                                      "rhesus", "rabbit", "dog"))
  expect_named(res, c("species_score", "site_detail", "gene_detail",
                       "recommendation"),
               ignore.order = TRUE)
})

test_that("species_match weights each unique gene once", {
  sites <- tibble::tibble(
    siRNA_name = rep("s1", 11),
    gene_name = c(rep("G1", 10), "G2"),
    cancer_gene_v2 = TRUE,
    AE_gene = FALSE,
    immune_gene_immport = FALSE
  )
  ortholog <- tibble::tibble(
    human_symbol = c("G1", "G2"),
    mouse_ortholog_type = c("one2one", "none"),
    mouse_identity_pct = c(95, NA_real_)
  )

  out <- species_match(sites, candidates = "mouse",
                       ortholog_table = ortholog)

  expect_equal(out$species_score$mouse, 0.5)
  expect_equal(nrow(out$gene_detail), 2)
  expect_equal(nrow(out$site_detail), 11)
  expect_match(out$recommendation, "2 critical genes")
})
