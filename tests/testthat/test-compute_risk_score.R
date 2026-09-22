test_that("core_columns and enriched_columns are character", {
  expect_type(core_columns(), "character")
  expect_length(core_columns(), 10)
  expect_true(all(core_columns() %in% enriched_columns()))
  expect_length(enriched_columns(), 13)
})

test_that("default_weights is four equal count dimensions summing to 1", {
  w <- default_weights()
  expect_equal(sum(w), 1, tolerance = 1e-8)
  expect_named(w)
  expect_length(w, 4)
  expect_setequal(names(w), c("cds_count", "seed_region_count",
                              "critical_gene_count", "high_risk_critical"))
  expect_true(all(abs(w - 0.25) < 1e-8))
})

test_that("compute_risk_score follows the documented formula", {
  rs <- compute_risk_score(
    region     = c("CDS", "CDS",   "3'UTR", "3'UTR", "5'UTR", "ncRNA"),
    match_type = c("full_complementarity", "partial_match",
                   "full_complementarity", "partial_match",
                   "unknown",               "unknown")
  )
  expect_equal(rs, c(7L, 5L, 5L, 3L, 1L, 1L))
  expect_true(all(rs <= 10L))
  expect_true(all(rs >= 1L))
})

test_that("compute_risk_score caps at 10", {
  rs <- compute_risk_score(rep("CDS", 3),
                            rep("full_complementarity", 3))
  expect_true(all(rs <= 10L))
})

test_that("compute_risk_score vector length matches input", {
  rs <- compute_risk_score(c("CDS"), c("partial_match"))
  expect_length(rs, 1)
})
