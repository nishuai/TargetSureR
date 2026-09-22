# Regression tests for the built-in reference_set.
#
# The first two guard a real bug: build_reference_set.R once called
# rank_sirna(enriched) without reference=, which defaults to "builtin" and
# normalised the cohort against a PREVIOUSLY saved reference_set.rda. The
# shipped composite_score then sat on a foreign scale (median 0.583) while
# rank_sirna(reference = "builtin") placed user siRNAs on the cohort's own
# scale (median 0.223), making every user-vs-reference comparison wrong.

test_that("shipped composite_score is on the same scale as reference='builtin'", {
  dims <- names(default_weights())
  m    <- reference_set$metrics
  th   <- reference_set$tier_thresholds
  w    <- default_weights()[dims]
  w    <- w / sum(w)

  normed <- vapply(dims, function(d) {
    mn <- th$min[th$dimension == d]
    mx <- th$max[th$dimension == d]
    if (length(mn) == 0 || mx == mn) return(rep(0, nrow(m)))
    v <- (m[[d]] - mn) / (mx - mn)
    pmin(pmax(v, 0), 1)
  }, numeric(nrow(m)))

  recomputed <- as.numeric(normed %*% as.numeric(w))
  expect_equal(recomputed, as.numeric(m$composite_score), tolerance = 1e-10)
})

test_that("tier_thresholds min/max are those of the stored metrics", {
  th <- reference_set$tier_thresholds
  m  <- reference_set$metrics
  for (d in names(default_weights())) {
    expect_equal(th$min[th$dimension == d], min(m[[d]]), info = d)
    expect_equal(th$max[th$dimension == d], max(m[[d]]), info = d)
  }
})

test_that("reference_set is internally consistent", {
  expect_equal(reference_set$n_siRNAs, 94L)
  expect_equal(nrow(reference_set$metadata), 94L)
  expect_equal(nrow(reference_set$metrics), 94L)
  expect_setequal(reference_set$metadata$siRNA_name,
                  reference_set$metrics$siRNA_name)
  # all guides share the 19 nt convention
  expect_true(all(reference_set$metadata$guide_len == 19L))
  expect_true(all(nchar(reference_set$metadata$guide) == 19L))
  # rank is a permutation of 1..94 ordered by ascending composite
  expect_setequal(reference_set$metrics$rank, seq_len(94))
  expect_false(is.unsorted(
    reference_set$metrics$composite_score[order(reference_set$metrics$rank)]))
})
