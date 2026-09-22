test_that("load_gene_lists returns three character vectors", {
  gl <- load_gene_lists()
  expect_named(gl, c("AE_gene", "cancer_gene_v2", "immune_gene_immport"),
               ignore.order = TRUE)
  for (nm in names(gl)) {
    expect_type(gl[[nm]], "character")
    expect_true(length(gl[[nm]]) > 0)
    expect_true(all(gl[[nm]] == toupper(gl[[nm]])))
  }
})

test_that("load_gtex_example has gene_symbol and tissue columns", {
  gx <- load_gtex_example()
  expect_true("gene_symbol" %in% names(gx))
  tissues <- setdiff(names(gx), "gene_symbol")
  expect_true(length(tissues) >= 1)
  for (t in tissues) expect_type(gx[[t]], "double")
})
