# (user library resolved by R itself; no .libPaths() override)

# Build then check
options(repos = c(CRAN = "https://cloud.r-project.org"))

cat("=== devtools::check() ===\n")
res <- devtools::check(
  pkg = ".",
  document = FALSE,        # we already documented
  manual   = FALSE,        # avoid pdflatex requirement
  vignettes = FALSE,       # speed; do separately later
  cran     = TRUE,
  error_on = "never",      # return result instead of stopping
  args     = c("--no-manual", "--as-cran"),
  build_args = "--no-manual"
)

cat("\n=== Summary ===\n")
print(res)
