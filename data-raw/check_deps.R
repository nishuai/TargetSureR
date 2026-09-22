cat(R.version.string, "\n\n")

repos <- getOption("repos")
if (is.null(repos) || repos[["CRAN"]] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
cat("CRAN repo:", getOption("repos")[["CRAN"]], "\n\n")

pkgs <- c(
  "BiocManager", "remotes",
  "dplyr", "tibble", "tidyr", "readr", "purrr", "rlang", "ggplot2", "scales",
  "testthat", "knitr", "rmarkdown",
  "SeedMatchR",
  # Bioconductor:
  "Biostrings", "ensembldb", "AnnotationHub", "AnnotationFilter",
  "GenomicRanges", "GenomicFeatures", "BiocGenerics", "IRanges", "S4Vectors",
  "EnsDb.Hsapiens.v86"
)

cat(sprintf("%-25s %s\n", "package", "status"))
cat(strrep("-", 50), "\n", sep="")
for (p in pkgs) {
  has <- requireNamespace(p, quietly = TRUE)
  ver <- if (has) as.character(utils::packageVersion(p)) else "-"
  cat(sprintf("%-25s %-12s %s\n", p, if (has) "INSTALLED" else "MISSING", ver))
}

cat("\nBiocManager version:\n")
if (requireNamespace("BiocManager", quietly = TRUE)) {
  print(BiocManager::version())
} else {
  cat("BiocManager not installed yet\n")
}
