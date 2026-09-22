# Install all dependencies needed for both Mode A and Mode B real-runs.
# Run with: <R-4.6>/Rscript --vanilla data-raw/install_deps.R

options(repos = c(CRAN = "https://cloud.r-project.org"),
        Ncpus = max(1, parallel::detectCores() - 1),
        timeout = 1200)

# Prefer user library; create and place first on .libPaths
user_lib <- Sys.getenv("R_LIBS_USER")
if (!nzchar(user_lib)) {
  user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library",
                         paste0(R.version$major, ".", strsplit(R.version$minor, "\\.")[[1]][1]))
}
user_lib <- path.expand(user_lib)
dir.create(user_lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))

cat("R:", R.version.string, "\n")
cat("Install target .libPaths()[1]:", .libPaths()[1], "\n")
cat("make:", Sys.which("make"), "\n")
cat("gcc :", Sys.which("gcc"),  "\n\n")

# 1) BiocManager + remotes + CRAN dependencies ---------------------------

cran_pkgs <- c("BiocManager", "remotes",
               "dplyr", "tibble", "tidyr", "readr", "purrr", "rlang",
               "ggplot2", "scales",
               "testthat", "knitr", "rmarkdown")

to_install <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1),
                                 quietly = TRUE)]
if (length(to_install)) {
  cat("Installing CRAN packages:\n  ", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install)
} else cat("CRAN packages already present.\n")

# 2) Bioconductor packages ----------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  stop("BiocManager not installed even after the CRAN step.")
}
cat("\nBiocManager version: ", as.character(BiocManager::version()), "\n")

bioc_pkgs <- c("Biostrings", "ensembldb", "AnnotationHub",
               "AnnotationFilter", "GenomicRanges", "GenomicFeatures",
               "BiocGenerics", "IRanges", "S4Vectors",
               "EnsDb.Hsapiens.v86")

bioc_to_install <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(bioc_to_install)) {
  cat("Installing Bioconductor packages:\n  ",
      paste(bioc_to_install, collapse = ", "), "\n")
  BiocManager::install(bioc_to_install, ask = FALSE, update = FALSE)
} else cat("Bioconductor packages already present.\n")

# 3) SeedMatchR (CRAN) --------------------------------------------------

if (!requireNamespace("SeedMatchR", quietly = TRUE)) {
  cat("\nInstalling SeedMatchR from CRAN...\n")
  install.packages("SeedMatchR")
}

# 4) Final report -------------------------------------------------------

cat("\n=== Final status ===\n")
all_pkgs <- c(cran_pkgs, bioc_pkgs, "SeedMatchR")
for (p in all_pkgs) {
  has <- requireNamespace(p, quietly = TRUE)
  ver <- if (has) as.character(utils::packageVersion(p)) else "-"
  cat(sprintf("  %-25s %-12s %s\n", p, if (has) "OK" else "MISSING", ver))
}
