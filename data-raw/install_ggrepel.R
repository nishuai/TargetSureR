options(repos = c(CRAN = "https://cloud.r-project.org"))
# (user library resolved by R itself; no .libPaths() override)
install.packages("ggrepel", quiet = TRUE)
cat("ggrepel installed?", requireNamespace("ggrepel", quietly = TRUE), "\n")
