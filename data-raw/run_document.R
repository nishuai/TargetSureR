# (user library resolved by R itself; no .libPaths() override)

# Make sure devtools/roxygen2 are present
need <- setdiff(c("devtools", "roxygen2"),
                rownames(installed.packages()))
if (length(need)) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  install.packages(need)
}

cat("R:", R.version.string, "\n")
cat("devtools :", as.character(packageVersion("devtools")), "\n")
cat("roxygen2 :", as.character(packageVersion("roxygen2")), "\n\n")

# Document with roxygen2
cat("=== devtools::document() ===\n")
devtools::document()

cat("\n=== Generated man/ files ===\n")
print(list.files("man"))
