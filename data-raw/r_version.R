cat(R.version.string, "\n")
cat("R version (numeric): ", R.version$major, ".", R.version$minor, "\n", sep = "")
cat("Repos:\n")
print(getOption("repos", default = c(CRAN = "https://cloud.r-project.org")))

# What is the matching Bioc release?
# R 4.1.x -> Bioc 3.13/3.14
# R 4.2.x -> Bioc 3.15/3.16
# R 4.3.x -> Bioc 3.17/3.18
# R 4.4.x -> Bioc 3.19/3.20
# R 4.5.x -> Bioc 3.21
