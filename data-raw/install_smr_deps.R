options(repos = c(CRAN = "https://cloud.r-project.org"),
        timeout = 1200,
        Ncpus = 4)
# (user library resolved by R itself; no .libPaths() override)

log <- function(...) {
  cat("[", format(Sys.time(), "%H:%M:%S"), "] ", sep = "")
  cat(..., "\n", sep = "")
  flush.console()
}

log("install target: ", .libPaths()[1])
log("make:           ", Sys.which("make"))

pkgs <- c("caTools", "gtools", "twosamples", "gt")
for (p in pkgs) {
  log("--- ", p, " ---")
  if (requireNamespace(p, quietly = TRUE)) {
    log("  already installed: ", utils::packageVersion(p))
    next
  }
  res <- tryCatch({
    install.packages(p, quiet = FALSE)
    "OK"
  }, error = function(e) paste("ERROR:", conditionMessage(e)),
     warning = function(w) paste("WARN:", conditionMessage(w)))
  ok <- requireNamespace(p, quietly = TRUE)
  log("  result=", res, " | installed=", ok,
      if (ok) paste0(" v", utils::packageVersion(p)) else "")
}

log("=== final ===")
for (p in c(pkgs, "SeedMatchR")) {
  has <- requireNamespace(p, quietly = TRUE)
  log("  ", sprintf("%-12s %s", p, if (has) "OK" else "MISSING"))
}
