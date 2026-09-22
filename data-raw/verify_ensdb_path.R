# Verify the new ensdb_path parameter works from both
# annotate_from_position() and annotate_sites().

# (user library resolved by R itself; no .libPaths() override)
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
for (f in list.files("data", pattern = "\\.rda$", full.names = TRUE)) {
  load(f, envir = .GlobalEnv)
}

v113 <- file.path(Sys.getenv("LOCALAPPDATA"),
                   "R/AnnotationHub/AH119325.sqlite")
stopifnot(file.exists(v113))

sample_input <- data.frame(
  siRNA_name = c("hsiR20", "hsiR20", "hsiR22"),
  transcript = c("ENST00000426385.4",  # PMP22 CDS
                 "ENST00000335508.11", # SF3B1 CDS
                 "ENST00000380152.7"), # BRCA2
  position   = c(359L, 3136L, 500L),
  stringsAsFactors = FALSE
)

cat("=== 1. annotate_from_position with ensdb_path ===\n")
r1 <- annotate_from_position(sample_input, ensdb_path = v113)
print(r1)

cat("\n=== 2. annotate_sites(..., ensdb_path=) passthrough ===\n")
r2 <- annotate_sites(sample_input, mode = "position", ensdb_path = v113)
print(r2)

stopifnot(identical(r1, r2))
cat("\nOK: both entry points produce identical output.\n")
