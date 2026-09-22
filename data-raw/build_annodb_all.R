# Build and cache the three feature-class annodbs used by
# build_reference_set.R. Thin wrapper around the exported build_annodb();
# see vignette("01_quickstart") step 3 for the user-facing version.
#
# Paths come from _paths.R; set TARGETSURER_DATA to where your EnsDb lives.
.a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
source(file.path(if (length(.a)) dirname(sub("^--file=", "", .a[1])) else "data-raw",
                 "_paths.R"))
load_pkg_source()

edb <- ensembldb::EnsDb(ensdb_path())
cat("EnsDb release:", ensembldb::ensemblVersion(edb), "\n")

for (feat in c("3UTR", "5UTR", "CDS")) {
  cache <- data_path(sprintf("hg38_%s_annodb.rds", feat))
  if (file.exists(cache)) {
    cat(sprintf("[skip] %s already cached\n", basename(cache)))
    next
  }
  ad <- build_annodb(feat, edb = edb)
  saveRDS(ad, cache, compress = "xz")
  cat(sprintf("  saved %s (%.1f MB)\n", basename(cache),
              file.info(cache)$size / 1024 / 1024))
}
