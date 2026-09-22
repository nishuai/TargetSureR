# Shared path resolution for the scripts in data-raw/.
#
# These scripts are development and provenance records, excluded from the
# installed package. They read inputs that are not redistributed with the
# package, so every location is resolved from an environment variable with a
# documented fallback rather than hard-coded.
#
# Set these before running any script in this directory:
#
#   TARGETSURER_HOME   package source root
#                      default: the parent of this file's directory
#   TARGETSURER_DATA   prepared reference data (EnsDb, annodb, GTEx, orthologs)
#                      default: ~/TargetSureR_data
#   TARGETSURER_INPUT  curated source tables (siRNA cohort, gene lists)
#                      default: <TARGETSURER_HOME>/../siRNA_prediction_model
#
# Example:
#   Sys.setenv(TARGETSURER_DATA = "~/TargetSureR_data")
#   source("data-raw/_paths.R")

.this_dir <- function() {
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a))
    return(dirname(normalizePath(sub("^--file=", "", a[1]), mustWork = FALSE)))
  for (i in rev(seq_len(sys.nframe()))) {
    f <- sys.frame(i)$ofile
    if (!is.null(f)) return(dirname(normalizePath(f, mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

.env_or <- function(var, default) {
  v <- Sys.getenv(var, unset = "")
  path.expand(if (nzchar(v)) v else default)
}

PKG_HOME <- .env_or("TARGETSURER_HOME", dirname(.this_dir()))
if (!dir.exists(file.path(PKG_HOME, "R")))
  stop("Package source not found at '", PKG_HOME,
       "'. Set TARGETSURER_HOME to the TargetSureR source root.", call. = FALSE)

DATA_DIR  <- .env_or("TARGETSURER_DATA",  "~/TargetSureR_data")
INPUT_DIR <- .env_or("TARGETSURER_INPUT",
                     file.path(dirname(PKG_HOME), "siRNA_prediction_model"))

pkg_path   <- function(...) file.path(PKG_HOME,  ...)
data_path  <- function(...) file.path(DATA_DIR,  ...)
input_path <- function(...) file.path(INPUT_DIR, ...)

# The EnsDb snapshot. AnnotationHub caches into a platform-specific directory,
# so prefer an explicit copy under TARGETSURER_DATA; fall back to the hub cache.
ensdb_path <- function(record = "AH119325") {
  p <- data_path("EnsDb.Hsapiens.v113.sqlite")
  if (file.exists(p)) return(p)
  hub <- tryCatch(tools::R_user_dir("AnnotationHub", "cache"),
                  error = function(e) NULL)
  if (!is.null(hub)) {
    q <- file.path(hub, paste0(record, ".sqlite"))
    if (file.exists(q)) return(q)
  }
  stop("EnsDb not found. Place it at '", p,
       "' or set TARGETSURER_DATA. See vignette('01_quickstart') step 2.",
       call. = FALSE)
}

# Load the package from source, mirroring how these scripts were run.
load_pkg_source <- function() {
  for (f in list.files(pkg_path("R"), pattern = "[.]R$", full.names = TRUE))
    source(f)
  for (f in list.files(pkg_path("data"), pattern = "[.]rda$", full.names = TRUE))
    load(f, envir = .GlobalEnv)
  invisible(TRUE)
}

invisible(NULL)
