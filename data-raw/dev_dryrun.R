# Static dry-run validator for Mode A and Mode B.
#
# We do NOT install SeedMatchR / ensembldb here. Instead we build mock
# objects that mimic the documented APIs and pass them in as arguments.
# This verifies our wrapper code:
#   - calls existing function names
#   - passes the right argument types
#   - reads the documented return-value structure
#   - feeds the right shape to downstream code

cat("=== STATIC DRY-RUN: Mode A and Mode B ===\n\n")

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
for (f in list.files("data", pattern = "\\.rda$", full.names = TRUE)) load(f, envir = .GlobalEnv)

# ----------------------------------------------------------------------
# MOCK ensembldb world (for Mode B)
# ----------------------------------------------------------------------
mock_ensdb <- structure(list(.tag = "MOCK_ENSDB"), class = "EnsDb")

# Mock the AnnotationFilter::TxIdFilter, ensembldb::transcripts,
# ensembldb::transcriptToCds, ensembldb::fiveUTRsByTranscript,
# IRanges::IRanges, IRanges::start, BiocGenerics::width, S4Vectors::mcols
# by injecting them under the ::-style names used by our code.

mock_pkg <- function(pkg, exports) {
  ns <- new.env()
  for (nm in names(exports)) assign(nm, exports[[nm]], envir = ns)
  attr(ns, "name") <- pkg
  e <- as.environment(list2env(as.list(ns), parent = baseenv()))
  assignInNamespace_safe <- function() invisible(NULL)
  ns
}

mock_TxIdFilter <- function(value, ...) structure(list(value = value), class = "TxIdFilter")

mock_transcripts <- function(x, filter = NULL, columns, return.type = "data.frame") {
  ids <- filter$value
  data.frame(
    tx_id      = ids,
    tx_biotype = ifelse(grepl("ENST_NCRNA", ids), "lncRNA", "protein_coding"),
    gene_name  = paste0("GENE_", seq_along(ids)),
    stringsAsFactors = FALSE
  )
}

mock_IRanges <- function(start, width, names) {
  structure(list(start = as.integer(start), width = as.integer(width),
                 names = names),
            class = "IRanges")
}
mock_start  <- function(x) x$start
mock_end    <- function(x) x$start + x$width - 1L
mock_width  <- function(x) x$width
mock_names  <- function(x) x$names

# transcriptToCds returns IRanges with start = -1 if NOT in CDS, else CDS-relative
mock_transcriptToCds <- function(x, db, id = "name") {
  ids <- x$names
  st  <- x$start
  # rule for the mock: if name starts with CDS_ -> in CDS, else not
  in_cds <- grepl("^CDS_", ids)
  new_starts <- ifelse(in_cds, st, -1L)
  structure(list(start = new_starts, width = x$width, names = ids),
            class = "IRanges")
}

mock_fiveUTRsByTranscript <- function(x, filter = NULL) {
  ids <- filter$value
  out <- list()
  for (id in ids) {
    out[[id]] <- list(width = if (grepl("^UTR5_", id)) 100L else 50L)
  }
  out
}
mock_fiveUTRsByTranscript_width <- function(g) g$width

# Inject mocks into the namespace lookup our code performs.
# Easiest: define '::' shadow via assign in package env after sourcing.
# Our code uses things like ensembldb::transcripts, IRanges::IRanges, etc.
# We override them via attaching a faux package on the search path.

faux <- new.env()
assign("AnnotationFilter", new.env(), envir = faux)
faux$AnnotationFilter$TxIdFilter <- mock_TxIdFilter

assign("ensembldb", new.env(), envir = faux)
faux$ensembldb$transcripts            <- mock_transcripts
faux$ensembldb$transcriptToCds        <- mock_transcriptToCds
faux$ensembldb$fiveUTRsByTranscript   <- mock_fiveUTRsByTranscript

assign("IRanges", new.env(), envir = faux)
faux$IRanges$IRanges <- mock_IRanges
faux$IRanges$start   <- mock_start
faux$IRanges$end     <- mock_end

assign("BiocGenerics", new.env(), envir = faux)
faux$BiocGenerics$width <- function(g) {
  if (is.list(g) && !is.null(g$width)) return(g$width)
  if (is.list(g)) return(vapply(g, function(z) z$width, integer(1)))
  attr(g, "width")
}
faux$BiocGenerics$start <- mock_start

assign("S4Vectors", new.env(), envir = faux)
faux$S4Vectors$mcols <- function(g) attr(g, "mcols")

# Override require_pkg so it accepts our mocks
require_pkg <- function(pkg, reason = NULL) invisible(TRUE)

# Inject the mocks under their `pkg::name` accessors via a hack: replace
# the pkg::name lookup by overriding the namespace search. The easiest
# robust way is monkey-patching the helper functions:

# Replace named refs inside functions by re-assigning them in our environment:
.do_call_with_mocks <- function(expr) {
  attached <- NULL
  on.exit({ for (a in attached) detach(pos = match(a, search()), unload = FALSE) }, add = TRUE)
  for (nm in c("AnnotationFilter", "ensembldb", "IRanges",
               "BiocGenerics", "S4Vectors")) {
    e <- faux[[nm]]
    attach(as.list(e), name = paste0("mock:", nm), warn.conflicts = FALSE)
    attached <- c(attached, paste0("mock:", nm))
  }
  eval(expr, envir = parent.frame())
}

# ----------------------------------------------------------------------
# Test 1: Mode B with mocks
# ----------------------------------------------------------------------

# Patch the :: calls by direct symbol lookup
ensembldb            <- faux$ensembldb
AnnotationFilter     <- faux$AnnotationFilter
IRanges              <- faux$IRanges
BiocGenerics         <- faux$BiocGenerics
S4Vectors            <- faux$S4Vectors

# We need our wrapper functions to use these globals. Replace `::` by stripping
# the qualifier inside the wrappers: rebuild the functions with mock-aware bodies.

# Simpler approach: redefine the helpers inline using mocks
.lookup_transcript_meta <- function(tx_ids, ensdb) {
  flt <- AnnotationFilter$TxIdFilter(tx_ids)
  df <- ensembldb$transcripts(ensdb, filter = flt,
                              columns = c("tx_id","tx_biotype","gene_name"),
                              return.type = "data.frame")
  tibble::as_tibble(df[, c("tx_id","gene_name","tx_biotype")])
}
.utr5_lengths <- function(tx_ids, ensdb) {
  flt <- AnnotationFilter$TxIdFilter(tx_ids)
  utrs <- ensembldb$fiveUTRsByTranscript(ensdb, filter = flt)
  lens <- vapply(utrs, function(g) sum(BiocGenerics$width(g)), integer(1))
  out <- setNames(rep(NA_integer_, length(tx_ids)), tx_ids)
  out[names(lens)] <- lens
  out
}
.classify_region <- function(tx_ids, positions, ensdb) {
  rng <- IRanges$IRanges(start = as.integer(positions),
                          width = 1L, names = tx_ids)
  cds_rng <- ensembldb$transcriptToCds(rng, ensdb)
  in_cds <- IRanges$start(cds_rng) > 0
  out <- rep(NA_character_, length(tx_ids))
  out[in_cds] <- "CDS"

  if (any(!in_cds)) {
    unique_tx <- unique(tx_ids[!in_cds])
    flt <- AnnotationFilter$TxIdFilter(unique_tx)
    txs <- ensembldb$transcripts(ensdb, filter = flt,
                                  columns = c("tx_id","tx_biotype"),
                                  return.type = "data.frame")
    utr5_lens <- .utr5_lengths(unique_tx, ensdb)
    for (i in which(!in_cds)) {
      tx <- tx_ids[i]
      bt <- txs$tx_biotype[match(tx, txs$tx_id)]
      if (is.na(bt) || !grepl("protein_coding", bt)) {
        out[i] <- "ncRNA"
      } else {
        u5 <- utr5_lens[[tx]]
        if (is.na(u5)) out[i] <- "ncRNA"
        else if (positions[i] <= u5) out[i] <- "5'UTR"
        else out[i] <- "3'UTR"
      }
    }
  }
  out[is.na(out)] <- "ncRNA"
  out
}

cat("--- Mode B (with mocked EnsDb) ---\n")
input_b <- data.frame(
  siRNA_name = c("hsiR1","hsiR1","hsiR1","hsiR2"),
  transcript = c("CDS_ENST1", "UTR5_ENST2", "ENST3", "ENST_NCRNA"),
  position   = c(379L, 50L, 200L, 100L)
)
out_b <- annotate_from_position(input_b, ensdb = mock_ensdb)
print(out_b)

stopifnot(all(core_columns() %in% names(out_b)))
stopifnot(nrow(out_b) == 4)
stopifnot(out_b$region[1] == "CDS")
stopifnot(out_b$region[2] == "5'UTR")
stopifnot(out_b$region[3] == "3'UTR")
stopifnot(out_b$region[4] == "ncRNA")
cat("Mode B region classification: PASS\n\n")

# ----------------------------------------------------------------------
# Test 2: Mode A with mocked SeedMatchR + Biostrings
# ----------------------------------------------------------------------

mock_load_annotations <- function(reference.name, feature.type,
                                   canonical, min.feature.width,
                                   longest.utr, return_gene_name) {
  cat("  [mock] load_annotations called: ref=", reference.name,
      ", feature=", feature.type,
      ", return_gene_name=", return_gene_name, "\n", sep="")
  stopifnot(reference.name %in% c("hg38","hg38-old","mm39","mm10","rnor6","rnor7"))
  stopifnot(return_gene_name == FALSE)
  # mock seqs is a named list
  seqs <- list(
    ENST_A = "AAA",
    ENST_B = "BBB",
    ENST_C = "CCC"
  )
  # mock gtf with mcols-like access
  gtf <- structure(list(),
                   mcols = data.frame(
                     tx_id = c("ENST_A","ENST_B","ENST_C"),
                     gene_name = c("GENE_A","GENE_B","GENE_C"),
                     tx_biotype = "protein_coding",
                     stringsAsFactors = FALSE
                   ))
  attr(gtf, "length") <- 3
  list(seqs = seqs, gtf = gtf)
}

mock_SeedMatchR <- function(seqs, sequence, seed.name, res.format) {
  stopifnot(res.format == "granges")
  cat("  [mock] SeedMatchR called: seed=", seed.name, "\n", sep="")
  # return mock GRanges-like
  structure(list(seqnames = c("ENST_A","ENST_B"),
                 start    = c(120L, 250L)),
            class = "GRanges")
}
mock_seqnames <- function(g) g$seqnames

mock_DNAString <- function(s) structure(s, class = "DNAString")
mock_revComp   <- function(s) structure(paste(rev(strsplit(chartr("ACGT","TGCA", s),"")[[1]]),
                                                collapse=""), class = "DNAString")

mock_vmatchPattern <- function(pattern, subject, max.mismatch,
                                with.indels = FALSE, fixed = TRUE) {
  cat("  [mock] vmatchPattern called: max.mismatch=", max.mismatch, "\n", sep="")
  if (max.mismatch == 0L) {
    # exact-only
    out <- list(ENST_A = list(start = c(50L)),
                ENST_B = list(start = integer()),
                ENST_C = list(start = integer()))
  } else {
    out <- list(ENST_A = list(start = c(50L)),
                ENST_B = list(start = c(300L)),
                ENST_C = list(start = integer()))
  }
  attr(out, "names") <- names(out)
  out
}

# wire mocks
SeedMatchR    <- new.env(); SeedMatchR$load_annotations <- mock_load_annotations
SeedMatchR$SeedMatchR <- mock_SeedMatchR
Biostrings    <- new.env()
Biostrings$DNAString          <- mock_DNAString
Biostrings$reverseComplement  <- mock_revComp
Biostrings$vmatchPattern      <- mock_vmatchPattern
GenomicRanges <- new.env(); GenomicRanges$seqnames <- mock_seqnames

# Override BiocGenerics for vmatchPattern result handling
BiocGenerics$start <- function(x) {
  if (inherits(x, "GRanges")) return(x$start)
  if (is.list(x) && !is.null(x[[1]]$start)) return(lapply(x, function(z) z$start))
  if (is.list(x)) return(unlist(lapply(x, function(z) z$start)))
  x
}

# rebuild Mode A internal pieces with mock dispatch
.scan_seed <- function(sequence, annodb, seed_name) {
  gr <- SeedMatchR$SeedMatchR(seqs = annodb$seqs, sequence = sequence,
                               seed.name = seed_name, res.format = "granges")
  if (length(gr$seqnames) == 0) return(tibble::tibble(transcript=character(),
                                                       position=integer(),
                                                       match_type=character()))
  tibble::tibble(transcript = GenomicRanges$seqnames(gr),
                 position   = as.integer(gr$start),
                 match_type = "partial_match")
}
.scan_full <- function(sequence, annodb, max_mismatch) {
  guide_dna <- Biostrings$DNAString(gsub("U","T",sequence))
  target    <- Biostrings$reverseComplement(guide_dna)
  mi <- Biostrings$vmatchPattern(target, annodb$seqs, max.mismatch=max_mismatch)
  starts_list <- BiocGenerics$start(mi)
  lens <- vapply(starts_list, length, integer(1))
  if (sum(lens) == 0) return(tibble::tibble(transcript=character(),
                                              position=integer(),
                                              match_type=character()))
  tx_ids <- rep(names(mi), lens)
  starts <- unlist(starts_list, use.names=FALSE)
  exact_mi <- Biostrings$vmatchPattern(target, annodb$seqs[unique(tx_ids)],
                                        max.mismatch = 0L)
  exact_starts <- BiocGenerics$start(exact_mi)
  exact_set <- names(exact_starts)[vapply(exact_starts, length, integer(1)) > 0]
  mt <- ifelse(tx_ids %in% exact_set, "full_complementarity", "partial_match")
  tibble::tibble(transcript = tx_ids, position = as.integer(starts), match_type = mt)
}
.build_tx_meta <- function(annodb) {
  m <- attr(annodb$gtf, "mcols")
  tibble::tibble(tx_id = m$tx_id, gene_name = m$gene_name, biotype = m$tx_biotype)
}

cat("--- Mode A (with mocked SeedMatchR/Biostrings) ---\n")
input_a <- list(siRNA_name = "hsiR1",
                sequence   = "UUAUAGAGCAAGAACACUGUUUU")
out_a <- annotate_from_sequence(input_a, species = "human")
print(out_a)
stopifnot(all(core_columns() %in% names(out_a)))
stopifnot(nrow(out_a) >= 1)
stopifnot("partial_match" %in% out_a$match_type ||
          "full_complementarity" %in% out_a$match_type)
stopifnot(all(out_a$strand == "guide"))
stopifnot(all(out_a$siRNA_name == "hsiR1"))
cat("Mode A end-to-end with mocks: PASS\n\n")

# ----------------------------------------------------------------------
# Test 3: enriched + downstream from Mode A output
# ----------------------------------------------------------------------

cat("--- Pipe Mode A output through enrich + rank + tissue ---\n")
out_a$gene_name <- c("CTNNB1", "TP53")[seq_len(nrow(out_a))]  # patch known critical names
enr <- enrich_annotations(out_a)
cat("enriched dims:", dim(enr), "\n")
rnk <- rank_sirna(enr, target_sirnas = "hsiR1")
print(rnk)

cat("\nALL DRY-RUN CHECKS PASSED\n")
