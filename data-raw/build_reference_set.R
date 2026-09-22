# Build the built-in reference_set from input_100sirna/human_siRNAs_filtered.tsv
# Workflow:
#   1. Load 94 known siRNAs, normalize sequences -> 19 nt guide RNA
#   2. Run Mode A on each (using cached annodb)
#   3. Aggregate into reference_set: metadata + per-siRNA metrics +
#      tier_thresholds + species/tissue summaries
#   4. Save to data/reference_set.rda

# Paths come from _paths.R; see that file for the environment variables.
.a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
source(file.path(if (length(.a)) dirname(sub("^--file=", "", .a[1])) else "data-raw",
                 "_paths.R"))
load_pkg_source()
setwd(PKG_HOME)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(Biostrings)
})

annodb_3utr <- readRDS(data_path("hg38_3UTR_annodb.rds"))
annodb_5utr <- readRDS(data_path("hg38_5UTR_annodb.rds"))
annodb_cds  <- readRDS(data_path("hg38_CDS_annodb.rds"))
cat(sprintf("annodb 3UTR: %d | 5UTR: %d | CDS: %d transcripts\n",
            length(annodb_3utr$seqs), length(annodb_5utr$seqs),
            length(annodb_cds$seqs)))

# ---- 1. Load + normalize 94 siRNAs ------------------------------------
src <- input_path("input_100sirna", "human_siRNAs_filtered.tsv")
stopifnot(file.exists(src))
raw <- read.delim(src, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("Loaded %d siRNA entries\n", nrow(raw)))

# Guide normalisation now lives in the package itself (R/normalize_guide.R,
# exported). Prefer the reported antisense strand; otherwise derive the guide
# by reverse-complementing the sense strand.
has_anti  <- !is.na(raw$Antisense) & nzchar(raw$Antisense)
raw$guide <- ifelse(has_anti,
                    normalize_guide(raw$Antisense, input_strand = "guide"),
                    normalize_guide(raw$Sense,     input_strand = "sense"))
raw <- raw[!is.na(raw$guide) & nchar(raw$guide) >= 17, ]
raw$siRNA_name <- paste0("ref_", sprintf("%04d", raw$MIT_ID))
cat(sprintf("After cleaning: %d siRNAs (all %d nt)
",
            nrow(raw), unique(nchar(raw$guide))[1]))
cat(sprintf("Length dist: "))
print(table(nchar(raw$guide)))

# ---- 2. Run Mode A x N (one-by-one, share annodb) --------------------
out_dir <- "data-raw/reference_set_build"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

t0 <- Sys.time()
all_sites <- list()
feat_dbs <- list(`3UTR` = annodb_3utr, `5UTR` = annodb_5utr, `cds` = annodb_cds)

for (i in seq_len(nrow(raw))) {
  s <- raw$siRNA_name[i]
  cat(sprintf("[%d/%d] %s (gene=%s, len=%d)... ",
              i, nrow(raw), s, raw$Gene[i], nchar(raw$guide[i])))
  flush.console()
  per_feat <- list()
  ok <- TRUE
  for (feat in names(feat_dbs)) {
    h <- tryCatch(
      annotate_sites(
        input  = list(siRNA_name = s, sequence = raw$guide[i]),
        mode   = "sequence", annodb = feat_dbs[[feat]],
        feature_type = feat,
        seed_name = "mer7m8", max_mismatch_full = 3L
      ),
      error = function(e) { cat(sprintf("FAILED %s: %s\n",
                                          feat, conditionMessage(e))); NULL }
    )
    if (!is.null(h) && nrow(h) > 0) per_feat[[feat]] <- h
  }
  hits <- if (length(per_feat) > 0) dplyr::bind_rows(per_feat) else NULL
  if (!is.null(hits) && nrow(hits) > 0) {
    cat(sprintf("%d hits (3UTR=%d, 5UTR=%d, CDS=%d)\n",
                nrow(hits),
                if (is.null(per_feat$`3UTR`)) 0 else nrow(per_feat$`3UTR`),
                if (is.null(per_feat$`5UTR`)) 0 else nrow(per_feat$`5UTR`),
                if (is.null(per_feat$cds))    0 else nrow(per_feat$cds)))
    all_sites[[s]] <- hits
  } else {
    cat("0 hits\n")
  }
}
cat(sprintf("\nTotal time: %.1f min\n",
            as.numeric(Sys.time() - t0, units = "mins")))

sites_all <- bind_rows(all_sites)
cat(sprintf("Combined: %d sites across %d siRNAs\n",
            nrow(sites_all), length(unique(sites_all$siRNA_name))))

# ---- 3. enrich + downstream summaries --------------------------------
gtex <- readr::read_tsv(data_path("gtex_median_tpm.tsv"), show_col_types = FALSE)
ortholog <- readr::read_csv(data_path("ortholog_biomart.csv"),
                            show_col_types = FALSE)
enriched <- enrich_annotations(sites_all, expression_matrix = gtex)
# reference = "self": normalise within this cohort. Must NOT be "builtin"
# (the default), which would normalise the cohort against a previously
# saved reference_set.rda and put composite_score on a foreign scale.
ranking  <- rank_sirna(enriched, reference = "self")
tissue   <- tissue_safety(enriched, primary_tissue = "Nerve_Tibial")
species  <- species_match(enriched,
              candidates = c("mouse","rat","cyno","rhesus","rabbit","dog"),
              ortholog_table = ortholog)

# ---- 4. extract reference_set ----------------------------------------
metadata <- tibble(
  siRNA_name  = raw$siRNA_name,
  guide       = raw$guide,
  guide_len   = nchar(raw$guide),
  target_gene = raw$Gene,
  refseq      = raw$RefSeq,
  source      = "MIT siRNA database (filtered)"
)

dims <- names(default_weights())

# Per-siRNA gene-level summary
gene_summary <- enriched |>
  group_by(siRNA_name) |>
  summarise(
    n_critical_genes = n_distinct(gene_name[cancer_gene_v2 |
                                              AE_gene |
                                              immune_gene_immport]),
    .groups = "drop"
  )

# Per-composite-dimension thresholds (for normalizing future siRNAs)
tier_thresholds <- lapply(dims, function(d) {
  v <- ranking[[d]]
  c(min = min(v, na.rm = TRUE),
    max = max(v, na.rm = TRUE),
    q25 = quantile(v, 0.25, na.rm = TRUE, names = FALSE),
    q50 = quantile(v, 0.50, na.rm = TRUE, names = FALSE),
    q75 = quantile(v, 0.75, na.rm = TRUE, names = FALSE))
})
names(tier_thresholds) <- dims
tier_thresholds <- as.data.frame(do.call(rbind, tier_thresholds))
tier_thresholds$dimension <- rownames(tier_thresholds)
tier_thresholds <- tier_thresholds[, c("dimension","min","max","q25","q50","q75")]

reference_set <- list(
  metadata          = metadata,
  metrics           = ranking,                    # 94 x 8 raw counts + composite
  tier_thresholds   = tier_thresholds,
  species_score     = species$species_score,      # 94 x 6 species
  gene_summary      = gene_summary,
  risk_distribution = enriched |>
    dplyr::count(.data$siRNA_name, .data$risk_score, name = "n") |>
    tidyr::complete(.data$siRNA_name,
                     risk_score = sort(unique(enriched$risk_score)),
                     fill = list(n = 0L)),
  built_at          = as.character(Sys.time()),
  pkg_version       = as.character(packageVersion("TargetSureR")),
  n_siRNAs          = nrow(metadata)
)

cat("\n=== reference_set summary ===\n")
str(reference_set, max.level = 1)
cat(sprintf("Total size in memory: %.1f KB\n",
            object.size(reference_set) / 1024))

# ---- 5. save ----------------------------------------------------------
save(reference_set, file = "data/reference_set.rda", compress = "xz")
cat(sprintf("Saved: data/reference_set.rda (%.1f KB on disk)\n",
            file.info("data/reference_set.rda")$size / 1024))

# Also save raw enriched to data-raw (NOT shipped) for ad-hoc audit
write.csv(enriched, file.path(out_dir, "all_enriched.csv"), row.names = FALSE)
saveRDS(list(sites_all = sites_all, enriched = enriched,
              ranking = ranking, tissue = tissue, species = species),
        file.path(out_dir, "build_intermediates.rds"), compress = "xz")
cat(sprintf("Intermediates saved to %s/\n", out_dir))
