# End-to-end TargetSureR demo
# Input  : real 7 hsiR off-target sites from v2.tsv (956 rows)
# Engine : EnsDb v113 (Oct 2024) + Ensembl REST API fallback
# Output : demo_output/  (CSVs + 6 PDF/PNG figures)
# Requires: network access (REST API), v113 sqlite locally cached

# (user library resolved by R itself; no .libPaths() override)

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
for (f in list.files("data", pattern = "\\.rda$", full.names = TRUE)) {
  load(f, envir = .GlobalEnv)
}

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(scales)
  library(ensembldb); library(AnnotationFilter)
})

out_root  <- "demo_output"
fig_dir   <- file.path(out_root, "figures")
csv_dir   <- file.path(out_root, "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)

hr <- function(x) cat("\n\n", strrep("=", 70), "\n", x, "\n",
                      strrep("=", 70), "\n\n", sep = "")

save_plot <- function(p, name, w = 8, h = 6) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  cat("  saved:", name, ".pdf + .png\n", sep = "")
}

theme_tsr <- function() {
  theme_classic(base_size = 13) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(color = "#555555"),
          panel.grid.major.y = element_line(color = "#eeeeee"),
          plot.background = element_rect(fill = "white", color = NA))
}

# ========================================================================
# Step 0. Load inputs
# ========================================================================
hr("STEP 0. Loading inputs")

v2_path <- file.path(Sys.getenv("TARGETSURER_INPUT", unset = "../siRNA_prediction_model"),
                     "siRNA_offtarget_annotated_final_v2.tsv")
v2 <- read.delim(v2_path, sep = "\t", stringsAsFactors = FALSE)
v2$position <- as.integer(v2$position)
cat("v2.tsv: ", nrow(v2), "rows,",
    n_distinct(v2$siRNA_name), "siRNAs,",
    n_distinct(v2$transcript), "transcripts\n")

minimal <- v2[, c("siRNA_name", "transcript", "position")]
cat("minimal input:", nrow(minimal), "rows x", ncol(minimal), "cols\n")

v113 <- file.path(Sys.getenv("LOCALAPPDATA"),
                   "R/AnnotationHub/AH119325.sqlite")
stopifnot(file.exists(v113))
cat("EnsDb v113:", v113, " (", round(file.info(v113)$size/1024/1024, 1),
    "MB)\n", sep = "")

# ========================================================================
# Step 1. annotate_sites (Mode B)
# ========================================================================
hr("STEP 1. annotate_sites() - Mode B + REST fallback")
t0 <- Sys.time()
sites <- annotate_sites(minimal, mode = "position", ensdb_path = v113,
                        rest_fallback = TRUE, with_lookup_source = TRUE)
cat(sprintf("done in %.2f s\n", as.numeric(Sys.time()-t0, units="secs")))
cat("output: ", nrow(sites), "rows x", ncol(sites), "cols\n")
cat("lookup_source distribution:\n"); print(table(sites$lookup_source))
cat(sprintf("unresolved (NA gene_name): %d / %d (%.1f%%)\n",
            sum(is.na(sites$gene_name)), nrow(sites),
            100*mean(is.na(sites$gene_name))))
# Re-introduce match_type from v2.tsv (we have it from their pipeline)
sites$match_type <- v2$match_type
sites$risk_score <- compute_risk_score(sites$region, sites$match_type)
cat("match_type distribution:\n"); print(table(sites$match_type))
cat("region distribution:\n"); print(table(sites$region))
write.csv(sites, file.path(csv_dir, "01_sites_core.csv"), row.names = FALSE)

# ========================================================================
# Step 2. enrich_annotations
# ========================================================================
hr("STEP 2. enrich_annotations()")
t0 <- Sys.time()
sites_for_enrich <- sites[, core_columns()]
enriched <- enrich_annotations(sites_for_enrich)
cat(sprintf("done in %.2f s\n", as.numeric(Sys.time()-t0, units="secs")))
cat("output: ", nrow(enriched), "rows x", ncol(enriched), "cols\n")
write.csv(enriched, file.path(csv_dir, "02_sites_enriched.csv"), row.names = FALSE)

# ========================================================================
# Step 3. rank_sirna
# ========================================================================
hr("STEP 3. rank_sirna()")
ranking <- rank_sirna(enriched,
                      target_sirnas = c("hsiR20","hsiR21","hsiR22",
                                         "hsiR23","hsiR25","hsiR30","hsiR31"))
print(ranking[, c("siRNA_name","composite_score","rank","risk_tier",
                   "high_risk_critical","critical_gene_count",
                   "full_comp_count","cds_count","total_offtargets")])
write.csv(ranking, file.path(csv_dir, "03_sirna_ranking.csv"), row.names = FALSE)

# ========================================================================
# Step 4. tissue_safety
# ========================================================================
hr("STEP 4. tissue_safety()")
tissue <- tissue_safety(enriched, primary_tissue = "Nerve_Tibial")
cat("top 10 genes by bio_attention:\n")
print(head(tissue[, c("gene_name","n_hits","primary_tpm","max_risk_score",
                        "bio_attention","is_key_gene")], 10))
write.csv(tissue, file.path(csv_dir, "04_gene_tissue_safety.csv"), row.names = FALSE)

# ========================================================================
# Step 5. species_match
# ========================================================================
hr("STEP 5. species_match()")
species <- species_match(sites_for_enrich,
                         candidates = c("mouse","rat","cyno",
                                        "rhesus","rabbit","dog"))
cat("species score matrix:\n"); print(species$species_score)
cat("\nrecommendation:\n")
cat(paste0("  ", species$recommendation, collapse = "\n"), "\n")
write.csv(species$species_score, file.path(csv_dir, "05_species_score.csv"),
          row.names = FALSE)
write.csv(species$site_detail,  file.path(csv_dir, "05_species_site_detail.csv"),
          row.names = FALSE)

# ========================================================================
# Step 6. Figures
# ========================================================================
hr("STEP 6. Generating figures")

# -- Figure 1: siRNA ranking bar chart
tier_pal <- rev(RColorBrewer::brewer.pal(4, "RdYlGn"))
names(tier_pal) <- c("Low", "Medium", "High", "Critical")
p1 <- ranking |>
  mutate(siRNA_name = factor(siRNA_name, levels = siRNA_name[order(composite_score)])) |>
  ggplot(aes(x = siRNA_name, y = composite_score, fill = risk_tier)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", composite_score)),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = tier_pal, drop = FALSE) +
  labs(title = "siRNA off-target risk ranking",
       subtitle = sprintf("Composite score across %d dimensions (lower = safer)",
                           length(default_weights())),
       x = NULL, y = "Composite risk score",
       fill = "Risk tier") +
  theme_tsr() +
  theme(axis.text.x = element_text(angle = 0))
save_plot(p1, "01_sirna_ranking_bar", w = 9, h = 5)

# -- Figure 2: risk_score distribution per siRNA
p2 <- enriched |>
  count(siRNA_name, risk_score) |>
  ggplot(aes(x = siRNA_name, y = n, fill = factor(risk_score))) +
  geom_col(position = position_stack(reverse = TRUE)) +
  scale_fill_brewer(palette = "OrRd", direction = 1,
                     name = "risk_score\n(CDS=5, +full=+2)") +
  labs(title = "Off-target site count by risk_score",
       subtitle = "Higher risk_score = more dangerous site",
       x = NULL, y = "Off-target sites") +
  theme_tsr()
save_plot(p2, "02_risk_distribution", w = 9, h = 5)

# -- Figure 3: gene bio-attention scatter (expr x risk)
scatter_pal <- RColorBrewer::brewer.pal(3, "Set1")[c(2, 1)]
top_genes <- head(tissue, 30)
p3 <- top_genes |>
  ggplot(aes(x = score_expr, y = score_risk,
             size = n_hits, color = is_key_gene)) +
  geom_point(alpha = 0.75) +
  ggrepel::geom_text_repel(
    aes(label = gene_name), size = 3.2,
    max.overlaps = 30) +
  scale_color_manual(values = c(`FALSE` = scatter_pal[1],
                                 `TRUE`  = scatter_pal[2]),
                      labels = c("Non-critical","Critical"),
                      name = "Gene category") +
  scale_size_continuous(range = c(2, 10), name = "# hits") +
  labs(title = "Off-target gene attention",
       subtitle = "Top 30 genes by bio_attention; x = Nerve_Tibial expression, y = max risk_score",
       x = "Normalized expression (Nerve_Tibial)",
       y = "Normalized max risk_score") +
  theme_tsr()
save_plot(p3, "03_gene_attention_scatter", w = 9, h = 6.5)

# -- Figure 4: species recommendation heatmap
sp_pal <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlGn"))(100)
sp_long <- species$species_score |>
  pivot_longer(-siRNA_name, names_to = "species", values_to = "score")
p4 <- sp_long |>
  ggplot(aes(x = species, y = siRNA_name, fill = score)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", score)), size = 3.5) +
  scale_fill_gradientn(colours = sp_pal, limits = c(0, 1),
                        name = "Conservation\nscore") +
  labs(title = "Animal model recommendation",
       subtitle = "Layer-1 ortholog conservation on critical off-target genes",
       x = NULL, y = NULL) +
  theme_tsr() +
  theme(panel.grid = element_blank())
save_plot(p4, "04_species_heatmap", w = 9, h = 5)

# -- Figure 5: tissue TPM heatmap for top critical-hit genes
tissues_of_interest <- c("Nerve_Tibial", "Brain_Cortex", "Liver",
                          "Heart_Left_Ventricle", "Muscle_Skeletal",
                          "Kidney_Cortex", "Lung", "Whole_Blood",
                          "Spleen", "Adipose_Subcutaneous")
top_gene_names <- head(tissue$gene_name, 20)
heat_df <- enriched |>
  dplyr::filter(gene_name %in% top_gene_names) |>
  group_by(gene_name) |>
  summarise(across(all_of(intersect(tissues_of_interest, names(enriched))),
                    ~ mean(.x, na.rm = TRUE))) |>
  pivot_longer(-gene_name, names_to = "tissue", values_to = "tpm")
p5 <- heat_df |>
  mutate(gene_name = factor(gene_name, levels = rev(top_gene_names))) |>
  ggplot(aes(x = tissue, y = gene_name, fill = log1p(tpm))) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_distiller(palette = "YlGnBu", direction = 1,
                        name = "log1p(TPM)") +
  labs(title = "Tissue expression of top off-target genes",
       subtitle = "Top 20 genes by bio_attention",
       x = NULL, y = NULL) +
  theme_tsr() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        panel.grid = element_blank())
save_plot(p5, "05_tissue_heatmap", w = 10, h = 6)

# -- Figure 6: siRNA x critical-gene bubble chart
critical_hits <- enriched |>
  dplyr::filter(cancer_gene_v2 | AE_gene | immune_gene_immport) |>
  count(siRNA_name, gene_name, name = "n_hits") |>
  left_join(tissue[, c("gene_name", "bio_attention")], by = "gene_name")
p6 <- critical_hits |>
  mutate(gene_name = factor(gene_name,
                             levels = unique(gene_name[order(bio_attention,
                                                              decreasing = FALSE)]))) |>
  ggplot(aes(x = siRNA_name, y = gene_name,
             size = n_hits, color = bio_attention)) +
  geom_point() +
  scale_size_continuous(range = c(2, 10), name = "# hits") +
  scale_color_distiller(palette = "YlOrRd", direction = 1,
                         name = "bio_attention") +
  labs(title = "siRNA x critical-gene off-target map",
       subtitle = "Only sites flagged as AE / cancer / immune",
       x = NULL, y = NULL) +
  theme_tsr()
save_plot(p6, "06_sirna_gene_bubble", w = 9, h = 7)

# ========================================================================
# Step 7. generate_report
# ========================================================================
hr("STEP 7. generate_report() bundle")
paths <- generate_report(enriched,
                          output_dir      = file.path(out_root, "report_bundle"),
                          target_sirnas   = c("hsiR22"),
                          primary_tissue  = "Nerve_Tibial",
                          species_result  = species)
cat("\nbundle written:\n")
for (p in unlist(paths)) cat("  ", p, "\n")

# ========================================================================
# Summary
# ========================================================================
hr("DONE")
cat("All outputs under:", normalizePath(out_root), "\n\n")

cat("== tables ==\n")
for (f in list.files(csv_dir, full.names = TRUE))
  cat(sprintf("  %-45s %6.1f KB\n", basename(f), file.info(f)$size/1024))

cat("\n== figures ==\n")
for (f in list.files(fig_dir, pattern = "\\.pdf$", full.names = TRUE))
  cat(sprintf("  %-45s %6.1f KB\n", basename(f), file.info(f)$size/1024))

cat("\n== report bundle ==\n")
for (f in list.files(file.path(out_root, "report_bundle"), full.names = TRUE))
  cat(sprintf("  %-45s %6.1f KB\n", basename(f), file.info(f)$size/1024))
