library(here)
library(DESeq2)
library(pheatmap)
library(tidyverse)

fig_dir <- here("results", "figures")
tab_dir <- here("results", "tables")

alpha <- 0.05
lfc_cut <- 1
n_top <- 20

dds <- readRDS(here("data", "processed", "dds.rds"))

de <- read_csv(
  file.path(tab_dir, "DE_Primary_vs_Normal.csv"),
  col_types = cols(entrez = col_character())
)

dds_np <- dds[, dds$tissue %in% c("Normal", "Primary")]
dds_np$tissue <- droplevels(dds_np$tissue)

table(dds_np$tissue)

vsd_np <- vst(dds_np, blind = TRUE)
mat <- assay(vsd_np)

dim(mat)

top_up <- de |>
  filter(padj < alpha, log2FoldChange >= lfc_cut, !is.na(symbol)) |>
  slice_min(padj, n = n_top, with_ties = FALSE)

top_down <- de |>
  filter(padj < alpha, log2FoldChange <= -lfc_cut, !is.na(symbol)) |>
  slice_min(padj, n = n_top, with_ties = FALSE)

top_genes <- bind_rows(top_up, top_down)

top_genes$symbol

stopifnot(all(top_genes$entrez %in% rownames(mat)))

mat_top <- mat[top_genes$entrez, ]

mat_z <- t(scale(t(mat_top)))
mat_z[mat_z > 2.5] <- 2.5
mat_z[mat_z < -2.5] <- -2.5

annot_col <- data.frame(
  tissue = dds_np$tissue,
  row.names = colnames(dds_np)
)

annot_colors <- list(
  tissue = c("Normal" = "#1B9E77", Primary = "#D95F02")
)

p_heat <- pheatmap(
  mat_z, 
  scale = "none",
  breaks = seq(-2.5, 2.5, length.out = 101),
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  labels_row = top_genes$symbol,
  annotation_col = annot_col,
  annotation_colors = annot_colors,
  show_colnames = FALSE,
  border_color = NA,
  main = "Top DE genes: Primary tumor vs normal colon (VST, row z-score)"
)

png(file.path(fig_dir, "heatmap_top_genes_primary_vs_normal.png"),
    width = 2400, height = 2600, res = 300)
grid::grid.draw(p_heat$gtable)
dev.off()

colnames(mat_z)[p_heat$tree_col$order]
