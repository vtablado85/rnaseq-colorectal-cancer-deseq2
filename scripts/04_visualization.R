library(here)
library(ggrepel)
library(tidyverse)

fig_dir <- here("results", "figures")
tab_dir <- here("results", "tables")

alpha <- 0.05
lfc_cut <- 1

de <- read_csv (
  file.path(tab_dir, "DE_Primary_vs_Normal.csv"),
  col_types = cols(entrez = col_character())
)

glimpse(de)

de <- de |>
  mutate(
    status = case_when(
      padj < alpha & log2FoldChange >=  lfc_cut ~ "Up",
      padj < alpha & log2FoldChange <= -lfc_cut ~ "Down",
      TRUE ~ "NS"
    ),
    status = factor (status, levels = c("Down", "NS", "Up")),
    neg_log10_padj = -log10(padj)
  ) |>
  arrange(desc(status == "NS"))

table(de$status)

n_lab <- 10

top_up <- de |>
  filter(status == "Up", !is.na(symbol)) |>
  slice_min(padj, n = n_lab, with_ties = FALSE)

top_down <- de |>
  filter(status == "Down", !is.na(symbol)) |>
  slice_min(padj, n = n_lab, with_ties = FALSE)

top_genes <- bind_rows(top_up, top_down)
top_genes$symbol

n_up <- sum(de$status == "Up")
n_down <- sum(de$status == "Down")

p_volcano <- ggplot(de, aes(x = log2FoldChange, y = neg_log10_padj, colour = status)) + 
  geom_point(alpha = 0.5, size = 1) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", colour = "grey40") + 
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", colour = "grey40") +
  geom_text_repel(
    data = top_genes,  aes(label = symbol),
    colour = "black", size = 3, max.overlaps = Inf,
    box.padding = 0.4, segment.colour = "grey50"
  ) + 
  scale_color_manual(values = c(Down = "#2166AC", NS = "grey75", Up = "#B2182B")) + 
  labs(
    title =  "Primary tumor vs normal colon",
    subtitle = paste0(n_up, " up and ", n_down, " down (padj < ", alpha,
                          ", |log2FC| >= ", lfc_cut, ")"),
    x = "Log2 fold change (apeglm shrunken)",
    y = "-log10(adjusted p-value)",
    colour = NULL
  ) + 
  theme_bw()

p_volcano

ggsave(file.path(fig_dir, "volcano_primary_vs_normal.png"), p_volcano,
       width = 8, height = 6, dpi = 300)

p_ma <- ggplot(de, aes(x = baseMean, y = log2FoldChange, colour = status)) + 
  geom_point(alpha = 0.5, size = 1) + 
  geom_hline(yintercept = 0, colour = "black") + 
  geom_hline(yintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", colour = "grey40") + 
  scale_x_log10(labels = scales::label_comma()) +
  scale_color_manual(values = c(Down = "#2166AC", NS = "grey75", Up = "#B2182B")) + 
  labs(
    title = "MA plot: Primary tumor vs normal colon",
    subtitle = paste0(n_up, " up and ", n_down, " down (padj < ", alpha, 
                          ", |log2FC| >= ", lfc_cut, ")"),
    x = "Mean of normalized counts (log10 scale)",
    y = "Log2 fold change (apeglm shrunken)",
    colour = NULL
  ) + 
  theme_bw()

p_ma

ggsave(file.path(fig_dir, "ma_primary_vs_normal.png"), p_ma,
       width = 8, height = 6, dpi = 300)

