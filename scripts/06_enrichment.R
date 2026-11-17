library(here)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(tidyverse)

fig_dir <- here("results", "figures")
tab_dir <- here("results", "tables")

alpha <- 0.05
lfc_cut <- 1

de <- read_csv(
  file.path(tab_dir, "DE_Primary_vs_Normal.csv"),
  col_types = cols(entrez = col_character())
)

sig <- de |>
  filter(padj < alpha, abs(log2FoldChange) >= lfc_cut)

genes_up <- sig |>
  filter(log2FoldChange > 0) |> pull(entrez)
genes_down <- sig |>
  filter(log2FoldChange < 0) |> pull(entrez)
universe <- de$entrez

length(genes_up)
length(genes_down)
length(universe)

go_up <- enrichGO(
  gene = genes_up,
  universe = universe,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff =  0.05,
  readable = TRUE
)

go_down <- enrichGO(
  gene = genes_down,
  universe = universe,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff =  0.05,
  readable = TRUE
)

nrow(as.data.frame(go_up))
nrow(as.data.frame(go_down))

go_up_s <- clusterProfiler::simplify(go_up, cutoff = 0.7, by = "p.adjust", select_fun = min)
go_down_s <- clusterProfiler::simplify(go_down, cutoff = 0.7, by = "p.adjust", select_fun = min)

write_csv(as.data.frame(go_up_s), file.path(tab_dir, "GO_BP_up_Primary_vs_Normal.csv"))
write_csv(as.data.frame(go_down_s), file.path(tab_dir, "GO_BP_down_Primary_vs_Normal.csv"))

as.data.frame(go_up_s) |>
  select(Description, Count, p.adjust) |>
  head(10)

as.data.frame(go_down_s) |>
  select(Description, Count, p.adjust) |>
  head(10)

p_up <- dotplot(go_up_s, showCategory = 15, label_format = 40) +
  labs(title = "GO Biological Process : Genes up in primary tumor")

p_down <- dotplot(go_down_s, showCategory = 15, label_format = 40) +
  labs(title = "GO Biological Process : Genes down in primary tumor")

p_up
p_down

ggsave(file.path(fig_dir, "go_bp_up_primary_vs_normal.png"), p_up, width = 8, height = 9, dpi = 300)
ggsave(file.path(fig_dir, "go_bp_down_primary_vs_normal.png"), p_down, width = 8, height = 9, dpi = 300)
