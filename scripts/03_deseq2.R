library(here)
library(DESeq2)
library(apeglm)
library(tidyverse)

fig_dir <- here("results", "figures")
tab_dir <- here("results", "tables")

alpha <- 0.05
lfc_cut <- 1

dds <- readRDS(here("data", "processed", "dds.rds"))
design(dds)
levels(dds$tissue)

dds <- DESeq(dds)
resultsNames(dds)

png(file.path(fig_dir, "dispersion_estimates.png"), width = 1800, height = 1400, res = 300)
plotDispEsts(dds)
dev.off()

get_res <- function(dds_obj, coef_name) {
  res <- results(dds_obj, name = coef_name, alpha = alpha)
  lfcShrink(dds_obj, coef = coef_name, res = res, type = "apeglm")
}

res_pri_nor <- get_res(dds, "tissue_Primary_vs_Normal")
res_met_nor <- get_res(dds, "tissue_Metastasis_vs_Normal")

dds_p <- readRDS(here("data", "processed", "dds.rds"))
dds_p$tissue <- relevel(dds_p$tissue, ref = "Primary")
dds_p <- DESeq(dds_p)
resultsNames(dds_p)

res_met_pri <- get_res(dds_p, "tissue_Metastasis_vs_Primary")

add_symbols <- function(res){
  as.data.frame(res) |>
    rownames_to_column("entrez") |>
    mutate(symbol = unname(AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = entrez, column = "SYMBOL",
      keytype = "ENTREZID", multiVals = "first"))) |>
    arrange(padj)
}

de_pri_nor <- add_symbols(res_pri_nor)
de_met_nor <- add_symbols(res_met_nor)
de_met_pri <- add_symbols(res_met_pri)

write_csv(de_pri_nor, file.path(tab_dir, "DE_Primary_vs_Normal.csv"))
write_csv(de_met_nor, file.path(tab_dir, "DE_Metastasis_vs_Normal.csv"))
write_csv(de_met_pri, file.path(tab_dir, "DE_Metastasis_vs_Primary.csv"))

count_de <- function(df) {
  df |>
    filter(!is.na(padj), padj < alpha, abs(log2FoldChange) >= lfc_cut) |>
    summarise(up = sum(log2FoldChange > 0), down = sum(log2FoldChange < 0))
}

bind_rows(
  Primary_vs_Normal = count_de(de_pri_nor),
  Metastasis_vs_Normal = count_de(de_met_nor),
  Metastasis_vs_Primary = count_de(de_met_pri),
  .id = "contrast"
  )

head(de_met_nor, 15)
de_met_nor |> filter(symbol %in% c("ALB", "APOA1", "TTR", "CDX2", "KRT20"))


