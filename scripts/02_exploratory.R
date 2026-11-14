# 02_exploratory.R
# Exploratory analysis of GSE50760: quality control, gene filtering, PCA and
# sample distances.
# Input : raw count matrix + sample_metadata.csv (created by 01_download.R)
# Output: QC table, figures and a filtered DESeq2 object (dds.rds)

# Load packages
library(here)
library(tidyverse)
library(DESeq2)
library(pheatmap)

# Paths: input files and output folders (here() starts at the project root)
counts_path <- here("data", "raw", "GSE50760_raw_counts_GRCh38.p13_NCBI.tsv.gz")
meta_path <- here("data", "processed", "sample_metadata.csv")
fig_dir <- here("results", "figures")
tab_dir <- here("results", "tables")

# Create the output folders if they do not exist yet
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load the data ---------------------------------------------------------
# Read the raw counts (GeneID as text, one numeric column per sample)
count_df <- read_tsv(
  counts_path,
  col_types = cols(GeneID = col_character(), .default = col_double())
)

# Convert to a matrix: genes in rows, samples in columns
count_mat <- count_df |>
  column_to_rownames("GeneID") |>
  as.matrix()

# Safety checks: no missing values and only whole numbers (raw counts)
stopifnot(
  "Counts contain NA values" = !anyNA(count_mat),
  "Counts must be integers" = all(count_mat == round(count_mat))
)
# DESeq2 needs integer counts
storage.mode(count_mat)   <-  "integer"

# Read the sample table; the tissue order sets "Normal" as the reference
meta <- read_csv(meta_path, show_col_types =  FALSE) |> 
  mutate(
    tissue = factor(tissue, levels = c("Normal", "Primary", "Metastasis")),
    patient = factor(patient)
  )

# Put the samples in the same order as the count columns, and check it
meta <- meta[match(colnames(count_mat), meta$gsm), ]
stopifnot(
  "Sample order differs between counts and metadata" = identical(meta$gsm, colnames(count_mat))
)

# Print the size of the data (genes x samples)
message("Loaded ", nrow(count_mat), " genes x ", ncol(count_mat), " samples")

# 2. Quality control -------------------------------------------------------
# Total counts (library size) and number of detected genes per sample
qc <- meta |> 
  mutate(
    library_sizess = colSums(count_mat),
    genes_detected = colSums(count_mat > 0)
  )
# Save the QC table
write_csv(qc, file.path(tab_dir, "qc_summary.csv"))

# Show the library sizes in millions (they should be similar across samples)
message("Library size (millions of assigned counts:")
print(summary(qc$library_sizess / 1e6))

# Bar plot of library size per sample, coloured by tissue
p_lib <- ggplot(qc, aes(x = reorder(gsm, library_sizess),
                        y = library_sizess / 1e6, fill = tissue)) + 
  geom_col() + 
  labs(x = "Sample", y = "Assigned counts (millions)", 
       title = "Library size per sample") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Save the plot
ggsave(file.path(fig_dir, "library_sizes.png"), p_lib, 
       width = 8, height = 4, dpi = 300)

# 3. DESeq2 object and gene filtering -------------------------------------
# Sample information as a data.frame, with sample IDs as row names
col_data <- meta |>
  column_to_rownames("gsm") |>
  as.data.frame()

# Create the DESeq2 object. Design: tissue effect within each patient (paired)
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = col_data,
  design = ~ patient + tissue
)

# Keep genes with at least 10 counts in at least 18 samples (removes noise)
n_before <- nrow(dds)
keep <- rowSums(counts(dds) >= 10) >= 18
dds <- dds[keep, ]
message("Genes after filtering: ", nrow(dds), " (from ", n_before, ")")

# 4. Transformation and PCA ------------------------------------------------
# Variance stabilising transformation: makes counts usable for PCA and distances
vsd <- vst(dds, blind = TRUE)

# PCA of the 500 most variable genes (we keep the table to plot it ourselves)
pca_df <- plotPCA(vsd, intgroup = c("tissue", "patient"), returnData = TRUE)
# Percentage of variance explained by each principal component
pct_var <- round(100 * attr(pca_df, "percentVar"))

# PCA coloured by tissue; grey lines join samples from the same patient
# How to read it: do the samples separate by tissue, or group by patient?
p_tissue <- ggplot(pca_df, aes(PC1, PC2, color = tissue)) + 
  geom_line(aes(group = patient), color = "grey80", linewidth = 0.4) + 
  geom_point(size = 3) + 
  labs(
    x = paste0("PC1:", pct_var[1], "% variance"),
    y = paste0("PC2:", pct_var[2], "% variance"),
    title = "PCA of VST-transformed counts (top 500 variable genes)",
    subtitle = "Grey lines join samples from the same patient"
  ) + 
  theme_bw()

# Same PCA coloured by patient (shape shows the tissue)
p_patient <- ggplot(pca_df, aes(PC1, PC2, color = patient, shape = tissue)) + 
  geom_point(size = 3) + 
  guides(color = "none") +
  labs(
    x = paste0("PC1:", pct_var[1], "% variance"),
    y = paste0("PC2:", pct_var[2], "% variance"),
    title = "Same PCA coloured by patient"
  ) + 
  theme_bw()

# Save both PCA plots
ggsave(file.path(fig_dir, "pca_by_tissue.png"), p_tissue,
       width = 7, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "pca_by_patient.png"), p_patient, 
       width = 7, height = 5, dpi = 300)

# 5. Sample distance heatmap -----------------------------------------------
# Distance between every pair of samples (small = similar expression)
sample_dist <- dist(t(assay(vsd)))
# Tissue and patient will appear as colour bars
annotation_df <- as.data.frame(colData(vsd)[, c("tissue", "patient")])

# Heatmap with clustering: similar samples are placed together
pheatmap(
  as.matrix(sample_dist),
  clustering_distance_rows = sample_dist,
  clustering_distance_cols = sample_dist,
  annotation_col = annotation_df,
  show_colnames = FALSE,
  show_rownames = FALSE,
  filename = file.path(fig_dir, "sample_distance_heatmap.png"),
  width = 10, height = 8
)

# 6. Save the result -------------------------------------------------------
# Save the filtered DESeq2 object for the next script
saveRDS(dds, here("data", "processed", "dds.rds"))
message("Done. Figures in results/figures, QC table in results/tables")

