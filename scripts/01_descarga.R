library(here)
library(tidyverse)
library(GEOquery)
library(readr)


gse_id <- "GSE50760"
counts_file <- "GSE50760_raw_counts_GRCh38.p13_NCBI.tsv.gz"
counts_url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=", gse_id, "&format=file&file=", counts_file)

raw_dir  <- here("data", "raw")
proc_dir <- here("data", "processed")
counts_path <- file.path(raw_dir, counts_file)

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Raw counts ------------------------------------------------------------
if (!file.exists(counts_path)) {
  message("Dowloading raw counts from NCBI GEO")
  options(timeout = max(300, getOption("timeout")))
  download.file(counts_url, destfile = counts_path, mode = "wb")
} else {
  message("Raw counts already present, skipping download: ", counts_path)
}

#Sanity chech: The first column must be GeneID, the rest GSM accessions
counts_header <- names(read_tsv(counts_path, n_max = 0, show_col_types = FALSE))
if (counts_header[1] != "GeneID") {
  stop("Unexpected counts file (first column is '", counts_header[1],
       "'). Delete ", counts_path, "and download it again.")
}


#2. Sample metadata from GEO    ------------------------------------------------------------
gse <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE)
pheno <- Biobase::pData(gse[[1]])

sample_meta <- pheno |> 
  trsnsmute(
    gsm = geo_accesion,
    title = title,
    label = str_to_lower(paste(title, source_name_ch1)),
    patient = str_extract(title, "AMC_[0-9]+"),
    tissue = case_when(
      str_detect(label, "normal")  ~ "Normal",
      str_detect(label, "metst")   ~ "Metastasis",
      str_detect(label, "primary") ~ "Primary",
      TRUE                         ~ NA_character_
    )
  )

if (anyNA(sample_meta$tissue) || anyNA(sample_meta$patient)) {
  stop ("Could not parse patient/tissue for some samples. Labels found: ",
        paste(unique(sample_meta$label), collapse = " | "))
}

sample_meta <- sample_meta |> 
  mutate(tissue = factor(tissue, levels = c("Normal", "Primary", "Metastasis")))  |>
  select(gsm, title, patient, tissue)  |>
  arrange(match(gsm, gsm_in_counts))

#3. Consistency checks  ------------------------------------------------------------
stopifnot(
  "Samples in metadata and counts differ" = setequal(sample_meta$gsm, gsm_in_counts),
  "Expected 18 patients" = n_distinct(sample_meta$patient),
  "Each patient must have 3 tissues" = all(table(sample_meta$patient, sample_meta$tissue) == 1)
)

write_csv(sample_meta, file.path(proc_dir, "sample_metadata.csv"))

message("Samples per tissue: ")
print(count.fields(sample_meta, tissue))
message("Done. Metadata written to data/processed/sample_metadata.csv")


                       