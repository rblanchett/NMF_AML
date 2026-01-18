# -------------------------
# Load libraries
# -------------------------
library(dplyr)
library(data.table)
library(ggplot2)
library(survival)
library(survminer)
library(pheatmap)
library(stringr)

# -------------------------
# Set working directory
# -------------------------
setwd("C:/Users/reida/Documents/Aim_4/TARGET_AML/")

#========================
# load clinical data
#=========================
tim_data <- read.csv("../AML_test/target_covs.eligible.csv")  # clinical data
beta_mat <- readRDS("machine_learning_method/NMF_G1_vs_G2_beta_filtered.rds")                       # your chrX beta matrix
#colnames(beta_mat) <- sub("_Dx$", "", colnames(beta_mat))  # clean sample names

#========================
# load saved nmf data
#=========================

nmf_fit <- readRDS("machine_learning_method/nmf_fit_20251231.rds")
nmf_df <- readRDS("machine_learning_method/nmf_df_20251231.rds")
groups_nmf <- read.csv("machine_learning_method/AML_EPIC_NMF_groups_20251231.csv")

top100 <- read.csv("machine_learning_method/Top100_module_specific_CpGs_updated_with_male.csv")

#========================
# load GDC clinical data
#========================
# ==============================
# 0️⃣ Load libraries
# ==============================
library(dplyr)
library(data.table)
library(DESeq2)
library(stringr)
library(tibble)

# ==============================
# 1️⃣ Read sample sheets
# ==============================
gdc_samples <- read.table("girls_RNA_seq/gdc_sample_sheet.2026-01-15.tsv", 
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)

head(groups_nmf)

# ==============================
# 2️⃣ Keep only primary tumors
# ==============================
primary_tumor_samples <- gdc_samples %>%
  filter(Tumor.Descriptor == "Primary") %>%
  select(File.Name, Case.ID, Sample.ID)

# ==============================
# 3️⃣ Match to groups_nmf samples
# ==============================
# Keep only rows where Case.ID contains any of the Sample IDs
sample_file_map <- primary_tumor_samples %>%
  rowwise() %>%
  filter(any(sapply(groups_nmf$Sample, function(s) str_detect(Case.ID, s)))) %>%
  ungroup() %>%
  # Map the Sample and NMF_Group
  rowwise() %>%
  mutate(
    Sample = groups_nmf$Sample[sapply(groups_nmf$Sample, function(s) str_detect(Case.ID, s))][1],
    NMF_Group = groups_nmf$NMF_Group[groups_nmf$Sample == Sample]
  ) %>%
  ungroup()

# ==============================
# 4️⃣ Deduplicate (keep first primary tumor per patient)
# ==============================
sample_file_map_unique <- sample_file_map %>%
  group_by(Sample) %>%
  dplyr::slice(1) %>%
  ungroup()

# ==============================
# 5️⃣ List all STAR count files
# ==============================
root_dir <- "girls_RNA_seq/"
all_files <- list.files(path = root_dir, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)

# Map full paths to our selected samples
sample_file_map_unique <- sample_file_map_unique %>%
  rowwise() %>%
  mutate(
    full_path = list.files(root_dir, pattern = File.Name, recursive = TRUE, full.names = TRUE)[1]
  ) %>%
  ungroup()

# ==============================
# 6️⃣ Load count matrices
# ==============================
count_list <- lapply(sample_file_map_unique$full_path, function(f){
  df <- fread(f, data.table = FALSE)
  # Use gene_id as rownames (usually unique)
  rownames(df) <- df$gene_id
  # Pick the counts column (unstranded)
  df_counts <- df[, "unstranded", drop = FALSE]
  return(df_counts)
})

count_matrix <- do.call(cbind, count_list)
colnames(count_matrix) <- sample_file_map_unique$Sample


# ==============================
# 7️⃣ Prepare DESeq2 metadata
# ==============================
col_data <- sample_file_map_unique %>%
  select(Sample, NMF_Group) %>%
  mutate(NMF_Group = factor(NMF_Group)) %>%
  tibble::column_to_rownames("Sample")

# ==============================
# 8️⃣ Run DESeq2
# ==============================
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = col_data,
                              design = ~ NMF_Group)

dds <- DESeq(dds)
res <- results(dds, contrast = c("NMF_Group", "1", "2"))
res <- lfcShrink(dds, coef = 2, type = "apeglm")

res_df <- as.data.frame(res) %>% tibble::rownames_to_column("gene_id")

# Map back to gene_name
gene_map <- fread(sample_file_map_unique$full_path[1], data.table = FALSE) %>%
  dplyr::select(gene_id, gene_name)

res_df <- res_df %>% left_join(gene_map, by = "gene_id")
#write.csv(res_df, "machine_learning_method/DESeq2_results_NMF1_vs_NMF2_primary_tumor.csv", row.names = FALSE)

exp_met <- res_df[res_df$gene_name %in% top100$gene_symbol,]

exp_met_filt <- exp_met %>%
  dplyr::filter(padj < 0.05)# %>%
  #dplyr::filter(abs(log2FoldChange) >1)

top100$gene_name <- top100$gene_symbol

gen_table <- exp_met_filt %>%
  left_join(top100, by = "gene_name") 

gen_table <- gen_table %>%
  dplyr::select(gene_id, log2FoldChange, padj, gene_symbol, Probe_ID, seqnames, context, delta_beta, adj.P.Val, methylation_direction,
                male_g1_adj, male_g1_dir, male_g2_adj, male_g2_dir)


write.csv(gen_table, "machine_learning_method/summary_DGE_DNAme_XX_XY_top100.csv")
# ==============================
# 9️⃣ Save results
# ==============================
#write.csv(res_df, "DESeq2_results_NMF1_vs_NMF2_primary_tumor_only.csv", row.names = FALSE)
