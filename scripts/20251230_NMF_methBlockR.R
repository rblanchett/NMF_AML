# ----------------------------
# 1️⃣ Libraries
# ----------------------------
library(dplyr)
library(data.table)
library(GenomicRanges)
library(ggplot2)
library(pheatmap)
library(methBlockR)
library(survival)
library(survminer)
library(matrixStats)
library(NMF)

# ----------------------------
# 2️⃣ Set working directory
# ----------------------------
setwd("C:/Users/reida/Documents/Aim_4/AML_test")

# ----------------------------
# 3️⃣ Load EPIC probe info
# ----------------------------
probes_EPIC <- fread("EPIC.hg38.manifest.tsv.gz", data.table = FALSE)

# Remove Y chromosome probes
probes_EPIC <- probes_EPIC[probes_EPIC$CpG_chrm != "chrY", ]

# ----------------------------
# 4️⃣ Load beta matrix
# ----------------------------
beta_mat <- fread("GSE190931_Norland_450k/GSE190931_Matrix_processed.csv.gz", data.table = FALSE)
rownames(beta_mat) <- beta_mat$ID_REF
beta_mat$ID_REF <- NULL

# Clean column names
colnames(beta_mat) <- sub("_Dx$", "", colnames(beta_mat))
beta_mat <- as.matrix(beta_mat)
mode(beta_mat) <- "numeric"

# ----------------------------
# 5️⃣ Load clinical data
# ----------------------------
tim_data <- read.csv("target_covs.eligible.csv", stringsAsFactors = FALSE)

# Create a new column for sample group:
# AML_Female, AML_Male, Control_Never, Control_Survived
tim_data <- tim_data %>%
  mutate(Group = case_when(
    Sex == "Female" & Diagnosis == "AML" ~ "AML_Female",
    Sex == "Male"   & Diagnosis == "AML" ~ "AML_Male",
    Diagnosis == "Control" & Status == "NeverLeukemia" ~ "Control_Never",
    Diagnosis == "Control" & Status == "Survived"     ~ "Control_Survived",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Group))

# Subset beta matrix to all selected samples
beta_mat_sub <- beta_mat[, colnames(beta_mat) %in% tim_data$USI]
stopifnot(ncol(beta_mat_sub) == nrow(tim_data))


# ----------------------------
# 6️⃣ Generate methBlocks
# ----------------------------
# Input: CpG positions from EPIC manifest + beta matrix
cpg_gr <- GRanges(
  seqnames = probes_EPIC$CpG_chrm,
  ranges   = IRanges(start = probes_EPIC$CpG_beg, end = probes_EPIC$CpG_end),
  Probe_ID = probes_EPIC$Probe_ID
)

# Align beta matrix to GRanges
beta_mat_sub <- beta_mat_sub[rownames(beta_mat_sub) %in% cpg_gr$Probe_ID, ]
cpg_gr <- cpg_gr[cpg_gr$Probe_ID %in% rownames(beta_mat_sub)]

# Confirm alignment
stopifnot(all(cpg_gr$Probe_ID == rownames(beta_mat_sub)))

# Create methBlocks
blocks <- createMethBlocks(
  betaMatrix = beta_mat_sub,
  cpgInfo    = cpg_gr,
  minCpGs    = 3,      # minimum CpGs per block
  method     = "average" # mean beta per block
)

# Save intermediate
saveRDS(blocks, "methBlocks_allSamples.rds")

# ----------------------------
# 7️⃣ Convert blocks to matrix for NMF
# ----------------------------
block_mat <- assay(blocks)  # block x sample matrix

# Optional: filter out low-variance blocks if performance is slow
block_var <- rowVars(block_mat, na.rm = TRUE)
block_mat <- block_mat[block_var > 0, ]

# ----------------------------
# 8️⃣ Run NMF
# ----------------------------
set.seed(123)
nmf_res <- nmf(block_mat, rank = 2:6, method = "brunet", nrun = 50)
plot(nmf_res)  # consensus plot to pick best K

# Pick K
best_k <- 4  # for example, based on consensus matrix

nmf_fit <- nmf(block_mat, rank = best_k, method = "brunet", nrun = 100)
H <- coef(nmf_fit)
nmf_group <- apply(H, 2, which.max)

nmf_df <- data.frame(
  Sample = colnames(block_mat),
  NMF_Group = factor(nmf_group, levels = 1:best_k)
)

tim_filt <- tim_data %>%
  left_join(nmf_df, by = c("USI" = "Sample"))

# Heatmap of block-level methylation
pheatmap(block_mat[1:500, ], # top 500 blocks for speed
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = tim_filt %>% select(NMF_Group, Group))

# Loadings visualization
W <- basis(nmf_fit)  # block loadings
W_df <- data.frame(
  Block = rownames(W),
  W
)

surv_df <- tim_filt %>%
  dplyr::select(USI, Group, NMF_Group, os.time = OS.time..days., os.evnt = OS.event.ID) %>%
  mutate(
    os.time = as.numeric(os.time),
    os_log  = ifelse(os.evnt == "Dead", 1, 0)
  ) %>%
  filter(!is.na(os.time) & !is.na(os_log))

fit <- survfit(Surv(os.time, os_log) ~ NMF_Group, data = surv_df)
ggsurvplot(fit, pval = TRUE, risk.table = TRUE,
           title = "Survival by NMF Block Clusters")

