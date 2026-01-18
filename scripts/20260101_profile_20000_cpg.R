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
setwd("C:/Users/reida/Documents/Aim_4/AML_test")

# -------------------------
# 1️⃣ Load chrX probes
# -------------------------
probes_EPIC <- fread("EPIC.hg38.manifest.tsv.gz", data.table = FALSE)
probes_EPIC_x <- probes_EPIC %>%
  filter(CpG_chrm == "chrX") %>%
  pull(Probe_ID)

# -------------------------
# 2️⃣ Load beta matrix (only needed columns)
# -------------------------
all_cols <- names(fread("GSE190931_Norland_450k/GSE190931_Matrix_processed.csv.gz", nrows = 0))
x_dx_cols <- grep("_Dx$", all_cols, value = TRUE)
x_dx_cols <- c("ID_REF", x_dx_cols)

df <- fread(
  "GSE190931_Norland_450k/GSE190931_Matrix_processed.csv.gz",
  select = x_dx_cols,
  data.table = FALSE
)

#df <- df[df$ID_REF %in% probes_EPIC_x, ]
rownames(df) <- df$ID_REF
df$ID_REF <- NULL
colnames(df) <- sub("_Dx$", "", colnames(df))

beta_mat <- as.matrix(df)

# -------------------------
# 0️⃣ Load clinical and beta data
# -------------------------
setwd("C:/Users/reida/Documents/Aim_4/AML_test")

tim_data <- read.csv("../AML_test/target_covs.eligible.csv")  # clinical data
beta_mat <- as.matrix(df)                        # your chrX beta matrix
colnames(beta_mat) <- sub("_Dx$", "", colnames(beta_mat))  # clean sample names

all_females <- tim_data %>% filter(Sex == "Female")

# beta_mat: CpG x Sample
# annotation: EPIC annotation with chr column
library(matrixStats)
library(NMF)

# ------------------------------------
# 1. Subset to females only
# ------------------------------------
beta_mat_xx <- beta_mat[, colnames(beta_mat) %in% all_females$USI]

stopifnot(ncol(beta_mat_xx) > 0)
stopifnot(nrow(beta_mat_xx) > 0)

beta_mat_xx <- as.matrix(beta_mat_xx)
mode(beta_mat_xx) <- "numeric"

# ------------------------------------
# 2. Handle missing values SAFELY
#    (do NOT use complete.cases here)
# ------------------------------------

# Fraction of missing values per CpG
na_frac <- rowMeans(is.na(beta_mat_xx))

# Keep CpGs with <=10% missing values
beta_clean <- beta_mat_xx[na_frac <= 0.1, ]

stopifnot(nrow(beta_clean) > 0)

# Compute variance
cpg_var <- rowVars(beta_clean, na.rm = TRUE)

# Remove CpGs with non-finite or zero variance
keep <- is.finite(cpg_var) & cpg_var > 0
beta_var <- beta_clean[keep, ]

# Impute remaining NA values using CpG-wise mean
na_idx <- is.na(beta_var)

if (any(na_idx)) {
  row_means <- rowMeans(beta_var, na.rm = TRUE)
  beta_var[na_idx] <- row_means[row(beta_var)[na_idx]]
}

stopifnot(all(is.finite(beta_var)))

cpg_var  <- cpg_var[keep]

stopifnot(
  nrow(beta_var) > 0,
  all(is.finite(cpg_var)),
  all(cpg_var > 0)
)

# Select top variable CpGs
top_n <- 20000
stopifnot(length(cpg_var) >= top_n)

top_idx <- order(cpg_var, decreasing = TRUE)[seq_len(top_n)]
beta_filt <- beta_var[top_idx, ]

stopifnot(
  all(is.finite(beta_filt)),
  all(rowVars(beta_filt) > 0)
)



sorted_vars <- sort(cpg_var, decreasing = TRUE)

plot(
  sorted_vars,
  type = "l",
  log = "y",
  xlab = "CpGs ranked by variance",
  ylab = "Variance (log scale)",
  main = "Ranked CpG Variance on chrX"
)

abline(v = 20000, col = "red", lty = 2)

cum_var <- cumsum(sorted_vars) / sum(sorted_vars)

plot(
  cum_var,
  type = "l",
  xlab = "Number of CpGs included",
  ylab = "Cumulative proportion of variance",
  main = "Cumulative Methylation Variance Explained"
)

abline(v = 20000, col = "red", lty = 2)

feature_sizes <- c(10000, 20000, 40000)

for (n in feature_sizes) {
  idx <- order(cpg_var, decreasing = TRUE)[seq_len(n)]
  nmf_res <- nmf(beta_var[idx, ], rank = 4, nrun = 10)
  print(paste("Features:", n, 
              "Cophenetic:", round(cophcor(nmf_res), 3)))
}

