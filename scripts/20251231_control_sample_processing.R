setwd("C://Users/reida/Documents/Aim_4/TARGET_AML/machine_learning_method/")

#never_clin <- read.table(gzfile("GSE124413_series_matrix.txt.gz"), sep = "\t")

library(GEOquery)

#------------------------------------
# Pediatric control sample processing, bone marrow never had AML
#------------------------------------

gse_never <- getGEO("GSE124413", GSEMatrix = TRUE)

pd_never <- pData(gse_never[[1]])

aml_idx_never <- grepl("RO", pd_never$title)
table(aml_idx_never)
pd_aml_never <- pd_never[aml_idx_never, ]
nrow(pd_aml_never)

gsm_aml_never <- pd_aml_never$geo_accession
head(gsm_aml_never)
idat_urls <- c(
  pd_aml_never$supplementary_file,
  pd_aml_never$supplementary_file.1
)

length(idat_urls)  # should be 2 × number of AML samples
head(idat_urls)

dir.create("AML_IDATs_marrow_never", showWarnings = FALSE)

idat_urls_https <- sub("^ftp://", "https://", idat_urls)


for (u in idat_urls_https) {
  dest <- file.path("AML_IDATs_marrow_never", basename(u))
  download.file(u, destfile = dest, mode = "wb", quiet = FALSE)
}

# ==========================================
# 1️⃣ Load libraries
# ==========================================
if (!requireNamespace("Sesame", quietly = TRUE)) BiocManager::install("Sesame")
library(sesame)

if (!requireNamespace("sesameData", quietly = TRUE)) BiocManager::install("sesameData")
library(sesameData)


# ==========================================
# 2️⃣ Define IDAT folder
# ==========================================
idat_dir <- "AML_IDATs_marrow_never"  # replace with your directory
idat_files <- list.files(idat_dir, pattern = ".idat$", full.names = TRUE)
length(idat_files)  # sanity check

#idat_files <- list.files("AML_IDATs_marrow_never", pattern = "idat.gz$", full.names = TRUE)
#sapply(idat_files, function(f) {
#  R.utils::gunzip(f, remove=FALSE)  # keeps original .gz
#})

# ==========================================
# 3️⃣ Read IDAT pairs
# ==========================================
files <- list.files(idat_dir, pattern = "\\.idat$", full.names = TRUE)

df <- data.frame(
  file = files,
  prefix = sub("_(Red|Grn)\\.idat$", "", basename(files)),
  type = ifelse(grepl("_Red\\.idat$", files), "Red", "Grn"),
  stringsAsFactors = FALSE
)

# Count Red/Grn per prefix
pair_table <- table(df$prefix, df$type)

# Find broken samples
bad_prefixes <- rownames(pair_table)[
  pair_table[, "Red"] != 1 | pair_table[, "Grn"] != 1
]

length(bad_prefixes)
bad_prefixes[1:10]

prefixes <- unique(sub("_(Red|Grn)\\.idat$", "", df$file))
length(prefixes)

sesame_list <- lapply(prefixes, function(p) {
  readIDATpair(p)
})

names(sesame_list) <- basename(prefixes)


# ==========================================
# 4️⃣ Apply noob background correction
# ==========================================
# SeSAMe implements the 'noob' method internally
# The default offset = 15 is added automatically
sesame_list <- lapply(sesame_list, noob)

# ==========================================
# 5️⃣ Apply non-linear dye bias correction
# ==========================================
sesame_list <- lapply(sesame_list, dyeBiasCorrTypeINorm)

# ==========================================
# 6️⃣ SNP and non-detection masking
# ==========================================
# Mask probes overlapping SNPs or with poor detection p-values
sesame_prequality <- sesame_list
sesame_list <- lapply(sesame_list, qualityMask)
# ==========================================
# 7️⃣ Compute beta values
# ==========================================
# Beta = M / (M + U + offset)
# SeSAMe computes this automatically; we extract
beta_list <- lapply(sesame_list, function(sa) getBetas(sa))

# Combine into one matrix: rows = probes, columns = samples
beta_mat_never <- do.call(cbind, beta_list)


beta_list_nomask <-lapply(sesame_prequality, function(sa) getBetas(sa))
beta_mat_nomask <- do.call(cbind, beta_list_nomask)
sex_inference <- inferSex(beta_mat_nomask)


sum(rownames(beta_mat_never) %in% probes_EPIC_x)
sum(rownames(beta_mat_never) %in% probes_EPIC$Probe_ID[probes_EPIC$CpG_chrm == "chrY"])

library(sesame)
anno <- sesameDataGet("EPIC.hg38")

x_ids <- probes_EPIC$Probe_ID[probes_EPIC$CpG_chrm == "chrX"]
probes_EPIC_y <- probes_EPIC %>%
  filter(CpG_chrm == "chrY") %>%
  pull(Probe_ID)

x_non_na <- colSums(!is.na(beta_mat_never[probes_EPIC_x, , drop = FALSE]))
y_non_na <- colSums(!is.na(beta_mat_never[probes_EPIC_y, , drop = FALSE]))

summary(x_non_na)
summary(y_non_na)

which(y_non_na < 2)

#“Sex was inferred using summary statistics of chromosome Y DNA methylation (standard deviation and mean beta), 
#as kernel density–based sex classifiers are unstable in mixed-sex cohorts due to the near-constant chrY signal 
#observed in female samples.”


y_sd <- apply(
  beta_mat_never[probes_EPIC_y, , drop = FALSE],
  2,
  sd,
  na.rm = TRUE
)

sex <- ifelse(y_sd > 0.2, "Male", "Female")
sex

# ==========================================
# ✅ 9) Result
# ==========================================
# beta_mat_never is your final methylation beta matrix
# Rows = high-quality probes
# Columns = samples
# Values between 0 (unmethylated) and 1 (fully methylated)

dim(beta_mat_never)
head(beta_mat_never[,1:5])

write.table(beta_mat_never, "../machine_learning_method/AML_marrow_never_beta_matrix.tsv")
write.csv(sex, "../machine_learning_method/AML_marrow_never_sex_ids.csv")
#------------------------------------
# Adult control sample processing, blood of former pAML
#------------------------------------
gse_adult <- getGEO("GSE169156", GSEMatrix = TRUE)

#df_adult <- read.table(gzfile("GSE169156_SJLIFE_IlluminaEPIC_2052samples_GEO_03162021_processed.txt.gz"))


gsm_ids <- rownames(pData(gse_adult[[1]]))

# After filtering gsm_ids to AML survivors:
pd <- pData(gse_adult[[1]])

aml_idx <- grepl("SJAML", pd$title)
table(aml_idx)
pd_aml <- pd[aml_idx, ]
nrow(pd_aml)
gsm_aml <- pd_aml$geo_accession
head(gsm_aml)
idat_urls <- c(
  pd_aml$supplementary_file,
  pd_aml$supplementary_file.1
)

length(idat_urls)  # should be 2 × number of AML samples
head(idat_urls)

idat_urls_truncated <- idat_urls[80:100]

dir.create("AML_IDATs", showWarnings = FALSE)

idat_urls_https <- sub("^ftp://", "https://", idat_urls_truncated)


for (u in idat_urls_https) {
  dest <- file.path("AML_IDATs", basename(u))
  download.file(u, destfile = dest, mode = "wb", quiet = FALSE)
}

# ==========================================
# 1️⃣ Load libraries
# ==========================================
if (!requireNamespace("Sesame", quietly = TRUE)) BiocManager::install("Sesame")
library(sesame)



# ==========================================
# 2️⃣ Define IDAT folder
# ==========================================
idat_files_gz <- list.files("../machine_learning_method/AML_IDATs/", pattern = "idat.gz$", full.names = TRUE)
sapply(idat_files_gz, function(f) {
 R.utils::gunzip(f, remove=FALSE)  # keeps original .gz
})

idat_dir <- "../machine_learning_method/AML_IDATs/"  # replace with your directory
idat_files <- list.files(idat_dir, pattern = "idat$", full.names = TRUE)
length(idat_files)  # sanity check


df <- data.frame(
  file = idat_files,
  prefix = sub("_(Red|Grn)\\.idat$", "", basename(idat_files)),
  type = ifelse(grepl("_Red\\.idat$", idat_files), "Red", "Grn"),
  stringsAsFactors = FALSE
)

# Count Red/Grn per prefix
pair_table <- table(df$prefix, df$type)

# Find broken samples
bad_prefixes <- rownames(pair_table)[
  pair_table[, "Red"] != 1 | pair_table[, "Grn"] != 1
]

length(bad_prefixes)
bad_prefixes[1:10]

prefixes <- unique(sub("_(Red|Grn)\\.idat$", "", df$file))
length(prefixes)

# ==========================================
# 3️⃣ Read IDATs using SeSAMe
# ==========================================
# This creates a SesameArray object for each sample

sesame_list <- lapply(prefixes, function(p) {
  readIDATpair(p)
})

names(sesame_list) <- basename(prefixes)


# ==========================================
# 4️⃣ Apply noob background correction
# ==========================================
# SeSAMe implements the 'noob' method internally
# The default offset = 15 is added automatically
sesame_list <- lapply(sesame_list, noob)

# ==========================================
# 5️⃣ Apply non-linear dye bias correction
# ==========================================
sesame_list <- lapply(sesame_list, dyeBiasCorrTypeINorm)


# ==========================================
# 6️⃣ SNP and non-detection masking
# ==========================================
# Mask probes overlapping SNPs or with poor detection p-values
sesame_list <- lapply(sesame_list, qualityMask)

# ==========================================
# 7️⃣ Compute beta values
# ==========================================
# Beta = M / (M + U + offset)
# SeSAMe computes this automatically; we extract
beta_list <- lapply(sesame_list, function(sa) getBetas(sa))

# Combine into one matrix: rows = probes, columns = samples
beta_mat <- do.call(cbind, beta_list)


write.table(beta_mat, "../machine_learning_method/AML_adults_beta_matrix.tsv")



