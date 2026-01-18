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
# 3️⃣ Load clinical data
# -------------------------
# -------------------------
# Required libraries
# -------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(survminer)
library(pheatmap)

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

# Impute remaining NAs using CpG-wise mean
#beta_clean[is.na(beta_clean)] <- rowMeans(beta_clean, na.rm = TRUE)

# ------------------------------------
# 3. Variance filtering (robust)
#    EPIC best practice: 10k–20k CpGs
# ------------------------------------
# Citation:
# 

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

# ------------------------------------
# 5. Minimal NMF test (debug step)
# ------------------------------------
set.seed(191)

#nmf_test <- nmf(
#  beta_filt,
#  rank   = 2,
#  method = "brunet",
#  nrun   = 1
#)

# ------------------------------------
# 6. Rank selection (full run)
# ------------------------------------
#nmf_res <- nmf(
#  beta_filt,
#  rank   = 2:5,
#  method = "brunet",
#  nrun   = 50
#)


#plot(nmf_res)


############
### Final K based on plot
best_k <- 2

set.seed(191)
#nmf_fit <- nmf(beta_filt,
#               rank = best_k,
#               method = "brunet",
#               nrun = 100)

#saveRDS(nmf_fit, "machine_learning_method/nmf_fit_20251231.rds")

H <- coef(nmf_fit)

#nmf_group <- apply(H, 2, which.max)
#saveRDS(nmf_group, "machine_learning_method/nmf_group_20251231.rds")

consensusmap(nmf_fit)
heatmap(H, scale = "row")

#nmf_df <- data.frame(
#  Sample = colnames(beta_filt),
#  NMF_Group = factor(nmf_group)
#)

#saveRDS(nmf_df, "machine_learning_method/nmf_df_20251231.rds")

probes_EPIC_x <- probes_EPIC %>%
  filter(CpG_chrm == "chrX") %>%
  pull(Probe_ID)

x_cpgs <- probes_EPIC_x

mean_x <- colMeans(beta_mat_xx[rownames(beta_mat_xx) %in% x_cpgs, , drop = FALSE], na.rm = TRUE)

boxplot(mean_x ~ nmf_df$NMF_Group,
        ylab = "Mean chrX beta",
        xlab = "NMF group")

wilcox.test(mean_x ~ nmf_df$NMF_Group)

W <- basis(nmf_fit)

W_df <- data.frame(
  CpG = rownames(W),
  Component1 = W[,1],
  Component2 = W[,2],
  chr = probes_EPIC$CpG_chrm[match(rownames(W), probes_EPIC$Probe_ID)]
)

aggregate(Component1 ~ chr, data = W_df, mean)
aggregate(Component2 ~ chr, data = W_df, mean)

library(survival)
library(survminer)

surv_df <- tim_data %>%
  filter(USI %in% nmf_df$Sample) %>%
  dplyr::select(USI, os.time = OS.time..days., os.evnt = OS.event.ID)

# Convert columns to correct types
surv_df <- surv_df %>%
  mutate(
    os.time = as.numeric(as.character(os.time)),       # numeric survival time
    os_log  = ifelse(as.character(os.evnt) == "Dead", 1, 0)  # numeric event: 1=Dead, 0=Alive
  )

# Remove samples with missing values
surv_df <- surv_df %>% filter(!is.na(os.time) & !is.na(os_log))

message("Number of samples with beta & survival info for AML: ", nrow(surv_df))

# 2️⃣ Align metric to filtered clinical samples
nmf_df_2 <- nmf_df[surv_df$USI,]

# 3️⃣ Build survival object
surv_obj <- Surv(time = surv_df$os.time, event = surv_df$os_log)


surv_df$NMF_Group <- nmf_df_2$NMF_Group
# 4️⃣ Fit Cox model
#fit_cox <- coxph(surv_obj ~ nmf_df)
#message("Cox model for ", group_label)

#print(summary(fit_cox))


#surv_df <- merge(tim_data, nmf_df, by = "USI")

fit <- survfit(Surv(time = surv_df$os.time, event = surv_df$os_log) ~ NMF_Group,
               data = surv_df)

ggsurvplot(fit, pval = TRUE, risk.table = TRUE, title = "AML (XX): Kaplan-Meier by Genome-wide Beta NMF",
           xlab = "Time (days)",
           ylab = "Survival Probability",
           legend.title = "NMF Group")

coxph(Surv(time = surv_df$os.time, event = surv_df$os_log) ~ NMF_Group,
      data = surv_df)
#==============================================

#Call: 
#coxph(formula = Surv(time = surv_df$os.time, event = surv_df$os_log) ~ NMF_Group, data = surv_df) 
#coef exp(coef) se(coef) z p 
#NMF_Group2 -0.7340 0.4800 0.1812 -4.051 5.1e-05 
#Likelihood ratio test=17.29 on 1 df, p=3.213e-05 
#n= 400, number of events= 134

#write.csv(nmf_df_2, "../TARGET_AML/machine_learning_method/AML_EPIC_NMF_groups_20251231.csv", row.names = FALSE)

####################
### Loadings
####################

rownames(W) <- rownames(beta_filt)
colnames(W) <- paste0("Module_", seq_len(ncol(W)))

top_cpgs <- lapply(colnames(W), function(mod) {
  ord <- order(W[, mod], decreasing = TRUE)
  data.frame(
    CpG = rownames(W)[ord],
    Loading = W[ord, mod],
    Module = mod
  )
})
names(top_cpgs) <- colnames(W)

head(top_cpgs$Module_1, 20)
head(top_cpgs$Module_2, 20)

loading_diff <- W[, 1] - W[, 2]

diff_df <- data.frame(
  CpG = rownames(W),
  Loading_1 = W[, 1],
  Loading_2 = W[, 2],
  Diff = loading_diff
)

diff_df <- diff_df[order(abs(diff_df$Diff), decreasing = TRUE), ]
head(diff_df, 50)

library(ggplot2)
library(tidyr)
library(dplyr)

W_long <- as.data.frame(W) %>%
  mutate(CpG = rownames(W)) %>%
  pivot_longer(
    cols = starts_with("Module"),
    names_to = "Module",
    values_to = "Loading"
  )

ggplot(W_long, aes(x = Loading)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "black") +
  facet_wrap(~ Module, scales = "free_y") +
  theme_classic() +
  labs(
    title = "Distribution of NMF CpG loadings",
    x = "Loading weight",
    y = "CpG count"
  )

diff_df <- data.frame(
  CpG = rownames(W),
  Module1 = W[, 1],
  Module2 = W[, 2],
  Diff = W[, 1] - W[, 2]
)

ggplot(diff_df, aes(x = Diff)) +
  geom_histogram(bins = 200, fill = "darkred", color = "black") +
  theme_classic() +
  labs(
    title = "CpG loading difference (Module 1 – Module 2)",
    x = "Loading difference",
    y = "CpG count"
  )

library(GenomicRanges)

# Ensure CpGs line up
stopifnot(all(rownames(W) %in% probes_EPIC$Probe_ID))

# Subset probe dictionary to CpGs used in NMF
gr_loadings <- probes_EPIC[rownames(W) %in% probes_EPIC$Probe_ID,]

library(GenomicRanges)

# Remove rows with missing coordinates
probes_clean <- probes_EPIC[!is.na(probes_EPIC$CpG_beg) & !is.na(probes_EPIC$CpG_end) & !is.na(probes_EPIC$CpG_chrm), ]

# Now convert to GRanges
gr_probes <- GRanges(
  seqnames = probes_clean$CpG_chrm,
  ranges   = IRanges(start = probes_clean$CpG_beg,
                     end   = probes_clean$CpG_end),
  strand   = "*",
  Probe_ID = probes_clean$Probe_ID
)

# Assign names for easy subsetting
names(gr_probes) <- gr_probes$Probe_ID

# Quick check
gr_probes[1:5]

# Ensure that all CpGs in W exist in gr_probes
common_cpgs <- intersect(rownames(W), gr_probes$Probe_ID)

# Subset both to matching CpGs
W_sub <- W[common_cpgs, , drop = FALSE]
gr_sub <- gr_probes[match(common_cpgs, gr_probes$Probe_ID)]

# Attach the loadings as metadata
mcols(gr_sub) <- cbind(mcols(gr_sub), as.data.frame(W_sub))

# Quick check
gr_sub[1:5]

# Convert GRanges metadata to a normal data.frame
df_loadings <- as.data.frame(mcols(gr_sub))

# Add seqnames as a regular column
df_loadings$chr <- as.character(seqnames(gr_sub))

# Now aggregate safely
agg <- aggregate(cbind(Module_1, Module_2) ~ chr, data = df_loadings, mean)
agg

library(GenomicRanges)
library(ggplot2)
library(dplyr)

# ----------------------------
# 1️⃣ Convert GRanges to data.frame
# ----------------------------
df_loadings <- as.data.frame(mcols(gr_sub))
df_loadings$chr <- as.character(seqnames(gr_sub))
df_loadings$start <- start(gr_sub)  # genomic coordinate
df_loadings$Probe_ID <- gr_sub$Probe_ID

# ----------------------------
# 2️⃣ Identify top CpGs for each component
# ----------------------------
top_n <- 1000  # top 1000 CpGs per component

top_comp1 <- df_loadings %>% arrange(desc(Module_1)) %>% slice_head(n = top_n)
top_comp2 <- df_loadings %>% arrange(desc(Module_2)) %>% slice_head(n = top_n)

# ----------------------------
# 3️⃣ Genome-wide scatter plot
# ----------------------------
ggplot(df_loadings, aes(x = start, y = Module_1, color = chr)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_point(data = top_comp1, aes(x = start, y = Module_1), color = "red", size = 1.5) +
  facet_wrap(~chr, scales = "free_x") +
  theme_bw() +
  labs(
    title = "Genome-wide distribution of NMF loadings (Component 1)",
    x = "Genomic coordinate",
    y = "Module 1 loading"
  )

# Repeat for Component 2
ggplot(df_loadings, aes(x = start, y = Module_2, color = chr)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_point(data = top_comp2, aes(x = start, y = Module_2), color = "blue", size = 1.5) +
  facet_wrap(~chr, scales = "free_x") +
  theme_bw() +
  labs(
    title = "Genome-wide distribution of NMF loadings (Component 2)",
    x = "Genomic coordinate",
    y = "Module 2 loading"
  )

# ----------------------------
# 4️⃣ Zoom in on chrX
# ----------------------------
df_X <- df_loadings %>% filter(chr == "chrX")
top_comp1_X <- top_comp1 %>% filter(chr == "chrX")
top_comp2_X <- top_comp2 %>% filter(chr == "chrX")

ggplot(df_X, aes(x = start, y = Module_1)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_point(data = top_comp1_X, aes(x = start, y = Module_1), color = "red", size = 2) +
  theme_bw() +
  labs(
    title = "Component 1 loadings on chrX",
    x = "Genomic coordinate",
    y = "Module 1 loading"
  )

ggplot(df_X, aes(x = start, y = Module_2)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_point(data = top_comp2_X, aes(x = start, y = Module_2), color = "blue", size = 2) +
  theme_bw() +
  labs(
    title = "Component 2 loadings on chrX",
    x = "Genomic coordinate",
    y = "Module 2 loading"
  )

# ----------------------------
# 5️⃣ Optional: chromosome-level summary
# ----------------------------
df_summary <- df_loadings %>%
  group_by(chr) %>%
  summarise(
    mean_module1 = mean(Module_1),
    mean_module2 = mean(Module_2),
    max_module1  = max(Module_1),
    max_module2  = max(Module_2),
    .groups = "drop"
  )

df_summary

###################
### X chromosome from df_loading
###################


##################
#### Survival with male samples
#################

library(survival)
library(survminer)


all_males <- tim_data %>% filter(Sex == "Male")
all_males <- all_males[all_males$USI %in% colnames(beta_mat),]
all_males$Comparison <- "Male"
nmf_df_2$Comparison <- ifelse(nmf_df_2$NMF_Group == "1", "NMF_G1", "NMF_G2")
all_females <- tim_data %>% filter(Sex == "Female")
all_females <- all_females[all_females$USI %in% colnames(beta_mat),]
nmf_df_3 <- nmf_df_2[all_females$USI,]
all_females$Comparison <- nmf_df_3$Comparison

tim_comb2 <- rbind.data.frame(all_males, all_females)

#write.csv(tim_comb, "machine_learning_method/AML_ML_clinical_with_groups.csv")

surv_df <- tim_comb2 %>%
  #filter(USI %in% nmf_df$Sample) %>%
  dplyr::select(USI, os.time = OS.time..days., os.evnt = OS.event.ID)

# Convert columns to correct types
surv_df <- surv_df %>%
  mutate(
    os.time = as.numeric(as.character(os.time)),       # numeric survival time
    os_log  = ifelse(as.character(os.evnt) == "Dead", 1, 0)  # numeric event: 1=Dead, 0=Alive
  )

# Remove samples with missing values
surv_df <- surv_df %>% filter(!is.na(os.time) & !is.na(os_log))
surv_df <- surv_df[surv_df$USI %in% colnames(beta_mat),]


message("Number of samples with beta & survival info for AML: ", nrow(surv_df))


# 3️⃣ Build survival object
surv_obj <- Surv(time = surv_df$os.time, event = surv_df$os_log)
surv_df$Comparison <- tim_comb2$Comparison


# 4️⃣ Fit Cox model
#fit_cox <- coxph(surv_obj ~ nmf_df)
#message("Cox model for ", group_label)

#print(summary(fit_cox))


#surv_df <- merge(tim_data, nmf_df, by = "USI")

fit <- survfit(Surv(time = surv_df$os.time, event = surv_df$os_log) ~ Comparison,
               data = surv_df)

ggsurvplot(fit, palette = c("green3", "salmon", "skyblue"), pval = TRUE, risk.table = TRUE, title = "AML Males and Females: Kaplan-Meier by Genome-wide Beta",
           xlab = "Time (days)",
           ylab = "Survival Probability",
           legend.title = "Group", legend.labs = c("Male", "Female NMF Group 1", "Female NMF Group 2"))

coxph(Surv(time = surv_df$os.time, event = surv_df$os_log) ~ Comparison,
      data = surv_df)

#Call:
#coxph(formula = Surv(time = surv_df$os.time, event = surv_df$os_log) ~ 
#        Comparison, data = surv_df)

#====================================

#coef exp(coef) se(coef)      z       p
#ComparisonNMF_G1  0.2496    1.2835   0.1327  1.882 0.05989
#ComparisonNMF_G2 -0.4908    0.6121   0.1655 -2.966 0.00302

#Likelihood ratio test=18.04  on 2 df, p=0.000121
#n= 854, number of events= 298 
