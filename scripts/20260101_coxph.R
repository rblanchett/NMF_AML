library(dplyr)
library(data.table)
library(ggplot2)
library(survival)
library(survminer)
library(pheatmap)
library(stringr)

setwd("C:/Users/reida/Documents/Aim_4/TARGET_AML/")

nmf_fit <- readRDS("machine_learning_method/nmf_fit_20251231.rds")
nmf_groups <- readRDS("machine_learning_method/nmf_group_20251231.rds")
nmf_df <- readRDS("machine_learning_method/nmf_df_20251231.rds")

covariate <- read.csv("machine_learning_method/XX_NMF_covariattes_all_20251231.csv")

# -------------------------
# 1️⃣ Load probes
# -------------------------
probes_EPIC <- fread("../AML_test/EPIC.hg38.manifest.tsv.gz", data.table = FALSE)
probes_EPIC_x <- probes_EPIC %>%
  filter(CpG_chrm == "chrX") %>%
  pull(Probe_ID)

# -------------------------
# 2️⃣ Load beta matrix (only needed columns)
# -------------------------
all_cols <- names(fread("../AML_test/GSE190931_Norland_450k/GSE190931_Matrix_processed.csv.gz", nrows = 0))
x_dx_cols <- grep("_Dx$", all_cols, value = TRUE)
x_dx_cols <- c("ID_REF", x_dx_cols)

df <- fread(
  "../AML_test/GSE190931_Norland_450k/GSE190931_Matrix_processed.csv.gz",
  select = x_dx_cols,
  data.table = FALSE
)

#df <- df[df$ID_REF %in% probes_EPIC_x, ]
rownames(df) <- df$ID_REF
df$ID_REF <- NULL
colnames(df) <- sub("_Dx$", "", colnames(df))

beta_mat <- as.matrix(df)

tim_data <- read.csv("../AML_test/target_covs.eligible.csv")  # clinical data
beta_mat <- as.matrix(df)                        # your chrX beta matrix
colnames(beta_mat) <- sub("_Dx$", "", colnames(beta_mat))  # clean sample names

all_females <- tim_data %>% filter(Sex == "Female")


# ------------------------------------
# 1. Subset to females only
# ------------------------------------
beta_mat_xx <- beta_mat[, colnames(beta_mat) %in% all_females$USI]

stopifnot(ncol(beta_mat_xx) > 0)
stopifnot(nrow(beta_mat_xx) > 0)

beta_mat_xx <- as.matrix(beta_mat_xx)
mode(beta_mat_xx) <- "numeric"

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


groups_nmf <- read.csv("machine_learning_method/AML_EPIC_NMF_groups_20251231.csv")
# -----------------------------
# Build modeling dataframe
# -----------------------------
cox_df <- surv_df %>%
  left_join(
    groups_nmf %>% 
      dplyr::select(Sample, NMF_Group),
    by = c("USI" = "Sample")
  ) %>%
  left_join(
    all_females %>%
      dplyr::select(
        USI,
        Age = Age.Category,
        #WBC = WBC..x10.3.MicroLiter..levels,
        Risk = Risk.group.classification..by.protocol.1031
      ),
    by = "USI"
  ) %>%
  filter(complete.cases(.))

cox_base <- coxph(
  Surv(os.time, os_log) ~ NMF_Group + Age + Risk,
  data = cox_df
)

summary(cox_base)

ph_test <- cox.zph(cox_base)
print(ph_test)
ggcoxzph(ph_test)

cox_strat <- coxph(
  Surv(os.time, os_log) ~ NMF_Group + Age + strata(Risk),
  data = cox_df
)

summary(cox_strat)

ggadjustedcurves(
  cox_strat,
  data = cox_df,
  variable = "NMF_Group",
  method = "average"
)

library(ggplot2)
library(survminer)

hr_text <- "Stratified Cox model:\nHR = 0.54 (95% CI 0.36–0.81)\np = 0.003"

p_adj <- ggadjustedcurves(
  cox_strat,
  data = cox_df,
  variable = "NMF_Group",
  method = "average",
  conf.int = TRUE,
  legend.title = "NMF methylation subtype",
  legend.labs = c("Subtype A", "Subtype B")
)
p_adj

p_adj <- p_adj +
  labs(
    title = "Adjusted Overall Survival by NMF-Defined Methylation Subtype",
    subtitle = "Survival curves adjusted for age and stratified by protocol risk",
    x = "Time since diagnosis (days)",
    y = "Adjusted survival probability"
  ) +
  scale_color_manual(
    values = c("salmon", "skyblue")
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11))
      

p_adj





