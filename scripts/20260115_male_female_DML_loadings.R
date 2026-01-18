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

# -------------------------
# 1️⃣ Load chrX probes
# -------------------------
probes_EPIC <- fread("../AML_test/EPIC.hg38.manifest.tsv.gz", data.table = FALSE)

# -------------------------
# 2️⃣ Load beta matrix (only needed columns)
# -------------------------
all_cols <- names(fread("../AML_test/GSE190931_Norland_450k/GSE190931_Matrix_processed.csv.gz", nrows = 0))
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

setwd("C:/Users/reida/Documents/Aim_4/AML_test")
#========================
# load clinical data
#=========================
tim_data <- read.csv("../AML_test/target_covs.eligible.csv")  # clinical data
beta_mat <- as.matrix(df)                        # your chrX beta matrix
colnames(beta_mat) <- sub("_Dx$", "", colnames(beta_mat))  # clean sample names

#========================
# load saved nmf data
#=========================

nmf_fit <- readRDS("machine_learning_method/nmf_fit_20251231.rds")
nmf_df <- readRDS("machine_learning_method/nmf_df_20251231.rds")
groups_nmf <- read.csv("machine_learning_method/AML_EPIC_NMF_groups_20251231.csv")

#========================
# DML analysis against male 
#========================
cpgs_filt <- rownames(beta_filt)

beta_mat_filt <- beta_mat[rownames(beta_mat) %in% cpgs_filt, ]

library(sesame)
library(SummarizedExperiment)

# beta_mat: CpGs x samples
se <- SummarizedExperiment(
  assays = list(beta = beta_mat_filt)
)

sample_ids <- colnames(beta_mat_filt)

sample_annot <- data.frame(
  Sample = sample_ids,
  stringsAsFactors = FALSE
)

sample_annot <- sample_annot %>%
  dplyr::left_join(groups_nmf, by = "Sample") %>%
  dplyr::mutate(
    Label = ifelse(is.na(Comparison), "Male", Comparison)
  )

colData(se)$group <- factor(sample_annot$Label)

smry <- DML(se, ~ group)

dml_male_vs_g1 <- lapply(names(smry), function(cpg) {
  fit <- smry[[cpg]]
  
  if("groupNMF_G1" %in% rownames(fit$coefficients)) {
    slope <- fit$coefficients["groupNMF_G1", "Estimate"]
    pval  <- fit$coefficients["groupNMF_G1", "Pr(>|t|)"]
  } else {
    slope <- NA
    pval <- NA
  }
  
  data.frame(
    CpG = cpg,
    slope = slope,
    pval  = pval,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows() %>%
  mutate(
    direction = ifelse(slope > 0, "Hyper in NMF_G1", "Hypo in NMF_G1"),
    adj.P.Val = p.adjust(pval, method = "BH"),   # <-- FDR correction
    sig = adj.P.Val < 0.05 
  )



dml_male_vs_g2 <- lapply(names(smry), function(cpg) {
  fit <- smry[[cpg]]
  
  if("groupNMF_G2" %in% rownames(fit$coefficients)) {
    slope <- fit$coefficients["groupNMF_G2", "Estimate"]
    pval  <- fit$coefficients["groupNMF_G2", "Pr(>|t|)"]
  } else {
    slope <- NA
    pval <- NA
  }
  
  data.frame(
    CpG = cpg,
    slope = slope,
    pval  = pval,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows() %>%
  mutate(
    direction = ifelse(slope > 0, "Hyper in NMF_G2", "Hypo in NMF_G2"),
    adj.P.Val = p.adjust(pval, method = "BH"),   # <-- FDR correction
    sig = adj.P.Val < 0.05                        # significance based on FDR
  )

#========================
# DML analysis G1 v G2
#========================


library(sesame)
library(SummarizedExperiment)

# beta_mat: CpGs x samples
se <- SummarizedExperiment(
  assays = list(beta = beta_filt)
)

sample_ids <- colnames(beta_filt)

sample_annot <- data.frame(
  Sample = sample_ids,
  stringsAsFactors = FALSE
)

sample_annot <- sample_annot %>%
  dplyr::left_join(groups_nmf, by = "Sample") %>%
  dplyr::mutate(
    Label = ifelse(is.na(Comparison), "Male", Comparison)
  )

colData(se)$group <- factor(sample_annot$Label)

smry <- DML(se, ~ group)
dmContrasts(smry) 

dml_g1_vs_g2 <- lapply(names(smry), function(cpg) {
  fit <- smry[[cpg]]
  
  if("groupNMF_G2" %in% rownames(fit$coefficients)) {
    slope <- fit$coefficients["groupNMF_G2", "Estimate"]
    pval  <- fit$coefficients["groupNMF_G2", "Pr(>|t|)"]
  } else {
    slope <- NA
    pval <- NA
  }
  
  data.frame(
    CpG = cpg,
    slope = slope,
    pval  = pval,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows() %>%
  mutate(
    direction = ifelse(slope > 0, "Hyper in NMF_G2", "Hypo in NMF_G2"),
    adj.P.Val = p.adjust(pval, method = "BH"),   # <-- FDR correction
    sig = adj.P.Val < 0.05 
  )


top100 <- read.csv("machine_learning_method/Top100_module_specific_CpGs_annotated.csv")

topcpg <- top100$Probe_ID

dml_g1_vs_g2_top <- dml_g1_vs_g2[dml_g1_vs_g2$CpG %in% topcpg,]
dml_male_vs_g1_top <- dml_male_vs_g1[dml_male_vs_g1$CpG %in% topcpg,]
dml_male_vs_g2_top <- dml_male_vs_g2[dml_male_vs_g2$CpG %in% topcpg,]

colnames(dml_g1_vs_g2_top)[5] <- "sesame.adj.p"

top100 <- top100 %>%
  dplyr::left_join(
    dml_g1_vs_g2_top %>% dplyr::select(CpG, sesame.adj.p),
    by = c("Probe_ID" = "CpG")
  ) %>%
  dplyr::rename(
    dml_g1_g2_adj = sesame.adj.p
  )

# Male vs G1
top100 <- top100 %>%
  dplyr::left_join(
    dml_male_vs_g1_top %>%
      dplyr::select(CpG, slope, adj.P.Val, direction) %>%
      dplyr::rename(
        male_g1_slope = slope,
        male_g1_adj = adj.P.Val,
        male_g1_dir = direction
      ),
    by = c("Probe_ID" = "CpG")
  )

# Male vs G2
top100 <- top100 %>%
  dplyr::left_join(
    dml_male_vs_g2_top %>%
      dplyr::select(CpG, slope, adj.P.Val, direction) %>%
      dplyr::rename(
        male_g2_slope = slope,
        male_g2_adj = adj.P.Val,
        male_g2_dir = direction
      ),
    by = c("Probe_ID" = "CpG")
  )

library(ComplexUpset)
library(ggplot2)
library(dplyr)

# Build the UpSet dataframe
upset_df <- top100 %>%
  transmute(
    Probe_ID,
    sig_g1_g2   = dml_g1_g2_adj < 0.05,
    sig_male_g1 = male_g1_adj < 0.05,
    sig_male_g2 = male_g2_adj < 0.05,
    dir_g1_Hyper = male_g1_dir == "Hyper in NMF_G1",
    dir_g1_Hypo  = male_g1_dir == "Hypo in NMF_G1",
    dir_g2_Hyper = male_g2_dir == "Hyper in NMF_G2",
    dir_g2_Hypo  = male_g2_dir == "Hypo in NMF_G2"
  )

# Fix theme issues with ggplot2 ≥3.4
custom_theme <- ggplot2::theme_bw() +
  theme(
    axis.title.x = element_text(),
    axis.title.y = element_text()
  )

# UpSet plot
ComplexUpset::upset(
  upset_df,
  intersect = c(
    "sig_g1_g2", "sig_male_g1", "sig_male_g2",
    "dir_g1_Hyper", "dir_g1_Hypo",
    "dir_g2_Hyper", "dir_g2_Hypo"
  ),
  name = "CpG Intersections",
  width_ratio = 0.1,
  min_size = 1,
  themes = list(
    overall_sizes = custom_theme,
    intersections_matrix = custom_theme
  )
)

library(dplyr)
library(pheatmap)

# Make matrix: genes x comparisons, values = direction as numeric
heatmap_mat <- top100 %>%
  dplyr::select(gene_symbol, delta_beta, 
                g1_g2_adj = dml_g1_g2_adj,
                male_g1_slope, male_g2_slope,
                male_g1_adj, male_g2_adj) %>%
  dplyr::mutate(
    g1_g2 = ifelse(g1_g2_adj < 0.05,
                   ifelse(delta_beta > 0, 1, -1),
                   0),
    male_g1 = ifelse(male_g1_adj < 0.05,
                     ifelse(male_g1_slope > 0, 1, -1),
                     0),
    male_g2 = ifelse(male_g2_adj < 0.05,
                     ifelse(male_g2_slope > 0, 1, -1),
                     0)
  ) %>%
  dplyr::select(gene_symbol, g1_g2, male_g1, male_g2)

# Fix duplicates by adding _1 to duplicates
gene_names <- heatmap_mat$gene_symbol
dup_idx <- duplicated(gene_names)
gene_names[dup_idx] <- paste0(gene_names[dup_idx], "_1")

# Assign row names directly
rownames(heatmap_mat) <- gene_names
heatmap_mat <- heatmap_mat %>% dplyr::select(-gene_symbol)  # remove original column

# Plot
pheatmap(
  as.matrix(heatmap_mat),
  color = c("blue", "white", "red"), # Hypo, NS, Hyper
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  legend_breaks = c(-1, 0, 1),
  legend_labels = c("Hypo", "NS", "Hyper"),
  main = "Top 100 CpGs: methylation direction across comparisons"
)


library(dplyr)
library(tidyr)
library(ggplot2)
# Prepare dot plot data with unambiguous colors
# Prepare dot plot data with unambiguous colors
dot_df <- top100 %>%
  # Fix duplicates
  dplyr::mutate(
    gene_label = ifelse(duplicated(gene_symbol),
                        paste0(gene_symbol, "_1"),
                        gene_symbol)
  ) %>%
  # Select needed columns
  dplyr::select(gene_label, delta_beta,
                dml_g1_g2_adj, male_g1_slope, male_g1_adj,
                male_g2_slope, male_g2_adj) %>%
  # Pivot to long format
  tidyr::pivot_longer(
    cols = c(dml_g1_g2_adj, male_g1_adj, male_g2_adj),
    names_to = "comparison",
    values_to = "adj_p"
  ) %>%
  # Assign slope and direction per comparison
  dplyr::mutate(
    slope = case_when(
      comparison == "dml_g1_g2_adj" ~ delta_beta,
      comparison == "male_g1_adj" ~ male_g1_slope,
      comparison == "male_g2_adj" ~ male_g2_slope
    ),
    direction = case_when(
      comparison == "dml_g1_g2_adj" & adj_p < 0.05 & delta_beta > 0 ~ "Hyper in G2",
      comparison == "dml_g1_g2_adj" & adj_p < 0.05 & delta_beta < 0 ~ "Hypo in G2",
      comparison == "male_g1_adj" & adj_p < 0.05 & male_g1_slope > 0 ~ "Hyper in G1",
      comparison == "male_g1_adj" & adj_p < 0.05 & male_g1_slope < 0 ~ "Hypo in G1",
      comparison == "male_g2_adj" & adj_p < 0.05 & male_g2_slope > 0 ~ "Hypo in Male",
      comparison == "male_g2_adj" & adj_p < 0.05 & male_g2_slope < 0 ~ "Hyper in Male",
      TRUE ~ "NS"
    ),
    comparison = recode(comparison,
                        "dml_g1_g2_adj" = "G1 vs G2",
                        "male_g1_adj"   = "Male vs G1",
                        "male_g2_adj"   = "Male vs G2")
  )

# Define unambiguous colors
direction_colors <- c(
  "Hyper in G2"        = "salmon",
  "Hypo in G2"         = "skyblue",
  "Hyper in G1" = "red",
  "Hypo in G1"  = "blue",
  "Hyper in Male" = "darkred",
  "Hypo in Male"  = "darkblue",
  "NS"                  = "grey"
)

# Plot
ggplot(dot_df, aes(x = gene_label, y = comparison, color = direction)) +
  geom_point(size = 4) +
  scale_color_manual(values = direction_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 10)
  ) +
  labs(
    x = "Gene",
    y = "Comparison",
    color = "Methylation direction",
    title = "Top 100 CpGs: methylation direction across comparisons"
  )


#write.csv(top100, "machine_learning_method/Top100_module_specific_CpGs_updated_with_male.csv")

top100_general <- top100 %>%
  #dplyr::filter(promoter == "TRUE") %>%
  dplyr::filter(dml_g1_g2_adj < 0.05 & male_g1_adj < 0.05 & male_g2_adj < 0.05) %>%
  dplyr::select(gene_symbol, seqnames, methylation_direction, male_g1_dir, male_g2_dir)



library(dplyr)
library(tidyr)
library(ggplot2)

# Filter to CpGs significant in all three comparisons
top100_sig <- top100 %>%
  filter(dml_g1_g2_adj < 0.05 & male_g1_adj < 0.05 & male_g2_adj < 0.05) %>%
  select(gene_symbol, delta_beta,
         dml_g1_g2_adj, male_g1_slope, male_g1_adj,
         male_g2_slope, male_g2_adj,
         methylation_direction, male_g1_dir, male_g2_dir)

# Prepare dot plot data
dot_df <- top100_sig %>%
  # Fix duplicates in gene names
  mutate(gene_label = ifelse(duplicated(gene_symbol),
                             paste0(gene_symbol, "_1"),
                             gene_symbol)) %>%
  # Pivot to long format for plotting
  pivot_longer(
    cols = c(dml_g1_g2_adj, male_g1_adj, male_g2_adj),
    names_to = "comparison",
    values_to = "adj_p"
  ) %>%
  mutate(
    slope = case_when(
      comparison == "dml_g1_g2_adj" ~ delta_beta,
      comparison == "male_g1_adj" ~ male_g1_slope,
      comparison == "male_g2_adj" ~ male_g2_slope
    ),
    direction = case_when(
      comparison == "dml_g1_g2_adj" & delta_beta > 0 ~ "Hyper in G2",
      comparison == "dml_g1_g2_adj" & delta_beta < 0 ~ "Hypo in G2",
      comparison == "male_g1_adj" & male_g1_slope > 0 ~ "Hyper in Male vs G1",
      comparison == "male_g1_adj" & male_g1_slope < 0 ~ "Hypo in Male vs G1",
      comparison == "male_g2_adj" & male_g2_slope > 0 ~ "Hyper in Male vs G2",
      comparison == "male_g2_adj" & male_g2_slope < 0 ~ "Hypo in Male vs G2"
    ),
    comparison = recode(comparison,
                        "dml_g1_g2_adj" = "G1 vs G2",
                        "male_g1_adj"   = "Male vs G1",
                        "male_g2_adj"   = "Male vs G2")
  )

# Define colors for unambiguous hyper/hypo
direction_colors <- c(
  "Hyper in G2"         = "red",
  "Hypo in G2"          = "blue",
  "Hyper in Male vs G1" = "orange",
  "Hypo in Male vs G1"  = "cyan",
  "Hyper in Male vs G2" = "darkred",
  "Hypo in Male vs G2"  = "darkblue"
)

# Plot
ggplot(dot_df, aes(x = gene_label, y = comparison, color = direction)) +
  geom_point(size = 4) +
  scale_color_manual(values = direction_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 10)
  ) +
  labs(
    x = "Gene",
    y = "Comparison",
    color = "Methylation direction",
    title = "Top CpGs significant in all three comparisons"
  )



library(ComplexHeatmap)
library(circlize)
library(dplyr)

# ----------------------------
# 1️⃣ Prepare beta matrix
# ----------------------------
heat_mat <- beta_mat_filt[top100$Probe_ID, ]

# Make rownames unique to avoid duplicates
gene_labels_unique <- make.unique(top100$gene_symbol, sep = "_")
rownames(heat_mat) <- gene_labels_unique

# Scale rows for better visualization
heat_mat_scaled <- t(scale(t(heat_mat)))

# ----------------------------
# 2️⃣ Column annotations (samples)
# ----------------------------
col_annot <- data.frame(
  Group = ifelse(colnames(heat_mat_scaled) %in% groups_nmf$Sample,
                 groups_nmf$Comparison[match(colnames(heat_mat_scaled), groups_nmf$Sample)],
                 "Male")
)
rownames(col_annot) <- colnames(heat_mat_scaled)

col_colors <- list(
  Group = c(NMF_G1 = "salmon", NMF_G2 = "skyblue", Male = "grey50")
)

ha_col <- HeatmapAnnotation(
  df = col_annot,
  col = col_colors,
  annotation_name_side = "left"
)

# ----------------------------
# 3️⃣ Row annotations (CpGs)
# ----------------------------
row_annot <- data.frame(
  Context = ifelse(top100$promoter, "Promoter", "Gene body"),
  G1_vs_G2 = ifelse(top100$dml_g1_g2_adj < 0.05,
                    ifelse(top100$delta_beta > 0, "Hyper", "Hypo"),
                    "NS"),
  Male_vs_G1 = ifelse(top100$male_g1_adj < 0.05,
                      ifelse(top100$male_g1_slope > 0, "Hyper", "Hypo"),
                      "NS"),
  Male_vs_G2 = ifelse(top100$male_g2_adj < 0.05,
                      ifelse(top100$male_g2_slope > 0, "Hyper", "Hypo"),
                      "NS")
)
rownames(row_annot) <- gene_labels_unique

row_colors <- list(
  Context = c(Promoter = "darkgreen", `Gene body` = "brown"),
  G1_vs_G2 = c(Hyper = "red", Hypo = "blue", NS = "grey90"),
  Male_vs_G1 = c(Hyper = "red", Hypo = "blue", NS = "grey90"),
  Male_vs_G2 = c(Hyper = "red", Hypo = "blue", NS = "grey90")
)

ha_row <- rowAnnotation(
  df = row_annot,
  col = row_colors,
  show_legend = TRUE
)

# ----------------------------
# 4️⃣ Draw heatmap
# ----------------------------
Heatmap(
  heat_mat_scaled,
  name = "Beta (row-scaled)",
  top_annotation = ha_col,
  right_annotation = ha_row,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_split = col_annot$Group,  # split by NMF group / Male
  heatmap_legend_param = list(title = "Beta")
)

library(ComplexHeatmap)
library(circlize)
library(dplyr)

# ----------------------------
# 1️⃣ Prepare beta matrix
# ----------------------------
heat_mat <- beta_mat_filt[top100$Probe_ID, ]

# Make rownames unique to avoid duplicates
gene_labels_unique <- make.unique(top100$gene_symbol, sep = "_")
rownames(heat_mat) <- gene_labels_unique

# Scale rows for better visualization
heat_mat_scaled <- t(scale(t(heat_mat)))

# ----------------------------
# 2️⃣ Column annotations (samples)
# ----------------------------
col_annot <- data.frame(
  Group = ifelse(colnames(heat_mat_scaled) %in% groups_nmf$Sample,
                 groups_nmf$Comparison[match(colnames(heat_mat_scaled), groups_nmf$Sample)],
                 "Male")
)
rownames(col_annot) <- colnames(heat_mat_scaled)

col_colors <- list(
  Group = c(NMF_G1 = "salmon", NMF_G2 = "skyblue", Male = "grey50")
)

ha_col <- HeatmapAnnotation(
  df = col_annot,
  col = col_colors,
  annotation_name_side = "left"
)

# ----------------------------
# 3️⃣ Row annotations (CpGs)
# ----------------------------
row_annot <- data.frame(
  Context = ifelse(top100$promoter, "Promoter", "Gene body"),
  G1_vs_G2 = ifelse(top100$dml_g1_g2_adj < 0.05,
                    ifelse(top100$delta_beta > 0, "Hyper", "Hypo"),
                    "NS"),
  Male_vs_G1 = ifelse(top100$male_g1_adj < 0.05,
                      ifelse(top100$male_g1_slope > 0, "Hyper", "Hypo"),
                      "NS"),
  Male_vs_G2 = ifelse(top100$male_g2_adj < 0.05,
                      ifelse(top100$male_g2_slope > 0, "Hyper", "Hypo"),
                      "NS")
)
rownames(row_annot) <- gene_labels_unique

row_colors <- list(
  Context = c(Promoter = "darkgreen", `Gene body` = "darkorange"),
  G1_vs_G2 = c(Hyper = "red", Hypo = "blue", NS = "grey90"),
  Male_vs_G1 = c(Hyper = "red", Hypo = "blue", NS = "grey90"),
  Male_vs_G2 = c(Hyper = "red", Hypo = "blue", NS = "grey90")
)

ha_row <- rowAnnotation(
  df = row_annot,
  col = row_colors,
  show_legend = TRUE
)

# ----------------------------
# 4️⃣ Draw heatmap
# ----------------------------
Heatmap(
  heat_mat_scaled,
  name = "Beta (row-scaled)",
  top_annotation = ha_col,
  right_annotation = ha_row,
  cluster_rows = TRUE,               # hierarchical clustering within each row_split
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_split = col_annot$Group,    # split columns by NMF group / Male
  row_split = row_annot$Context,     # split rows by Promoter / Gene body
  heatmap_legend_param = list(title = "Beta")
)



