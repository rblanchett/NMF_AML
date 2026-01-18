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
# 5. Read in SAVED nmf data
# ------------------------------------

nmf_fit <- readRDS("machine_learning_method/nmf_fit_20251231.rds")
nmf_groups <- readRDS("machine_learning_method/nmf_group_20251231.rds")
nmf_df <- readRDS("machine_learning_method/nmf_df_20251231.rds")


H <- coef(nmf_fit)
consensusmap(nmf_fit)
heatmap(H, scale = "row")

W <- basis(nmf_fit)

W_df <- data.frame(
  CpG = rownames(W),
  Component1 = W[,1],
  Component2 = W[,2],
  chr = probes_EPIC$CpG_chrm[match(rownames(W), probes_EPIC$Probe_ID)]
)

aggregate(Component1 ~ chr, data = W_df, mean)
aggregate(Component2 ~ chr, data = W_df, mean)

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
top_n <- 100  # top 100 CpGs per component

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
    title = "Genome-wide distribution of NMF loadings (Group 1)",
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
    title = "Genome-wide distribution of NMF loadings (Group 2)",
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
### differences in loadings plots
###################

ggplot(diff_df, aes(x = Diff)) +
  geom_histogram(bins = 200, fill = "darkred", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(
    title = "CpG loading difference (Module 1 – Module 2)",
    x = "Loading difference",
    y = "CpG count"
  )

ggplot(df_loadings, aes(x = start, y = Module_1 - Module_2)) +
  geom_point(alpha = 0.25, size = 0.5) +
  facet_wrap(~ chr, scales = "free_x") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Genome-wide CpG loading difference (Module 1 − Module 2)",
    x = "Genomic coordinate",
    y = "Loading difference"
  )

top_diff <- diff_df %>%
  slice_max(order_by = abs(Diff), n = 100)

stopifnot("CpG" %in% colnames(top_diff))
stopifnot("Probe_ID" %in% colnames(df_loadings))

top_diff_plot <- top_diff %>%
  left_join(
    df_loadings %>%
      dplyr::select(Probe_ID, start, chr),
    by = c("CpG" = "Probe_ID")
  )

# Sanity check
stopifnot(!any(is.na(top_diff_plot$start)))

ggplot(df_loadings, aes(x = start, y = Module_1 - Module_2)) +
  geom_point(alpha = 0.2, size = 0.5) +
  geom_point(
    data = top_diff_plot,
    aes(x = start, y = Diff),
    color = "red",
    size = 1.5
  ) +
  facet_wrap(~ chr, scales = "free_x") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Top 100 module-specific CpGs by loading difference",
    x = "Genomic coordinate",
    y = "Loading difference (Module 1 − Module 2)"
  )


############
## Annotate these CpGs
############

library(GenomicRanges)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

setdiff(gr_sub$Probe_ID, rownames(W_sub))   # in GRanges but not in W_sub
setdiff(rownames(W_sub), gr_sub$Probe_ID)   # in W_sub but not in GRanges
common_cpgs <- intersect(gr_sub$Probe_ID, rownames(W_sub))
W_sub <- W_sub[common_cpgs, , drop = FALSE]
gr_cpg <- gr_sub[match(common_cpgs, gr_sub$Probe_ID)]
stopifnot(gr_cpg$Probe_ID == rownames(W_sub))  # should now pass

stopifnot(gr_cpg$Probe_ID == rownames(W_sub))

gr_cpg$Diff <- gr_cpg$Module_1 - gr_cpg$Module_2

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr <- genes(txdb)

idx <- nearest(gr_cpg, genes_gr)

gr_cpg$gene_id <- genes_gr$gene_id[idx]

gr_cpg$gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = gr_cpg$gene_id,
  keytype = "ENTREZID",
  column = "SYMBOL"
)

promoters_gr <- promoters(genes_gr, upstream = 2000, downstream = 500)

gr_cpg$promoter <- overlapsAny(gr_cpg, promoters_gr)
gr_cpg$gene_body <- overlapsAny(gr_cpg, genes_gr)

gr_cpg$context <- ifelse(
  gr_cpg$promoter, "Promoter",
  ifelse(gr_cpg$gene_body, "Gene body", "Intergenic")
)

gr_cpg$island <- annot_sub$Relation_to_Island[
  match(gr_cpg$Probe_ID, annot_sub$Probe_ID)
]

# Make sure the sample order matches
group_vec <- nmf_df$NMF_Group[match(colnames(beta_mat_xx), rownames(nmf_df))]

stopifnot(length(group_vec) == ncol(beta_mat_xx))


mean_Group1 <- rowMeans(beta_mat_xx[, group_vec == 1, drop = FALSE], na.rm = TRUE)
mean_Group2 <- rowMeans(beta_mat_xx[, group_vec == 2, drop = FALSE], na.rm = TRUE)

beta_means <- data.frame(
  CpG = rownames(beta_mat_xx),
  mean_Group1 = mean_Group1,
  mean_Group2 = mean_Group2,
  delta_beta = mean_Group2 - mean_Group1
)

beta_means$delta_beta <- beta_means$mean_Group2 - beta_means$mean_Group1

gr_cpg$delta_beta <- beta_means$delta_beta[
  match(gr_cpg$Probe_ID, beta_means$CpG)
]

gr_cpg$methylation_direction <- ifelse(
  gr_cpg$delta_beta > 0,
  "Hyper in Group2",
  "Hypo in Group2"
)

df <- as.data.frame(gr_cpg)

ggplot(df, aes(x = start, y = Diff, color = methylation_direction)) +
  geom_point(alpha = 0.4, size = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ seqnames, scales = "free_x") +
  theme_bw() +
  labs(
    title = "Module-specific CpGs with methylation direction overlay",
    x = "Genomic position (hg38)",
    y = "Module 1 − Module 2 loading"
  )

top_diff <- df %>%
  arrange(desc(abs(Diff))) %>%
  slice_head(n = 50)

ggplot(df, aes(start, Diff)) +
  geom_point(alpha = 0.3) +
  geom_point(data = top_diff, color = "red", size = 1.2) +
  geom_text(
    data = top_diff,
    aes(label = gene_symbol),
    size = 3,
    vjust = -0.5
  ) +
  facet_wrap(~ seqnames, scales = "free_x") +
  theme_bw()
df %>%
  count(context, island, methylation_direction) %>%
  ggplot(aes(context, n, fill = methylation_direction)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(title = "Genomic context of module-specific CpGs")


group <- factor(nmf_groups)
group2 <- ifelse(group  == "1", "Group1", "Group2")
group2 <- as.factor(group2)

stopifnot(nlevels(group2) == 2)

library(limma)

design <- model.matrix(~ 0 + group2)
colnames(design) <- levels(group2)
design
fit <- lmFit(beta_filt, design)

contrast <- makeContrasts(
  Group2_vs_Group1 = Group2 - Group1,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

dml <- topTable(
  fit2,
  coef = "Group2_vs_Group1",
  number = Inf,
  adjust.method = "BH"
)

dml$CpG <- rownames(dml)

dml$direction <- ifelse(
  dml$logFC > 0,
  "Hyper in Group2",
  "Hypo in Group2"
)

dml$sig <- dml$adj.P.Val < 0.05

mcols(gr_cpg)$logFC <- dml$logFC[
  match(gr_cpg$Probe_ID, dml$CpG)
]

mcols(gr_cpg)$adj.P.Val <- dml$adj.P.Val[
  match(gr_cpg$Probe_ID, dml$CpG)
]

#######

gr_cpg$direction <- ifelse(
  gr_cpg$logFC > 0,
  "Hyper in Group2",
  "Hypo in Group2"
)

gr_cpg$sig <- gr_cpg$adj.P.Val < 0.05

df <- as.data.frame(gr_cpg)
df_sig <- df %>%
  filter(adj.P.Val < 0.05)

df_sig$seqnames <- factor(df_sig$seqnames, levels = chr_levels)
# Define chromosome order: 1-22 + X
# Define chromosome order: 1-22 + X
#chr_levels <- c(paste0("chr", 1:22), "X")
# Remove 'chr' prefix if present
#df_sig$chr_clean <- sub("chr", "", df_sig$seqnames)
#df_sig$chr_clean <- factor(df_sig$chr_clean, levels = chr_levels)

ggplot(df_sig,
       aes(x = start,
           y = Diff,
           color = delta_beta)) +
  geom_point(size = 0.7) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Δβ (Group2 − Group1)"
  ) +
  facet_wrap(~ seqnames, scales = "free_x") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Module-specific CpGs with significant hyper/hypomethylation",
    x = "Genomic position (hg38)",
    y = "Module 1 − Module 2 loading"
  )


########################

library(dplyr)
library(ggplot2)

# Ensure gr_cpg has:
# Diff = Module1 - Module2
# logFC, adj.P.Val, direction, sig
gr_cpg$Diff <- gr_cpg$Module_1 - gr_cpg$Module_2

# Convert to data.frame for plotting
df_plot <- as.data.frame(gr_cpg)

# Identify top 100 per chromosome based on absolute Diff
top_per_chr <- df_plot %>%
  group_by(seqnames) %>%
  slice_max(order_by = abs(Diff), n = 100, with_ties = FALSE) %>%
  ungroup()

# Define methylation color based on logFC and significance
top_per_chr$color <- ifelse(
  top_per_chr$sig & top_per_chr$logFC > 0, "Hypermethylated in Group2",
  ifelse(top_per_chr$sig & top_per_chr$logFC < 0, "Hypomethylated in Group2", "Not significant")
)

# Define chromosome order factor
chr_levels <- c(paste0("chr", 1:22), "chrX")
top_per_chr$seqnames <- factor(top_per_chr$seqnames, levels = chr_levels)

# Plot
ggplot(top_per_chr, aes(x = start, y = Diff, color = color)) +
  geom_point(alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c(
      "Hypermethylated in Group2" = "red",
      "Hypomethylated in Group2" = "blue",
      "Not significant" = "gray"
    )
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ seqnames, scales = "free_x") +
  theme_bw() +
  labs(
    title = "Top 100 module-specific CpGs per chromosome",
    subtitle = "Colored by statistically significant hyper/hypomethylation in Group 2",
    x = "Genomic coordinate",
    y = "Module 1 − Module 2 loading",
    color = "Methylation direction"
  )

#################################################

library(dplyr)

# Compute loading difference
gr_cpg$Diff <- gr_cpg$Module_1 - gr_cpg$Module_2

# Convert to data.frame
df <- as.data.frame(gr_cpg)

# Make sure methylation direction is clear (using statistically significant DMLs)
df$methylation_direction <- ifelse(
  df$sig & df$logFC > 0, "Hypermethylated in Group2",
  ifelse(df$sig & df$logFC < 0, "Hypomethylated in Group2", "Not significant")
)

# Select top 100 CpGs by absolute loading difference
top100_df <- df %>%
  arrange(desc(abs(Diff))) %>%
  slice_head(n = 100) %>%
  dplyr::select(
    Probe_ID,
    seqnames, start, end,
    gene_symbol,
    context, promoter, gene_body, 
    Module_1, Module_2, Diff,
    delta_beta, logFC, adj.P.Val, methylation_direction
  )

# Optional: save table
#write.csv(top100_df, "machine_learning_method/Top100_module_specific_CpGs_annotated.csv", row.names = FALSE)

top100_filt <- top100_df %>%
  filter(adj.P.Val < 0.05) %>%
  filter(promoter == TRUE)



library(ggplot2)

# Order chromosomes
library(ggrepel)

# Ensure chromosomes are ordered
chr_levels <- c(paste0("chr", 1:22), "chrX")
top100_df$seqnames <- factor(top100_df$seqnames, levels = chr_levels)

# Label only top 20 by absolute Diff
top_labels <- top100_df %>%
  arrange(desc(abs(Diff))) %>%
  slice_head(n = 20)

ggplot(top100_df, aes(x = start, y = Diff, color = methylation_direction, shape = context)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(
    data = top_labels,
    aes(label = gene_symbol),
    size = 3,
    max.overlaps = 50,
    nudge_y = 0.05
  ) +
  scale_color_manual(
    values = c(
      "Hypermethylated in Group2" = "red",
      "Hypomethylated in Group2" = "blue",
      "Not significant" = "gray"
    )
  ) +
  facet_wrap(~ seqnames, scales = "free_x") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_blank(),      # remove x-axis numbers
    axis.ticks.x = element_blank()      # remove x-axis ticks
  ) +
  labs(
    title = "Top 100 module-specific CpGs",
    subtitle = "Annotated with methylation direction and genomic context",
    x = NULL,
    y = "Module 1 − Module 2 loading",
    color = "Methylation direction",
    shape = "Genomic context"
  )
