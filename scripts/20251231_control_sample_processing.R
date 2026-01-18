setwd("C://Users/reida/Documents/Aim_4/TARGET_AML/machine_learning_method/")

#never_clin <- read.table(gzfile("GSE124413_series_matrix.txt.gz"), sep = "\t")

library(GEOquery)

gse_never <- getGEO("GSE124413", GSEMatrix = TRUE)
gse_adult <- getGEO("GSE169156", GSEMatrix = TRUE)

df_never <- read.csv(gzfile("GSE124413_matrix_processed.csv.gz"))
df_adult <- read.rable(gzfile("GSE169156_SJLIFE_IlluminaEPIC_2052samples_GEO_03162021_processed.txt.gz"))



