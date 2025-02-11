#plot heatmap of MAGMA cell type scores for MAGMA enriched genes 
#170 GWAS genes with MAGMA Bonferonni p-value < 2.71e-06
#pearson correlation complete - full table highlight program and cell_origin

install.packages("pheatmap")
library(pheatmap)
library(dplyr)
library(tidyverse)

#n=170 genes (177 genes, remove 7 genes with all NA values for gene program scores)
data_170 <- read.delim("data/Heatmap_celltype_values_170.txt", header = T, row.names = "GeneID")
my_sample_col_170 <- read.delim("data/Heatmap_celltype_annotation_full.txt", header = T, row.names = "Cell_type")
ann_colors <- list(Program = c(H = 'green', LS = 'darkviolet', DP = 'deeppink2'),Cell_origin = c(immune='gold', "non-immune"='cyan2'))

#drop additional genes for plotting purposes, any gene with total summary of gene programs < 1.0
data_146 <- read.delim("data/Heatmap_celltype_values_146.txt", header = T, row.names = "GeneID")
head(data_146)
#drop extra cols
data_146 <- data_146 %>% select(-X,-Novel)
data_146 <- data_146 %>% select(-X.1,-SUM.L...DP)
head(data_146)

#correlation
pheatmap(data_146,
         annotation_col = my_sample_col_170,
         color = colorRampPalette(c('darkblue', 'black', 'yellow'))(100),
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation",
         clustering_method = "complete",
         annotation_colors = ann_colors,
         show_rownames = FALSE,
         annotation_names_col= FALSE,
         treeheight_row = 0,
         cutree_cols = 6,
         cellwidth = 10,
         main = "Cell Type Enrichment (146 genes) - Gene Program Scores - Pearson")

ggsave("extfig.celltype_cluster_heatmap.pdf", width=5, height=5, limitsize = FALSE) 
