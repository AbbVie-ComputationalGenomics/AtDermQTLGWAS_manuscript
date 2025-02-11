#!/usr/bin/env /usr/apps/general/shpc/modules/shpc-registry/R/4.2.2/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
print(args[1])
# test if there is at least one argument: if not, return an error
if (length(args)==0 || !args[1]%in%c('arm','leg','axilla','back','face','scalp','acral')) {
  stop("At least one argument - arm,leg,scalp,back,face,scalp,acral - must be supplied", call.=FALSE)
}

library(Seurat)
library(SeuratData) # Not installed?
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(reshape2)
library(tidyverse)
library(RColorBrewer)
tissue <- args[1] 
print(tissue)

obj <- readRDS('data/Umich_epiderm_atlas_seurat_v3.RDS') # read Seurat object

table(obj@meta.data[c('Body_site','Tissue','bodysite')]) # Display sample sites

# Subset cells for a bodysite
Idents(obj)= obj$bodysite # Idents: Get, set, and manipulate an object's identity classes. Will be defined by "bodysite" column
obj = subset(obj, idents = tissue) # Subset by bodysite=$tissue entries.
Idents(obj)= obj$Tissue
obj = subset(obj, idents = "skin") # Subset skin cells, discard epi and dermis cells
filename=paste0('/scratch/users/olivamx2/AD/data/Umich_epiderm_atlas_seurat_v3.',tissue,'.RDS')
table(obj@meta.data[c('Body_site','Tissue','bodysite')]) # Display sample sites
saveRDS(file = filename, obj)

quit(save="N")

# If I split the object, I need to normalize again.
Diseases <- FindVariableFeatures(scalp) # Normalize, 
Diseases <- ScaleData(scalp)
Diseases <- RunPCA(scalp, assay = "GeoMx", verbose = FALSE)
Diseases <- FindNeighbors(Diseases, reduction = "pca", dims = seq_len(30))
Diseases <- FindClusters(Diseases, verbose = FALSE)
Diseases <- RunUMAP(Diseases, reduction = "pca", dims = seq_len(30))
DimPlot(Diseases, reduction = "umap", label = TRUE, group.by = "segmentclass")

#####################


Idents(Diseases) <- Diseases$segmentclass
de_genes <- FindMarkers(Diseases, ident.1 = "AD_plaque", ident.2 = "AD_plaque_close", logfc.threshold = -Inf)
de_genes <- de_genes[order(abs(de_genes$avg_log2FC), decreasing = TRUE),]


