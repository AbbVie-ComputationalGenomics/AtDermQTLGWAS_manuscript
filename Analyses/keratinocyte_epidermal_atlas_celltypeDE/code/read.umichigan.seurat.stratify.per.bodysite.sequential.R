library(Seurat)
library(SeuratData) # Not installed?
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(reshape2)
library(tidyverse)
library(RColorBrewer)
#library(SCArray.sat)

obj <- readRDS('data/Umich_epiderm_atlas_seurat_v3.RDS') # read Seurat object

table(obj@meta.data[c('Body_site','Tissue','bodysite')]) # Display sample sites

# Subset acral cells, from palm + sole skin (discard epi and dermis cells)
Idents(obj)= obj$bodysite # Idents: Get, set, and manipulate an object's identity classes. Will be defined by "bodysite" column
acral = subset(obj, idents = "acral") # Subset by bodysite='acral' entries.
Idents(acral)= acral$Tissue
acral = subset(acral, idents = "skin") # Subset skin cells, discard epi and dermis cells
saveRDS(file = '/scratch/users/olivamx2/AD/data/Umich_epiderm_atlas_seurat_v3.acral.RDS',acral)
acral <- NULL

# Subset arm cells, from arm skin (include  Body_site=unspecified cells and Body_site=arm cells)
Idents(obj)= obj$bodysite
arm = subset(obj, idents = "arm") # Subset by bodysite='arm' entries.
Idents(arm)= arm$Tissue
arm = subset(arm, idents = "skin") # Subset skin cells, discard epi and dermis cells
saveRDS(file = '/scratch/users/olivamx2/AD/data/Umich_epiderm_atlas_seurat_v3.arm.RDS',arm)
arm <- NULL

# Subset axilla skin cells
Idents(obj)= obj$bodysite
axilla = subset(obj, idents = "axilla") # Subset by bodysite='acral' entries.
saveRDS(file = '/scratch/users/olivamx2/AD/data/Umich_epiderm_atlas_seurat_v3.axilla.RDS',axilla)
axilla <- NULL

# Subset back skin cells
back = subset(obj, idents = "back") # Subset by bodysite='back' entries.
saveRDS(file = '/scratch/users/olivamx2/AD/data/Umich_epiderm_atlas_seurat_v3.back.RDS',back)
back <- NULL

# Subset face skin cells
face = subset(obj, idents = "face") # Subset by bodysite='face' entries.
saveRDS(file = '/scratch/users/olivamx2/AD/data/Umich_epiderm_atlas_seurat_v3.face.RDS',face)
face <- NULL

# Subset leg cells, from thigh + hip + unspecified skin (discard epi and dermis cells)
leg = subset(obj, idents = "leg") # Subset by bodysite='leg' entries.
Idents(leg)= leg$Tissue
leg = subset(leg, idents = "skin") # Subset skin cells, discard epi and dermis cells
saveRDS(file = '/scratch/users/olivamx2/AD/data/Umich_epiderm_atlas_seurat_v3.leg.RDS',leg)
leg <- NULL

# Subset scalp skin cells
Idents(obj)= obj$bodysite
scalp = subset(obj, idents = "scalp") # Subset by bodysite='scalp' entries.
saveRDS(file = '/scratch/users/olivamx2/AD/data/Umich_epiderm_atlas_seurat_v3.scalp.RDS',scalp)
scalp <- NULL

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


