#!/usr/bin/env /usr/apps/general/shpc/modules/shpc-registry/R/4.2.2/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
print(args[1])
# test if there is at least one argument: if not, return an error
if (length(args)==0 || !args[1]%in%c('arm','leg','axilla','back','face','scalp','acral')) {
  stop("At least one argument - arm,leg,scalp,back,face,scalp,acral - must be supplied", call.=FALSE)
}

library(Seurat)
tissue <- args[1] 
keratinocyte_subtype <- args[2] 
print(tissue)
print(keratinocyte_subtype)

cs <- read.table('../../data/other/coloc_score.tsv',header=T)
file <- paste0('/scratch/users/olivamx2/AD/data/Umich_epiderm_atlas_seurat_v3.',tissue,'.RDS')
obj <- readRDS(file)

obj@meta.data$celltype <- gsub(' ','_',obj@meta.data$celltype)
#keratinocyte_subtypes <- unique(obj@meta.data$celltype)[grep('Keratinocytes',unique(obj@meta.data$celltype))]
obj@meta.data$KeratinocyteBinary <- unlist(lapply(obj@meta.data$celltype,function(x) if(!grepl('Keratinocytes',x)) {print("Others")} else  {print(x)}))

table(obj@meta.data$KeratinocyteBinary)
Idents(obj)= obj$KeratinocyteBinary
#for (keratinocyte_subtype in keratinocyte_subtypes) {
  KeratinocyteDE.markers <- FindMarkers(obj, ident.1 = keratinocyte_subtype, ident.2 = "Others", verbose = T, logfc.threshold = 0, min.pct = 0.1)
  KeratinocyteDE.markers <- KeratinocyteDE.markers[order(abs(KeratinocyteDE.markers$avg_log2FC),decreasing=T),]
  write.table(file=paste0(tissue,'.',keratinocyte_subtype,'.Markers.tsv'),KeratinocyteDE.markers,quote=F)
#}

quit(save="no")
