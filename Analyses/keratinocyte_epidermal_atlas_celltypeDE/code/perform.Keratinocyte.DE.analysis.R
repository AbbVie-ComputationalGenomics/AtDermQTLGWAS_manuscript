#!/usr/bin/env /usr/apps/general/shpc/modules/shpc-registry/R/4.2.2/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
print(args[1])
# test if there is at least one argument: if not, return an error
if (length(args)==0 || !args[1]%in%c('arm','leg','axilla','back','face','scalp','acral')) {
  stop("At least one argument - arm,leg,scalp,back,face,scalp,acral - must be supplied", call.=FALSE)
}

library(Seurat)
tissue <- args[1] 
print(tissue)

cs <- read.table('../../data/other/coloc_score.tsv',header=T)
file=paste0('/scratch/users/olivamx2/AD/data/Umich_epiderm_atlas_seurat_v3.',tissue,'.RDS')
obj <- readRDS(file)
obj@meta.data$KeratinocyteBinary <- gsub('TRUE','Keratinocyte',grepl('Keratinocyte',obj@meta.data$celltype))
obj@meta.data$KeratinocyteBinary <- gsub('FALSE','Other',obj@meta.data$KeratinocyteBinary)
table(obj@meta.data$KeratinocyteBinary)
Idents(obj)= obj$KeratinocyteBinary
KeratinocyteDE.markers <- FindMarkers(obj, ident.1 = "Keratinocyte", ident.2 = "Other", verbose = T, logfc.threshold = 0, min.pct = 0.1)
KeratinocyteDE.markers <- KeratinocyteDE.markers[order(abs(KeratinocyteDE.markers$avg_log2FC),decreasing=T),]
write.table(file=paste0(tissue,'.Keratinocyte.Markers.tsv'),KeratinocyteDE.markers,quote=F)
pdf(paste0(tissue,'.Keratinocyte.DE.coloc.weak.strong.pdf'))
plot(FeaturePlot(obj, features = subset(cs,colocScore%in%'Weak' & Include%in%'Y' ,select = 'CandidateGene')[,1]))
plot(FeaturePlot(obj, features = subset(cs,colocScore%in%'Strong' & Include%in%'Y' ,select = 'CandidateGene')[,1]))
dev.off()

quit(save="no")


