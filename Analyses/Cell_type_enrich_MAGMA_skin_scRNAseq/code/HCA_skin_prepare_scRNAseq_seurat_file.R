# step4. Prepare HCA skin portal scRNAseq seurat file
wget https://zenodo.org/record/4569496/files/submission_210120.h5ad?download=1
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-8142/E-MTAB-8142.sdrf.txt

setwd("/projects/abv/users/olivamx2/Projects/AD/reproduce_fanying/magma/data/hca_skin_portal")

#import H5AD object as seurat object
renv::init()
renv::install("reticulate")
renv::use_python()
#pkgs <- c("reticulate","ggplot2","Seurat","scater")
#renv::install(pkgs)
py_pkgs <- c("scanpy")
reticulate::py_install("scanpy")
renv::install("BiocManager")
renv::install("bioc::BiocGenerics")
renv::install("bioc::S4Vectors")
renv::install("bioc::SummarizedExperiment")
renv::install("bioc::DelayedMatrixStats")
renv::install("bioc::beachmat")
renv::install("bioc::BiocNeighbors")
renv::install("bioc::scater")

suppressPackageStartupMessages({
  library("reticulate")
  library("ggplot2")
  library("scater")
  library("Seurat")
})

sc <- import("scanpy")
adata <- sc$read_h5ad("submission_210120.h5ad")
adata
counts <- t(adata$X)

head(adata$obs$sample_id)
head(adata$var$'feature_types-SKN8090524')

#identical(adata$var$'gene_ids-SKN8090524', adata$var$'gene_ids-4820STDY7389014')
#gene_names <- adata$var$'gene_ids-SKN8090524'
colnames(counts) <- adata$obs_names$to_list()
rownames(counts) <- adata$var_names$to_list()
saveRDS(counts, "counts_adult.RDS")

#26 groups * 17369 cells in each group. The matrix is too big for file transformation
desired_length <- 26
group_num <- 17369 
mat_list <- vector(mode="list", length= desired_length)
for(i in c(1:desired_length)) {
  print(i)
  test <- counts[, (group_num*(i-1)+1):(group_num*i)]
  mat <- Matrix::Matrix(as.matrix(test), sparse=T)
  mat_list[[i]] <- mat
}

counts <- do.call(cbind, mat_list)
seurat <- CreateSeuratObject(counts)
seurat <- AddMetaData(seurat, adata$obs)
Idents(seurat) <- "full_clustering"
saveRDS(seurat, "AD_HCA_adult_seurat.RDS")

# remember add meta data for patient information
# remember to do data normalization
#seurat <- readRDS("AD_HCA_adult_seurat.RDS")
patient_info <- read.table("E-MTAB-8142.sdrf.txt", header=T, stringsAsFactors=F, sep="\t")
patient_ID <- patient_info$Characteristics.individual.[match(seurat@meta.data$sample_id, patient_info$Source.Name)]
seurat@meta.data$patient_ID <- patient_ID
seurat <- NormalizeData(seurat)
saveRDS(seurat, "AD_HCA_adult_seurat_normalized.RDS")

#seurat <- readRDS("AD_HCA_adult_seurat_normalized.RDS")
#head(Idents(seurat))
#seurat_healthy <- subset(seurat, Status=="Healthy")
#seurat_eczema_lesion <- subset(seurat, Status=="Eczema" & Site=="lesion")
#seurat_healthy_eczema <- subset(seurat, Status%in%c("Healthy", "Eczema"))

#Cell type level 1 and 2 clustering
BiocManager::install("edgeR")
BiocManager::install("scran")
library(edgeR)
library(magrittr)
library(tibble)
library(scran)
library(Libra)

input <- readRDS("AD_HCA_adult_seurat_normalized.RDS")
head(input@meta.data)

# Clubbing cells together:
cell_type_level1 <- input@meta.data$full_clustering
cell_type_level1 <- as.character(cell_type_level1)
unique(cell_type_level1)

cell_type_level1[which(cell_type_level1%in%c("Differentiated_KC",
                                             "Undifferentiated_KC",
                                             "Differentiated_KC*",
                                             "Proliferating_KC"))] <- "Keratinocytes"

cell_type_level1[which(cell_type_level1%in%c("LC_4",
                                             "LC_3",
                                             "LC_2",
                                             "LC_1"))] <- "Langerhans_cell"

cell_type_level1[which(cell_type_level1%in%c("Tc17_Th17",
                                             "Tc_IL13_IL22",
                                             "Th",
                                             "Tc",
                                             "Treg"))] <- "T_cells"

cell_type_level1[which(cell_type_level1%in%c("MigDC",
                                             "moDC_2",
                                             "DC1",
                                             "moDC_3",
                                             "moDC_1",
                                             "DC2"))] <- "Dendritic_cells"

cell_type_level1[which(cell_type_level1%in%c("ILC1_NK",
                                             "ILC2",
                                             "NK",
                                             "ILC1_3"))] <- "ILC_NK_cells"

cell_type_level1[which(cell_type_level1%in%c("LE2",
                                             "LE1"))] <- "Lymphatic_endothelium"

cell_type_level1[which(cell_type_level1%in%c("VE1",
                                             "VE2",
                                             "VE3"))] <- "Vascular_endothelium"

cell_type_level1[which(cell_type_level1%in%c("F1",
                                             "F2",
                                             "F3"))] <- "Fibroblasts"

cell_type_level1[which(cell_type_level1%in%c("Mono_mac",
                                             "Macro_1",
                                             "Macro_2",
                                             "Inf_mac"))] <- "Macrophages"

cell_type_level1[which(cell_type_level1%in%c("Schwann_1" ,
                                             "Schwann_2"))] <- "Schwann_cells"

cell_type_level1[which(cell_type_level1%in%c("Pericyte_1" ,
                                             "Pericyte_2"))] <- "Pericytes" 

cell_type_level1[which(cell_type_level1%in%c("Mast_cell"))] <- "Mast_cells" 

cell_type_level1[which(cell_type_level1%in%c("Plasma"))] <- "Plasma_cells" 

cell_type_level1[which(cell_type_level1%in%c("Melanocyte"))] <- "Melanocytes" 

input@meta.data$celltype_level1 <- cell_type_level1


# Re-name cells in the orignial cluster
cell_type_level2 <- input@meta.data$full_clustering
cell_type_level2 <- as.character(cell_type_level2)
unique(cell_type_level2)

cell_type_level2[which(cell_type_level2%in%c("F3"))] <- "Fibroblast_3"
cell_type_level2[which(cell_type_level2%in%c("F2"))] <- "Fibroblast_2"
cell_type_level2[which(cell_type_level2%in%c("F1"))] <- "Fibroblast_1"

cell_type_level2[which(cell_type_level2%in%c("LC_2"))] <- "Langerhans_cell_2"
cell_type_level2[which(cell_type_level2%in%c("LC_3"))] <- "Langerhans_cell_3"
cell_type_level2[which(cell_type_level2%in%c("LC_1"))] <- "Langerhans_cell_1"
cell_type_level2[which(cell_type_level2%in%c("LC_4"))] <- "Langerhans_cell_4"

cell_type_level2[which(cell_type_level2%in%c("VE1"))] <- "Vascular_endothelium_1"
cell_type_level2[which(cell_type_level2%in%c("VE2"))] <- "Vascular_endothelium_2"
cell_type_level2[which(cell_type_level2%in%c("VE3"))] <- "Vascular_endothelium_3"

cell_type_level2[which(cell_type_level2%in%c("LE2"))] <- "Lymphatic_endothelium_2"
cell_type_level2[which(cell_type_level2%in%c("LE1"))] <- "Lymphatic_endothelium_1"

cell_type_level2[which(cell_type_level2%in%c("Macro_1"))] <- "Macrophage_1"
cell_type_level2[which(cell_type_level2%in%c("Macro_2"))] <- "Macrophage_2"
cell_type_level2[which(cell_type_level2%in%c("Inf_mac"))] <- "Inflam_Macrophage"
cell_type_level2[which(cell_type_level2%in%c("Mono_mac"))] <- "Mono_Macrophage"

unique(cell_type_level2)

input@meta.data$celltype_level2 <- cell_type_level2

saveRDS(input, "AD_HCA_adult_seurat_celltypelevel1.RDS")

#input <- readRDS("AD_HCA_adult_seurat_celltypelevel1.RDS")

str(input@meta.data)
input@meta.data$orig.ident <- as.character(input@meta.data$orig.ident)
input@meta.data$sample_id <- as.character(input@meta.data$sample_id)
input@meta.data$Status <- as.character(input@meta.data$Status)
input@meta.data$Site <- as.character(input@meta.data$Site)
input@meta.data$Tissue <- as.character(input@meta.data$Tissue)
input@meta.data$Enrichment <- as.character(input@meta.data$Enrichment)
input@meta.data$Location <- as.character(input@meta.data$Location)
input@meta.data$Sex <- as.character(input@meta.data$Sex)
input@meta.data$stage <- as.character(input@meta.data$stage)
input@meta.data$full_clustering <- as.character(input@meta.data$full_clustering)

str(input@meta.data)

unique(input@meta.data$patient_ID)

input_sub <- subset(input, Status=="Healthy"|(Status=="Eczema" & Site=="lesion"))
table(input_sub@meta.data$Status, input_sub@meta.data$Site)
#        lesion non_lesion
#Eczema   63512          0
#Healthy      0     195739

saveRDS(input_sub, "HCA_AD_lesion_Healthy_seurat_celltypelevel1.RDS")