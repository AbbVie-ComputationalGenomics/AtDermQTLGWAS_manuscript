#step5. Change Seurat object to AnnData

setwd("/projects/abv/users/olivamx2/Projects/AD/reproduce_fanying/magma/data/hca_skin_portal")
library(Seurat)
devtools::install_github("satijalab/seurat-data", force = TRUE)
library(SeuratData)
dyn.load("/usr/local/hdf5/hdf5-1.10.6-openmpi-3.1.5/lib/libhdf5_hl.so.100")
library(hdf5r)
library(SeuratDisk)

input <- readRDS("HCA_AD_lesion_Healthy_seurat_celltypelevel1.RDS")

head(input@meta.data)
#orig.ident nCount_RNA nFeature_RNA  sample_id mad_prd Status   Site    Tissue
#AAACCTGAGAGGGCTT-1-SKN8090528 SeuratProject       2210          967 SKN8090528   FALSE Eczema lesion Epidermis
#AAACGGGAGAAGGACA-1-SKN8090528 SeuratProject       2637         1080 SKN8090528   FALSE Eczema lesion Epidermis

input <- subset(input, (Status=="Healthy")|(Status=="Eczema" & Site=="lesion"))
table(input@meta.data$Status, input@meta.data$Site)
#        lesion non_lesion
#Eczema   63512          0
#Healthy      0     195739


#option1: Lesion and H (need to manually remove celltypes with fewer than 3 cells in either one of the group)
which(table(input$Status, input$celltype_level1)<3)
[1] 21

# for level1
input_sub <- subset(input, celltype_level1!="Plasma_cells")
unique(input_sub$Status)
#[1] "Eczema"  "Healthy"
table(input_sub$Status, input_sub$celltype_level1)
#       Dendritic_cells Fibroblasts ILC_NK_cells Keratinocytes Langerhans_cell Lymphatic_endothelium Macrophages
#Eczema             5745       10236          541         11475            2419                   789        4165
#Healthy           15898       17456         6772         53844           12545                  4924       13890

#Mast_cells Melanocytes Pericytes Schwann_cells T_cells Vascular_endothelium
#Eczema          94         596      3901            53   18159                 5338
#Healthy        557        3682      5006           269   35040                25785

which(table(input_sub$Status, input_sub$celltype_level1)<3)
sample_name <- "HCA_AD_Lesion_vs_Healthy_level1"

# Change all factor to characters
str(input_sub@meta.data)

# remove the scale data
input_sub_slim <- DietSeurat(input_sub)
VariableFeatures(input_sub_slim) <- NULL


SaveH5Seurat(input_sub_slim, filename = paste0(sample_name, ".h5Seurat"))
Convert(paste0(sample_name, ".h5Seurat"), 
        dest = "h5ad")
# count is saved to adata.raw.X
# log normalized data is saved to adata.X


#### for level2 cell types
input_sub <- subset(input, celltype_level2%in%c("Plasma", "ILC1_NK", "Tc_IL13_IL22"), invert=T)
unique(input_sub$Status)
unique(input_sub$celltype_level2)
table(input_sub$Status, input_sub$celltype_level2)
which(table(input_sub$Status, input_sub$celltype_level2)<3)
sample_name <- "HCA_AD_Lesion_vs_Healthy_level2"
# Change all factor to characters
str(input_sub@meta.data)
# remove the scale data
input_sub_slim <- DietSeurat(input_sub)
VariableFeatures(input_sub_slim) <- NULL
SaveH5Seurat(input_sub_slim, filename = paste0(sample_name, ".h5Seurat"))
Convert(paste0(sample_name, ".h5Seurat"), 
        dest = "h5ad")
# count is saved to adata.raw.X
# log normalized data is saved to adata.X

####option2: Lesion only:
input_sub <- subset(input, Status=="Eczema")
unique(input_sub$Status)
sample_name <- "HCA_AD_Lesion_level1"
# Change all factor to characters
str(input_sub@meta.data)
# remove the scale data
input_sub_slim <- DietSeurat(input_sub)
VariableFeatures(input_sub_slim) <- NULL
SaveH5Seurat(input_sub_slim, filename = paste0(sample_name, ".h5Seurat"))
Convert(paste0(sample_name, ".h5Seurat"), 
        dest = "h5ad")
# count is saved to adata.raw.X
# log normalized data is saved to adata.X

####option3: Healthy only:
input_sub <- subset(input, Status=="Healthy")
unique(input_sub$Status)
sample_name <- "HCA_Healthy_level1"
# Change all factor to characters
str(input_sub@meta.data)
# remove the scale data
input_sub_slim <- DietSeurat(input_sub)
VariableFeatures(input_sub_slim) <- NULL
SaveH5Seurat(input_sub_slim, filename = paste0(sample_name, ".h5Seurat"))
Convert(paste0(sample_name, ".h5Seurat"), 
        dest = "h5ad")
# count is saved to adata.raw.X
# log normalized data is saved to adata.X
