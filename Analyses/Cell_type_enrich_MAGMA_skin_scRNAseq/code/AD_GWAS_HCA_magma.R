#Code to run magma-sc cell type enrichment analysis on AD GWAS results and scRNAseq HCA skin datasets

#############################################################################################
#step 1. prepare AD GWAS results
#############################################################################################
library(data.table)
library(tidyverse)
input_GWAS <- fread("/LDSC_Immune_ATAC/data/AD_GWAS_EUR_meta_LDSC.txt")
head(input_GWAS)
rsID CHR   POS REF ALT     PVAL      BETA         Z      N
1:  rs367896724   1 10177   A  AC 0.232566 -0.025689 -1.193723 221190
2:  rs540431307   1 10235   T  TA 0.968801 -0.033948 -0.039095 221190
dim(input_GWAS)
[1] 58440810        9

unique(input_GWAS$CHR)
#[1]  1  2 12  7 10  9 42 11  3 19  4 22  5  6  8 20 13 15 16 17 18 21

#Prepare the input SNPs 
output <- input_GWAS[ , c("rsID", "CHR", "POS")]
head(output)
colnames(output) <- c("SNP", "CHR", "BP")
output <- na.omit(output)
write.table(output, "data/AD_EUR_SNP_LOC_MAGMA.txt", quote=F, col.names=F, row.names=F, sep="\t")

# Prepare the input SNP with pvalues
output <- input_GWAS[, c("rsID", "PVAL", "N")]
colnames(output) <- c("SNP", "P", "N")
head(output)
output <- na.omit(output)
write.table(output, "data/AD_EUR_pvalue_MAGMA.txt", quote=F, col.names=T, row.names=F, sep="\t")

#############################################################################################
#step 2. install magma v1.10
#############################################################################################
wget https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.10_source.zip
wget https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.zip
wget https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip

unzip magma_v1.10_source.zip
make
./magma
#No arguments specified. Please consult manual for usage instructions.
#Exiting MAGMA. Goodbye.

unzip NCBI37.3.zip
unzip g1000_eur.zip

#############################################################################################
#step 3. Prepare GWAS magma gene annotation
#############################################################################################
projects/abv/users/olivamx2/Projects/AD/reproduce_fanying/magma/src
# for gene_body:
./magma --annotate --snp-loc /data/AD_EUR_SNP_LOC_MAGMA.txt --gene-loc /data/NCBI37.3.gene.loc --out /data/AD_EUR_NCBI37.3

#output: 
AD_UKBB_EAGLE_NCBI37.3.genes.annot 

#############################################################################################
#step 4. Run GWAS magma gene analysis
#############################################################################################
# models: mean
#gene_mean
#Test of mean SNP association
#Uses sum of squared SNP Z-statistics as test statistic
#Uses simulation-based empirical p-values as backstop for genes where numerical integration fails
#snp=mean is more attuned to the main SNP association, and skews towards associations in areas of higher LD in a gene
./magma --gene-model snp-wise=mean --bfile /data/EUR/g1000_eur --gene-annot /data/AD_EUR_NCBI37.3.genes.annot --pval /projects/abv/users/olivamx2/Projects/AD/reproduce_fanying/magma/AD_EUR_pvalue_MAGMA.txt ncol=3 --out /data/AD_EUR_NCBI37.3_mean

#output files
AD_EUR_NCBI37.3_mean

#############################################################################################
# step4. Prepare HCA skin portal scRNAseq seurat file
#############################################################################################
#download scRNAseqdata HCA portal - Reynolds et al (2021) Science "Developmental cell programs are co-opted in inflammatory skin disease"
wget https://zenodo.org/record/4569496/files/submission_210120.h5ad?download=1
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-8422/E-MTAB-8422.sdrf.txt

#code: HCA_skin_prepare_scRNAseq_seurat_file.R

#input:/data/hca_skin_portal/
submission_210120.h5ad

#output:/data/hca_skin_portal/
HCA_AD_lesion_Healthy_seurat_celltypelevel1.RDS


################################################
#step5. Change Seurat object to AnnData
################################################
#Code:
seurat_to_anndata.R

#input:/data/hca_skin_portal/
HCA_AD_lesion_Healthy_seurat_celltypelevel1.RDS

#Output: /data/hca_skin_portal/
HCA_AD_Lesion_level1.h5ad
HCA_Healthy_level1.h5ad
HCA_AD_Lesion_vs_Healthy_level1.h5ad
HCA_AD_Lesion_vs_Healthy_level2.h5ad

#############################################################################################
# step6. Run the wilcoxon rank test to get zscore, pvalue, log2FC
#calculate cell type programs and disease progression programs
#############################################################################################
#code:
/projects/abv/users/olivamx2/Projects/AD/reproduce_fanying/magma/src/AD_HCA_scanpy_CellType.ipynb
/projects/abv/users/olivamx2/Projects/AD/reproduce_fanying/magma/src/AD_HCA_scanpy_DiseaseProgression.ipynb

library(rmarkdown)
convert_ipynb(AD_HCA_scanpy_CellType.ipynb, AD_HCA_scanpy_CellType.Rmd)
convert_ipynb(AD_HCA_scanpy_DiseaseProgression.ipynb, AD_HCA_scanpy_DiseaseProgression.Rmd)

#converted ipynb to RMD file
/projects/abv/users/olivamx2/Projects/AD/reproduce_fanying/magma/src/AD_HCA_scanpy_CellType.Rmd
/projects/abv/users/olivamx2/Projects/AD/reproduce_fanying/magma/src/AD_HCA_scanpy_DiseaseProgression.Rmd

#input: 
HCA_AD_Lesion_level1.h5ad
HCA_Healthy_level1.h5ad
HCA_AD_Lesion_vs_Healthy_level1.h5ad
HCA_AD_Lesion_vs_Healthy_level2.h5ad

#Output:/data/AD_HCA_scanpy_output
42_clusters_AD_Lesion_CT
42_clusters_AD_Lesion_CT 
42_clusters_HC_CT   
42_clusters_HC_CT
42_clusters_AD_Lesion_vs_HC_DP
42_clusters_AD_Lesion_vs_HC_DP

skin_logfold.csv
skin_pval.csv
skin_score.csv

#############################################################################################
# step7. Prepare scRNAseq MAGMA signatures
#############################################################################################
#code: 
CellType_AD_HCA_continuous_input_v2.R
DiseaseProgression_AD_HCA_continuous_input_v2.R

#input:/data/AD_HCA_scanpy_output/*/
skin_score.csv

#Output: /data/markers_continuous_values/
AD_HCA_skin_H_42_clusters_CellType_signature_MAGMA.txt
AD_HCA_skin_H_42_clusters_CellType_signature_MAGMA.txt 
AD_HCA_skin_LS_42_clusters_CellType_signature_MAGMA.txt 
AD_HCA_skin_LS_42_clusters_CellType_signature_MAGMA.txt
AD_HCA_skin_LS_vs_H_42_clusters_DiseaseProgression_signature_MAGMA.txt
AD_HCA_skin_LS_vs_H_42_clusters_DiseaseProgression_signature_MAGMA.txt

#############################################################################################
# step8. Run MAGMA on cluster
#############################################################################################
#run using mean MAGMA model for GWAS results

input_name <- "H_42_clusters_CellType"
input_name <- "H_42_clusters_CellType"
input_name <- "LS_42_clusters_CellType"
input_name <- "LS_42_clusters_CellType"
input_name <- "LS_vs_H_42_clusters_DiseaseProgression"
input_name <- "LS_vs_H_42_clusters_DiseaseProgression"

setwd("/data")
## Genebody:
sample_name <- input_name
#mean
system(paste0("./magma --gene-results /data/AD_EUR_NCBI37.3_mean.genes.raw ",
              "--gene-covar /data/markers_continuous_values/AD_HCA_skin_", sample_name,"_signature_MAGMA.txt",
              " --model direction=positive --out ",
              "/data/continuous_output/HCA_",sample_name,"_mean"))


#############################################################################################
# step8. Make barplots
#############################################################################################

#code: MAGMA_HCA_continuous_scanpy_barplot_v1.R 

#output:
/projects/abv/users/olivamx2/Projects/AD/reproduce_fanying/magma/results/continuous_output/*.pdf
/projects/abv/users/olivamx2/Projects/AD/reproduce_fanying/magma/results/continuous_output/*.pdf


#############################
#step9. Additional barplot code for publication
##############################

#Plot cell-type enrichment results
####14 and 42 clusters

######################################
#14 clusters
#combine 14 cluster plots of LS, H and DP
library(data.table)
library(tidyverse)
setwd("/projects/abv/users/olivamx2/Projects/AD/reproduce_fanying/magma/results/continuous_output")
LS_14 <- read.table("HCA_LS_14_clusters_CellType_mean.gsa.out", stringsAsFactors = F, header = T)
H_14 <- read.table("HCA_H_14_clusters_CellType_mean.gsa.out", stringsAsFactors = F, header = T)
DP_14 <- read.table("HCA_LS_vs_H_14_clusters_DiseaseProgression_mean.gsa.out", stringsAsFactors = F, header = T)

#for DP strip the "Disease_" from cell type
name_temp <- DP_14$VARIABLE
name_temp <- as.vector(sapply(name_temp, function(x) unlist(strsplit(x, "Disease_", fixed = T))[2]))
DP_14$VARIABLE <- name_temp

#log transform P
LS_14$log_p <- -log10(LS_14$P)
H_14$log_p <- -log10(H_14$P)
DP_14$log_p <- -log10(DP_14$P)

#create FDR col
pvalue_fdr <- p.adjust(LS_14$P , method = "fdr", n = length(LS_14$P))
pvalue_fdr_sorted <- pvalue_fdr[order(pvalue_fdr, decreasing = F)]
LS_14$fdr <- pvalue_fdr
head(LS_14)

pvalue_fdr <- p.adjust(H_14$P , method = "fdr",  n = length(H_14$P))
pvalue_fdr_sorted <- pvalue_fdr[order(pvalue_fdr, decreasing = F)]
H_14$fdr <- pvalue_fdr
head(H_14)

pvalue_fdr <- p.adjust(DP_14$P , method = "fdr", n = length(DP_14$P))
pvalue_fdr_sorted <- pvalue_fdr[order(pvalue_fdr, decreasing = F)]
DP_14$fdr <- pvalue_fdr
head(DP_14)


#combine LS, H and DP in single file
library(dplyr)
all_14 <- left_join(H_14, LS_14, by = "VARIABLE")  
head(all_14)
all_14 <- left_join(all_14, DP_14, by = "VARIABLE")  
head(all_14)

colnames(all_14)[9] <- "H_FDR"
colnames(all_14)[17] <- "LS_FDR"
colnames(all_14)[25] <- "DP_FDR"
colnames(all_14)[8] <- "H_log_p"
colnames(all_14)[16] <- "LS_log_p"
colnames(all_14)[24] <- "DP_log_p"

#drop extra columns
all_14_filter <- all_14 %>% select('VARIABLE', 'H_log_p', 'H_FDR', 'LS_log_p', 'LS_FDR', 'DP_log_p', 'DP_FDR')
head(all_14_filter)
all_14_filter$VARIABLE <- as.character(all_14_filter$VARIABLE)
all_14_filter$sample_name <- all_14_filter$VARIABLE

#create just log_p file
all_14_filter2 <- all_14_filter %>% select('VARIABLE', 'H_log_p', 'LS_log_p', 'DP_log_p')
write.table(all_14_filter2, "14_clusters_logp_plots.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
head(all_14_filter2)

library(ggplot2)

dfm1 <- pivot_longer(all_14_filter2, -VARIABLE, names_to="variable", values_to="value")
dfm1$variable <- factor(dfm1$variable, levels = c('DP_log_p','LS_log_p','H_log_p'))
levels(dfm1$variable)

bp <- ggplot(dfm1,aes(x = VARIABLE,y = value)) +
  geom_bar(aes(fill = variable),stat = "identity",position = "dodge", width = 0.7) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=0), legend.title=element_blank()) + 
  labs(title= "AD GWAS Cell-Type Enrichment - 14 skin cell clusters", x="Cell_Type", y="-log10(pvalue)") +
  coord_flip() + 
  geom_hline(yintercept = (-log10(0.05)), linetype="dashed", color="grey", size=0.3) +
  geom_hline(yintercept=-(log10(0.02368141)), linetype="dashed", color = "red",size=0.3)

bp
bp + scale_fill_grey(breaks = rev(levels(dfm1$variable)), labels=c("Healthy", "Lesional", "Disease Progression"))
ggsave(filename = "HCA_cell_type_enrichment_14_clusters.tiff", device = "tiff", dpi=300, width = 180, height = 200, units = "mm")

######################################
#42 clusters
#combine 42 cluster plots of LS, H and DP
library(data.table)
library(tidyverse)
setwd("/projects/abv/users/olivamx2/Projects/AD/reproduce_fanying/magma/results/continuous_output")
LS_42 <- read.table("HCA_LS_42_clusters_CellType_mean.gsa.out", stringsAsFactors = F, header = T)
H_42 <- read.table("HCA_H_42_clusters_CellType_mean.gsa.out", stringsAsFactors = F, header = T)
DP_42 <- read.table("HCA_LS_vs_H_42_clusters_DiseaseProgression_mean.gsa.out", stringsAsFactors = F, header = T)

#log trasnform P
LS_42$log_p <- -log10(LS_42$P)
H_42$log_p <- -log10(H_42$P)
DP_42$log_p <- -log10(DP_42$P)

#create FDR col
pvalue_fdr <- p.adjust(LS_42$P , method = "fdr", n = length(LS_42$P))
pvalue_fdr_sorted <- pvalue_fdr[order(pvalue_fdr, decreasing = F)]
LS_42$fdr <- pvalue_fdr
head(LS_42)

pvalue_fdr <- p.adjust(H_42$P , method = "fdr", n = length(H_42$P))
pvalue_fdr_sorted <- pvalue_fdr[order(pvalue_fdr, decreasing = F)]
H_42$fdr <- pvalue_fdr
head(H_42)

pvalue_fdr <- p.adjust(DP_42$P , method = "fdr", n = length(DP_42$P))
pvalue_fdr_sorted <- pvalue_fdr[order(pvalue_fdr, decreasing = F)]
DP_42$fdr <- pvalue_fdr
head(DP_42)

#combine LS, H and DP in single file
library(dplyr)
all_42 <- left_join(H_42, LS_42, by = "VARIABLE")  
all_42 <- left_join(all_42, DP_42, by = "VARIABLE")  
head(all_42)

colnames(all_42)[9] <- "H_FDR"
colnames(all_42)[17] <- "LS_FDR"
colnames(all_42)[26] <- "DP_FDR"
colnames(all_42)[8] <- "H_log_p"
colnames(all_42)[16] <- "LS_log_p"
colnames(all_42)[25] <- "DP_log_p"
head(all_42)

#drop extra columns
all_42_filter <- all_42 %>% select('VARIABLE', 'H_log_p', 'H_FDR', 'LS_log_p', 'LS_FDR', 'DP_log_p', 'DP_FDR')
head(all_42_filter)
all_42_filter$VARIABLE <- as.character(all_42_filter$VARIABLE)
all_42_filter$sample_name <- all_42_filter$VARIABLE

#create just log_p file
all_42_filter2 <- all_42_filter %>% select('VARIABLE', 'H_log_p', 'LS_log_p', 'DP_log_p')
write.table(all_42_filter2, "42_clusters_logp_plots.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
head(all_42_filter2)

library(ggplot2)
dfm1 <- pivot_longer(all_42_filter2, -VARIABLE, names_to="variable", values_to="value")
write.table(dfm1, "42_clusters_logp_plot_to_edit.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)

#manually edits cell types to match Reynolds 2021 cluster names
#https://www.science.org/doi/10.1126/science.aba6500?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed#sec-3

dfm2 <- fread("42_clusters_logp_plot_edited.txt") 
head(dfm2)

dfm2$variable <- factor(dfm2$variable, levels = c('DP_log_p','LS_log_p','H_log_p'))
levels(dfm2$variable)

bp <- ggplot(dfm2,aes(x = VARIABLE,y = value)) +
  geom_bar(aes(fill = variable),stat = "identity",position = "dodge", width = 0.7) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=0), legend.title=element_blank(), plot.title = element_text(hjust=0.5)) + 
  labs(title= "AD GWAS Cell-Type Enrichment - 42 skin cell clusters", x="Cell_Type", y="-log10(pvalue)") +
  coord_flip() + 
  geom_hline(yintercept = (-log10(0.05)), linetype="dashed", color="grey", size=0.3) +
  geom_hline(yintercept=-(log10(0.02368421)), linetype="dashed", color = "red",size=0.3) +
  theme(plot.margin = margin(t = 10, b = 10, l = 10, r = 10))

bp

bp + scale_fill_grey(breaks = rev(levels(dfm2$variable)), labels=c("Healthy", "Lesional", "Disease Progression"))
ggsave(filename = "HCA_cell_type_enrichment_42_clusters.tiff", device = "tiff", dpi=300, width = 180, height = 200, units = "mm")


######################################
#42 clusters - filtered for significance in 14 cluster plot
#Plot cell-type enrichment filtered 42 clusters - significant 14 clusters
dfm3 <- fread("42_clusters_logp_plot_filtered_edited.txt") 
head(dfm3)

dfm3$variable <- factor(dfm3$variable, levels = c('DP_log_p','LS_log_p','H_log_p'))

levels(dfm3$variable)

bp <- ggplot(dfm3,aes(x = VARIABLE,y = value)) +
  geom_bar(aes(fill = variable),stat = "identity",position = "dodge", width = 0.7) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=0), legend.title=element_blank()) + 
  labs(title= "AD GWAS Cell-Type Enrichment - 42 skin cell clusters filtered", x="Cell_Type", y="-log10(pvalue)") +
  coord_flip() + 
  geom_hline(yintercept = (-log10(0.05)), linetype="dashed", color="grey", size=0.2) +
  geom_hline(yintercept=-(log10(0.02368421)), linetype="dashed", color = "red",size=0.2)

bp
bp + scale_fill_grey(breaks = rev(levels(dfm3$variable)), labels=c("Healthy", "Lesional", "Disease Progression"))
ggsave(filename = "HCA_cell_type_enrichment_42_clusters_filtered.tiff", device = "tiff", dpi=300, width = 180, height = 200, units = "mm")


#############################
#END
##############################