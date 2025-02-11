#fig LDSC enrichment - Immune Cells barplot, remove progenitors

library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggpubr)
library(dplyr)

# step1. Prepare input data
sample_dir <- "data_/"
setwd(sample_dir)
input_name <- "Immune_Cells_Enrichment_LDSC"

meta_info <- read.csv("GSE118189_GSE74912_sample_list_modified.csv", header = T)
unique(meta_info$cell_group)
table(meta_info$cell_group)
meta_info$cell_group[35] <- "B"

meta_info <- meta_info[order(meta_info$status), ]
meta_info <- meta_info[order(meta_info$Name), ]
meta_info <- meta_info[order(meta_info$cell_group), ]
meta_info$cell_name <- as.vector(sapply(meta_info$Name, function(x) unlist(strsplit(x, "-", fixed = T))[1]))
dim(meta_info)
#[1] 58  5
head(meta_info)

#drop progenitor
meta_info <- meta_info %>% filter(!cell_group == 'Progenitor')

meta_info_sub <- meta_info[-which(meta_info$source=="GSE74912_Corces" & meta_info$cell_group!="Progenitor"),]
dim(meta_info_sub)
[1] 45  5

#meta_info <- meta_info[-c(1, 8), ]
input_ldsc <- read.table("AD_ATAC_EURmeta.cell_type_results.txt", stringsAsFactors = F, header = T)

sample_input_func <- function(input_dir) {
  input <- read.table(input_dir, stringsAsFactors = F, header = T)
  input$sample <- as.vector(sapply(input$Name, function(x) unlist(strsplit(x, "__", fixed = T))[1]))
  input_reorder <- input[match(meta_info_sub$Name, input$sample),]
  input_reorder$cell_group <- meta_info_sub$cell_group
  input_reorder$status <- meta_info_sub$status
  input_reorder$logPvalue <- -log10(input_reorder$Coefficient_P_value)
  input_reorder$FDR <- p.adjust(input_reorder$Coefficient_P_value, method = "fdr", n = length(input_reorder$Coefficient_P_value))
  return(input_reorder)
}

input_AD <- sample_input_func("AD_ATAC_EURmeta.cell_type_results.txt")

pvalue_mat <- cbind(input_AD$logPvalue)
rownames(pvalue_mat) <-input_AD$sample
colnames(pvalue_mat) <- c("AD")
head(pvalue_mat)

FDR_mat <- cbind(input_AD$FDR)
rownames(FDR_mat) <-input_AD$sample
colnames(FDR_mat) <- c("AD")
head(FDR_mat)

####################################################################################################
# Make barplots
####################################################################################################

input_disease <- "AD"

input_mat <- as.matrix(pvalue_mat)
head(meta_info_sub)
head(input_mat)

tail(meta_info_sub)
tail(input_mat)

dim(meta_info_sub)
dim(input_mat)

input_mat_df <- reshape2::melt(input_mat, id.vars='-log10pvalue', variable.name="AD")
meta_info_sub_reorder <- meta_info_sub[match(input_mat_df$Var1, meta_info_sub$Name),]
input_mat_df$cell_group <- meta_info_sub_reorder$cell_group
input_mat_df$status  <- meta_info_sub_reorder$status
input_mat_df$cell_name <- meta_info_sub_reorder$cell_name

input_mat_df$Var2 <- as.character(input_mat_df$Var2)
unique(input_mat_df$Var2)
# AD

head(input_mat_df)
input_mat_df$general_type <- "PBMC"
#input_mat_df$general_type[input_mat_df$cell_group=="Progenitor"] <- "Progenitors"
input_mat_df$general_type[input_mat_df$status=="U" & input_mat_df$general_type=="PBMC"] <- "PBMC_U"
input_mat_df$general_type[input_mat_df$status=="S" & input_mat_df$general_type=="PBMC"] <- "PBMC_S"


colnames(input_mat_df)[c(1,2)] <- c("sample_name", "disease")
input_mat_df_sub <- input_mat_df[which(input_mat_df$disease==input_disease),]
rownames(input_mat_df_sub) <- c(1:dim(input_mat_df_sub)[1])
input_mat_df_sub$cell_name <- factor(input_mat_df_sub$cell_name, levels=rev(unique(input_mat_df_sub$cell_name)))

# calculate the pvalue threhold
raw_pvalue_mat <- cbind(input_AD$Coefficient_P_value)
colnames(raw_pvalue_mat) <- c("AD")
rownames(raw_pvalue_mat) <- rownames(input_AD)

pvalue_mat_select <- as.data.frame(raw_pvalue_mat[, input_disease])
colnames(pvalue_mat_select) <- "p_value"
pvalue_mat_fdr <- p.adjust(pvalue_mat_select$p_value , method = "fdr", n = length(pvalue_mat_select$p_value))
pvalue_mat_fdr_sorted <- pvalue_mat_fdr[order(pvalue_mat_fdr, decreasing = F)]
rank_thre <- max(rank(pvalue_mat_fdr_sorted[which(pvalue_mat_fdr_sorted < 0.05)]))
pvalue_thre <- rank_thre/length(pvalue_mat_fdr_sorted) * 0.05

############################################
# make the barplot:
pdf(paste0(input_name,"_",input_disease,"_barplot.pdf"), height = 5, width = 4)
p <- ggplot(input_mat_df_sub, aes(fill=status, y=value, x=cell_name)) +
  geom_bar(position = position_dodge(preserve = "single"), stat="identity", width = 0.6)+
  theme_bw()+
  scale_fill_grey()+
  scale_y_continuous(breaks=seq(0,9,by=2)) +
  theme(axis.text.x = element_text(angle=0),text = element_text(size=8))+
  labs(title="AD GWAS Cell-Type Enrichment - Immune Cell ATACseq",
       x="Cell_Type", y="-log10(pvalue)") +
  theme(plot.title = element_text(hjust=0.5))
p+coord_flip() +
  geom_hline(yintercept=-(log10(pvalue_thre)), linetype="dashed", color = "red",linewidth=0.3)
dev.off()

ggplot(input_mat_df_sub, aes(fill=status, y=value, x=cell_name)) +
  geom_bar(position = position_dodge(preserve = "single"), stat="identity", width = 0.6)+
  theme_bw()+
  scale_fill_grey()+
  scale_y_continuous(breaks=seq(0,12,by=2)) +
  theme(axis.text.x = element_text(angle=0),text = element_text(size=8),axis.text=element_text(size=8))+
  labs(title="AD GWAS Cell-Type Enrichment - Immune Cell ATACseq",
       x="Cell_Type", y="-log10(pvalue)") +
  theme(plot.title = element_text(hjust=0.5)) +
  coord_flip() +
  geom_hline(yintercept=-(log10(pvalue_thre)), linetype="dashed", color = "red",linewidth=0.3)

ggsave("fig.LDSC_cell_type_enrichment_barplot.pdf", width=5, height=5, limitsize = FALSE) 
