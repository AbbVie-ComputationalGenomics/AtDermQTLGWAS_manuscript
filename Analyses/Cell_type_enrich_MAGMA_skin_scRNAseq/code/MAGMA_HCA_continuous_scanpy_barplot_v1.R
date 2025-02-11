library(ggplot2)

sample_dir <- "/projects/abv/users/olivamx2/Projects/AD/reproduce_fanying/magma/results/continuous_output/"
setwd(sample_dir)

# barplot code for continuous results: Cell Type Program; Disease Progression Program

##############################################################
input_dir <- "/projects/abv/users/olivamx2/Projects/AD/reproduce_fanying/magma/results/continuous_output/"
input <- read.table("HCA_H_14_clusters_CellType_mean.gsa.out", stringsAsFactors = F, header = T)
head(input)
output_name <- "HCA_H_14_clusters_CellType"


input$sample_name <- input$VARIABLE
input$log_p <- -log10(input$P)

pvalue_fdr <- p.adjust(input$P , method = "fdr", 
                           n = length(input$P))
pvalue_fdr_sorted <- pvalue_fdr[order(pvalue_fdr, decreasing = F)]
input$fdr <- pvalue_fdr


if(min(pvalue_fdr) < 0.05) {
  rank_thre <- max(rank(pvalue_fdr[which(pvalue_fdr_sorted < 0.05)]))
  pvalue_thre <- rank_thre/length(pvalue_fdr_sorted) * 0.05
} else {
  print("No results pass FDR")
  pvalue_thre <- 0.05
}

pvalue_thre

input$VARIABLE <- as.character(input$VARIABLE)

if(length(pvalue_fdr[which(pvalue_fdr< 0.05)]) > 0) {
  pdf(paste0(output_name, ".pdf"), height = 5, width = 4)
  p <- ggplot(input, aes(y=log_p, x=sample_name)) +
    geom_bar(position = position_dodge(preserve = "single"), stat="identity", width = 0.6)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=0))+
    labs(title=output_name,
         x="Cell_Type", y="-log10(pvalue)")
  p1 <- p+coord_flip() +
    geom_hline(yintercept = (-log10(0.05)), linetype="dashed", color="grey", size=0.2) +
    geom_hline(yintercept=-(log10(pvalue_thre)), linetype="dashed", color = "red",size=0.2)
  print(p1)
  dev.off()
} else {
  pdf(paste0(output_name, ".pdf"), height = 5, width = 5)
  p <- ggplot(input, aes(y=log_p, x=sample_name)) +
    geom_bar(position = position_dodge(preserve = "single"), stat="identity", width = 0.6)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=0))+
     labs(title=output_name,
         x="Cell_Type", y="-log10(pvalue)")
  p1 <- p+coord_flip() +
    geom_hline(yintercept = (-log10(0.05)), linetype="dashed", color="grey", size=0.2)
  print(p1)
  dev.off()
}


print(length(pvalue_fdr[which(pvalue_fdr< 0.05)]))


#########
#Disease Progression
input_dir <- "/projects/abv/users/olivamx2/Projects/AD/reproduce_fanying/magma/results/continuous_output/"
input <- read.table("HCA_LS_vs_H_42_clusters_DiseaseProgression_mean.gsa.out", stringsAsFactors = F, header = T)
head(input)
output_name <- "LS_vs_H_42_DiseaseProgression"

name_temp <- input$VARIABLE
name_temp <- as.vector(sapply(name_temp, function(x) unlist(strsplit(x, "Disease_", fixed = T))[2]))
input$VARIABLE <- name_temp



input$sample_name <- input$VARIABLE
input$log_p <- -log10(input$P)

pvalue_fdr <- p.adjust(input$P , method = "fdr", 
                       n = length(input$P))
pvalue_fdr_sorted <- pvalue_fdr[order(pvalue_fdr, decreasing = F)]
input$fdr <- pvalue_fdr


if(min(pvalue_fdr) < 0.05) {
  rank_thre <- max(rank(pvalue_fdr[which(pvalue_fdr_sorted < 0.05)]))
  pvalue_thre <- rank_thre/length(pvalue_fdr_sorted) * 0.05
} else {
  print("No results pass FDR")
  pvalue_thre <- 0.05
}

input$VARIABLE <- as.character(input$VARIABLE)

if(length(pvalue_fdr[which(pvalue_fdr< 0.05)]) > 0) {
  pdf(paste0(output_name, ".pdf"), height = 5, width = 5)
  p <- ggplot(input, aes(y=log_p, x=sample_name)) +
    geom_bar(position = position_dodge(preserve = "single"), stat="identity", width = 0.6)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=0))+
    labs(title=output_name,
         x="Cell_Type", y="-log10(pvalue)")
  p1 <- p+coord_flip() +
    geom_hline(yintercept = (-log10(0.05)), linetype="dashed", color="grey", size=0.2) +
    geom_hline(yintercept=-(log10(pvalue_thre)), linetype="dashed", color = "red",size=0.2)
  print(p1)
  dev.off()
} else {
  pdf(paste0(output_name, ".pdf"), height = 5, width = 5)
  p <- ggplot(input, aes(y=log_p, x=sample_name)) +
    geom_bar(position = position_dodge(preserve = "single"), stat="identity", width = 0.6)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=0))+
    labs(title=output_name,
         x="Cell_Type", y="-log10(pvalue)")
  p1 <- p+coord_flip() +
    geom_hline(yintercept = (-log10(0.05)), linetype="dashed", color="grey", size=0.2)
  print(p1)
  dev.off()
}


print(length(pvalue_fdr[which(pvalue_fdr< 0.05)]))

