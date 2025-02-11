# Data obtained from UMichigan from correspondance with mksarkar <mksarkar@med.umich.edu>
library(ComplexHeatmap)
library(reshape2)
library(ggplot2)

data1 <- read.table('data/1st_siRNA_batch_of_P43-14_For_AbbVie.KO_eff.tsv',header=T)
data2 <- read.table('data/2nd_siRNA_batch_of_P43-23_For_AbbVie.KO_eff.tsv',header=T)
data <- rbind(data1,data2)
data$Batch <- paste0('Batch',data$Batch)

data <- subset(data,Include%in%'Y')# exclude SPINK7, IL22RA
data$siRNA <- factor(data$siRNA)
data$Treatment <- factor(data$Treatment,levels=c('NoTx','IL-22','IL-13','IL-13+IL-22'))

pdf('fig.keratinocyte_siRNA.KO_eff.pdf',width = 11)
ggplot(data, aes(x=Treatment, y=Expression)) + 
  geom_boxplot(aes(fill=siRNA)) + 
  facet_wrap(CandidateGene~Batch,scales = 'free')+
  theme_classic()
dev.off()


