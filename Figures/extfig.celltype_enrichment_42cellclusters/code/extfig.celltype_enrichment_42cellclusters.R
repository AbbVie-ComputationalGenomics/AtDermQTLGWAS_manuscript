#extended figure cell type enrichment 42 cell-type clusters

library(data.table)
library(tidyr)
library(ggplot2)
all_42_filter <- fread("data_/42_clusters_logp_plot.txt") 
head(all_42_filter)

all_42_filter$variable <- factor(all_42_filter$variable, levels = c('DP_log_p','LS_log_p','H_log_p'))

levels(all_42_filter$variable)

ggplot(all_42_filter,aes(x = VARIABLE,y = value)) +
  geom_bar(aes(fill = variable),stat = "identity",position = "dodge", width = 0.7) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=0), legend.title=element_blank(),text = element_text(size=8),axis.text=element_text(size=8)) +
  scale_y_continuous(breaks=seq(0,12,by=2)) +
  labs(title= "AD GWAS Cell-Type Enrichment - 42 skin cell clusters", x="Cell_Type", y="-log10(pvalue)") +
  coord_flip() +
  geom_hline(yintercept = (-log10(0.05)), linetype="dashed", color="grey", size=0.2) +
  geom_hline(yintercept=-(log10(0.02368421)), color = "red",size=0.2) +
  scale_fill_grey(breaks = rev(levels(all_42_filter$variable)), labels=c("Healthy", "Lesional", "Disease Progression"))

ggsave("extfig.celltype_enrichment_42_clusters.pdf", width=5, height=5, limitsize = FALSE)                  

