#figure cell type enrichment 14 cell-type clusters

library(data.table)
library(ggplot2)
library(tidyverse)
all_14_filter <- fread("data_/14_clusters_logp_plots.txt") 
head(all_14_filter)

dfm1 <- pivot_longer(all_14_filter, -VARIABLE, names_to="variable", values_to="value")

dfm1$variable <- factor(dfm1$variable, levels = c('DP_log_p','LS_log_p','H_log_p'))

levels(dfm1$variable)

ggplot(dfm1,aes(x = VARIABLE,y = value)) +
  geom_bar(aes(fill = variable),stat = "identity",position = "dodge", width = 0.7) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=0), legend.title=element_blank(),text = element_text(size=8),axis.text=element_text(size=8)) +
  scale_y_continuous(breaks=seq(0,12,by=2)) +
  labs(title= "AD GWAS Cell-Type Enrichment - 14 skin cell clusters", x="Cell_Type", y="-log10(pvalue)") +
  coord_flip() +
  geom_hline(yintercept = (-log10(0.05)), linetype="dashed", color="grey", size=0.2) +
  geom_hline(yintercept=-(log10(0.02368421)), color = "red",size=0.2) +
  scale_fill_grey(breaks = rev(levels(dfm1$variable)), labels=c("Healthy", "Lesional", "Disease Progression"))

ggsave("fig.MAGMA_cell_type_enrichment_14_clusters.pdf", width=5, height=5, limitsize = FALSE)         