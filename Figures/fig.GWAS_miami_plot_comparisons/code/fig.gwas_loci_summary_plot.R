#fig2: GWAS loci summary plot (ancestry and known/novel atopic march)

library(data.table)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggh4x)
library(ggrepel)

#################

#0=known
#1=novel
#NA= not reported for allergy and asthma
known_novel <- fread("data/AD_GWAS_loci_plot_known_novel.txt")
head(known_novel)
known_novel_long <- known_novel %>% 
  pivot_longer(
    cols = `AD/eczema`:`Asthma`)
head(known_novel_long)

known_novel_long <- known_novel_long %>% na.omit(value)

known_novel_plot <- ggplot(known_novel_long) +
  geom_tile(aes(y=weave_factors(Locus_number,GWAS_Locus), x=name, fill=value, alpha = value),color = "white",) +
  scale_fill_gradient(low = "steelblue",
                      high = "steelblue") +
  theme(plot.margin = unit(c(0, 0, 0, 0),"pt"),
        legend.position = "none", 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(size = 9,angle = 90,hjust = 1,vjust=0.5),
        axis.text.y=element_text(size = 7,hjust = 0.2),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_line(colour = "steelblue"),
        ggh4x.axis.nestline.y = element_line(linewidth = 1.4),
        ggh4x.axis.nesttext.y = element_text(hjust = 0.2)) +
  scale_y_discrete(position = "right",guide =guide_axis_nested(n.dodge=1), breaks = interaction(known_novel_long$Locus_number,known_novel_long$GWAS_Locus)[c(11,48,49,75,85,100,101,107,114,135,148,164,181,210,212)]) +
  coord_equal(0.4)



###################
#ancestry stacked tile plot
ancestry <- fread("data/AD_GWAS_loci_plot_ancestry.txt")
head(ancestry)

ancestry_long <- ancestry %>% 
  pivot_longer(
    cols = `Multi-ancestry`:`AFR`)
head(ancestry_long)

anc_order <- c("Multi-ancestry","EUR","ASN","AFR")


ancestry_plot <- ggplot(ancestry_long, aes(y=weave_factors(Locus_number,Chromosome),x=factor(name,level=anc_order), fill=name, alpha = value)) +
  geom_tile() +
  theme(plot.margin = unit(c(0, 0, 0, 0),"pt"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(size = 9, angle = 90, hjust = 1,vjust=0.5),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size = 7.5),
        legend.position = "none",
        panel.border = element_blank())+
  ylab("GWAS Loci") +
  scale_y_discrete(guide = 'axis_nested', breaks = interaction(ancestry_long$Locus_number,ancestry_long$Chromosome)[c(1,48,49,96,97,120,121,124,125,156,157,184,185,200,201,220,221,232,233,252,253,296,297,320,321,324,325,336,337,344,345,352,353,372,373,376,377,388,389,400,401)]) +
  coord_equal(0.4) 



##############
#combine plots
library(patchwork)

#2panel vertical
ancestry_plot + known_novel_plot  + plot_annotation(title = "AD GWAS Loci")
ggsave(filename = "fig.GWAS_compare_ancestry_novel_plot.pdf", device = "pdf", width = 5, height = 8.5, units = "in")
          