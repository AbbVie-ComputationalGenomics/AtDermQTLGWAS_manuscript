#extended figure 4: prioritized genes at AD GWAS loci with >=2 lines of evidence
#tables: V2G, nearest gene, DEPICT, MAGMA, Coloc all, Coloc Trans signal, coding var: CS+LD+lit, OMIM, Tsoi GSE121212 Healthy vs ADL (FDR < 0.05, logfc > [1])

library(data.table)
library(tidyverse)
library(dplyr)
setwd("data_/")

#STEP1: create HC priority gene table
#combine all files 
# read file path
all_paths <- list.files(pattern = "*priority.txt",
                        full.names = TRUE)
# read file content
all_content <- all_paths %>% lapply(read.table,
                                    header = TRUE,
                                    sep = "\t",
                                    encoding = "UTF-8")
# read file name
all_filenames <- all_paths %>% basename() %>% as.list()
# combine file content list and file name list
all_lists <- mapply(c, all_content, all_filenames, SIMPLIFY = FALSE)
# unlist all lists and change column name
all_result <- rbindlist(all_lists, fill = T)
head(all_result)

library(reshape2)
combined <- dcast(all_result, GWAS_locus + Gene ~ V1)
head(combined)

library(tidyr)
combined <- combined %>% drop_na()

#rename columns
col_names = list('GWAS_locus','Gene','cis_espQTL','trans_pQTL','Coding_CS', 'Coding_LD', 'Coding_Lit','DEPICT','MAGMA','Nearest_gene','Novel_loci','OMIM','DEG_AD','V2G')
names(combined) <- c(col_names)
combined <- select(combined,'GWAS_locus','Gene','Novel_loci','Nearest_gene','V2G','DEPICT','MAGMA','cis_espQTL','trans_pQTL','Coding_CS', 'Coding_LD', 'Coding_Lit','OMIM','DEG_AD')
head(combined)

#annotate with cytoband and lead SNP ID
cyto <- fread("Cytoband_table_AD_GWAS.txt")
combined_annotate <- inner_join(cyto, combined, by="GWAS_locus")
head(combined_annotate)
write.table(combined_annotate, "combined_annotate.txt", sep = '\t', row.names = F, col.names = T, quote = F)
n_distinct(combined_annotate$Gene)
[1] 478
# genes = 478 unique


############################
#STEP2: plot HC genes
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggh4x)

head(combined_annotate)
combined_annotate <- combined_annotate[order(combined_annotate$GWAS_locus),]
head(combined_annotate)

combined_long <- combined_annotate %>% 
  pivot_longer(
    cols = `Novel_loci`:`DEG_AD`)
head(combined_long)

#any value >0 = 1
combined_long <- combined_long %>%
  mutate(across(.cols = c(6), .fns = function(x) ifelse(x > 0, 1, 0)))

#change cytoband "." to "_"
combined_long$Cytoband <- gsub('[.]', '_', combined_long$Cytoband)
head(combined_long)

#drop ENSG genes
combined_long <- combined_long %>% filter(!str_detect(Gene, 'ENSG'))
#drop RP1* genes
combined_long <- combined_long %>% filter(!str_detect(Gene, 'RP1'))
#drop LOC genes
combined_long <- combined_long %>% filter(!str_detect(Gene, 'LOC'))
#drop RNU genes
combined_long <- combined_long %>% filter(!str_detect(Gene, 'RNU'))
#drop MIR
combined_long <- combined_long %>% filter(!str_detect(Gene, 'MIR'))
#drop AC0* genes
combined_long <- combined_long %>% filter(!str_detect(Gene, 'AC0'))

n_distinct(combined_long$Gene)
[1] 441

#remove duplicate genes and "-", ETS1_ and SATB1_ in 2 different clumps
combined_long <- fread("data_/combined_long.txt")
combined_annotate <- fread("data_/combined_annotate.txt")

head(combined_long)
head(combined_annotate)

n_distinct(combined_long$Gene)
[1] 443
n_distinct(combined_annotate$Gene)
[1] 480
#note did not remove specific genes from this file

#set y-axis order
head(combined_long)
combined_long <- combined_long[order(combined_long$GWAS_locus),]

#set evidence order
evidence_order <- c('Novel_loci','Nearest_gene','V2G','DEPICT','MAGMA','cis_espQTL','trans_pQTL','Coding_CS', 'Coding_LD', 'Coding_Lit','OMIM','DEG_AD')
level=evidence_order

#filter for genes with summed evidence >1
combined_long_filter <- combined_long %>% group_by(Gene) %>% summarize(value = sum(value)) %>% filter(value >1)
head(combined_long_filter)
head(combined_long)
combined_long_filter_plot <- combined_long %>% filter(Gene %in% combined_long_filter$Gene)
head(combined_long_filter_plot)

combined_long_filter_plot$evidence_group <- ifelse(combined_long_filter_plot$name %in% c("Novel_loci"), "Novel",
                                                   ifelse(combined_long_filter_plot$name %in% c("Nearest_gene", "V2G","DEPICT","MAGMA"), "Variant-to-Gene",
                                                          ifelse(combined_long_filter_plot$name %in% c("cis_espQTL","trans_pQTL"), "QTL",
                                                                 ifelse(combined_long_filter_plot$name %in% c("Coding_CS","Coding_LD","Coding_Lit"), "Coding Variant",
                                                                        ifelse(combined_long_filter_plot$name %in% c("OMIM","DEG_AD"), "AD phenotype","")))))

head(combined_long_filter_plot)
n_distinct(combined_long_filter_plot$Gene)
[1] 239

#plot
library(RColorBrewer)
library(viridis)

plot2 <- ggplot(combined_long_filter_plot, aes(x=factor(name,level=evidence_order), y=weave_factors(Gene,Cytoband,GWAS_locus), color = evidence_group, alpha = value)) +
  geom_tile(aes(width=1, height=1),color= "black", fill= "white", size = 0.1) +
  geom_point(aes(fill=evidence_group),shape = 21, size = 2,color="grey50",stroke=0.1) +
  scale_fill_viridis_d(direction = -1) +
  scale_alpha_continuous(range = c(-0.5, 1)) +
  labs(y="Genes and Loci",x="Evidence") + 
  labs() + 
  theme_bw() +
  theme(legend.position = "none",axis.ticks=element_line(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.background=element_blank(), panel.border=element_blank(),plot.title = element_text(color="black", size=11,hjust = 0),axis.text.x=element_text(angle=90, size = 9, hjust=0.95,vjust=0.2),axis.text.y=element_text(size = 9)) + 
  scale_size(range=c(1,2)) +
  scale_y_discrete(guide = 'axis_nested',limits=rev)


#arrange histogram plot to same order as geom_tile
dens1 <- combined_long_filter_plot %>% group_by(Gene) %>% summarize(value = sum(value))
head(dens1)
head(combined_annotate)
combo_dens <- inner_join(combined_annotate, dens1, by="Gene")
head(combo_dens)
combo_dens <- combo_dens[order(combo_dens$GWAS_locus),]
head(combo_dens)
tail(combo_dens)
#lock gene order
combo_dens$Gene <- factor(combo_dens$Gene, levels = combo_dens$Gene)

dens2_plot <- ggplot(combo_dens, aes(x=Gene,y = value)) + 
  scale_x_discrete(limits=rev) +
  scale_y_discrete(limits = factor(1:8)) +
  geom_bar(stat = "identity", position = "dodge",fill="grey80",color="black") +
  coord_flip() +
  theme(axis.text.y=element_blank(),axis.title.x=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())

library(patchwork)

plot2 + dens2_plot + plot_spacer() + plot_annotation(title = "ALL Loci - Gene Prioritization (sum evidence >1)") 
ggsave("extfig_prioritized_genes_evidence>1_plot.pdf", width=25, height=25, limitsize = FALSE)
