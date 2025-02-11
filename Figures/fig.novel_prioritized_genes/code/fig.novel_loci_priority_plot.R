#figure 5 novel loci priority plot
#create HC priority gene table
#tables: V2G, nearest gene, DEPICT, MAGMA, Coloc all, Coloc Trans signal, coding var: CS+LD+lit, OMIM, Tsoi GSE121212 Healthy vs ADL (FDR < 0.05, logfc > [1])

library(data.table)
library(tidyverse)
library(dplyr)

setwd("data_/")

#STEP1: create HC priority gene table
#tables: V2G, nearest gene, DEPICT, MAGMA, Coloc all, Coloc Trans signal, coding var: CS+LD+lit, OMIM, Tsoi GSE121212 Healthy vs ADL (FDR < 0.05, logfc > [1])
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

head(combined_long)
head(combined_annotate)

#select novel loci only
head(combined_long)
novel_loci <- c("5","21","22","34","38","44","45","48","51","62","68","75","83","93","94")
combined_long_novel <- combined_long %>% filter(GWAS_locus %in% novel_loci)
head(combined_long_novel)

#drop novel rows
combined_long_novel <- combined_long_novel %>% filter(!str_detect(name, 'Novel_loci'))

evidence_order_novel <- c('Nearest_gene','V2G','DEPICT','MAGMA','cis_espQTL','trans_pQTL','Coding_CS', 'Coding_LD', 'Coding_Lit','OMIM','DEG_AD')
level=evidence_order_novel

#drop genes without any evidence
combined_long_novel <- combined_long_novel %>% filter(!str_detect(Gene, 'FLVCR2'))

#group evidence by type
combined_long_novel$evidence_group <- ifelse(combined_long_novel$name %in% c("Nearest_gene", "V2G","DEPICT","MAGMA"), "Variant-to-Gene",
                                             ifelse(combined_long_novel$name %in% c("cis_espQTL","trans_pQTL"), "QTL",
                                                    ifelse(combined_long_novel$name %in% c("Coding_CS","Coding_LD","Coding_Lit"), "Coding Variant",
                                                           ifelse(combined_long_novel$name %in% c("OMIM","DEG_AD"), "AD phenotype",""))))

head(combined_long_novel)

plot1 <- ggplot(combined_long_novel, aes(x=factor(name,level=evidence_order_novel), y=weave_factors(Gene,Lead_SNP,Cytoband), color = evidence_group, alpha = value)) +
  geom_tile(aes(width=1, height=1),color= "black", fill= "white", size = 0.1) +
  geom_point(aes(fill=evidence_group),shape = 21,size = 2.5,color="grey10",stroke=0.4) +
  scale_alpha_continuous(range = c(-0.5, 1)) +
  scale_fill_viridis_d(direction = -1) +
  theme_bw() +
  theme(legend.position = "none",axis.title.x =element_blank(),axis.title.y =element_blank(),axis.ticks=element_line(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.background=element_blank(), panel.border=element_blank(),plot.title = element_blank(),axis.text.x=element_text(angle=90, size = 8, hjust=0.95,vjust=0.2),axis.text.y=element_text(size = 8),text = element_text(size=8),axis.text=element_text(size=8)) + 
  scale_size(range=c(1,2)) +
  scale_y_discrete(guide = 'axis_nested',limits=rev)

gene_order <- c('CD52','CEP85','SH3BGRL3','UBXN11','ZNF593','GLS','NAB1','PPIL3','STRADB','TRAK2','ANKRD55','IL6ST','ITK','IL22RA2','SCAF8','TIAM2','IL6','STEAP1B','CHD7','ANK3','ANK3-DT','CDK1','C11orf80','CTSF','GSTP1','NDUFV1','RAD52','WNK1','BATF','PGS1','SOCS3','ALPK2')

dens5 <- combined_long_novel %>% group_by(Gene) %>% summarize(value = sum(value)) %>% 
  ggplot(aes(x =factor(Gene,level=gene_order),y = value)) + 
  scale_x_discrete(limits=rev) +
  geom_bar(stat = "identity", position = "dodge",fill="grey80",color="black") +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y=element_blank(),axis.title.x=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.line.x = element_line(size = .2),axis.text.x=element_text(size = 8))

library(patchwork)

plot1 + dens5 + plot_spacer()
ggsave("fig.novel_loci_priority_plot.pdf", width=6, height=5, limitsize = FALSE)

