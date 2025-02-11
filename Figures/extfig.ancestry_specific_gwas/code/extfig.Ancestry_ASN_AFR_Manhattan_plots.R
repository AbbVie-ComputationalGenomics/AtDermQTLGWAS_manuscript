#extfig Ancestry-specific Manhattan plots

library(dplyr)
library(ggplot2)
library(data.table)
requireNamespace("RColorBrewer", quietly = TRUE)
library(ggrepel)
library(gridExtra)

##########
#ASN meta-analysis Manhattan plot
#filter and format summary stats

ASN <- fread("data_/AD_GWAS_ASN_meta_FE_hg38.txt.gz")
head(ASN)
dim(ASN)
[1] 27070036       16
ASN <- filter(ASN, PVAL < 0.3)
head(ASN)
dim(ASN)
[1] 8033254      16

colnames(ASN)[1] <- "chr"
colnames(ASN)[2] <- "pos"
colnames(ASN)[6] <- "rsid"
colnames(ASN)[7] <- "pval"
colnames(ASN)[8] <- "beta"
colnames(ASN)[9] <- "se"
head(ASN)

#select columns
ASN <- ASN %>% select(rsid,chr,pos,beta,se,pval)
head(ASN)

#remove rows with no rsID
ASN <- ASN %>% filter(rsid != '.')
head(ASN)
dim(ASN)
[1] 7703910       6

#assign relative position to each SNP along genome
ASN_prep <- ASN %>% group_by(chr) %>% summarise(chrlength = max(pos)) 
ASN_prep <- ASN_prep %>% mutate(cumulativechrlength = cumsum(as.numeric(chrlength)) - chrlength) 
ASN <- ASN_prep %>% select(-chrlength) %>% left_join(ASN, ASN_prep, by = "chr")
ASN <- ASN %>%  arrange(chr, pos) %>% mutate(rel_pos = pos + cumulativechrlength) %>% select(-cumulativechrlength)
head(ASN)
tail(ASN)

#provide the chromosome position labels relative to the entire genome.
axis_data <- ASN %>%
  dplyr::mutate(chr = chr) %>%
  dplyr::group_by(chr) %>%
  dplyr::summarize(chr_center = (max(rel_pos) + min(rel_pos)) / 2)

# To make it easier for later, create function-named chromosome and logged p-value columns
ASN <- ASN %>%
  dplyr::mutate(logged_p = -log10(ASN$pval))

#dplyr::rename(chr = as.name(chr))
maxp <- ceiling(max(ASN$logged_p, na.rm = TRUE))
#36

#plot ASN manhattan
ggplot() +
  geom_point(data = ASN, 
             aes(x = rel_pos, y = logged_p, color = as.factor(chr)), 
             size = 0.25) +
  scale_color_manual(values = rep(c("black", "grey"), nrow(axis_data))) +
  scale_x_continuous(labels = axis_data$chr, 
                     breaks = axis_data$chr_center, 
                     expand = expansion(mult = 0.01),
                     guide = guide_axis(check.overlap =TRUE)) +
  scale_y_continuous(limits = c(0, maxp), 
                     expand = expansion(mult = c(0.02, 0))) +
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed", 
             size = 0.3) +
  labs(x = "Chromosome", y = bquote(atop('-log'[10]*'(p)'))) + 
  theme_classic() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.margin = margin(t = 10, b = 0, l = 10, r = 10))

ggsave(filename = "ASN_GWAS_meta_Manhattan_plot.png", device = "png", dpi=300, width = 180, height = 115, units = "mm")


########################
#######################
#Make AFR meta-analysis manhattan plot

AFR <- fread("data_/AD_GWAS_AFR_meta_FE_hg38.txt.gz")
head(AFR)
dim(AFR)
[1] 35441158       16
AFR <- filter(AFR, PVAL < 0.3)
head(AFR)
dim(AFR)
[1] 10235466       16
colnames(AFR)[1] <- "chr"
colnames(AFR)[2] <- "pos"
colnames(AFR)[6] <- "rsid"
colnames(AFR)[7] <- "pval"
colnames(AFR)[8] <- "beta"
colnames(AFR)[9] <- "se"
head(AFR)

#select columns
AFR <- AFR %>% select(rsid,chr,pos,beta,se,pval)
head(AFR)

#remove rows with no rsID
AFR <- AFR %>% filter(rsid != '.')
head(AFR)
dim(AFR)
[1] 10153289        6


#assign relative position to each SNP along genome
AFR_prep <- AFR %>% group_by(chr) %>% summarise(chrlength = max(pos)) 
AFR_prep <- AFR_prep %>% mutate(cumulativechrlength = cumsum(as.numeric(chrlength)) - chrlength) 
AFR <- AFR_prep %>% select(-chrlength) %>% left_join(AFR, AFR_prep, by = "chr")
AFR <- AFR %>%  arrange(chr, pos) %>% mutate(rel_pos = pos + cumulativechrlength) %>% select(-cumulativechrlength)
head(AFR)
tail(AFR)

#provide the chromosome position labels relative to the entire genome.
axis_data <- AFR %>%
  dplyr::mutate(chr = chr) %>%
  dplyr::group_by(chr) %>%
  dplyr::summarize(chr_center = (max(rel_pos) + min(rel_pos)) / 2)

# To make it easier for later, create function-named chromosome and logged p-value columns
AFR <- AFR %>%
  dplyr::mutate(logged_p = -log10(AFR$pval))

#dplyr::rename(chr = as.name(chr))
maxp <- ceiling(max(AFR$logged_p, na.rm = TRUE))
#9
head(AFR)

#plot AFR manhattan
ggplot() +
  geom_point(data = AFR, 
             aes(x = rel_pos, y = logged_p, color = as.factor(chr)), 
             size = 0.25) +
  scale_color_manual(values = rep(c("black", "grey"), nrow(axis_data))) +
  scale_x_continuous(labels = axis_data$chr, 
                     breaks = axis_data$chr_center, 
                     expand = expansion(mult = 0.01),
                     guide = guide_axis(check.overlap =TRUE)) +
  scale_y_continuous(limits = c(0, maxp), 
                     expand = expansion(mult = c(0.02, 0))) +
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed", 
             size = 0.3) +
  labs(x = "Chromosome", y = bquote(atop('-log'[10]*'(p)'))) + 
  theme_classic() +
  theme(legend.position = "none",
        plot.margin = margin(t = 10, b = 0, l = 10, r = 10))

ggsave(filename = "AFR_GWAS_meta_Manhattan_plot.png", device = "png", dpi=300, width = 180, height = 115, units = "mm")