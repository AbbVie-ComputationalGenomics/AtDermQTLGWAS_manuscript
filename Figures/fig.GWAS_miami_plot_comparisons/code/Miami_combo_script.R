#Run R script using slurm: sbatch Run_miami.slurm
library(miamiplot)
library(dplyr)
library(ggplot2)
library(data.table)
requireNamespace("RColorBrewer", quietly = TRUE)

setwd("/projects/abv/GRC/rileybm/projects/AD_metagwas_paper/Figures/fig.GWAS_miami_plot_comparisons/data")

#Get Miami_combo results
miami_combo <- fread("miami_combo_filter_12_EUR_hg38")
miami_combo_filter <- filter(miami_combo, pval < 0.1)
miami_combo_filter <- miami_combo_filter %>% filter(rsid != '.')

# Highlight novel loci green.
studyA_highlight_rsid <- fread("MultiAncestry_rsID_highlights_final_plot.txt")

# Highlight novel loci green (overlap multi) and magenta (EUR novel).
studyB_highlight_rsid <- fread("EUR_rsID_highlights_plot.txt")

#add to plot
p <- ggmiami(data = miami_combo_filter, split_by = "study", split_at = "A", p = "pval", 
             upper_ylab = "Multi-Ancestry", lower_ylab = "European-only",
             suggestive_line = NULL,
             upper_highlight = studyA_highlight_rsid$rsid, 
             upper_highlight_col = "rsid", 
             upper_highlight_color = studyA_highlight_rsid$color,
             lower_highlight = studyB_highlight_rsid$rsid, 
             lower_highlight_col = "rsid", 
             lower_highlight_color = studyB_highlight_rsid$color)


#ggsave(p, filename = "MiamiPlot_Multi_EUR_updated.png", device = "png", dpi=300, width = 180, height = 115, units = "mm")
ggsave(p, filename = "MiamiPlot_Multi_EUR_updated.tiff", device = "tiff", dpi=300, width = 180, height = 115, units = "mm")