library(reshape2)
library(ComplexHeatmap)
library(ggforestplot)

cs <- read.table('../../data/other/coloc_score.tsv',header=T)
colocscore <- subset(cs,Include%in%'Y')# exclude SPINK7, IL22RA

keratinocyte_de_stats_all <- t(data.frame(rep('NA',8)))
colnames(keratinocyte_de_stats_all) <- c('gene','p_val','avg_log2FC','pct.1','pct.2','p_val_adj','body_site','keratinocyte_subtype')
for (body_site in c('acral','arm','axilla','back','face','leg','scalp')) {
  for (keratinocyte_subtype in c('Cycling_Keratinocytes','Keratinized_Keratinocytes','Basal_Keratinocytes','Differentiated_Keratinocytes')) {
    keratinocyte_de_stats_body_site <- read.table(file=paste0('./data/',body_site,'.',keratinocyte_subtype,'.Markers.tsv'),header=T,row.names = NULL)
    colnames(keratinocyte_de_stats_body_site)[1] <- 'gene'
    keratinocyte_de_stats_body_site$body_site <- body_site
    keratinocyte_de_stats_body_site$keratinocyte_subtype <- keratinocyte_subtype
    keratinocyte_de_stats_all <- rbind(keratinocyte_de_stats_all,keratinocyte_de_stats_body_site)
  }
}
keratinocyte_de_stats_all <- keratinocyte_de_stats_all[-1,]
keratinocyte_de_stats_candidates <- keratinocyte_de_stats_all[keratinocyte_de_stats_all$gene%in%colocscore$CandidateGene,]
keratinocyte_de_stats_candidates_m <- dcast(keratinocyte_de_stats_candidates[,c('gene','body_site','avg_log2FC','keratinocyte_subtype')],formula = gene~body_site+keratinocyte_subtype,value.var = 'avg_log2FC')
genes <-  keratinocyte_de_stats_candidates_m[,1]
keratinocyte_de_stats_candidates_m <- keratinocyte_de_stats_candidates_m[,-1]
keratinocyte_de_stats_candidates_m <- data.frame(apply(keratinocyte_de_stats_candidates_m, 2, function(x) as.numeric(as.character(x))))
rownames(keratinocyte_de_stats_candidates_m) <- genes; 

# Generate coloc annotation for genes
annotation_row <- as.data.frame(colocscore[,2])
rownames(annotation_row) <- colocscore[,1]
colnames(annotation_row) <- 'colocScore'

paletteLength <- 500
myColor <- colorRampPalette(c("blue", "skyblue1","white", "orange","red"))(paletteLength)
# use floor and ceiling to deal with even/odd length pallettelengths
top_val <- max(abs(c(min(keratinocyte_de_stats_candidates_m,na.rm = T),max(keratinocyte_de_stats_candidates_m,na.rm = T))))
myBreaks <- c(seq(-top_val, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(top_val/paletteLength, top_val, length.out=floor(paletteLength/2)))

avg_log2FC <- as.matrix(keratinocyte_de_stats_candidates_m)
avg_log2FC[is.na(avg_log2FC)]<-0

# Valuable resource to select colors
# https://emilhvitfeldt.github.io/r-color-palettes/discrete.html

## Do NOT USE p-values of t-test in Seurat, they are inflated 'cause cells are treated as biological replicates.
## To do a proper DE test, pseudo-bulk estimates should be generated
## https://www.biostars.org/p/9473703/

pdf('fig.keratinocyte_DE_epidermal_atlas.pdf',width = 8)
pheatmap(avg_log2FC, 
         name='avg_log2FC',
         #annotation_col = annotation_col,
         annotation_row = annotation_row, 
         row_split = annotation_row$colocScore,
         cluster_cols=FALSE,
         color=myColor,
         breaks=myBreaks,
         main = "Keratinocyte vs other cells",
         annotation_colors = list(colocScore = c(Strong = "#03A62CFF", Weak = "#98FB98"))
)
dev.off()