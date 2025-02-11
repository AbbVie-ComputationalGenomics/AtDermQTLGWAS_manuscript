library(reshape2)
library(ComplexHeatmap)
library(ggforestplot)

cs <- read.table('../../data/other/coloc_score.tsv',header=T)
colocscore <- subset(cs,Include%in%'Y')# exclude SPINK7, IL22RA

keratinocyte_de_stats_all <- t(data.frame(rep('NA',8)))
colnames(keratinocyte_de_stats_all) <- c('gene','p_val','avg_log2FC','pct.1','pct.2','p_val_adj','body_site','keratinocyte_subtype')
#for (body_site in c('acral','arm','axilla','back','face','leg','scalp')) {
for (body_site in c('arm','axilla','back','face','leg','scalp')) { # acral is an outlier.
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


avg_log2FC <- as.matrix(keratinocyte_de_stats_candidates_m)
avg_log2FC[is.na(avg_log2FC)]<-0

# Valuable resource to select colors
# https://emilhvitfeldt.github.io/r-color-palettes/discrete.html

## Do NOT USE p-values of t-test in Seurat, they are inflated 'cause cells are treated as biological replicates.
## To do a proper DE test, pseudo-bulk estimates should be generated
## https://www.biostars.org/p/9473703/
## We will use the average fold change bet. keratinocyte and non-keratinocytes to represent keratinocyte gene specificity

genes <- c("RORA","CEBPA","ANK3","AQP3","RGS14","WNK1","RTF1","SCAMP3","NAB1")
dat <- data.frame(t(rep(NA,4)))
colnames(dat) <- c('Gene','Mean_avg_log2FC','Mean_avg_log2FC_sd','keratinocyte_subtype')
for (keratinocyte_subtype in c('Cycling_Keratinocytes','Keratinized_Keratinocytes','Basal_Keratinocytes','Differentiated_Keratinocytes')) {
  avg_log2FC_sub <- avg_log2FC[,colnames(avg_log2FC)[grep(pattern = keratinocyte_subtype, colnames(avg_log2FC))]]
  dats <- as.data.frame(cbind(rownames(avg_log2FC_sub),apply(avg_log2FC_sub,1,mean),apply(avg_log2FC_sub,1,sd)))
  colnames(dats) <- c('Gene','Mean_avg_log2FC','Mean_avg_log2FC_sd')
  dats <- dats[genes,]
  dats$keratinocyte_subtype <- keratinocyte_subtype
  dat <- rbind(dat,dats)
}
dat <- dat[-1,]
#dat <- cbind(rownames(avg_log2FC),apply(avg_log2FC,1,mean),apply(avg_log2FC,1,sd))
#colnames(dat)<-c('Gene','Mean_avg_log2FC','SD_Mean_avg_log2FC_cross_body_sites')
colocscore <- colocscore[order(colocscore$colocScore),]
dat <- merge(dat,colocscore,by.x=1,by.y=1,sort = F)
dat$Mean_avg_log2FC <- as.numeric(dat$Mean_avg_log2FC)
dat$Mean_avg_log2FC_sd <- as.numeric(dat$Mean_avg_log2FC_sd)
dat$Gene <- factor(dat$Gene,levels=rev(unique(dat$Gene)))
dat$keratinocyte_subtype<-gsub('_Keratinocytes','',dat$keratinocyte_subtype)

pdf('fig.keratinocyte_DE_epidermal_atlas.forestplot.pdf',width = 4, height = 8.75) # 102W x 222.25H (mm)
ggplot(dat,aes(x=Mean_avg_log2FC,y=Gene,
               shape=keratinocyte_subtype,
               group=keratinocyte_subtype,
               color=colocScore,
               fill=colocScore))+
  geom_vline(xintercept=0)+
  geom_stripes(aes(y=Gene),inherit.aes = FALSE)+
  geom_linerange(aes(xmin=Mean_avg_log2FC-1.96*Mean_avg_log2FC_sd,
                     xmax=Mean_avg_log2FC+1.96*Mean_avg_log2FC_sd,
                     group=keratinocyte_subtype),
                 position=position_dodge(width=0.6),
                 size = 0.5,
                 color='black')+
  scale_shape_manual(
    values = c(0,1,2,5))+
  geom_point(position =position_dodge(width=0.6),size = 2.75,colour = "black")+
  geom_point(position =position_dodge(width=0.6),size=2)+
  scale_shape_manual(
    values = c(15,17,18,19))+
  scale_color_manual(values = c("#03A62CFF", "#98FB98")) +
  theme_forest()+
  ggtitle('Keratinocyte\nspecificity')+
  #labs(tag='A') +
  #theme(plot.tag = element_text(face='bold'))+ 
  guides(fill=FALSE) + 
  xlab("Log2(FC)")+
  guides(shape=guide_legend(title="Keratinocyte\nsubtype"))+
  guides(color=guide_legend(title="Colocalization\nscore"))+
  theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5),legend.position="left",legend.title.align = 0.5)
dev.off()

