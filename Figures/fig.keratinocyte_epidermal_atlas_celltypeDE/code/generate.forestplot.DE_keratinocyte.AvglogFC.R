library(reshape2)
library(ComplexHeatmap)
library(ggforestplot)

cs <- read.table('../../data/other/coloc_score.tsv',header=T)
colocscore <- cs

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
keratinocyte_de_stats_candidates <- keratinocyte_de_stats_all[keratinocyte_de_stats_all$gene%in%colocscore$Gene,]
keratinocyte_de_stats_candidates_m <- reshape2::dcast(keratinocyte_de_stats_candidates[,c('gene','body_site','avg_log2FC','keratinocyte_subtype')],formula = gene~body_site+keratinocyte_subtype,value.var = 'avg_log2FC')
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
dat <- merge(dat,colocscore,by.x=1,by.y=1,sort = F)
dat$Mean_avg_log2FC <- as.numeric(dat$Mean_avg_log2FC)
dat$Mean_avg_log2FC_sd <- as.numeric(dat$Mean_avg_log2FC_sd)
dat$keratinocyte_subtype<-gsub('_Keratinocytes','',dat$keratinocyte_subtype)

prioritized_genes <- subset(colocscore,Prioritized%in%'Y','Gene')[,1]
datp <- subset(dat,Gene%in%prioritized_genes)
datp$Gene <- factor(datp$Gene,levels=rev(prioritized_genes))
datnp <- subset(dat,!Gene%in%prioritized_genes)


pdf('fig.keratinocyte_DE_epidermal_atlas.forestplot.wo_legend.pdf',width = 1.5, height = 6) 
ggplot(datp,aes(x=Mean_avg_log2FC,y=Gene,
               shape=keratinocyte_subtype,
               group=keratinocyte_subtype))+
  geom_vline(xintercept=0, linetype=2)+
  geom_stripes(aes(y=Gene),inherit.aes = FALSE)+
  geom_linerange(aes(xmin=Mean_avg_log2FC-1.96*Mean_avg_log2FC_sd,
                     xmax=Mean_avg_log2FC+1.96*Mean_avg_log2FC_sd,
                     group=keratinocyte_subtype),
                 position=position_dodge(width=0.6),
                 size = 0.5,
                 color='black')+
  scale_shape_manual(
    values = c(0,1,2,5))+
  geom_point(fill='limegreen',position =position_dodge(width=0.6),size = 2.75,colour = "black")+
  scale_shape_manual(
    values = c(22,23,24,25))+
  theme_forest()+
  ggtitle('Keratinocyte\nspecificity')+
  #labs(tag='A') +
  #theme(plot.tag = element_text(face='bold'))+ 
  guides(fill=FALSE) + 
  xlab("Log(FC)")+
  guides(shape=guide_legend(title="Keratinocyte\nsubtype"))+
  guides(color=guide_legend(title="Colocalization\nscore"))+
  theme(axis.title.y = element_blank(), axis.text.x = element_text(size=8), axis.text.y = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="none",legend.title.align = 0.5) +
  theme(plot.title = element_text(size=12))
dev.off()

pdf('fig.keratinocyte_DE_epidermal_atlas.forestplot.w_legend.pdf',width = 3.74, height = 6) # 102W x 222.25H (mm)
ggplot(datp,aes(x=Mean_avg_log2FC,y=Gene,
               shape=keratinocyte_subtype,
               group=keratinocyte_subtype))+
  geom_vline(xintercept=0, linetype=2)+
  geom_stripes(aes(y=Gene),inherit.aes = FALSE)+
  geom_linerange(aes(xmin=Mean_avg_log2FC-1.96*Mean_avg_log2FC_sd,
                     xmax=Mean_avg_log2FC+1.96*Mean_avg_log2FC_sd,
                     group=keratinocyte_subtype),
                 position=position_dodge(width=0.6),
                 size = 0.5,
                 color='black')+
  scale_shape_manual(
    values = c(0,1,2,5))+
  geom_point(fill='limegreen',position =position_dodge(width=0.6),size = 2.75,colour = "black")+
  scale_shape_manual(
    values = c(22,23,24,25))+
  theme_forest()+
  ggtitle('Keratinocyte\nspecificity')+
  #labs(tag='A') +
  #theme(plot.tag = element_text(face='bold'))+ 
  guides(fill=FALSE) + 
  xlab("Log(FC)")+
  guides(shape=guide_legend(title="Keratinocyte\nSpecificity\n(Subtype)"))+
  guides(color=guide_legend(title="Colocalization\nscore"))+
  theme(axis.title.y = element_blank(), axis.text.x = element_text(size=8), plot.title = element_text(hjust = 0.5),legend.position="left",legend.title.align = 0.5) +
  theme(plot.title = element_text(size=12))
dev.off()

nonprioritized_genes <- subset(colocscore,Prioritized%in%'N','Gene')[,1]
datp <- merge(nonprioritized_genes,dat,all.x=T,by.x=1,by.y=1)
colnames(datp)[1] <- 'Gene' #missinggenes <- c("GRID2IP","IL2RB") 
datp$Gene <- factor(datp$Gene,levels=rev(nonprioritized_genes))
datp$keratinocyte_subtype[is.na(datp$keratinocyte_subtype)] <- datp$keratinocyte_subtype[1]

pdf('fig.keratinocyte_DE_epidermal_atlas.forestplot.wo_legend.nonprioritized.pdf',width = 1.5, height = 8) 
ggplot(datp,aes(x=Mean_avg_log2FC,y=Gene,
                shape=keratinocyte_subtype,
                group=keratinocyte_subtype))+
  geom_vline(xintercept=0, linetype=2)+
  geom_stripes(aes(y=Gene),inherit.aes = FALSE)+
  geom_linerange(aes(xmin=Mean_avg_log2FC-1.96*Mean_avg_log2FC_sd,
                     xmax=Mean_avg_log2FC+1.96*Mean_avg_log2FC_sd,
                     group=keratinocyte_subtype),
                 position=position_dodge(width=0.6),
                 size = 0.5,
                 color='black')+
  scale_shape_manual(
    values = c(0,1,2,5))+
  geom_point(fill='limegreen',position =position_dodge(width=0.6),size = 2.75,colour = "black")+
  scale_shape_manual(
    values = c(22,23,24,25))+
  theme_forest()+
  ggtitle('Keratinocyte\nspecificity')+
  #labs(tag='A') +
  #theme(plot.tag = element_text(face='bold'))+ 
  guides(fill=FALSE) + 
  xlab("Log(FC)")+
  guides(shape=guide_legend(title="Keratinocyte\nsubtype"))+
  guides(color=guide_legend(title="Colocalization\nscore"))+
  theme(axis.title.y = element_blank(), axis.text.x = element_text(size=8), axis.text.y = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="none",legend.title.align = 0.5) +
  theme(plot.title = element_text(size=12))
dev.off()

pdf('fig.keratinocyte_DE_epidermal_atlas.forestplot.w_legend.nonprioritized.pdf',width = 3.75, height = 8) # 102W x 222.25H (mm)
ggplot(datp,aes(x=Mean_avg_log2FC,y=Gene,
                shape=keratinocyte_subtype,
                group=keratinocyte_subtype))+
  geom_vline(xintercept=0, linetype=2)+
  geom_stripes(aes(y=Gene),inherit.aes = FALSE)+
  geom_linerange(aes(xmin=Mean_avg_log2FC-1.96*Mean_avg_log2FC_sd,
                     xmax=Mean_avg_log2FC+1.96*Mean_avg_log2FC_sd,
                     group=keratinocyte_subtype),
                 position=position_dodge(width=0.6),
                 size = 0.5,
                 color='black')+
  scale_shape_manual(
    values = c(0,1,2,5))+
  geom_point(fill='limegreen',position =position_dodge(width=0.6),size = 2.75,colour = "black")+
  scale_shape_manual(
    values = c(22,23,24,25))+
  theme_forest()+
  ggtitle('Keratinocyte\nspecificity')+
  #labs(tag='A') +
  #theme(plot.tag = element_text(face='bold'))+ 
  guides(fill=FALSE) + 
  xlab("Log(FC)")+
  guides(shape=guide_legend(title="Keratinocyte\nSpecificity\n(Subtype)"))+
  guides(color=guide_legend(title="Colocalization\nscore"))+
  theme(axis.title.y = element_blank(), axis.text.x = element_text(size=8), plot.title = element_text(hjust = 0.5),legend.position="left",legend.title.align = 0.5) +
  theme(plot.title = element_text(size=12))
dev.off()

# Compare abs log2FC between prioritized and non prioritized genes
wilcox.test(abs(datp$Mean_avg_log2FC),abs(datnp$Mean_avg_log2FC))

