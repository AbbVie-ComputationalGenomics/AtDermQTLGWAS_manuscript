library(ggplot2)
library(ggforestplot)

fpkm <- read.csv('data/FPKM_Epidermal_Raft_Data.csv',header=T,check.names = T,row.names = 1)
dat <- read.csv('data/DE.raft.Differentiation.csv')

gene_names <- fpkm$symbol
gene_ids <- rownames(fpkm)
gene_dict <- cbind(gene_ids,gene_names)
rownames(gene_dict) <- gene_dict[,'gene_ids']
have_na <- rownames(fpkm)[apply(fpkm,1,function(x) is.na(x))[2,]] # some genes have NA values, eliminate them
fpkm <- fpkm[!rownames(fpkm)%in%have_na,]
log2fpkm <- log2(fpkm[,-1]+0.1)
colocscore <- read.table('data/coloc_score.tsv',header = T)
markers <- c('KRT1','KRT10','FLG','LOR') # keratinocyte diff. markers
markers <- cbind(markers,'Diff. Marker')
colnames(markers) <- c('Gene','Prioritized')
colocscore <- rbind(colocscore,markers)
rownames(colocscore) <- colocscore$Gene
selected <- gene_dict[gene_dict[,2]%in%colocscore$Gene,] 
selected <- selected[!duplicated(selected[,2]),] # SCAMP3 is duplicated!
#selected <- selected[-18,] # MAP3K14 is missing!
rownames(selected) <- selected[,1]
#colocscore<- merge(colocscore,selected,by.x=1,by.y=2)
#rownames(colocscore) <- colocscore$gene_ids

## FPKM trajectories
log2fpkm <- log2fpkm[c(rownames(selected)),]
rownames(log2fpkm) <- c(selected[,2])
log2fpkm<-reshape2::melt(cbind(rownames(log2fpkm),log2fpkm),id.vars = 1)
replicates <- unlist(lapply(str_split(log2fpkm$variable,'\\.'),function(x) x[2]))
# Model 0="Subconfluent", 1="D0confluent", 3="D3confluent", 4="D3raft", 7-="D6raft", 10="D9raft", 13="D12raft"   
timepoints <- as.numeric(str_extract(str_extract(unlist(lapply(str_split(gsub('D3confluent','D2confluent',gsub('Subconfluent','D-1Subconfluent',log2fpkm$variable)),'\\.'),function(x) x[3])),'D\\d+'),'\\d+'))+1
timepoints[is.na(timepoints)] <- 0
states <- str_replace(unlist(lapply(str_split(log2fpkm$variable,'\\.'),function(x) x[3])),'D\\d+','')
log2fpkm <- cbind(states,timepoints,replicates,log2fpkm)
colnames(log2fpkm) <- c('State','Timepoint','Replicate','Gene','Id','Expression')
log2fpkm$Prioritized <- colocscore[log2fpkm$Gene,'Prioritized']

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
df2 <- data_summary(log2fpkm, varname="Expression", 
                    groupnames=c("Gene", "Timepoint","State","Prioritized"))
orderedgenes <- c(rownames(subset(colocscore,Prioritized%in%'Diff. Marker')),rownames(subset(colocscore,Prioritized%in%'Y')),rownames(subset(colocscore,Prioritized%in%'N')))
orderedgenes <- orderedgenes[!orderedgenes%in%'MAP3K14']
df2$Gene <- factor(df2$Gene,levels = orderedgenes)

# Default bar plot
p2<- ggplot(df2, aes(x=factor(Timepoint), y=Expression)) +
  geom_line(aes(group = Gene)) +
  geom_errorbar(aes(ymin=Expression-sd, ymax=Expression+sd), width=.2,
                position=position_dodge(.9)) +
  geom_point(aes(fill=Prioritized),shape=21,size=2) +
  facet_wrap(.~Gene, scales = 'free_y') +
  scale_fill_manual(values = c("yellow","beige","green")) +
  theme_bw() +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=8)) +
  ylab("log2(fpkm+0.1)") + xlab("Timepoint") +
  scale_x_discrete(guide = guide_axis(angle = 45),labels=c("0" = "Subconfluent", 
                            "1" = "Confluent D0",
                            "3" = "Confluent D3",
                            "4" = "Raft D3",
                            "7" = "Raft D6",
                            "10" = "Raft D9",
                            "13" = "Raft D12")) +
  theme(axis.text.x = element_text(size=8), plot.title = element_text(hjust = 0.5),legend.position="bottom",legend.title.align = 0.5)
  
pdf('fig.keratinocyte_differentiation_raft.trajectories.pdf',width = 8,height = 5.5)
p2
dev.off()
                            
## DE forestplots (FIX, introduced a mistake)

selected <- subset(as.data.frame(selected), !gene_names%in%markers)
dat <- subset(dat,gene_ids%in%rownames(selected))
#dat <- merge(dat,colocscore,by.x=1,by.y=4,sort = F)
dat$SE <- (dat$CI.R-dat$CI.L)/3.92

prioritized_genes <- subset(colocscore,Prioritized%in%'Y','Gene')[,1]
datp <- subset(dat,gene_names%in%prioritized_genes)
datp$Gene <- factor(datp$gene_names,levels=rev(prioritized_genes))
datnp <- subset(dat,!gene_names%in%c(prioritized_genes,markers[,1]))


pdf('fig.keratinocyte_differentiation_raft.forestplot.pdf',width = 1.5, height = 6) # 127W x 222.25H (mm)
ggplot(datp,aes(x=logFC, y=Gene))+
  geom_vline(xintercept=0, linetype=2)+
  geom_stripes(aes(y=Gene),inherit.aes = FALSE)+
  geom_linerange(aes(xmin=logFC-1.96*SE,
                     xmax=logFC+1.96*SE),
                 position=position_dodge(width=0.6),
                 size = 0.5,
                 color='black')+
  geom_point(shape=21,fill='limegreen',position =position_dodge(width=0.6),size = 2.75,colour = "black")+
  theme_forest()+
  ggtitle('Keratinocyte\ndifferentiation')+
  #labs(tag='B') +
  #theme(plot.tag = element_text(face='bold'))+ 
  guides(fill=FALSE) + 
  xlab("Log(FC)")+
  theme(axis.title.y = element_blank(), axis.text.x = element_text(size=8), axis.text.y = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="none",legend.title.align = 0.5) +
  theme(plot.title = element_text(size=12))
dev.off()

pdf('fig.keratinocyte_differentiation_raft.forestplot.w_legend.pdf',width = 3.74, height = 6) # 127W x 222.25H (mm)
ggplot(datp,aes(x=logFC,y=Gene))+
  geom_vline(xintercept=0, linetype=2)+
  geom_stripes(aes(y=Gene),inherit.aes = FALSE)+
  geom_linerange(aes(xmin=logFC-1.96*SE,
                     xmax=logFC+1.96*SE),
                 position=position_dodge(width=0.6),
                 size = 0.5,
                 color='black')+
  geom_point(shape=21,fill='limegreen',position =position_dodge(width=0.6),size = 3,colour = "black")+
  theme_forest()+
  ggtitle('Keratinocyte\ndifferentiation')+
  #labs(tag='B') +
  #theme(plot.tag = element_text(face='bold'))+ 
  guides(fill=FALSE) + 
  xlab("Log(FC)")+
  theme(axis.title.y = element_blank(), axis.text.x = element_text(size=8), axis.text.y = element_text(size=8), plot.title = element_text(hjust = 0.5),legend.position="left",legend.title.align = 0.5) +
  theme(plot.title = element_text(size=12))
dev.off()

nonprioritized_genes <- subset(colocscore,Prioritized%in%'N','Gene')[,1]
datnp <- subset(dat,gene_names%in%nonprioritized_genes)
datnp[nrow(datp)+1,] <- NA
datnp$gene_names[length(datp$gene_names)] <- 'MAP3K14' 
datnp$Gene <- factor(datnp$gene_names,levels=rev(nonprioritized_genes))

pdf('fig.keratinocyte_differentiation_raft.forestplot.nonprioritized.pdf',width = 1.5, height = 8) # 127W x 222.25H (mm)
ggplot(datnp,aes(x=logFC, y=Gene))+
  geom_vline(xintercept=0, linetype=2)+
  geom_stripes(aes(y=Gene),inherit.aes = FALSE)+
  geom_linerange(aes(xmin=logFC-1.96*SE,
                     xmax=logFC+1.96*SE),
                 position=position_dodge(width=0.6),
                 size = 0.5,
                 color='black')+
  geom_point(shape=21,fill='limegreen',position =position_dodge(width=0.6),size = 2.75,colour = "black")+
  theme_forest()+
  ggtitle('Keratinocyte\ndifferentiation')+
  #labs(tag='B') +
  #theme(plot.tag = element_text(face='bold'))+ 
  guides(fill=FALSE) + 
  xlab("Log(FC)")+
  theme(axis.title.y = element_blank(), axis.text.x = element_text(size=8), axis.text.y = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="none",legend.title.align = 0.5) +
  theme(plot.title = element_text(size=12))
dev.off()

pdf('fig.keratinocyte_differentiation_raft.forestplot.w_legend.nonprioritized.pdf',width = 3, height = 8) # 127W x 222.25H (mm)
ggplot(datp,aes(x=logFC,y=Gene))+
  geom_vline(xintercept=0, linetype=2)+
  geom_stripes(aes(y=Gene),inherit.aes = FALSE)+
  geom_linerange(aes(xmin=logFC-1.96*SE,
                     xmax=logFC+1.96*SE),
                 position=position_dodge(width=0.6),
                 size = 0.5,
                 color='black')+
  geom_point(shape=21,fill='limegreen',position =position_dodge(width=0.6),size = 2.75,colour = "black")+
  theme_forest()+
  ggtitle('Keratinocyte\ndifferentiation')+
  #labs(tag='B') +
  #theme(plot.tag = element_text(face='bold'))+ 
  guides(fill=FALSE) + 
  xlab("Log(FC)")+
  theme(axis.title.y = element_blank(), axis.text.x = element_text(size=8), plot.title = element_text(hjust = 0.5),legend.position="left",legend.title.align = 0.5) +
  theme(plot.title = element_text(size=12))
dev.off()

# Compare abs(logfc) between prioritized 5/22 and non prioritized 17/22 genes
wilcox.test(abs(datp$logFC),abs(datnp$logFC))

