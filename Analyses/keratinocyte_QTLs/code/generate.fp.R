library(ggplot2)
library(ggforestplot)

dat <- read.csv('data/candidates.stats.csv',header=T)
dats <- cbind(tapply(dat$V13,list(dat$V1),mean),tapply(dat$V13,list(dat$V1),sd))
colocscore <- read.table('../keratinocyte_siRNA/data/coloc_score.tsv',header = T)
colocscore <- subset(colocscore,Include%in%'Y')# exclude SPINK7, IL22RA
rownames(colocscore) <- colocscore$CandidateGene
colnames(dats) <- c('MeanEffSize','MeanEffSizeSD')
colocscore <- colocscore[order(colocscore$colocScore),]
dats <- rbind(dats,c(0,0),c(0,0))
rownames(dats)[nrow(dats)]<-'CEBPA'
rownames(dats)[nrow(dats)-1]<-'NAB1'
dats <- merge(dats,colocscore,by.x=0,by.y=0,sort = F)

genes <- c("RORA","CEBPA","ANK3","AQP3","RGS14","WNK1","RTF1","SCAMP3","NAB1")
dats$CandidateGene <- factor(dats$CandidateGene,levels=rev(genes))
dats$Display='Y'
dats$Display[nrow(dats)]<-'N'
dats$Display[nrow(dats)-1]<-'N'
ggplot(dats,aes(x=MeanEffSize,
                y=CandidateGene,
                color=colocScore,
                fill=colocScore,
                alpha=Display))+
  geom_vline(xintercept=0, linetype=2)+
  geom_stripes(aes(y=CandidateGene),inherit.aes = FALSE)+
  geom_linerange(aes(xmin=MeanEffSize-MeanEffSizeSD*2,
                     xmax=MeanEffSize+MeanEffSizeSD*2),
                 position=position_dodge(width=0.6),
                 size = 0.5,
                 color='black')+
  geom_point(position =position_dodge(width=0.6),size = 3.75,colour = "black", shape=25)+
  geom_point(position =position_dodge(width=0.6),size=2.5, shape=25)+
  scale_color_manual(values = c("#03A62CFF", "#98FB98")) +
  scale_fill_manual(values = c("#03A62CFF", "#98FB98")) +
  scale_alpha_manual(values = c(0,1)) +
  theme_forest()+
  ggtitle('Keratinocyte\neQTL')+
  #labs(tag='B') +
  #theme(plot.tag = element_text(face='bold'))+ 
  guides(fill=FALSE) + 
  xlab("eQTL Eff Size")+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none", axis.title.y = element_blank())
