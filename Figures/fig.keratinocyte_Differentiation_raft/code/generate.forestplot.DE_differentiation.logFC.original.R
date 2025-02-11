library(ggplot2)

dat <- read.csv('data/DE.raft.Differentiation.csv')
colocscore <- read.table('data/coloc_score.tsv',header = T)
colocscore <- subset(colocscore,Include%in%'Y')# exclude SPINK7, IL22RA
gene_dict <- cbind(dat$gene_ids,dat$gene_names)
selected <- gene_dict[gene_dict[,2]%in%colocscore$CandidateGene,] # SCAMP3 is duplicated
selected <- selected[-8,] # SCAMP3 is repeated
rownames(selected) <- selected[,1]
colocscore<- merge(colocscore,selected,by.x=1,by.y=2)
rownames(colocscore) <- colocscore$gene_ids
colocscore <- colocscore[order(colocscore$colocScore),]

dat <- subset(dat,gene_ids%in%rownames(selected))
dat <- merge(dat,colocscore,by.x=1,by.y=4,sort = F)
dat$SE <- (dat$CI.R-dat$CI.L)/3.92
genes <- c("RORA","CEBPA","ANK3","AQP3","RGS14","WNK1","RTF1","SCAMP3","NAB1")
dat$gene_names <- factor(dat$gene_names,levels=rev(genes))
pdf('fig.keratinocyte_differentiation_raft.forestplot.pdf',width = 2, height = 8.75) # 127W x 222.25H (mm)
ggplot(dat,aes(x=logFC,
               y=gene_names,
               color=colocScore,
               fill=colocScore))+
  geom_vline(xintercept=0, linetype=2)+
  geom_stripes(aes(y=gene_names),inherit.aes = FALSE)+
  geom_linerange(aes(xmin=logFC-1.96*SE,
                     xmax=logFC+1.96*SE,
                     group=keratinocyte_subtype),
                 position=position_dodge(width=0.6),
                 size = 0.5,
                 color='black')+
  geom_point(position =position_dodge(width=0.6),size = 3.75,colour = "black", shape=25)+
  geom_point(position =position_dodge(width=0.6),size=2.5, shape=25)+
  scale_color_manual(values = c("#03A62CFF", "#98FB98")) +
  scale_fill_manual(values = c("#03A62CFF", "#98FB98")) +
  theme_forest()+
  ggtitle('Keratinocyte\ndifferentiation')+
  #labs(tag='B') +
  #theme(plot.tag = element_text(face='bold'))+ 
  guides(fill=FALSE) + 
  xlab("Log(FC)")+
  theme(plot.title = element_text(hjust = 0.5),legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank())
dev.off()
