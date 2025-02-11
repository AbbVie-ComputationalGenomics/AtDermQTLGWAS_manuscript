library(ggplot2)
library(ggforestplot)

dat <- read.csv('data/candidates.stats.csv',header=T)
dat2 <- read.table('data/candidates.tsv',header=T,sep='\t')
dat2$AD_riskAllele <- gsub('.*\\|','',dat2$Variant.LD.AD_RiskA)
dat <- merge(dat,dat2,by.x=1,by.y=1)
dat$Sign <- 1
dat$Sign[dat$alt != dat$AD_riskAllele] <- -1
dat$beta <- dat$beta*dat$Sign

dat$beta_se <- abs(dat$beta/qnorm(dat$pvalue/2))
missinggenes <- c('CLIP1','CEBPA','WNK1','NDUFA4L2') # No eQTL stats for those, for diff. reasons. Check Sup Table.
dat <- rbind(dat,merge(dat,unlist(lapply(missinggenes,function(x) rep(x,8))),by.x=1,by.y=1,sort = F,all.y=T))

colocscore <- read.table('../../data/other/coloc_score.tsv',header = T)
rownames(colocscore) <- colocscore$Gene

prioritized_genes <- rownames(subset(colocscore,Prioritized%in%'Y'))
datp <- subset(dat,Gene%in%prioritized_genes)
datp <- merge(prioritized_genes,dat,all.x=T,by.x=1,by.y=1)
colnames(datp)[1] <- 'Gene' 
datp$Gene <- factor(datp$Gene,levels=rev(prioritized_genes))
datp$Stimulus[is.na(datp$Stimulus)] <- datp$Stimulus[1]
betasp <- na.omit(abs(datp$beta))

pdf('fig.keratinocyte_eQTLs.forestplot.pdf',width = 1.5, height = 6)
ggplot(datp,aes(x=beta,
               y=Gene,
               #shape=Stimulus,
               group=Stimulus,
               color=Stimulus,
               fill=Stimulus))+
  geom_vline(xintercept=0, linetype = "longdash")+
  geom_stripes(aes(y=Gene),inherit.aes = FALSE)+
  geom_linerange(aes(xmin=beta-1.96*beta_se,
                     xmax=beta+1.96*beta_se,
                     group=Stimulus),
                 position=position_dodge(width=0.6),
                 size = 0.5,
                 color='black')+
  geom_point(shape=19,position =position_dodge(width=0.6),size = 2.75,colour = "black")+
  geom_point(shape=19,position =position_dodge(width=0.6),size=2)+
  scale_color_manual(values = c("limegreen", "green","red","blue","orange","violet","purple","yellow")) +
  scale_fill_manual(values = c("limegreen", "green","red","blue","orange","violet","purple","yellow")) +
  theme_forest()+
  ggtitle('Keratinocyte\neQTLs')+
  #labs(tag='B') +
  #theme(plot.tag = element_text(face='bold'))+ 
  guides(fill=FALSE) +
  xlab("Effect Size")+
  guides(color=guide_legend(title="Keratinocyte\neQTLs\n(Stimulus)"))+
  theme(axis.title.y = element_blank(), axis.text.x = element_text(size=8), axis.text.y = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="none",legend.title.align = 0.5) +
  theme(plot.title = element_text(size=12))
dev.off()


pdf('fig.keratinocyte_eQTLs.forestplot.w_legend.pdf',width = 3.74, height = 6)
ggplot(datp,aes(x=beta,
                y=Gene,
                #shape=Stimulus,
                group=Stimulus,
                color=Stimulus,
                fill=Stimulus))+
  geom_vline(xintercept=0, linetype = "longdash")+
  geom_stripes(aes(y=Gene),inherit.aes = FALSE)+
  geom_linerange(aes(xmin=beta-1.96*beta_se,
                     xmax=beta+1.96*beta_se,
                     group=Stimulus),
                 position=position_dodge(width=0.6),
                 size = 0.5,
                 color='black')+
  geom_point(shape=19,position =position_dodge(width=0.6),size = 2.75,colour = "black")+
  geom_point(shape=19,position =position_dodge(width=0.6),size=2)+
  scale_color_manual(values = c("limegreen", "green","red","blue","orange","violet","purple","yellow")) +
  scale_fill_manual(values = c("limegreen", "green","red","blue","orange","violet","purple","yellow")) +
  theme_forest()+
  ggtitle('Keratinocyte\neQTLs')+
  #labs(tag='B') +
  #theme(plot.tag = element_text(face='bold'))+ 
  guides(fill=FALSE) +
  xlab("Effect Size")+
  guides(color=guide_legend(title="Keratinocyte\neQTLs\n(Stimulus)"))+
  theme(axis.title.y = element_blank(), axis.text.x = element_text(size=8), plot.title = element_text(hjust = 0.5),legend.position="left",legend.title.align = 0.5) +
  theme(plot.title = element_text(size=12))
dev.off()


nonprioritized_genes <- rownames(subset(colocscore,Prioritized%in%'N'))
datp <- subset(dat,Gene%in%nonprioritized_genes)
datp <- merge(nonprioritized_genes,dat,all.x=T,by.x=1,by.y=1)
colnames(datp)[1] <- 'Gene' #missinggenes <- c("GRID2IP","IL2RB") 
datp$Gene <- factor(datp$Gene,levels=rev(nonprioritized_genes))
datp$Stimulus[is.na(datp$Stimulus)] <- datp$Stimulus[1]
betasnp <- na.omit(abs(datp$beta))


pdf('fig.keratinocyte_eQTLs.forestplot.nonprioritized.pdf',width = 1.5, height = 8)
ggplot(datp,aes(x=beta,
                y=Gene,
                #shape=Stimulus,
                group=Stimulus,
                color=Stimulus,
                fill=Stimulus))+
  geom_vline(xintercept=0, linetype = "longdash")+
  geom_stripes(aes(y=Gene),inherit.aes = FALSE)+
  geom_linerange(aes(xmin=beta-1.96*beta_se,
                     xmax=beta+1.96*beta_se,
                     group=Stimulus),
                 position=position_dodge(width=0.6),
                 size = 0.5,
                 color='black')+
  geom_point(shape=19,position =position_dodge(width=0.6),size = 2.75,colour = "black")+
  geom_point(shape=19,position =position_dodge(width=0.6),size=2)+
  scale_color_manual(values = c("limegreen", "green","red","blue","orange","violet","purple","yellow")) +
  scale_fill_manual(values = c("limegreen", "green","red","blue","orange","violet","purple","yellow")) +
  theme_forest()+
  ggtitle('Keratinocyte\neQTLs')+
  #labs(tag='B') +
  #theme(plot.tag = element_text(face='bold'))+ 
  guides(fill=FALSE) +
  xlab("Effect Size")+
  guides(color=guide_legend(title="Keratinocyte\neQTLs\n(Stimulus)"))+
  theme(axis.title.y = element_blank(), axis.text.x = element_text(size=8), axis.text.y = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="none",legend.title.align = 0.5) +
  theme(plot.title = element_text(size=12))
dev.off()


pdf('fig.keratinocyte_eQTLs.forestplot.w_legend.nonprioritized.pdf',width = 3.74, height = 8)
ggplot(datp,aes(x=beta,
                y=Gene,
                #shape=Stimulus,
                group=Stimulus,
                color=Stimulus,
                fill=Stimulus))+
  geom_vline(xintercept=0, linetype = "longdash")+
  geom_stripes(aes(y=Gene),inherit.aes = FALSE)+
  geom_linerange(aes(xmin=beta-1.96*beta_se,
                     xmax=beta+1.96*beta_se,
                     group=Stimulus),
                 position=position_dodge(width=0.6),
                 size = 0.5,
                 color='black')+
  geom_point(shape=19,position =position_dodge(width=0.6),size = 2.75,colour = "black")+
  geom_point(shape=19,position =position_dodge(width=0.6),size=2)+
  scale_color_manual(values = c("limegreen", "green","red","blue","orange","violet","purple","yellow")) +
  scale_fill_manual(values = c("limegreen", "green","red","blue","orange","violet","purple","yellow")) +
  theme_forest()+
  ggtitle('Keratinocyte\neQTLs')+
  #labs(tag='B') +
  #theme(plot.tag = element_text(face='bold'))+ 
  guides(fill=FALSE) +
  xlab("Effect Size")+
  guides(color=guide_legend(title="Keratinocyte\neQTLs\n(Stimulus)"))+
  theme(axis.title.y = element_blank(), axis.text.x = element_text(size=8), plot.title = element_text(hjust = 0.5),legend.position="left",legend.title.align = 0.5) +
  theme(plot.title = element_text(size=12))
dev.off()

stats <- rbind(unlist(dat[1:8,c('beta','beta_se','pvalue')]),unlist(dat[9:16,c('beta','beta_se','pvalue')]),unlist(dat[17:24,c('beta','beta_se','pvalue')]),unlist(dat[25:32,c('beta','beta_se','pvalue')]),unlist(dat[33:40,c('beta','beta_se','pvalue')]),unlist(dat[41:48,c('beta','beta_se','pvalue')]),unlist(dat[49:56,c('beta','beta_se','pvalue')]),unlist(dat[57:64,c('beta','beta_se','pvalue')]),unlist(dat[65:72,c('beta','beta_se','pvalue')]),unlist(dat[73:80,c('beta','beta_se','pvalue')]),unlist(dat[81:88,c('beta','beta_se','pvalue')]),unlist(dat[89:96,c('beta','beta_se','pvalue')]),unlist(dat[97:104,c('beta','beta_se','pvalue')]),unlist(dat[105:112,c('beta','beta_se','pvalue')]),unlist(dat[113:120,c('beta','beta_se','pvalue')]),unlist(dat[121:128,c('beta','beta_se','pvalue')]),unlist(dat[129:136,c('beta','beta_se','pvalue')]),unlist(dat[137:144,c('beta','beta_se','pvalue')]))
colnames(stats) <- c(paste('beta_',unique(dat$Stimulus),collpase='',sep = ''),paste('beta_se_',unique(dat$Stimulus),collpase='',sep = ''),paste('pavalue_',unique(dat$Stimulus),collpase='',sep = ''))
write.csv(file = 'data/candidates.stats.wide.csv', stats, quote=F)

# Test differences in effect sizes between prioritized and non prioritized genes
wilcox.test(betasp,betasnp)
