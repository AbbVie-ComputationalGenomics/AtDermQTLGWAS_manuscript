# Data obtained from UMichigan from emails with mksarkar <mksarkar@med.umich.edu>
library(ComplexHeatmap)
library(reshape2)
library(metafor)
library(stringr)
library(ggplot2)
source('code/geom_stripes.R')


filename1='data/1st_siRNA_batch_of_P43-14_For_AbbVie.tsv'
filename2='data/2nd_siRNA_batch_of_P43-23_For_AbbVie.tsv'
filename3='data/3rd_siRNA_batch_of_43-26_For_AbbVie.tsv'
filename4='data/4th_siRNA_batch_of_43-34_For_AbbVie.tsv'

#data <- na.omit(rbind(read.table(filename1,header=T),read.table(filename2,header=T)))
data <- na.omit(rbind(read.table(filename1,header=T),read.table(filename2,header=T),read.table(filename3,header=T),read.table(filename4,header=T)))
data <- subset(data,Include%in%'Y')# exclude SPINK7, IL22RA
data <- subset(data,!BiologicalProcess%in%'KeratinocyteDifferentiation')# exclude KeratinocyteDifferentiation
data$BiologicalProcess<-gsub('PositiveRegulation','',data$BiologicalProcess)

data$MarkerBiologicalProcess <- as.factor(data$MarkerBiologicalProcess)
missinggenes <- c('LIME1','MAP3K14')
data[nrow(data)+1,] <- NA
data$CandidateGene[length(data$CandidateGene)] <- 'LIME1' 
data[nrow(data)+1,] <- NA
data$CandidateGene[length(data$CandidateGene)] <- 'MAP3K14' 
data$CandidateGene <- as.factor(as.character(data$CandidateGene))
data$Treatment <- factor(data$Treatment)

colocscore <- read.table('data/coloc_score.tsv',header = T)

# For each candidate target gene, compare target KO vs NTC (non-target control) marker expression, for each marker (N=3) of each pathway (N=3) in each Tx (treatment) condition (N=4)
# Compare means with 'regular' t-tests. There are 3 replicates per experiment
# Generate standarized mean difference (SMD) effect size and corresponding sampling variance using f(x) 'escalc' from 'metafor' v.4.0
# SMD values will be used for score aggregation and meta-analysis purposes, to obtain one score per gene x pathway 
# https://www.metafor-project.org/doku.php/tips:assembling_data_smd

mean_tests <- t(as.data.frame(rep(NA,8)))
colnames(mean_tests) <- c("BiologicalProcess","MarkerBiologicalProcess","Treatment","CandidateGene","Batch","signed_log10P_ttest","yi","vi")
for (candidateGene in levels(data$CandidateGene)[!levels(data$CandidateGene)%in%'NTC']) {
  for (markerBiologicalProcess in levels(data$MarkerBiologicalProcess)) {
    for (treatment in levels(data$Treatment)) {
      for (batch in unique(data$Batch)) {
      subdata_ntc <- subset(data,CandidateGene%in%'NTC' & MarkerBiologicalProcess%in%markerBiologicalProcess & Treatment%in%treatment & Batch%in%batch)
      subdata <- subset(data,CandidateGene%in%candidateGene & MarkerBiologicalProcess%in%markerBiologicalProcess & Treatment%in%treatment & Batch%in%batch)
      if(nrow(subdata)<1) {next}
      signed_log10P_ttest <- -log10(t.test(subdata_ntc$Value,subdata$Value)$p.value)*sign(t.test(subdata$Value,subdata_ntc$Value)$statistic)
      #mean_tests <- rbind(mean_tests,cbind(unique(subdata[c('BiologicalProcess','MarkerBiologicalProcess','Treatment','CandidateGene','Batch')]),signed_log10P_ttest))
      smd_stats <- escalc(measure="SMD", m1i=mean(subdata$Value), sd1i=sd(subdata$Value), n1i=3,
             m2i=mean(subdata_ntc$Value), sd2i=sd(subdata_ntc$Value), n2i=3)
      mean_tests <- rbind(mean_tests,cbind(unique(subdata[c('BiologicalProcess','MarkerBiologicalProcess','Treatment','CandidateGene','Batch')]),signed_log10P_ttest,as.matrix(smd_stats)))
      }
    }
  }
}
mean_tests <- mean_tests[-1,]
#colnames(mean_tests) <- c("BiologicalProcess","MarkerBiologicalProcess","Treatment","CandidateGene","Batch","signed_log10P_ttest","yi","vi")

mean_tests$Treatment <- factor(mean_tests$Treatment,levels=c('NoTx','IL-22','IL-13','IL-13+IL-22'))
mean_tests_m <- reshape2::dcast(mean_tests[,c('CandidateGene','MarkerBiologicalProcess','Treatment','signed_log10P_ttest')],CandidateGene~MarkerBiologicalProcess+Treatment)
rownames(mean_tests_m) <- mean_tests_m[,1]; mean_tests_m <- mean_tests_m[,-1]

# Identify most extreme signed_log10P_ttest value and corresponding Tx, per marker-gene combination
m <- na.omit(reshape2::melt(tapply(mean_tests$signed_log10P_ttest,list(mean_tests$CandidateGene,mean_tests$MarkerBiologicalProcess,mean_tests$BiologicalProcess),function(x) max(abs(na.omit(x))))*tapply(mean_tests$signed_log10P_ttest,list(mean_tests$CandidateGene,mean_tests$MarkerBiologicalProcess,mean_tests$BiologicalProcess),function(x) sign(max(na.omit(x))-abs(min(na.omit(x)))))))
max_mean_tests <- merge(mean_tests,m,by.x=c('BiologicalProcess','MarkerBiologicalProcess','CandidateGene','signed_log10P_ttest'),by.y=c('Var3','Var2','Var1','value'))
max_mean_tests_m <- reshape2::dcast(max_mean_tests[,c('CandidateGene','MarkerBiologicalProcess','signed_log10P_ttest')],CandidateGene~MarkerBiologicalProcess)
rownames(max_mean_tests_m) <- max_mean_tests_m[,1]; max_mean_tests_m <- max_mean_tests_m[,-1]

# Select double cyto. Tx  signed_log10P_ttest value per marker-gene combination
doublec_mean_tests <- subset(mean_tests,Treatment%in%'IL-13+IL-22')
doublec_mean_tests_m <- reshape2::dcast(doublec_mean_tests[,c('CandidateGene','MarkerBiologicalProcess','signed_log10P_ttest')],CandidateGene~MarkerBiologicalProcess)
rownames(doublec_mean_tests_m) <- doublec_mean_tests_m[,1]; doublec_mean_tests_m <- doublec_mean_tests_m[,-1]

# Generate coloc annotation for genes
annotation_row <- as.data.frame(subset(colocscore,!Gene%in%missinggenes)[,2])
rownames(annotation_row) <- subset(colocscore,!Gene%in%missinggenes)[,1]
colnames(annotation_row) <- 'Prioritized'

# Generate annotation for conditions
MarkerBiologicalProcess_array <- unique(mean_tests[,c('MarkerBiologicalProcess','BiologicalProcess')])[,2]
names(MarkerBiologicalProcess_array) <- unique(mean_tests[,c('MarkerBiologicalProcess','BiologicalProcess')])[,1]
annotation_col <- cbind(unlist(lapply(str_split(colnames(mean_tests_m),'_'),function(x) x[1])),unlist(lapply(str_split(colnames(mean_tests_m),'_'),function(x) x[2])))
annotation_col <- cbind(annotation_col,MarkerBiologicalProcess_array[annotation_col[,1]])
rownames(annotation_col) <- colnames(mean_tests_m)
colnames(annotation_col) <- c('Marker','Tx','Pathway')
annotation_col <- as.data.frame(annotation_col)
annotation_col$Tx <- factor(annotation_col$Tx,levels=c('NoTx','IL-22','IL-13','IL-13+IL-22'))

paletteLength <- 500
myColor <- colorRampPalette(c("blue", "skyblue1","white", "orange","red"))(paletteLength)
# use floor and ceiling to deal with even/odd length pallettelengths
top_val <- max(abs(c(min(mean_tests_m,na.rm = T),max(mean_tests_m,na.rm = T))))
myBreaks <- c(seq(-top_val, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(top_val/paletteLength, top_val, length.out=floor(paletteLength/2)))
              
SignedLog10P <- as.data.frame(mean_tests_m)
#SignedLog10P[nrow(SignedLog10P)+1,] <- NA
#rownames(SignedLog10P)[nrow(SignedLog10P)] <- 'LIME1' 
#SignedLog10P[nrow(SignedLog10P)+1,] <- NA
#rownames(SignedLog10P)[nrow(SignedLog10P)] <- 'MAP3K14' 
SignedLog10P[is.na(SignedLog10P)]<-0
SignedLog10P <- SignedLog10P[order(rownames(SignedLog10P)),]
SignedLog10P <- as.matrix(SignedLog10P)


# Valuable resource to select colors
# https://emilhvitfeldt.github.io/r-color-palettes/discrete.html

pdf('fig.keratinocyte_siRNA.all.pdf',width = 8,height = 6)
p1<-pheatmap(SignedLog10P, 
             name='SignedLog10P',
             annotation_col = annotation_col,
             annotation_row = annotation_row, 
             row_split = annotation_row$Prioritized,
             column_split = annotation_col$Pathway,
             cluster_cols=FALSE,
             cluster_rows=FALSE, 
             color=myColor,
             breaks=myBreaks,
             cell_fun = function(j, i, x, y, width, height, fill) {
               if(SignedLog10P[i, j] > 2 || SignedLog10P[i, j] < -2)
                 grid.text(sprintf("%.1f", SignedLog10P[i, j]), x, y, gp = gpar(fontsize = 7, col='black'))
               else if((SignedLog10P[i, j] > 1.30103 && SignedLog10P[i, j] < 2)  || (SignedLog10P[i, j] < -1.30103 && SignedLog10P[i, j] > -2))
                 grid.text(sprintf("%.1f", SignedLog10P[i, j]), x, y, gp = gpar(fontsize = 7, col='darkgrey'))
             }, 
             #annotation_legend = FALSE, 
             #legend=FALSE, 
             fontsize = 8,
             annotation_colors = list(Prioritized = c(Y = "#03A62CFF", N = "#C5E1A5FF"),
                                      Pathway = c(`IL-13` = "#00FF00FF", `IL-22` = "#EE1289FF"),
                                      Tx = c(NoTx = "White", `IL-13` = "#00FF00FF", `IL-22` = "#EE1289FF", `IL-13+IL-22`="#A020F0FF"),
                                      Marker = c(CCL26 = "red", CISH = "blue", HSD3B1 = "yellow", S100A7 = "brown", S100A8 = "violet", S100A9 = "orange"))
)
print(p1)
dev.off()

eff <- read.table('../../data/keratinocyte_siRNA/KO_efficiency.tsv',header=T)
eff$id <- paste0(eff$CandidateGene,eff$siRNA,eff$Treatment)

# A: AD keratinocyte-linked gene candidate. 
# B: Silence RNA (siRNA) type, either targeting control or candidate gene
# C: Interleukin treatment applied to plated keratinocytes
# D: Expression of gene in A, as measured by qPCR
# E:  Experiment batch
#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables

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
df2 <- data_summary(eff, varname="Expression", 
                    groupnames=c("CandidateGene", "siRNA","Treatment"))

# Default bar plot
p2<- ggplot(df2, aes(x=Treatment, y=Expression, fill=siRNA)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Expression-sd, ymax=Expression+sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(.~CandidateGene, scales = 'free_y') +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=8))
pdf('fig.keratinocyte_siRNA.KO_eff.pdf',width = 8, height = 4.5)
print(p2)
dev.off()

pdf('fig.keratinocyte_siRNA.agg_Tx.all.pdf',width = 15)

annotation_col2 <- unique(annotation_col[,c('Marker','Pathway')])
rownames(annotation_col2) <- annotation_col2$Marker
annotation_col2$Marker <- NULL

SignedLog10P <- as.matrix(max_mean_tests_m)
SignedLog10P[is.na(SignedLog10P)]<-0
pheatmap(SignedLog10P, 
         name='SignedLog10P',
         annotation_col = annotation_col2,
         annotation_row = annotation_row, 
         row_split = annotation_row$colocScore,
         column_split = annotation_col2$Pathway,
         cluster_cols=FALSE,
         cluster_rows=FALSE,
         color=myColor,
         breaks=myBreaks,
         cell_fun = function(j, i, x, y, width, height, fill) {
           if(SignedLog10P[i, j] > 2 || SignedLog10P[i, j] < -2)
             grid.text(sprintf("%.1f", SignedLog10P[i, j]), x, y, gp = gpar(fontsize = 10, col='black'))
           else if((SignedLog10P[i, j] > 1.30103 && SignedLog10P[i, j] < 2)  || (SignedLog10P[i, j] < -1.30103 && SignedLog10P[i, j] > -2))
             grid.text(sprintf("%.1f", SignedLog10P[i, j]), x, y, gp = gpar(fontsize = 10, col='darkgrey'))
         },
         main = 'Aggregation:Max Value across Tx',
         annotation_colors = list(Prioritized = c(Y = "#03A62CFF", N = "#C5E1A5FF"),
                                  Pathway = c(`IL-13` = "#00FF00FF", `IL-22` = "#EE1289FF"),
                                  Tx = c(NoTx = "White", `IL-13` = "#00FF00FF", `IL-22` = "#EE1289FF", `IL-13+IL-22`="#A020F0FF"))
)

SignedLog10P <- as.matrix(doublec_mean_tests_m)
SignedLog10P[is.na(SignedLog10P)]<-0

pheatmap(SignedLog10P, 
         name='SignedLog10P',
         annotation_col = annotation_col2,
         annotation_row = annotation_row, 
         row_split = annotation_row$colocScore,
         column_split = annotation_col2$Pathway,
         cluster_cols=FALSE,
         cluster_rows=FALSE,
         color=myColor,
         breaks=myBreaks,
         cell_fun = function(j, i, x, y, width, height, fill) {
           if(SignedLog10P[i, j] > 2 || SignedLog10P[i, j] < -2)
             grid.text(sprintf("%.1f", SignedLog10P[i, j]), x, y, gp = gpar(fontsize = 10, col='black'))
           else if((SignedLog10P[i, j] > 1.30103 && SignedLog10P[i, j] < 2)  || (SignedLog10P[i, j] < -1.30103 && SignedLog10P[i, j] > -2))
             grid.text(sprintf("%.1f", SignedLog10P[i, j]), x, y, gp = gpar(fontsize = 10, col='darkgrey'))
         },
         main = 'Aggregation:IL-22+IL-13',
         annotation_colors = list(Prioritized = c(Y = "#03A62CFF", N = "#C5E1A5FF"),
                                  Pathway = c(`IL-13` = "#00FF00FF", `IL-22` = "#EE1289FF"),
                                  Tx = c(NoTx = "White", `IL-13` = "#00FF00FF", `IL-22` = "#EE1289FF", `IL-13+IL-22`="#A020F0FF"))
)
      
dev.off()

# Meta-analyze across markers per pathway
# https://www.metafor-project.org/doku.php/tips:assembling_data_smd
xmarker_max_metaanalysis <-as.data.frame(rep(NA,5))
xmarker_doublec_metaanalysis <-as.data.frame(rep(NA,5))
i=1
for (candidateGene in unique(max_mean_tests$CandidateGene)) {
  for (biologicalProcess in unique(max_mean_tests$BiologicalProcess)) {
    
    subdata_max <- subset(max_mean_tests,CandidateGene%in%candidateGene & BiologicalProcess%in%biologicalProcess)
    subdata_doublec <- subset(doublec_mean_tests,CandidateGene%in%candidateGene & BiologicalProcess%in%biologicalProcess)
    
    xmarker_max_metaanalysis_i <- rma(yi=subdata_max$yi, vi=subdata_max$vi)
    xmarker_max_metaanalysis <- cbind(xmarker_max_metaanalysis,unlist(xmarker_max_metaanalysis_i[c(2,3,5,6,7)]))
    colnames(xmarker_max_metaanalysis)[i+1] <- paste(candidateGene,biologicalProcess)
    
    xmarker_doublec_metaanalysis_i <- rma(yi=subdata_doublec$yi, vi=subdata_doublec$vi)
    xmarker_doublec_metaanalysis <- cbind(xmarker_doublec_metaanalysis,unlist(xmarker_doublec_metaanalysis_i[c(2,3,5,6,7)]))
    colnames(xmarker_doublec_metaanalysis)[i+1] <- paste(candidateGene,biologicalProcess)
    i=i+1
    
  }
}  
xmarker_max_metaanalysis <- xmarker_max_metaanalysis[,-1]
xmarker_doublec_metaanalysis <- xmarker_doublec_metaanalysis[,-1]
xmarker_max_metaanalysis<-as.data.frame(t(xmarker_max_metaanalysis))
xmarker_doublec_metaanalysis<-as.data.frame(t(xmarker_doublec_metaanalysis))
g <- unlist(lapply(str_split(rownames(xmarker_doublec_metaanalysis)," "),function(x) x[1]))
p <- unlist(lapply(str_split(rownames(xmarker_doublec_metaanalysis)," "),function(x) x[2]))
Signedlog10P<- sign(xmarker_doublec_metaanalysis$beta)*(-log10(xmarker_doublec_metaanalysis$pval))
xmarker_doublec_metaanalysis <- cbind(g,p,Signedlog10P,xmarker_doublec_metaanalysis)
Signedlog10P<- sign(xmarker_max_metaanalysis$beta)*(-log10(xmarker_max_metaanalysis$pval))
xmarker_max_metaanalysis <- cbind(g,p,Signedlog10P,xmarker_max_metaanalysis)
xmarker_doublec_metaanalysis_m <- reshape2::dcast(xmarker_doublec_metaanalysis[,c('g','p','Signedlog10P')],g~p)
xmarker_max_metaanalysis_m <- reshape2::dcast(xmarker_max_metaanalysis[,c('g','p','pval','Signedlog10P')],g~p)
rownames(xmarker_doublec_metaanalysis_m) <- xmarker_doublec_metaanalysis_m$g; xmarker_doublec_metaanalysis_m<-xmarker_doublec_metaanalysis_m[,-1]
rownames(xmarker_max_metaanalysis_m) <- xmarker_max_metaanalysis_m$g; xmarker_max_metaanalysis_m <- xmarker_max_metaanalysis_m[,-1]

pdf('fig.keratinocyte_siRNA.agg_Tx_agg_Marker.pdf')
top_val2 <- max(abs(c(min(xmarker_max_metaanalysis_m,na.rm = T),max(xmarker_max_metaanalysis_m,na.rm = T))))
myBreaks2 <- c(seq(-top_val2, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(top_val2/paletteLength, top_val2, length.out=floor(paletteLength/2)))
pheatmap(as.matrix(xmarker_max_metaanalysis_m[,c(1,2)]), 
         name='SignedLog10P',
         cluster_rows=FALSE,
         cluster_cols=FALSE,
         row_split = annotation_row$colocScore,
         display_numbers = TRUE,
         annotation_row = annotation_row,          
         color=myColor,
         breaks=myBreaks,
         main = 'Aggregation:Max Value across Tx',
         annotation_colors = list(Prioritized = c(Y = "#03A62CFF", N = "#C5E1A5FF"))
)
top_val2 <- max(abs(c(min(xmarker_doublec_metaanalysis_m,na.rm = T),max(xmarker_max_metaanalysis_m,na.rm = T))))
myBreaks2 <- c(seq(-top_val2, 0, length.out=ceiling(paletteLength/2) + 1), 
               seq(top_val2/paletteLength, top_val2, length.out=floor(paletteLength/2)))
pheatmap(as.matrix(xmarker_doublec_metaanalysis_m[,c(1,2)]), 
         name='SignedLog10P',
         cluster_rows=FALSE,
         cluster_cols=FALSE,
         row_split = annotation_row$colocScore,
         display_numbers = TRUE,
         annotation_row = annotation_row,          
         color=myColor,
         breaks=myBreaks,
         main = 'Aggregation:IL-22+IL-13',
         annotation_colors = list(Prioritized = c(Y = "#03A62CFF", N = "#C5E1A5FF"))
)
dev.off()


#xmarker_doublec_metaanalysis <- merge(xmarker_doublec_metaanalysis,colocscore,by.x='g',by.y=1,sort = F)
prioritized_genes <- subset(colocscore,Prioritized%in%'Y','Gene')[,1]
xmarker_doublec_metaanalysisp <- subset(xmarker_doublec_metaanalysis,g%in%prioritized_genes)
xmarker_doublec_metaanalysisp$g <- factor(xmarker_doublec_metaanalysisp$g,levels=rev(prioritized_genes))

pdf('fig.keratinocyte_siRNA.forestplot.pdf',width = 1.5, height = 6) 
ggplot(xmarker_doublec_metaanalysisp,
       aes(x=beta,
           y=g,
           shape=p,
           group=p))+
  geom_vline(xintercept=0, linetype=2)+
  geom_stripes(aes(y=g),inherit.aes = FALSE)+
  geom_linerange(aes(xmin=beta-1.96*se,
                     xmax=beta+1.96*se,
                     group=p),
                 position=position_dodge(width=0.6),
                 size = 0.5,
                 color='black')+
  geom_point(fill='limegreen',position =position_dodge(width=0.6),size = 2.5,colour = "limegreen",shape=21)+
  geom_point(fill='limegreen',position =position_dodge(width=0.6),size = 2.75,colour = "black")+
  scale_shape_manual(
    values = c(10,13))+
  theme_forest()+
  #ggtitle('Keratinocyte\ndifferentiation')+
  ggtitle('IL-22 and IL-13\npathways')+
  #labs(tag='C') +
  #theme(plot.tag = element_text(face='bold'))+ 
  guides(fill=FALSE) + 
  xlab("SMD")+
  guides(shape=guide_legend(title="Pathway"))+
  guides(color=guide_legend(title="Colocalization\nscore"))+
  theme(legend.title.align = 0.5)+
  theme(axis.title.y = element_blank())+
  guides(color="none")+
  #theme(legend.position="top") # Add legend in Affinity Designer
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size=8), legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank())+
  theme(plot.title = element_text(size=12))
dev.off()

pdf('fig.keratinocyte_siRNA.forestplot.w_legend.pdf',width = 3.74, height = 6) 
ggplot(xmarker_doublec_metaanalysisp,
       aes(x=beta,
           y=g,
           shape=p,
           group=p))+
  geom_vline(xintercept=0, linetype=2)+
  geom_stripes(aes(y=g),inherit.aes = FALSE)+
  geom_linerange(aes(xmin=beta-1.96*se,
                     xmax=beta+1.96*se,
                     group=p),
                 position=position_dodge(width=0.6),
                 size = 0.5,
                 color='black')+
  geom_point(fill='limegreen',position =position_dodge(width=0.6),size = 2.5,colour = "limegreen",shape=21)+
  geom_point(fill='limegreen',position =position_dodge(width=0.6),size = 2.75,colour = "black")+
  scale_shape_manual(
    values = c(10,13))+
  theme_forest()+
  #ggtitle('Keratinocyte\ndifferentiation')+
  ggtitle('IL-22 and IL-13\npathways')+
  #labs(tag='C') +
  #theme(plot.tag = element_text(face='bold'))+ 
  guides(fill=FALSE) + 
  xlab("SMD")+
  guides(shape=guide_legend(title="Pathway"))+
  guides(color=guide_legend(title="Colocalization\nscore"))+
  theme(legend.title.align = 0.5)+
  theme(axis.title.y = element_blank())+
  guides(shape=guide_legend(title="IL-22 and IL-13\nPathways\n(Pathway)"))+
  #theme(legend.position="top") # Add legend in Affinity Designer
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size=8), legend.position = "left", axis.title.y = element_blank()) +
  theme(plot.title = element_text(size=12))
dev.off()

#xmarker_doublec_metaanalysis <- merge(xmarker_doublec_metaanalysis,colocscore,by.x='g',by.y=1,sort = F)
nonprioritized_genes <- subset(colocscore,Prioritized%in%'N','Gene')[,1]
xmarker_doublec_metaanalysisp <- subset(xmarker_doublec_metaanalysis,g%in%nonprioritized_genes)
xmarker_doublec_metaanalysisp[nrow(xmarker_doublec_metaanalysisp)+1,] <- NA
xmarker_doublec_metaanalysisp$g[length(xmarker_doublec_metaanalysisp$g)] <- 'LIME1' 
xmarker_doublec_metaanalysisp[nrow(xmarker_doublec_metaanalysisp)+1,] <- NA
xmarker_doublec_metaanalysisp$g[length(xmarker_doublec_metaanalysisp$g)] <- 'MAP3K14'
xmarker_doublec_metaanalysisp$g <- factor(xmarker_doublec_metaanalysisp$g,levels=rev(nonprioritized_genes))

pdf('fig.keratinocyte_siRNA.forestplot.nonprioritized.pdf',width = 1.5, height = 8) 
ggplot(xmarker_doublec_metaanalysisp,
       aes(x=beta,
           y=g,
           shape=p,
           group=p))+
  geom_vline(xintercept=0, linetype=2)+
  geom_stripes(aes(y=g),inherit.aes = FALSE)+
  geom_linerange(aes(xmin=beta-1.96*se,
                     xmax=beta+1.96*se,
                     group=p),
                 position=position_dodge(width=0.6),
                 size = 0.5,
                 color='black')+
  geom_point(fill='limegreen',position =position_dodge(width=0.6),size = 2.5,colour = "limegreen",shape=21)+
  geom_point(fill='limegreen',position =position_dodge(width=0.6),size = 2.75,colour = "black")+
  scale_shape_manual(
    values = c(10,13))+
  theme_forest()+
  #ggtitle('Keratinocyte\ndifferentiation')+
  ggtitle('IL-22 and IL-13\npathways')+
  #labs(tag='C') +
  #theme(plot.tag = element_text(face='bold'))+ 
  guides(fill=FALSE) + 
  xlab("SMD")+
  guides(shape=guide_legend(title="Pathway"))+
  guides(color=guide_legend(title="Colocalization\nscore"))+
  theme(legend.title.align = 0.5)+
  theme(axis.title.y = element_blank())+
  guides(color="none")+
  #theme(legend.position="top") # Add legend in Affinity Designer
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size=8), legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank())+
  theme(plot.title = element_text(size=12))
dev.off()

pdf('fig.keratinocyte_siRNA.forestplot.w_legend.nonprioritized.pdf',width = 3.74, height = 8) 
ggplot(xmarker_doublec_metaanalysisp,
       aes(x=beta,
           y=g,
           shape=p,
           group=p))+
  geom_vline(xintercept=0, linetype=2)+
  geom_stripes(aes(y=g),inherit.aes = FALSE)+
  geom_linerange(aes(xmin=beta-1.96*se,
                     xmax=beta+1.96*se,
                     group=p),
                 position=position_dodge(width=0.6),
                 size = 0.5,
                 color='black')+
  geom_point(fill='limegreen',position =position_dodge(width=0.6),size = 2.5,colour = "limegreen",shape=21)+
  geom_point(fill='limegreen',position =position_dodge(width=0.6),size = 2.75,colour = "black")+
  scale_shape_manual(
    values = c(10,13))+
  theme_forest()+
  #ggtitle('Keratinocyte\ndifferentiation')+
  ggtitle('IL-22 and IL-13\npathways')+
  #labs(tag='C') +
  #theme(plot.tag = element_text(face='bold'))+ 
  guides(fill=FALSE) + 
  xlab("SMD")+
  guides(shape=guide_legend(title="Pathway"))+
  guides(color=guide_legend(title="Colocalization\nscore"))+
  theme(legend.title.align = 0.5)+
  theme(axis.title.y = element_blank())+
  guides(shape=guide_legend(title="IL-22 and IL-13\nPathways\n(Pathway)"))+
  #theme(legend.position="top") # Add legend in Affinity Designer
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size=8), legend.position = "left", axis.title.y = element_blank())+
  theme(plot.title = element_text(size=12))
dev.off()

# Compare abs(SMD) between prioritized and non-prioritized genes.
wilcox.test(abs(unlist(subset(xmarker_doublec_metaanalysis,g%in%subset(colocscore,Prioritized%in%'N')$Gene,'beta'))),abs(unlist(subset(xmarker_doublec_metaanalysis,g%in%subset(colocscore,Prioritized%in%'Y')$Gene,'beta'))))

