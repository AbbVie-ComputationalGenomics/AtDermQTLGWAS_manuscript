library(stringr)
library(limma)
library(qqman)
library(zFPKM)

# Ideally, DEG should be performed with counts: https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
# If FPKMs are the only thing available, then:
# https://support.bioconductor.org/p/56275/#56299

fpkm <- read.csv('data/FPKM_Epidermal_Raft_Data.csv',header=T,check.names = T,row.names = 1)
gene_names <- fpkm$symbol
gene_ids <- rownames(fpkm)
gene_dict <- cbind(gene_ids,gene_names)
rownames(gene_dict) <- gene_dict[,'gene_ids']
have_na <- rownames(fpkm)[apply(fpkm,1,function(x) is.na(x))[2,]] # some genes have NA values, eliminate them
fpkm <- fpkm[!rownames(fpkm)%in%have_na,]
log2fpkm <- log2(fpkm[,-1]+0.1)
zFPKM <- zFPKM(fpkm[,-1])
timepoints <- as.numeric(str_extract(str_extract(unlist(lapply(str_split(gsub('D3confluent','D2confluent',gsub('Subconfluent','D-1Subconfluent',colnames(log2fpkm))),'\\.'),function(x) x[3])),'D\\d+'),'\\d+'))+1
timepoints[is.na(timepoints)] <- 0
zFPKMmean <- t(apply(zFPKM,1,function(x) tapply(x,list(timepoints),mean)))
genes_active <- rownames(zFPKMmean)[apply(zFPKMmean>-3,1,sum)>2]
#genes_active <- rownames(zFPKM)[apply(zFPKM>-3,1,sum)>2]
log2fpkm <- log2fpkm[genes_active,] # using zFPKM > -3, in at least 2 timepoints, to select expressed genes.
gene_dict <- gene_dict[genes_active,]

replicates <- unlist(lapply(str_split(colnames(log2fpkm),'\\.'),function(x) x[2]))
states <- str_replace(unlist(lapply(str_split(colnames(log2fpkm),'\\.'),function(x) x[3])),'D\\d+','')
design <- cbind(as.factor(replicates),timepoints)

paste0('There are ',length(unique(timepoints)),' timepoints, ',length(unique(replicates)),' replicates and ',length(unique(states)),' states: ',paste0(unique(states),collapse=','))# https://support.bioconductor.org/p/56275/#56299
plot(density(as.matrix(log2fpkm))) # Eliminate lowly exp. genes? Important to fit correlation.

design <- model.matrix(~0+timepoints)
corfit <- duplicateCorrelation(log2fpkm,design,block=replicates)
fit <- lmFit(log2fpkm,design,block=replicates,correlation=corfit$consensus)
fit <- eBayes(fit)
tt <- topTable(fit, coef="timepoints",number = Inf,confint = T)
tt <- cbind(gene_dict[rownames(tt),],tt)
# https://support.bioconductor.org/p/37524/
write.csv(file = 'DE.raft.Differentiation.csv',row.names = F,quote = F,tt)

#colocscore <- read.table('../keratinocyte_siRNA/data/coloc_score.tsv',header = T)
#colocscore <- subset(colocscore,Include%in%'Y')# exclude SPINK7, IL22RA
#selected <- gene_dict[gene_dict[,2]%in%colocscore$CandidateGene,] # SCAMP3 is duplicated
#selected <- selected[-8,] # SCAMP3 is repeated
#colocscore<- merge(colocscore,selected,by.x=1,by.y=2)
#rownames(colocscore) <- colocscore$gene_ids

#pdf('Differentiation.candidate.genes.pdf')
#for (numrow in seq(1,nrow(selected))) {
#  gene_id <- selected[numrow,'gene_ids']
#  gene_name <- selected[numrow,'gene_names']
#  boxplot(unlist(log2fpkm[gene_id,])~timepoints,main=paste(gene_id,gene_name))
#}
#plotSA(fit, main="Final model: Mean-variance trend", ylab = "Sqrt( standard deviation )") # Looks VERY bad
#qq(tt$P.Value) # Extreme inflation!
#boxplot(abs(tt[selected[,'gene_ids'],'logFC'])~colocscore[selected[,'gene_ids'],'colocScore'],ylab='|logFC|',xlab='colocScore',main='Gene DE scores per coloc category')
#boxplot(sign(tt[selected[,'gene_ids'],'logFC'])*-log10(tt[selected[,'gene_ids'],'P.Value'])~colocscore[selected[,'gene_ids'],'colocScore'],ylab='signed -log10P',xlab='colocScore',main='Gene DE scores per coloc category')
#dev.off()
