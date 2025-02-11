#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE);

library(stringr)

## set the phenotype and chromosome from the command line or hard code it (for testing)
## eg for debugging: pheno = 'AtopicDerm_GWAS_meta_FE_hg38'; chr = 10; window = 2; pvalue = 5e-08
window <- args[1]; # e.g. 2
thres <- args[2];  # e.g. 5e-08
pheno <- args[3];  # e.g. 'AD_GWAS_meta_FE_DF9_hg38_forMeri'
chr <- args[4];  # e.g. 10

## read in the clumps file
clumps <- read.table(paste0('/ui/abv/olivamx2/scratch-global/AD/',pheno,'.chr',chr,'.UKBunrelatedLD.pvalue',thres,'.window',window,'Mb.clumped'),header=T)
clumps$SNP <- as.character(clumps$SNP)

## read in the genomic positions in the bim file
bim_data <- paste0('/ui/abv/olivamx2/scratch-global/UKB/ukb_imp_chr',chr,'_v3.15000_ukb_unrelated_samples.filtered_variants.info_score_dot3.maf_dotzerozeroone.bim')
bim <- read.table(bim_data, header=F)[,c(2,4)]
colnames(bim) <- c('SNP','BP')
bim$SNP <- as.character(bim$SNP)

## sort the data so that it is in genomic order
clumps <- clumps[order(clumps$BP),]
bim <- bim[order(bim$BP),]

### 1) Identify the genomic regions (loci) from the ld clumping
final_result = as.data.frame(rep(NA,4))
for (i in seq(1,nrow(clumps))) {
  ## annotate the index snp
  index_snp <- clumps[i,'SNP'] 
  posindex <- bim$BP[bim$SNP%in%index_snp]
  ## identify the first and last snp in the clump of this index snp, including the index snp itself
  snplist <- c(unlist(strsplit(as.character(clumps[i,'SP2']),'\\(1\\),*')),index_snp)

  ## positions of the first and last snps in the window
  snp1 <- subset(bim,SNP%in%snplist)[1,'SNP'] # start snp

  if ( snp1%in%'NONE' ) {
   ## if there are no SNPs in the clump, the clump is the snp itself (plus buffer, see below)
   snp1 = index_snp
   pos1 = posindex
   snp2 = index_snp
   pos2 = pos1
  } else {
  ## else, set the start and end of the cluster to the most upstream and downstream snp, respectively
    bim_sub <- subset(bim,SNP%in%snplist)
    snp1 = bim_sub[1,'SNP']
    pos1 = bim_sub[1,'BP']
    snp2 = bim_sub[nrow(bim_sub),'SNP']
    pos2 = bim_sub[nrow(bim_sub),'BP']
}
  ## add a tiny bit of buffer to around the positions
  pos1 = pos1 - 1000
  pos2 = pos2 + 1000

  ## store the results
  result = c(chr, index_snp, pos1, pos2)
  final_result <- cbind(final_result,result)
}
final_result <- as.data.frame(t(final_result))[-1,]
write.table(file=paste0('/ui/abv/olivamx2/scratch-global/AD/',pheno,'.chr',chr,'.UKBunrelatedLD.pvalue',thres,'.window',window,'Mb.clumped.raw_windows.txt'),final_result,sep=' ',quote=F,row.names=F)

### 2) Find the overlapping windows and merge them
colnames(final_result) <- c('chr','SNP','start','end')
final_result$start <- as.numeric(as.character(final_result$start));
final_result$end <- as.numeric(as.character(final_result$end));
final_result$SNP <- as.character(final_result$SNP)

## sort the intervals that are given to the function by their lower bound
final_result <- final_result[order(final_result$start),]

start <- final_result[1,'start']
end <- final_result[1,'end']

merged_final_result = final_result[1,]
for (i in seq(1,nrow(final_result))) {
  if(final_result[i,'start'] > start & final_result[i,'start'] < end) {
    merged_final_result[nrow(merged_final_result),'end'] <- final_result[i,'end']
    merged_final_result[nrow(merged_final_result),'SNP'] <- paste0(merged_final_result[nrow(merged_final_result),'SNP'],"|",final_result[i,'SNP'])
  } else {
    merged_final_result <- rbind(merged_final_result,final_result[i,])
  }
  start <- merged_final_result[nrow(merged_final_result),'start']
  end <- merged_final_result[nrow(merged_final_result),'end']
}
merged_final_result <- merged_final_result[-1,]
write.table(file=paste0('/ui/abv/olivamx2/scratch-global/AD/',pheno,'.chr',chr,'.UKBunrelatedLD.pvalue',thres,'.window',window,'Mb.clumped.windows.txt'),merged_final_result,sep=' ',quote=F,row.names=F)

## use the genomic windows to find the list of snps in each locus. write out the list to be used in gcta
for (i in seq(1,nrow(merged_final_result))) {
  locus_snps <- bim[bim$BP >= merged_final_result[i,'start'] & bim$BP <= merged_final_result[i,'end'],'SNP'];
  write.table(file=paste0('/ui/abv/olivamx2/scratch-global/AD/',pheno,'.chr',chr,'.UKBunrelatedLD.pvalue',thres,'.window',window,'Mb.clumped.loci',i,'.txt'),as.data.frame(locus_snps),quote=F,row.names=F)
}

quit(save='no')
