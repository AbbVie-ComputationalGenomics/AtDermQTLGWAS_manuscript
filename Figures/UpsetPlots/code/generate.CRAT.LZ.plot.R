f1 <- read.table('data/clump57.chr9.128720519.129216000.rs3124499.9.129106712.pQTL.ARIC_EA.SeqId_12637_7_ENSG00000095321_P43155_CRAT.CRAT.A.G.neg.rs1107329.CROSSANCESTRY.A.G.neg.txt',header=T)
f2 <- read.table('data/tmp.GCST90264396_buildGRCh37.clump57.chr9.128720519.129216000.txt',header=T)
f3 <- merge(f2[,c('marker_id','p_value')],f1,by.x=1,by.y=2)
write.table(file='data/LZ.data.txt',f3,sep='\t',quote=F)

library(locuscomparer)
library(cowplot)

pdf('lz.CRAT.pdf',width=3.5,height=6.5);
ll.top <- locuscompare('data/LZ.data.txt','data/LZ.data.txt',marker_col1='marker_id',pval_col1='eQTLp',marker_col2 = 'marker_id', pval_col2 = 'GWASp',snp = 'rs1107329',genome='hg38',combine=F)
ll.bot <- locuscompare('data/LZ.data.txt','data/LZ.data.txt',marker_col1='marker_id',pval_col1='eQTLp',marker_col2 = 'marker_id', pval_col2 = 'p_value',snp = 'rs1107329',genome='hg38',combine=F)
plot_grid(as_grob(ll.top[2]$locuszoom1),as_grob(ll.top[3]$locuszoom2),as_grob(ll.bot[3]$locuszoom2),nrow = 3)# 250 x 350
dev.off()  
