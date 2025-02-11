library(ggplot2)
library(ComplexUpset)
library(reshape2)
library(cowplot)

dat <- read.table('data/summary.per_colocalized_instance.November2023.full.txt',header = T)
loci_info <- read.table('data/AD_GWAS_loci.txt',sep='\t',header=T,skip = 1)

qtls_per_clumpm <- table(dat[,c('GWAS_clump','QTL_type')])
qtls_per_clumpm = qtls_per_clumpm >= 1
qtls_per_clumpm <- as.data.frame(qtls_per_clumpm)
#qtls_per_clump <- rbind(qtls_per_clump,c(FALSE,FALSE,TRUE,FALSE)) #Add fake clump to populate pQTL-only GWAS hits to align with other plots. Modify in Aff Des to make it zero.
#qtls_per_clump_upset <- upset(qtls_per_clump, colnames(qtls_per_clump), name='GWAS loci', width_ratio=0.1, keep_empty_groups=TRUE, min_size=0)


# Map novel loci and ancestry
#qtls_per_clumpm <- qtls_per_clump[-nrow(qtls_per_clump),] # remove the artificially added row
clump_indexs <- as.numeric(gsub('.chr.*','',gsub('clump','',rownames(qtls_per_clumpm)))) # clump index. 86/101 clumps have a QTL link
qtls_per_clumpm$clump <- rownames(qtls_per_clumpm)
rownames(qtls_per_clumpm) <- clump_indexs

qtls_per_clumpm <- cbind(qtls_per_clumpm,loci_info[as.character(clump_indexs),c('Known.Novel','ancestry_endpoint')])
colnames(qtls_per_clumpm)[6:7] <- c('Novelty','Ancestry')
qtls_per_clumpm$Ancestry <- gsub('CROSSANCESTRY','MULTI',qtls_per_clumpm$Ancestry)

pdf('fig4.a.pdf',width = 4.6, height = 3)
upset(data=qtls_per_clumpm, 
      intersect=colnames(qtls_per_clumpm)[1:4],
      name='GWAS loci', width_ratio=0.1, keep_empty_groups=TRUE, min_size=0, 
      base_annotations=list(
        'Intersection size'=intersection_size(
          counts=TRUE,
          mapping=aes(fill=Novelty, alpha=Ancestry, color='black')
          )
        + scale_alpha_manual(values = c(0.3,1)))
      ) #+ theme_bw(base_size = 8)
dev.off()

datsm <- subset(dat,!QTL_type%in%'mQTL') 
datsmeth <- subset(dat,QTL_type%in%'mQTL') 
methgenes <- unique(datsmeth$Gene)
methgenessplit <- unique(unlist(strsplit(datsmeth$Gene,';')))

qtls_per_gene <- table(datsm[,c('Gene','QTL_type')])
qtls_per_gene = qtls_per_gene >= 1
qtls_per_gene <- as.data.frame(qtls_per_gene)
qtls_per_gene <- cbind(qtls_per_gene,rownames(qtls_per_gene)%in%methgenessplit)
colnames(qtls_per_gene)[4] <- 'mQTL'
qtls_per_gene_upset <- upset(qtls_per_gene, colnames(qtls_per_gene), name='Genes', width_ratio=0.1, keep_empty_groups=TRUE, min_size=0)

pdf('fig4.b.pdf',width = 3.5, height = 3)
qtls_per_gene_upset
dev.off()

qtlsupport <- qtls_per_gene[,1]+qtls_per_gene[,3]
top_sup_genes <- rownames(subset(qtls_per_gene,pQTL & eQTL & mQTL & sQTL ))
qtls_per_top_sup_genes <- table(datsm[,c('Gene','QTL_type')])[top_sup_genes,]
#qtls_per_top_sup_genes <- as.data.frame(qtls_per_top_sup_genes)
methgenessplit <- unlist(strsplit(datsmeth$Gene,';'))
qtls_per_top_sup_genes <- cbind(qtls_per_top_sup_genes,table(methgenessplit[methgenessplit%in%top_sup_genes])[rownames(qtls_per_top_sup_genes)])
colnames(qtls_per_top_sup_genes)[4] <- 'mQTL'
qtls_per_top_sup_genes <- as.data.frame(qtls_per_top_sup_genes)
qtls_per_top_sup_genes$Total <- apply(qtls_per_top_sup_genes,1,sum)
qtls_per_top_sup_genesm <- melt(cbind(rownames(qtls_per_top_sup_genes),qtls_per_top_sup_genes))
colnames(qtls_per_top_sup_genesm) <- c('Gene','QTL','# QTL\nendpoints')
qtls_per_top_sup_genesm$Gene <- factor(qtls_per_top_sup_genesm$Gene,levels=rev(unique(qtls_per_top_sup_genesm$Gene)))

qtlsupport_bubble <- ggplot(subset(qtls_per_top_sup_genesm,QTL%in%c('eQTL','pQTL','sQTL','mQTL')), aes(y=Gene, size=`# QTL\nendpoints`, x = QTL)) + 
  geom_point(aes(fill = `# QTL\nendpoints`),colour='black',pch=21) + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(angle=30, vjust = 0.5, hjust=1)) +
  scale_size_continuous(name = "# QTL\nendpoints", 
                        range = c(0.5, 7),  
                        limits = c(0, 100), 
                        breaks = c(5,20,40,80)) + 
  scale_fill_binned(type = "viridis", breaks = c(5, 20, 30, 60)) + 
  #theme(legend.position="top") +
  guides(fill="none") + guides(size="none")

qtlsupport_bar <- ggplot(data=subset(qtls_per_top_sup_genesm,QTL%in%'Total'), aes(y=`# QTL\nendpoints`, x = Gene)) +
  geom_bar(position = 'dodge', stat='identity', colour = 'black', aes(fill=`# QTL\nendpoints`)) +
  geom_text(aes(label=`# QTL\nendpoints`), position=position_dodge(width=0.9), vjust = 0.5, hjust=-0.75) +
  coord_flip() + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(angle=45, vjust = 0.5, hjust=1.2),
        axis.text.y=element_blank()) +
  #theme(legend.position="top") +
  guides(fill="none") +
  scale_fill_binned(name = "# QTL endpoints", type = "viridis", breaks = c(5, 20, 30, 60),  limits = c(0, 120))

pdf('qtl_support_inlet.2.pdf',width = 2.5, height = 1.25)
plot_grid(qtlsupport_bubble, qtlsupport_bar,rel_widths = c(2,1))
dev.off()

# Create legend manually in AffinityDesigner




  