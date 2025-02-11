library(Seurat)

obj <- readRDS('data/Umich_epiderm_atlas_seurat_v3.RDS') # read Seurat object

cs <- read.table('../../data/other/coloc_score.tsv',header=T)
cs <- cs[order(cs$Prioritized,decreasing = T),]

plotsl <- list()
for (genen in seq(1,22)) {
  plotsl[[genen]] <- FeaturePlot(obj, features = cs[genen,1],raster=T) + 
    theme(plot.title = element_text(size=8), 
          axis.title = element_blank(),
          legend.text = element_text(size = 8),
          axis.text = element_text(size = 8)) + 
    scale_color_continuous(low = "lightgrey", high = "purple") +
    theme(legend.position = "none")
}
plotsl[[23]] <- FeaturePlot(obj, features = cs[22,1],raster=T) + 
  theme(plot.title = element_text(size=8), 
        axis.title = element_blank(),
        legend.text = element_text(size = 8),
        axis.text = element_text(size = 8)) + 
  scale_color_continuous(limits=c(0, 2), breaks=seq(0,2,by=1), low = "lightgrey", high = "purple") 
pdf('ExpUMAP.pdf',width = 8, height = 6.5)
plot_grid(plotsl[[1]],plotsl[[2]],plotsl[[3]],plotsl[[4]],plotsl[[5]],plotsl[[6]],plotsl[[7]],plotsl[[8]],plotsl[[9]],plotsl[[10]],plotsl[[11]],plotsl[[12]],plotsl[[13]],plotsl[[14]],plotsl[[15]],plotsl[[16]],plotsl[[17]],plotsl[[18]],plotsl[[19]],plotsl[[20]],plotsl[[21]],plotsl[[22]], plotsl[[23]], ncol=5,nrow=5, scale=1.2)
dev.off()