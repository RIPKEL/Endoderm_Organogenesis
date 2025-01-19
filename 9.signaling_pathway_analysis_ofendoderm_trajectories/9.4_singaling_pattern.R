#===============================================================================
#>>> 9.signaling pathways analysis 
#===============================================================================

#===============================================================================
#>>> 9.4 signaling pattern
#===============================================================================


#-- pattern of signaling
for(i.pathway in c("FRS-mediated FGFR1 signaling",
                   "PI-3K cascade:FGFR2",
                   "SHC-mediated cascade:FGFR3",
                   "Beta-catenin independent WNT signaling")){
  p = ggplot()+ geom_point(
    data = as.data.frame(cbind(src.27ss.integrated@reductions$umap@cell.embeddings,
                               seurat.plsda.pathway.endoderm.time.list$ss27$merge$Pathway@scale.data[
                                 i.pathway,])),
    mapping = aes(x=as.numeric(UMAP_1), y=as.numeric(UMAP_2), color=V3), size=3.5)+ 
    theme_void() +
    theme(# legend.position = "none",
      aspect.ratio = 1) +
    scale_color_viridis()
  
  png(filename = paste("update_result/Pattern.pathway.", gsub(" ","_",i.pathway), ".27SS.png", sep = ""),
      width = 1000,height = 1000,pointsize = 20)
  print(p)
  dev.off()  
}

#-- pathway associated genes：FGF, MAPK, WNT 
gi=read.csv("~/genome/mm10/mm10.gene.inf.merge.10X.v1.csv",stringsAsFactors = F)
gi$gene=gi$Symbol; rownames(gi)=gi$Symbol2
gi$Symbol = gsub("_","-",gi$Symbol)
expressed_ensembl = gi[match(rownames(
  src@assays$RNA@data[rowSums(src)>0,]), gi$Symbol),]$EnsemblGeneID
expressed_ensembl = MymoveStrange(expressed_ensembl)

marker.set.1 = intersect(reactome[reactome$Pathway%in%"FRS-mediated FGFR1 signaling",]$Ensembl.id, expressed_ensembl)
marker.set.2 = intersect(reactome[reactome$Pathway%in%"PI-3K cascade:FGFR2",]$Ensembl.id, expressed_ensembl)
marker.set.3 = intersect(reactome[reactome$Pathway%in%"SHC-mediated cascade:FGFR3",]$Ensembl.id, expressed_ensembl)

marker.set.4 = intersect(reactome[reactome$Pathway%in%"MAPK1/MAPK3 signaling",]$Ensembl.id, expressed_ensembl)
marker.set.5 = intersect(reactome[reactome$Pathway%in%"Beta-catenin independent WNT signaling",]$Ensembl.id, expressed_ensembl)

pdf("update_result/Pattern.MAPK.27ss.pdf", 7, 7)
for(feature in c(gi[gi$EnsemblGeneID%in%unique(c(marker.set.4, marker.set.5)),]$Symbol)){
  print(FeaturePlot(src.27ss.integrated, features = feature, reduction = "umap")+
          theme_void() +
          theme(legend.position = "none",
                aspect.ratio = 1))
}
dev.off()

pdf("update_result/Pattern.FGF.21ss.pdf", 7, 7)
for(feature in c(gi[gi$EnsemblGeneID%in%unique(c(marker.set.1, marker.set.2, marker.set.3)),]$Symbol)){
  print(FeaturePlot(src.21ss.integrated, features = feature, reduction = "umap")+
          theme_void() +
          theme(legend.position = "none",
                aspect.ratio = 1))
}
dev.off()

pdf("update_result/Pattern.MAPK.21ss.pdf", 7, 7)
for(feature in c(gi[gi$EnsemblGeneID%in%unique(c(marker.set.4)),]$Symbol)){
  print(FeaturePlot(src.21ss.integrated, features = feature, reduction = "umap")+
          theme_void() +
          theme(legend.position = "none",
                aspect.ratio = 1))
}
dev.off()

pdf("update_result/Pattern.WNT.21ss.pdf", 7, 7)
for(feature in c(gi[gi$EnsemblGeneID%in%unique(c(marker.set.5)),]$Symbol)){
  print(FeaturePlot(src.21ss.integrated, features = feature, reduction = "umap")+
          theme_void() +
          theme(legend.position = "none",
                aspect.ratio = 1))
}
dev.off()
#===============================================================================
