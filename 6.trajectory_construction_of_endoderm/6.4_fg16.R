#----------------------------------------------------------------------------
#>   Part1: FG1 AND FG6
#----------------------------------------------------------------------------

src.endoderm.fg16.ext.v1.1 = merge(src.endoderm.fg1.ext_v1.1,
                                   src.endoderm.fg6.ext_v1.1)

src.endoderm.fg16.ext.v1.1.selectgene = select_gene_fg16
pdf("figure.v08.07/try.pdf")
#-----------------
fg16.ext.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.fg16.ext_v1.1@assays$RNA@data[src.endoderm.fg16.ext.v1.1.selectgene,]),
            type = "row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.fg16.ext_v1.1$Time,colors.time),
              MyName2Col(src.endoderm.fg16.ext_v1.1$cluster.extract.v1.1,cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.fg16.ext_v1.1$cluster.v06.26.re,cluster.endoderm.color.v5)
            ),
            ColSideColorsSize = 1.5,
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            #Rowv = "none",
            return.tree = "row",
            graph = T)
#-------------------
dev.off()

fg16.ext.re.rowtree_1 = as.dendrogram(fg16.ext.re.rowtree)

src.endoderm.fg16.ext_v1.1.gene.fin = setdiff(
  src.endoderm.fg16.ext_v1.1.selectgene,
  c(labels(fg16.ext.re.rowtree_1[[2]][[1]]),
    labels(fg16.ext.re.rowtree_1[[2]][[2]][[2]][[1]])))

src.endoderm.fg16.ext_v1.1 = 
  seurat_mnn(seurat_name = "src.endoderm.fg16.ext_v1.1",
             seurat.selectgene = src.endoderm.fg16.ext_v1.1.gene.fin)
src.endoderm.fg16.ext_v1.1 = 
  seurat_int(seurat_name = "src.endoderm.fg16.ext_v1.1",
             seurat.selectgene = src.endoderm.fg16.ext_v1.1.gene.fin)

src.endoderm.fg16.ext_v1.1 = FindNeighbors(src.endoderm.fg16.ext_v1.1, dims=1:30)
src.endoderm.fg16.ext_v1.1 = FindClusters(src.endoderm.fg16.ext_v1.1, resolution = 2)

DimPlot(src.endoderm.fg16.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.2', reduction = "umap_integrated")
DimPlot(src.endoderm.fg16.ext_v1.1, group.by = 'Time', 
        reduction = "umap_mnn", cols = colors.time)
DimPlot(src.endoderm.fg16.ext_v1.1, group.by = 'Time', 
        reduction = "umap_integrated", cols = colors.time)
DimPlot(src.endoderm.fg16.ext_v1.1, group.by = 'cluster.extract.v1.1', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.fg16.ext_v1.1, group.by = 'cluster.v06.26.re', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)

#src.endoderm.fg16.ext_v1.1$cluster.v06.26 = src.endoderm.fg16.ext_v1.1$cluster.v06.26.re
src.endoderm.fg16.ext_v1.1$cluster.v06.26.re = src.endoderm.fg16.ext_v1.1$cluster.v06.26
src.endoderm.fg16.ext_v1.1@meta.data[
  src.endoderm.fg16.ext_v1.1$cluster.v06.26%in%c("FG.1-Pharynx.organ.2"),]$cluster.v06.26.re = "FG.1"
src.endoderm.fg16.ext_v1.1@meta.data[
  src.endoderm.fg16.ext_v1.1$cluster.v06.26%in%c("FG.1-Esophagus"),]$cluster.v06.26.re = "FG.6"
src.endoderm.fg16.ext_v1.1@meta.data[
  src.endoderm.fg16.ext_v1.1$RNA_snn_res.2%in%c(15),]$cluster.v06.26.re = "Pharynx.organ.2"
src.endoderm.fg16.ext_v1.1@meta.data[
  src.endoderm.fg16.ext_v1.1$RNA_snn_res.2%in%c(2,12),]$cluster.v06.26.re = "FG.1-Pharynx.organ.1"

a = FNN::knn(
  src.endoderm.fg16.ext_v1.1@reductions$umap_integrated@cell.embeddings[
    rownames(src.endoderm.fg16.ext_v1.1@meta.data[!src.endoderm.fg16.ext_v1.1$cluster.v06.26.re%in%NA,]),],
  src.endoderm.fg16.ext_v1.1@reductions$umap_integrated@cell.embeddings[
    rownames(src.endoderm.fg16.ext_v1.1@meta.data[src.endoderm.fg16.ext_v1.1$cluster.v06.26.re%in%NA,]),],
  src.endoderm.fg16.ext_v1.1@meta.data[!src.endoderm.fg16.ext_v1.1$cluster.v06.26.re%in%NA,]$cluster.v06.26.re, k = 10)
src.endoderm.fg16.ext_v1.1@meta.data[src.endoderm.fg16.ext_v1.1$cluster.v06.26.re%in%NA,]$cluster.v06.26.re = as.character(a)



#========================================================================================
src.endoderm.fg16.ext_v1.1.re = src.endoderm.fg16.ext_v1.1[,!(
  src.endoderm.fg16.ext_v1.1$Time%in%"ss15"&
    src.endoderm.fg16.ext_v1.1$cluster.v06.26.re%in%"FG.1-Pharynx.organ.1")]
src.endoderm.fg16.ext_v1.1.re = ScaleData(src.endoderm.fg16.ext_v1.1.re, split.by = "batch_phase",
                                          features = rownames(src.endoderm.fg16.ext_v1.1.re))
src.endoderm.fg16.ext_v1.1.re = 
  seurat_mnn(seurat_name = "src.endoderm.fg16.ext_v1.1.re",
             seurat.selectgene = src.endoderm.fg16.ext_v1.1.gene.fin)
src.endoderm.fg16.ext_v1.1.re = 
  seurat_int(seurat_name = "src.endoderm.fg16.ext_v1.1.re",
             seurat.selectgene = src.endoderm.fg16.ext_v1.1.gene.fin)

# FDL
src.endoderm.fg16.ext_v1.1.re = 
  seurat_fdl("src.endoderm.fg16.ext_v1.1.re", "pca","RNA")
src.endoderm.fg16.ext_v1.1.re = 
  seurat_fdl("src.endoderm.fg16.ext_v1.1.re", "mnn","RNA")
src.endoderm.fg16.ext_v1.1.re = 
  seurat_fdl("src.endoderm.fg16.ext_v1.1.re", "pca_integrated","integrated")

src.endoderm.fg16.ext_v1.1.re = FindNeighbors(src.endoderm.fg16.ext_v1.1.re, dims=1:30)
src.endoderm.fg16.ext_v1.1.re = FindClusters(src.endoderm.fg16.ext_v1.1.re, resolution = 2)

DimPlot(src.endoderm.fg16.ext_v1.1.re, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.2', reduction = "umap_integrated")


# Correct cell type
#-----------------------
src.endoderm.fg16.ext_v1.1.re$cluster.v06.26.re = 
  src.endoderm.fg16.ext_v1.1@meta.data[colnames(src.endoderm.fg16.ext_v1.1.re),]$cluster.v06.26.re
src.endoderm.fg16.ext_v1.1.re@meta.data[
  src.endoderm.fg16.ext_v1.1.re$cluster.v06.26.re%in%"FG.1-Pharynx.organ.1",]$cluster.v06.26.re = "FG.1"

src.endoderm.fg16.ext_v1.1.re@meta.data[src.endoderm.fg16.ext_v1.1.re$cluster.v06.26.re%in%"Esophagus"&
                                          src.endoderm.fg16.ext_v1.1.re$Time%in%'ss9',]$cluster.v06.26.re = "FG.6"
src.endoderm.fg16.ext_v1.1.re@meta.data[src.endoderm.fg16.ext_v1.1.re$cluster.v06.26.re%in%"Pharynx.organ.2"&
                                          src.endoderm.fg16.ext_v1.1.re$Time%in%'ss9',]$cluster.v06.26.re = "FG.1"
src.endoderm.fg16.ext_v1.1.re@meta.data[src.endoderm.fg16.ext_v1.1.re$cluster.v06.26.re%in%"FG.6"&
                                          src.endoderm.fg16.ext_v1.1.re$Time%in%c('ss27',"ss24"),]$cluster.v06.26.re = "Esophagus"
src.endoderm.fg16.ext_v1.1.re@meta.data[src.endoderm.fg16.ext_v1.1.re$cluster.v06.26.re%in%"FG.1"&
                                          src.endoderm.fg16.ext_v1.1.re$Time%in%c('ss27',"ss24"),]$cluster.v06.26.re = "Pharynx.organ.2"


# Raw
src.endoderm.fg16.ext_v1.1@meta.data[
  paste("ss15_", rownames(src.15ss.integrated.merge@meta.data[
    src.15ss.integrated.merge$cluster_snn_res.2%in%44,]),sep = ""),]$cluster.extract.v1.1 = "FG.1"
src.endoderm.fg16.ext_v1.1@meta.data[
  paste("ss15_", rownames(src.15ss.integrated.merge@meta.data[
    src.15ss.integrated.merge$cluster_snn_res.2%in%37,]),sep = ""),]$cluster.extract.v1.1 = "FG.6"
# Re
src.endoderm.fg16.ext_v1.1.re@meta.data[
  paste("ss15_", rownames(src.15ss.integrated.merge@meta.data[
    src.15ss.integrated.merge$cluster_snn_res.2%in%44,]),sep = ""),]$cluster.extract.v1.1 = "FG.1"
src.endoderm.fg16.ext_v1.1.re@meta.data[
  paste("ss15_", rownames(src.15ss.integrated.merge@meta.data[
    src.15ss.integrated.merge$cluster_snn_res.2%in%37,]),sep = ""),]$cluster.extract.v1.1 = "FG.6"

# src.endoderm.fg16.ext_v1.1.re@meta.data = 
#   src.endoderm.fg16.ext_v1.1.re@meta.data[!src.endoderm.fg16.ext_v1.1.re@meta.data$Time%in%NA,]
#-----------------------


DimPlot(src.endoderm.fg16.ext_v1.1.re,group.by = 'cluster.v06.26.re', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5,
        pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
DimPlot(src.endoderm.fg16.ext_v1.1.re,group.by = 'cluster.extract.v1.1', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5,
        pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
DimPlot(src.endoderm.fg16.ext_v1.1.re,group.by = 'Time', 
        reduction = "umap_integrated", cols = colors.time,
        pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)


# Graph
#--------------------------
graph.fg16.ext_v1.1.re = 
  seurat_GCN(seurat_name = "src.endoderm.fg16.ext_v1.1.re",
             src.endoderm.fg16.ext_v1.1.gene.fin,
             reduction = "pca_integrated", assay = "integrated")
graph.fg16_cluster = cluster_walktrap(graph.fg16.ext_v1.1.re, steps = 1)

pdf("figure.v08.07/organ_development_re/fg16_summary_graph.re.pdf",12,7)
# SNN
for(color in c("cluster.extract.v1.1","cluster.v06.26.re","Time")){
  for(edge.color in c(NA,"lightgray")){
    for(edge.color in c(NA,"lightgray")){
    set.seed(1)
    plot.igraph(graph.fg16.ext_v1.1.re,
                layout = layout_with_fr(graph.fg16.ext_v1.1.re),
                edge.color = edge.color,
                vertex.size=2.5, vertex.label=NA,
                vertex.label.color = "black", 
                vertex.frame.color = NA, vertex.frame.width = 0.5,
                #vertex.color =  colors.time[src.endoderm.fg16.ext_v1.1.re[["Time"]][graph.fg16_cluster$names,]],
                #vertex.color =  cluster.endoderm.color.v5[src.endoderm.fg16.ext_v1.1.re[["cluster.v06.26.re"]][graph.fg16_cluster$names,]],
                vertex.color =  colors_bg[src.endoderm.fg16.ext_v1.1.re[[color]][graph.fg16_cluster$names,]])
}}}
# Reduction
for(red in c("umap_integrated","umap_mnn",
             "fdl_pca",'fdl_mnn',"fdl_pca_integrated")){
  p1 = DimPlot(src.endoderm.fg16.ext_v1.1.re,group.by = 'cluster.v06.26.re', 
          reduction = red, cols = cluster.endoderm.color.v5,
          pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p2 = DimPlot(src.endoderm.fg16.ext_v1.1.re,group.by = 'cluster.extract.v1.1', 
          reduction = red, cols = cluster.endoderm.color.v5,
          pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p3 = DimPlot(src.endoderm.fg16.ext_v1.1.re,group.by = 'Time', 
          reduction = red , cols = colors.time,
          pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  print(p1); print(p2); print(p3)
}
dev.off()


pdf("figure.v08.07/organ_development_re/fg16_summary_proportion.re.pdf",12,7)
#-----------------------------
cell_type = c("FG.1","FG.6","FG.1-Pharynx.organ.1","Pharynx.organ.2","Esophagus")
meta_fg16 = src.endoderm.fg16.ext_v1.1.re@meta.data
meta_fg16$Time = factor(meta_fg16$Time, levels = names(colors.time))
meta_fg16$cluster.extract.v1.1 = factor(meta_fg16$cluster.extract.v1.1, levels = c("FG.1","FG.6"))
meta_fg16$cluster.v06.26.re = factor(meta_fg16$cluster.v06.26.re, levels = rev(cell_type))

meta_fg16_B0 = meta_fg16
meta_fg16_B1 = meta_fg16[meta_fg16$batch%in%1&!meta_fg16$Time%in%"ss9",]
meta_fg16_B2 = meta_fg16[meta_fg16$batch%in%2,]

meta_fg16_BO.fg1 = meta_fg16[meta_fg16$cluster.extract.v1.1%in%"FG.1",]
meta_fg16_B1.fg1 = meta_fg16[meta_fg16$batch%in%1&!meta_fg16$Time%in%"ss9"&
                               meta_fg16$cluster.extract.v1.1%in%"FG.1",]
meta_fg16_B2.fg1 = meta_fg16[meta_fg16$batch%in%2&
                               meta_fg16$cluster.extract.v1.1%in%"FG.1",]

meta_fg16_BO.fg6 = meta_fg16[meta_fg16$cluster.extract.v1.1%in%"FG.6",]
meta_fg16_B1.fg6 = meta_fg16[meta_fg16$batch%in%1&!meta_fg16$Time%in%"ss9"&
                               meta_fg16$cluster.extract.v1.1%in%"FG.6",]
meta_fg16_B2.fg6 = meta_fg16[meta_fg16$batch%in%2&
                               meta_fg16$cluster.extract.v1.1%in%"FG.6",]

for(data in c("meta_fg16_B0","meta_fg16_B1","meta_fg16_B2",
              "meta_fg16_BO.fg1","meta_fg16_B1.fg1","meta_fg16_B2.fg1",
              "meta_fg16_BO.fg6","meta_fg16_B1.fg6","meta_fg16_B2.fg6")){
  p = ggplot() +
    scale_fill_manual(values = cluster.endoderm.color.v5)+
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle=90, size = 13),
          axis.text.y = element_text(size = 13),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          #legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          plot.title = element_blank(),
          aspect.ratio=1/0.5) +
    guides(color = guide_legend(override.aes = list(size = 5))) 
  
  print(p + 
          geom_bar(data = get(data), 
                   mapping = aes(x = Time,
                                 group = cluster.extract.v1.1, fill = cluster.extract.v1.1),
                   stat = "count", position = 'fill') + 
          xlab("Time") + ylab("Cluster proportion of
 tracing code")+ ggtitle(data))
  
  print(p + 
          geom_bar(data = get(data), 
                   mapping = aes(x = Time,
                                 group = cluster.v06.26.re, fill = cluster.v06.26.re),
                   stat = "count", position = 'fill') + 
          xlab("Time") + ylab("Cluster proportion of
 cell type")+ ggtitle(data))
  
}
#-----------------------------
dev.off()

save(src.endoderm.fg16.ext_v1.1.re,
     file = "figure.v08.07/organ_development_re/src.endoderm.fg16.ext_v1.1.re.Rdata")

#--------------------------
#-- Marker all for FG  ----

src.endoderm.fg16.ext_v1.1.re = 
  SetIdent(src.endoderm.fg16.ext_v1.1.re,
           value = src.endoderm.fg16.ext_v1.1.re$cluster.v06.26.re)
marker.fg16.ext_v1.1.re = 
  FindAllMarkers(src.endoderm.fg16.ext_v1.1.re,assay = "RNA")

src.endoderm.fg16.ext_v1.1.re.markergene = marker.fg16.ext_v1.1.re[
  marker.fg16.ext_v1.1.re$avg_log2FC>0.25&
    marker.fg16.ext_v1.1.re$p_val_adj<0.1,"gene"]



fg16.umap.embedding = src.endoderm.fg16.ext_v1.1.re[["fdl_pca_integrated"]]@cell.embeddings[,1:3]
fg16.umap.embedding = src.endoderm.fg16.ext_v1.1.re[["pca_integrated"]]@cell.embeddings[,1:4]

# Cell order :: Princurve
#---------------------------
#  FG.1
cellorder.fg16_fg1 = 
  rownames(src.endoderm.fg16.ext_v1.1.re@meta.data[
    src.endoderm.fg16.ext_v1.1.re$cluster.v06.26.re%in%"FG.1",])
princurve.fg16_fg1 = princurve::principal_curve(
  x = fg16.umap.embedding[cellorder.fg16_fg1,],  smoother = "smooth.spline")
cellorder.fg16_fg1 =  names(
  princurve.fg16_fg1$lambda[cellorder.fg16_fg1][order(princurve.fg16_fg1$lambda[cellorder.fg16_fg1])])

# FG.1 - Pha.2
cellorder.fg16_pha2 = 
  rownames(src.endoderm.fg16.ext_v1.1.re@meta.data[
    src.endoderm.fg16.ext_v1.1.re$cluster.v06.26.re%in%"Pharynx.organ.2",])
princurve.fg16_pha2 = princurve::principal_curve(
  x = fg16.umap.embedding[cellorder.fg16_pha2,],  smoother = "smooth.spline")
cellorder.fg16_pha2 =  names(
  princurve.fg16_pha2$lambda[cellorder.fg16_pha2][order(princurve.fg16_pha2$lambda[cellorder.fg16_pha2])])

# FG.6
cellorder.fg16_fg6 = 
  rownames(src.endoderm.fg16.ext_v1.1.re@meta.data[
    src.endoderm.fg16.ext_v1.1.re$cluster.v06.26.re%in%"FG.6",])
princurve.fg16_fg6 = princurve::principal_curve(
  x = fg16.umap.embedding[cellorder.fg16_fg6,],  smoother = "smooth.spline")
cellorder.fg16_fg6 =  names(
  princurve.fg16_fg6$lambda[cellorder.fg16_fg6][order(princurve.fg16_fg6$lambda[cellorder.fg16_fg6])])

# FG.6 - Eso
cellorder.fg16_eso = 
  rownames(src.endoderm.fg16.ext_v1.1.re@meta.data[
    src.endoderm.fg16.ext_v1.1.re$cluster.v06.26.re%in%"Esophagus",])
princurve.fg16_eso = princurve::principal_curve(
  x = fg16.umap.embedding[cellorder.fg16_eso,],  smoother = "smooth.spline")
cellorder.fg16_eso =  names(
  princurve.fg16_eso$lambda[cellorder.fg16_eso][order(princurve.fg16_eso$lambda[cellorder.fg16_eso])])
#---------------------------

plot(c(1:length(cellorder.fg16_fg1)),
     c(1:length(cellorder.fg16_fg1)),
     col = colors.time[src.endoderm.fg16.ext_v1.1.re[,cellorder.fg16_fg1]$Time])
plot(c(1:length(cellorder.fg16_fg6)),
     c(1:length(cellorder.fg16_fg6)),
     col = colors.time[src.endoderm.fg16.ext_v1.1.re[,cellorder.fg16_fg6]$Time])
plot(c(1:length(cellorder.fg16_eso)),
     c(1:length(cellorder.fg16_eso)),
     col = colors.time[src.endoderm.fg16.ext_v1.1.re[,cellorder.fg16_eso]$Time])
plot(c(1:length(cellorder.fg16_pha2)),
     c(1:length(cellorder.fg16_pha2)),
     col = colors.time[src.endoderm.fg16.ext_v1.1.re[,cellorder.fg16_pha2]$Time])

cellorder.fg16_fg6_FPI = cellorder.fg16_fg6 # 1:3
cellorder.fg16_eso_FPI = cellorder.fg16_eso # 1:3
cellorder.fg16_fg1_PI = cellorder.fg16_fg1
cellorder.fg16_pha2_FPI = cellorder.fg16_pha2 # 1:3

cellorder.fg16 = c(cellorder.fg16_fg6,cellorder.fg16_eso)
#---------------------------------------------------------------------------------



#----------------------------------------------------------------------------
#>   Part2: FG1 AND FG6, seperated
#----------------------------------------------------------------------------
src.endoderm.fg1.ext_v1.1.re = FindVariableFeatures(src.endoderm.fg1.ext_v1.1.re, nfeatures = 2000)
src.endoderm.fg6.ext_v1.1.re = FindVariableFeatures(src.endoderm.fg6.ext_v1.1.re, nfeatures = 2000)

src.endoderm.fg1.ext_v1.1.re.selectgene = 
  Myfilter(as.matrix(src.endoderm.fg1.ext_v1.1.re@assays$RNA@data),
           gene = src.endoderm.fg1.ext_v1.1.re@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.endoderm.fg6.ext_v1.1.re.selectgene = 
  Myfilter(as.matrix(src.endoderm.fg6.ext_v1.1.re@assays$RNA@data),
           gene = src.endoderm.fg6.ext_v1.1.re@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)

#----=================-----
# FG.6 Filter gene
#----=================-----
src.endoderm.fg6.ext_v1.1.re = SetIdent(src.endoderm.fg6.ext_v1.1.re,
                                        value = src.endoderm.fg6.ext_v1.1.re$cluster.v06.26.re)
src.endoderm.fg6.ext_v1.1.re.markergene = 
  FindAllMarkers(src.endoderm.fg6.ext_v1.1.re, assay = "RNA")
src.endoderm.fg6.ext_v1.1.re.filtergene = 
  src.endoderm.fg6.ext_v1.1.re.markergene[
    src.endoderm.fg6.ext_v1.1.re.markergene$avg_log2FC>0.25&
      src.endoderm.fg6.ext_v1.1.re.markergene$p_val_adj<0.1,"gene"]

cellorder.fg16 = c(cellorder.fg16_fg6_FPI, rev(cellorder.fg16_eso_FPI))
selectgene = src.endoderm.fg6.ext_v1.1.re.filtergene_fin
cellorder.fg16 = intersect(cellorder.fg16,
  colnames(src.endoderm.fg6.ext_v1.1.re))

pdf("figure.v08.07/organ_development_re/fg6_heatmap_marker.re.pdf",6,9)
#-------------------
fg16.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.fg16.ext_v1.1.re@assays$RNA@data[selectgene, cellorder.fg16]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.fg16.ext_v1.1.re$Time[cellorder.fg16], colors.time),
              MyName2Col(src.endoderm.fg16.ext_v1.1.re$cluster.v06.26.re[cellorder.fg16], cluster.endoderm.color.v5),
              # MyName2Col(src.endoderm.fg16.ext_v1.1.re$cluster.extract.v1.1[cellorder.fg16], cluster.endoderm.color.v5)
              MyName2Col(rep("FG.6",length(cellorder.fg16)), cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            #return.tree = "row",
            labCol="none",
            #column_split = data.frame(src.endoderm.fg1.ext.re@meta.data[cell_fg1_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

# fg16.ext_v1.1.re.rowtree.1 = as.dendrogram(fg16.ext_v1.1.re.rowtree)
src.endoderm.fg6.ext_v1.1.re.filtergene_fin = c(
  labels(fg16.ext_v1.1.re.rowtree.1[[1]][[1]][[1]]),
  labels(fg16.ext_v1.1.re.rowtree.1[[1]][[2]]),
  labels(fg16.ext_v1.1.re.rowtree.1[[1]][[1]][[2]]),
  labels(fg16.ext_v1.1.re.rowtree.1[[2]][[1]]),
  labels(fg16.ext_v1.1.re.rowtree.1[[2]][[2]]))
names(src.endoderm.fg6.ext_v1.1.re.filtergene_fin) = c(
  rep(4, length(labels(fg16.ext_v1.1.re.rowtree.1[[1]]))),
  rep(7, length(labels(fg16.ext_v1.1.re.rowtree.1[[2]])))
)

save(src.endoderm.fg6.ext_v1.1.re.filtergene_fin,
     src.endoderm.fg6.ext_v1.1.re.selectgene,
     cellorder.fg16_fg6_FPI, cellorder.fg16_eso_FPI,
     file= "figure.v08.07/organ_development_re/fg6_heatmap_parameter.Rdata")
save(src.endoderm.fg6.ext_v1.1,
     src.endoderm.fg6.ext_v1.1.re,
     file = "figure.v08.07/organ_development_re/src.endoderm.fg6.ext_v1.1.re.Rdata")
#----=================-----



#----=================-----
# FG.1 Filter gene
#----=================-----

src.endoderm.fg1.ext_v1.1.re = SetIdent(src.endoderm.fg1.ext_v1.1.re,
                                        value = src.endoderm.fg1.ext_v1.1.re$cluster.v06.26.re)
src.endoderm.fg1.ext_v1.1.re.markergene = 
  FindAllMarkers(src.endoderm.fg1.ext_v1.1.re, assay = "RNA")
src.endoderm.fg1.ext_v1.1.re.filtergene = 
  src.endoderm.fg1.ext_v1.1.re.markergene[
    src.endoderm.fg1.ext_v1.1.re.markergene$avg_log2FC>0.25&
      src.endoderm.fg1.ext_v1.1.re.markergene$p_val_adj<0.1,"gene"]

# selectgene = unique(c(
#   src.endoderm.fg1.ext_v1.1.re.filtergene,
#   src.endoderm.fg1.ext_v1.1.re.selectgene))

cellorder.fg16 = c(cellorder.fg16_fg1_PI, 
                   cellorder.fg16_pha2_FPI)
selectgene = src.endoderm.fg1.ext_v1.1.re.filtergene_fin

cellorder.fg16 = intersect(
  cellorder.fg16,
  colnames(src.endoderm.fg1.ext_v1.1.re))

pdf("figure.v08.07/organ_development_re/fg1_heatmap_marker.re.pdf",6,9)
#-------------------
fg16.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.fg16.ext_v1.1.re@assays$RNA@data[selectgene, cellorder.fg16]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.fg16.ext_v1.1.re$Time[cellorder.fg16], colors.time),
              MyName2Col(src.endoderm.fg16.ext_v1.1.re$cluster.v06.26.re[cellorder.fg16], cluster.endoderm.color.v5),
              #MyName2Col(src.endoderm.fg16.ext_v1.1.re$cluster.extract.v1.1[cellorder.fg16], cluster.endoderm.color.v5)
              MyName2Col(rep("FG.1",length(cellorder.fg16)), cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            #return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.fg1.ext.re@meta.data[cell_fg1_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

fg16.ext_v1.1.re.rowtree.2 = as.dendrogram(fg16.ext_v1.1.re.rowtree)
src.endoderm.fg1.ext_v1.1.re.filtergene_fin = c(
  labels(fg16.ext_v1.1.re.rowtree.2[[2]])) # New

fg16.ext_v1.1.re.rowtree.3 = as.dendrogram(fg16.ext_v1.1.re.rowtree)
src.endoderm.fg1.ext_v1.1.re.filtergene_fin = c(
  labels(fg16.ext_v1.1.re.rowtree.3[[1]]),
  rev(labels(fg16.ext_v1.1.re.rowtree.3[[2]][[2]][[2]])),
  rev(labels(fg16.ext_v1.1.re.rowtree.3[[2]][[1]][[1]])))


fg16.ext_v1.1.re.rowtree.4 = as.dendrogram(fg16.ext_v1.1.re.rowtree)
src.endoderm.fg1.ext_v1.1.re.filtergene_fin = c(
  labels(fg16.ext_v1.1.re.rowtree.4[[1]]),
  labels(fg16.ext_v1.1.re.rowtree.4[[2]][[2]]),
  rev(labels(fg16.ext_v1.1.re.rowtree.4[[2]][[1]][[2]][[2]])),
  rev(labels(fg16.ext_v1.1.re.rowtree.4[[2]][[1]][[1]])))


fg16.ext_v1.1.re.rowtree.5 = as.dendrogram(fg16.ext_v1.1.re.rowtree)
src.endoderm.fg1.ext_v1.1.re.filtergene_fin = c(
  labels(fg16.ext_v1.1.re.rowtree.5[[1]][[1]]),
  labels(fg16.ext_v1.1.re.rowtree.5[[1]][[2]][[1]]),
  labels(fg16.ext_v1.1.re.rowtree.5[[2]][[2]][[2]]),
  labels(fg16.ext_v1.1.re.rowtree.5[[2]][[1]]),
  labels(fg16.ext_v1.1.re.rowtree.5[[2]][[2]][[1]]))

cell_Ext = c("S100a6","Cryab","Anxa5","Tagln","Myl9","Sfrp1","Ntng1","Gpc4")
src.endoderm.fg1.ext_v1.1.re.filtergene_fin = setdiff(
  src.endoderm.fg1.ext_v1.1.re.filtergene_fin, cell_Ext)

names(src.endoderm.fg1.ext_v1.1.re.filtergene_fin) = c(
  rep(4, length(labels(fg16.ext_v1.1.re.rowtree.5[[1]][[1]])) + 
        length(labels(fg16.ext_v1.1.re.rowtree.5[[1]][[2]][[1]]))),
  rep(7, length(labels(fg16.ext_v1.1.re.rowtree.5[[2]][[2]][[2]]))+
        length(labels(fg16.ext_v1.1.re.rowtree.5[[2]][[1]]))+
        length(labels(fg16.ext_v1.1.re.rowtree.5[[2]][[2]][[1]]))-
        length(cell_Ext)))


save(src.endoderm.fg1.ext_v1.1.re.filtergene_fin,
     src.endoderm.fg1.ext_v1.1.re.selectgene,
     cellorder.fg16_fg1_FPI, cellorder.fg16_eso_FPI, cellorder.fg16_fg1_PI, 
     file= "figure.v08.07/organ_development_re/fg1_heatmap_parameter.Rdata")
save(src.endoderm.fg1.ext_v1.1,
     src.endoderm.fg1.ext_v1.1.re,
     file = "figure.v08.07/organ_development_re/src.endoderm.fg1.ext_v1.1.re.Rdata")
#----=================-----




cellorder.fg16 = c(cellorder.fg16_fg6_FPI, rev(cellorder.fg16_eso_FPI),
                   cellorder.fg16_fg1_PI, cellorder.fg16_pha2_FPI)
selectgene = src.endoderm.fg16.ext_v1.1.re.markergene

selectgene = src.endoderm.fg16.ext_v1.1.re.filtergene_all
pdf("figure.v08.07/organ_development_re/fg16_heatmap_marker.re.pdf",6,9)
#-------------------
fg16.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.fg16.ext_v1.1.re@assays$RNA@data[selectgene, cellorder.fg16]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.fg16.ext_v1.1.re$Time[cellorder.fg16], colors.time),
              MyName2Col(src.endoderm.fg16.ext_v1.1.re$cluster.v06.26.re[cellorder.fg16], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.fg16.ext_v1.1.re$cluster.extract.v1.1[cellorder.fg16], cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            # return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.fg1.ext.re@meta.data[cell_fg1_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

fg16.ext_v1.1.re.rowtree.all = as.dendrogram(fg16.ext_v1.1.re.rowtree)

src.endoderm.fg16.ext_v1.1.re.filtergene_all = c(
 labels(fg16.ext_v1.1.re.rowtree.all[[2]][[2]]), 
 labels(fg16.ext_v1.1.re.rowtree.all[[2]][[1]]),
 labels(fg16.ext_v1.1.re.rowtree.all[[1]][[1]]),
 labels(fg16.ext_v1.1.re.rowtree.all[[1]][[2]])
)

names(src.endoderm.fg16.ext_v1.1.re.filtergene_all) = c(
  rep(4,length(
    labels(fg16.ext_v1.1.re.rowtree.all[[2]][[2]])
  )),
  rep(3,length(
    labels(fg16.ext_v1.1.re.rowtree.all[[2]][[1]])
  )),
  rep(9,length(
    labels(fg16.ext_v1.1.re.rowtree.all[[1]][[1]])
  )),
  rep(7,length(
    labels(fg16.ext_v1.1.re.rowtree.all[[1]][[2]])
  ))
)


#----------------------------------------------------------------------------
#>   Part3: FG1 AND FG6, renew (Finally !!)
#----------------------------------------------------------------------------
#------------------
#  FG.1
#------------------

src = src.fg1.integrated
src = NormalizeData(src, assay = "RNA", scale.factor = 10^5)
src = FindVariableFeatures(src, nfeatures = 2000)
src = ScaleData(src, features = rownames(src), assay = "RNA")
src.filtergene = 
  Myfilter(as.matrix(src@assays$RNA@data),
           gene = src@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)

src = SetIdent(src, value = src$cluster.v06.26.re)
marker_src = FindAllMarkers(src)
marker_src$pct.ratio = marker_src$pct.1 / marker_src$pct.2
markergene_src = unique(marker_src[(marker_src$pct.ratio)>2,]$gene)
coverage_rate(unique(c(src.filtergene, markergene_src)), 
              gene_src.fg1.tracing.re.rowtree.1)


# Basic::Marker+Filter
#----------------------------
# --- Marker
pdf("figure.v08.07/organ_development_re_v240115/try.src.pdf",9,7)
src.rowtree =
  MyHeatmap(as.matrix(src@assays$RNA@data[markergene_src,]),
            type = "row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            ColSideColors = cbind(
              MyName2Col(src$Time, colors.time),
              MyName2Col(src$cluster.extract.v1.1, cluster.endoderm.color.v5),
              MyName2Col(src$cluster.v06.26.re, cluster.endoderm.color.v5),
              MyName2Col(src$batch, colors.num)),
            ColSideColorsSize = 3,
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            #Rowv = "none",
            return.tree = "row",
            graph = T)
dev.off()
tree_src.rowtree = as.dendrogram(src.rowtree)
gene_src.rowtree = c(
  labels(tree_src.rowtree[[1]][[2]]),
  labels(tree_src.rowtree[[2]][[2]][[1]]))
coverage_rate(src.filtergene, gene_src.fg1.tracing.re.rowtree.1)

# --- Filter
pdf("figure.v08.07/organ_development_re_v240115/try.src.pdf",9,7)
src.rowtree.1 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    src.filtergene,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time, colors.time),
      MyName2Col(src$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re, cluster.endoderm.color.v5)),
    ColSideColorsSize = 3,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()
tree_src.rowtree.1 = as.dendrogram(src.rowtree.1)
gene_src.rowtree.1 = c(
  labels(tree_src.rowtree.1[[2]][[2]][[2]][[2]]),
  labels(tree_src.rowtree.1[[1]][[1]]))
coverage_rate(gene_src.rowtree.1, gene_src.fg1.tracing.re.rowtree.1)
#----------------------------

# Intersect
#----------------------------
pdf("figure.v08.07/organ_development_re_v240115/try.src.pdf",9,7)
src.coltree.2 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    intersect(unique(c(gene_src.rowtree.1, gene_src.rowtree)), 
              gene_src.fg1.tracing.re.rowtree.1),
    cellorder_src.fg1]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time[cellorder_src.fg1], 
                 colors.time),
      MyName2Col(src$cluster.extract.v1.1[cellorder_src.fg1],  
                 cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re[cellorder_src.fg1],
                 cluster.endoderm.color.v5)),
    ColSideColorsSize = 3,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    # Rowv = "none", Colv = "none",
    return.tree = "col",
    graph = T)
src.rowtree.2 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    intersect(unique(c(gene_src.rowtree.1, gene_src.rowtree)), 
              gene_src.fg1.tracing.re.rowtree.1),
    cellorder_src.fg1]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time[cellorder_src.fg1], 
                 colors.time),
      MyName2Col(src$tree[cellorder_src.fg1],  
                 color.temp ),
      MyName2Col(src$cluster.extract.v1.1[cellorder_src.fg1],  
                 cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re[cellorder_src.fg1],
                 cluster.endoderm.color.v5)),
    ColSideColorsSize = 3,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    # Rowv = "none", Colv = "none",
    return.tree = "row",
    graph = T)
dev.off()

tree_src.coltree.2 = as.dendrogram(src.coltree.2)

src$tree = NA
src@meta.data[labels(tree_src.coltree.2[[1]]),]$tree = 1
src@meta.data[labels(tree_src.coltree.2[[2]][[1]]),]$tree = 2
src@meta.data[labels(tree_src.coltree.2[[2]][[2]][[1]]),]$tree = 3
src@meta.data[labels(tree_src.coltree.2[[2]][[2]][[2]][[1]]),]$tree = 4
src@meta.data[labels(tree_src.coltree.2[[2]][[2]][[2]][[2]]),]$tree = 5

src$cluster.v06.26.re_correct = src$cluster.v06.26.re
src@meta.data[src$tree%in%c(1,3,4),]$cluster.v06.26.re_correct = "FG.1"
src@meta.data[src$tree%in%c(2,5),]$cluster.v06.26.re_correct = "Pharynx.organ.2"
DimPlot(src, reduction = "fdl_pca", group.by = "cluster.v06.26.re_correct")

tree_src.rowtree.2 = as.dendrogram(src.rowtree.2)
gene_src.rowtree.2 = c(labels(tree_src.rowtree.2[[1]]),
                       rev(labels(tree_src.rowtree.2[[2]])))
#----------------------------

# Setdiff
#----------------------------
pdf("figure.v08.07/organ_development_re_v240115/try.src.pdf",9,7)
src.rowtree.3 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    setdiff(unique(c(gene_src.rowtree.1, gene_src.rowtree)),
            gene_src.fg1.tracing.re.rowtree.1),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time, colors.time),
      MyName2Col(src$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re_correct, cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re, cluster.endoderm.color.v5)),
    ColSideColorsSize = 3,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()
tree_src.rowtree.3 = as.dendrogram(src.rowtree.3)
gene_src.rowtree.3 = c(
  labels(tree_src.rowtree.3[[2]][[2]][[2]]),
  labels(tree_src.rowtree.3[[1]][[2]][[2]][[1]]))

pdf("figure.v08.07/organ_development_re_v240115/try.src.pdf",9,7)
src.rowtree.4 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    gene_src.rowtree.3,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time, colors.time),
      MyName2Col(src$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re_correct, cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re, cluster.endoderm.color.v5)),
    ColSideColorsSize = 3,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()
tree_src.rowtree.4 = as.dendrogram(src.rowtree.4)
gene_src.rowtree.4 = c(
  labels(tree_src.rowtree.4[[2]][[1]][[1]]),
  labels(tree_src.rowtree.4[[2]][[1]][[2]][[1]]),
  labels(tree_src.rowtree.4[[1]]))
#-----------------------------

# Finally! 
#-----------------------------
pdf("figure.v08.07/organ_development_re_v240115/try.src.pdf",9,7)
src.rowtree.5 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    c(gene_src.rowtree.2,
      gene_src.rowtree.4),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time, colors.time),
      MyName2Col(src$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re_correct, cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re, cluster.endoderm.color.v5)),
    ColSideColorsSize = 3,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()

# tree_src.fg1.integrated.rowtree.5 = as.dendrogram(src.fg1.integrated.rowtree.5)
gene_src.fg1.integrated.rowtree.5 = c(
  labels(tree_src.fg1.integrated.rowtree.5[[1]]),
  labels(tree_src.fg1.integrated.rowtree.5[[2]]))
names(gene_src.fg1.integrated.rowtree.5) = c(
  rep(4, length(labels(tree_src.fg1.integrated.rowtree.5[[1]]))),
  rep(7, length(labels(tree_src.fg1.integrated.rowtree.5[[2]]))))

# pdf("figure.v08.07/organ_development_re_v240115/try.src.pdf",9,7)
pdf("figure.v08.07/organ_development_re_v240115/try.src.fg1.integrated.pdf",9,7)
src.rowtree.6 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    gene_src.fg1.integrated.rowtree.5,
    cellorder_src.fg1.integrated.fg1]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time[cellorder_src.fg1.integrated.fg1], 
                 colors.time),
      MyName2Col(src$cluster.extract.v1.1[cellorder_src.fg1.integrated.fg1],  
                 cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re_correct[cellorder_src.fg1.integrated.fg1],
                 cluster.endoderm.color.v5)),
    RowSideColors = t(cbind(
      MyName2Col(names(gene_src.fg1.integrated.rowtree.5),
                 colors.geneset))),
    ColSideColorsSize = 3,
    RowSideColorsSize = 1,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    Rowv = "none",
    Colv = "none",
    #return.tree = "row",
    margins = c(10,10), graph = T)
dev.off()

tf_ggene_src.fg1.integrated.rowtree.5 = intersect(
  gene_src.fg1.integrated.rowtree.5, gi[gi$TF%in%T,]$SymbolDedu)
pdf("figure.v08.07/organ_development_re_v240115/try.src.fg1.integrated_tf.pdf",9,7)
src.rowtree.6 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    tf_ggene_src.fg1.integrated.rowtree.5,
    cellorder_src.fg1.integrated.fg1]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time[cellorder_src.fg1.integrated.fg1], 
                 colors.time),
      MyName2Col(src$cluster.extract.v1.1[cellorder_src.fg1.integrated.fg1],  
                 cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re_correct[cellorder_src.fg1.integrated.fg1],
                 cluster.endoderm.color.v5)),
    ColSideColorsSize = 3,
    RowSideColorsSize = 1,
    labRow = tf_ggene_src.fg1.integrated.rowtree.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    Rowv = "none",
    Colv = "none",
    #return.tree = "row",
    margins = c(10,10), graph = T)
dev.off()

#-----------------------------


src.fg1.integrated$cluster.v06.26.re_correct = 
  src$cluster.v06.26.re_correct

for(i in c("markergene_src", "tree_src.rowtree", "gene_src.rowtree",
           "src.filtergene", "tree_src.rowtree.1", "gene_src.rowtree.1",
           "tree_src.coltree.2", "tree_src.rowtree.2", "gene_src.rowtree.2",
           "tree_src.rowtree.3", "gene_src.rowtree.3",
           "tree_src.rowtree.4", "gene_src.rowtree.4",
           "tree_src.rowtree.5", "gene_src.rowtree.5",
           "cellorder_src.fg1")){
  assign(gsub("src","src.fg1.integrated",i), get(i))
}

cellorder_src.fg1.integrated.fg1 = cellorder_src.fg1

save(src.fg1.integrated,
     file = "figure.v08.07/organ_development_re_v240115/src.fg1.integrated.Rdata")

save(markergene_src.fg1.integrated, tree_src.fg1.integrated.rowtree, gene_src.fg1.integrated.rowtree,
     src.fg1.integrated.filtergene, tree_src.fg1.integrated.rowtree.1, gene_src.fg1.integrated.rowtree.1,
     tree_src.fg1.integrated.coltree.2, tree_src.fg1.integrated.rowtree.2, gene_src.fg1.integrated.rowtree.2,
     tree_src.fg1.integrated.rowtree.3, gene_src.fg1.integrated.rowtree.3,
     tree_src.fg1.integrated.rowtree.4, gene_src.fg1.integrated.rowtree.4,
     tree_src.fg1.integrated.rowtree.5, gene_src.fg1.integrated.rowtree.5,
     cellorder_src.fg1.integrated.fg1,
     file = "~/Bioinformatic/project_20221224_endoderm.refine/figure_v6.05/figure.v08.07/organ_development_re_v240115/fg1_10x_parmeter.Rdata")




#================================================================================

#------------------
#  FG.6
#------------------
for(names in names(src.fg6.integrated@reductions)){
  print(DimPlot(src.fg6.integrated, reduction = names,
                group.by = "Time", cols = colors.time)+
          ggtitle(names))
}


src = src.fg6.integrated
src = NormalizeData(src, assay = "RNA", scale.factor = 10^5)
src = FindVariableFeatures(src, nfeatures = 2000)
src = ScaleData(src, features = rownames(src), assay = "RNA")

coord_src = src[["fdl_mnn"]]@cell.embeddings
cell_src = rownames(src@meta.data)
pcurve_src = princurve::principal_curve(
  x = coord_src[cell_src,], smoother = "smooth.spline")
src$lambda = pcurve_src$lambda[cell_src]
src$lambda = norm_range(src$lambda)
src$order = pcurve_src$ord[cell_src]
cellorder_src.fg6 = names(src$lambda[order(src$lambda)])
#cellorder_src.fg6 = names(src$order[order(src$order)])

src.filtergene = 
  Myfilter(as.matrix(src@assays$RNA@data),
           gene = src@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)

src = SetIdent(src, value = src$cluster.v06.26.re)
marker_src = FindAllMarkers(src)
marker_src$pct.ratio = marker_src$pct.1 / marker_src$pct.2
markergene_src = unique(marker_src[(marker_src$pct.ratio)>2,]$gene)
coverage_rate(unique(c(src.filtergene, markergene_src)), 
              gene_src.fg6.tracing.re.rowtree)


# Basic::Marker+Filter
#----------------------------
# --- Marker
pdf("figure.v08.07/organ_development_re_v240115/try.src.pdf",9,7)
src.rowtree =
  MyHeatmap(as.matrix(src@assays$RNA@data[markergene_src,]),
            type = "row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            ColSideColors = cbind(
              MyName2Col(src$Time, colors.time),
              MyName2Col(src$cluster.extract.v1.1, cluster.endoderm.color.v5),
              MyName2Col(src$cluster.v06.26.re, cluster.endoderm.color.v5),
              MyName2Col(src$batch, colors.num)),
            ColSideColorsSize = 3,
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            #Rowv = "none",
            return.tree = "row",
            graph = T)
dev.off()
tree_src.rowtree = as.dendrogram(src.rowtree)
gene_src.rowtree = c(
  labels(tree_src.rowtree[[1]]),
  labels(tree_src.rowtree[[2]][[1]][[1]]))
coverage_rate(src.filtergene, gene_src.fg6.tracing.re.rowtree)

# --- Filter
pdf("figure.v08.07/organ_development_re_v240115/try.src.pdf",9,7)
src.rowtree.1 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    src.filtergene,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time, colors.time),
      MyName2Col(src$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re, cluster.endoderm.color.v5)),
    ColSideColorsSize = 3,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()
tree_src.rowtree.1 = as.dendrogram(src.rowtree.1)
gene_src.rowtree.1 = c(
  labels(tree_src.rowtree.1[[1]]),
  labels(tree_src.rowtree.1[[2]][[2]][[2]][[2]][[2]]))
coverage_rate(gene_src.rowtree.1, gene_src.fg6.tracing.re.rowtree)
#----------------------------

# Intersect
#----------------------------
pdf("figure.v08.07/organ_development_re_v240115/try.src.pdf",9,7)
src.coltree.2 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    intersect(unique(c(gene_src.rowtree.1, gene_src.rowtree)), 
              gene_src.fg6.tracing.re.rowtree),
    cellorder_src.fg6]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time[cellorder_src.fg6], 
                 colors.time),
      MyName2Col(src$cluster.extract.v1.1[cellorder_src.fg6],  
                 cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re[cellorder_src.fg6],
                 cluster.endoderm.color.v5)),
    ColSideColorsSize = 3,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    # Rowv = "none", Colv = "none",
    return.tree = "col",
    graph = T)
src.rowtree.2 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    intersect(unique(c(gene_src.rowtree.1, gene_src.rowtree)), 
              gene_src.fg6.tracing.re.rowtree),
    cellorder_src.fg6]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time[cellorder_src.fg6], 
                 colors.time),
      MyName2Col(src$cluster.extract.v1.1[cellorder_src.fg6],  
                 cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re[cellorder_src.fg6],
                 cluster.endoderm.color.v5)),
    ColSideColorsSize = 3,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    # Rowv = "none", Colv = "none",
    return.tree = "row",
    graph = T)
dev.off()

tree_src.coltree.2 = as.dendrogram(src.coltree.2)
tree_src.rowtree.2 = as.dendrogram(src.rowtree.2)
gene_src.rowtree.2 = c(labels(tree_src.rowtree.2[[2]]),
                       rev(labels(tree_src.rowtree.2[[1]])))
#----------------------------

# Setdiff
#----------------------------
pdf("figure.v08.07/organ_development_re_v240115/try.src.pdf",9,7)
src.rowtree.3 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    setdiff(unique(c(gene_src.rowtree.1, gene_src.rowtree)),
            gene_src.fg6.tracing.re.rowtree),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time, colors.time),
      MyName2Col(src$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re, cluster.endoderm.color.v5)),
    ColSideColorsSize = 3,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()

tree_src.rowtree.3 = as.dendrogram(src.rowtree.3)
gene_src.rowtree.3 = c(
  labels(tree_src.rowtree.3[[2]][[1]][[2]]),
  labels(tree_src.rowtree.3[[2]][[2]]),
  labels(tree_src.rowtree.3[[1]][[1]][[2]][[2]]))

pdf("figure.v08.07/organ_development_re_v240115/try.src.pdf",9,7)
src.rowtree.4 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    gene_src.rowtree.3,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time, colors.time),
      MyName2Col(src$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re, cluster.endoderm.color.v5)),
    ColSideColorsSize = 3,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()

tree_src.rowtree.4 = as.dendrogram(src.rowtree.4)
gene_src.rowtree.4 = c(
  labels(tree_src.rowtree.4[[1]]),
  labels(tree_src.rowtree.4[[2]][[1]][[2]][[2]]))
#-----------------------------

# Finally! 
#-----------------------------
pdf("figure.v08.07/organ_development_re_v240115/try.src.pdf",9,7)
src.rowtree.5 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    c(gene_src.rowtree.2,
      gene_src.rowtree.4),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time, colors.time),
      MyName2Col(src$cluster.extract.v1.1, cluster.endoderm.color.v5),
      # MyName2Col(src$cluster.v06.26.re_correct, cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re, cluster.endoderm.color.v5)),
    ColSideColorsSize = 3,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()

tree_src.rowtree.5 = as.dendrogram(src.rowtree.5)
gene_src.rowtree.5 = c(
  labels(tree_src.rowtree.5[[1]][[2]][[1]][[2]][[1]]),
  labels(tree_src.rowtree.5[[1]][[2]][[1]][[1]]),
  labels(tree_src.rowtree.5[[1]][[1]]),
  setdiff(labels(tree_src.rowtree.5[[2]][[1]]),
          c(labels(tree_src.rowtree.5[[2]][[1]][[2]][[2]][[2]][[2]]),
            labels(tree_src.rowtree.5[[2]][[1]][[1]][[1]][[2]][[2]][[2]]),
            labels(tree_src.rowtree.5[[2]][[1]][[2]][[1]][[2]][[1]]),
            labels(tree_src.rowtree.5[[2]][[1]][[1]][[2]][[2]]))))
names(gene_src.rowtree.5) = c(
  rep(4, length(c(
    labels(tree_src.rowtree.5[[1]][[2]][[1]][[2]][[1]]),
    labels(tree_src.rowtree.5[[1]][[2]][[1]][[1]]),
    labels(tree_src.rowtree.5[[1]][[1]])))),
  rep(7, length(
    setdiff(labels(tree_src.rowtree.5[[2]][[1]]),
            c(labels(tree_src.rowtree.5[[2]][[1]][[2]][[2]][[2]][[2]]),
              labels(tree_src.rowtree.5[[2]][[1]][[1]][[1]][[2]][[2]][[2]]),
              labels(tree_src.rowtree.5[[2]][[1]][[2]][[1]][[2]][[1]]),
              labels(tree_src.rowtree.5[[2]][[1]][[1]][[2]][[2]]))))))
table(gene_src.rowtree.5)


pdf("figure.v08.07/organ_development_re_v240115/try.src.pdf",9,7)
# pdf("figure.v08.07/organ_development_re_v240115/try.src.fg6.integrated.pdf",9,7)
src.rowtree.6 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    unique(c(gene_src.rowtree.5,
             intersect(gene_src.fg6.tracing.re.rowtree[names(gene_src.fg6.tracing.re.rowtree)==7], 
                       gi[gi$TF%in%T,]$SymbolDedu))),
    cellorder_src.fg6]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time[cellorder_src.fg6], 
                 colors.time),
      MyName2Col(src$cluster.extract.v1.1[cellorder_src.fg6],  
                 cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re[cellorder_src.fg6],
                 cluster.endoderm.color.v5)),
    ColSideColorsSize = 3,
    RowSideColorsSize = 1,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    # Rowv = "none",
    # Colv = "none",
    return.tree = "row",
    margins = c(10,10), graph = T)
dev.off()
tree_src.rowtree.6 = as.dendrogram(src.rowtree.6)
gene_src.rowtree.6 = c(
  labels(tree_src.rowtree.6[[2]]),
  setdiff(labels(tree_src.rowtree.6[[1]]),
          c(labels(tree_src.rowtree.6[[1]][[2]][[2]][[2]][[2]]),
            labels(tree_src.rowtree.6[[1]][[2]][[1]][[2]][[1]]))))
names(gene_src.rowtree.6) = c(
  rep(4, length(labels(tree_src.rowtree.6[[2]]))),
  rep(7, length(
    setdiff(labels(tree_src.rowtree.6[[1]]),
            c(labels(tree_src.rowtree.6[[1]][[2]][[2]][[2]][[2]]),
              labels(tree_src.rowtree.6[[1]][[2]][[1]][[2]][[1]])))
  )))


pdf("figure.v08.07/organ_development_re_v240115/try.src.pdf",9,7)
# pdf("figure.v08.07/organ_development_re_v240115/try.src.fg6.integrated.pdf",9,7)
src.rowtree.7 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    gene_src.rowtree.6,
    cellorder_src.fg6]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time[cellorder_src.fg6], 
                 colors.time),
      MyName2Col(src$cluster.extract.v1.1[cellorder_src.fg6],  
                 cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re[cellorder_src.fg6],
                 cluster.endoderm.color.v5)),
    RowSideColors = t(cbind(
      MyName2Col(names(gene_src.rowtree.6),
                 colors.geneset))),
    ColSideColorsSize = 3,
    RowSideColorsSize = 1,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    Rowv = "none",
    Colv = "none",
    # return.tree = "row",
    margins = c(10,10), graph = T)
dev.off()

tf_gene_src.rowtree.6 = intersect(
  gene_src.rowtree.6, gi[gi$TF%in%T,]$SymbolDedu)
pdf("figure.v08.07/organ_development_re_v240115/try.src.fg6.integrated_tf.pdf",9,7)
src.rowtree.7 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    tf_gene_src.rowtree.6,
    cellorder_src.fg6]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time[cellorder_src.fg6], 
                 colors.time),
      MyName2Col(src$cluster.extract.v1.1[cellorder_src.fg6],  
                 cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re[cellorder_src.fg6],
                 cluster.endoderm.color.v5)),
    # RowSideColors = t(cbind(
    #   MyName2Col(names(gene_src.rowtree.5),
    #              colors.geneset))),
    ColSideColorsSize = 3,
    RowSideColorsSize = 1,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    labRow = tf_gene_src.rowtree.6,
    Rowv = "none",
    Colv = "none",
    #return.tree = "row",
    margins = c(10,10), graph = T)
dev.off()
coverage_rate(tf_gene_src.rowtree.6, 
              intersect(gene_src.fg6.tracing.re.rowtree, gi[gi$TF%in%T,]$SymbolDedu))
#-----------------------------

# Save and Rename
#-----------------------------
for(i in c("markergene_src", "tree_src.rowtree", "gene_src.rowtree",
           "src.filtergene", "tree_src.rowtree.1", "gene_src.rowtree.1",
           "tree_src.coltree.2", "tree_src.rowtree.2", "gene_src.rowtree.2",
           "tree_src.rowtree.3", "gene_src.rowtree.3",
           "tree_src.rowtree.4", "gene_src.rowtree.4",
           "tree_src.rowtree.5", "gene_src.rowtree.5",
           "tree_src.rowtree.6", "gene_src.rowtree.6",
           "tf_gene_src.rowtree.6", 
           "cellorder_src.fg6")){
  assign(gsub("src","src.fg6.integrated",i), get(i))
}

save(src.fg6.integrated,
     file = "figure.v08.07/organ_development_re_v240115/src.fg6.integrated.Rdata")
save(markergene_src.fg6.integrated, 
     tree_src.fg6.integrated.rowtree, gene_src.fg6.integrated.rowtree,
     src.fg6.integrated.filtergene, 
     tree_src.fg6.integrated.rowtree.1, gene_src.fg6.integrated.rowtree.1,
     tree_src.fg6.integrated.coltree.2, tree_src.fg6.integrated.rowtree.2, gene_src.fg6.integrated.rowtree.2,
     tree_src.fg6.integrated.rowtree.3, gene_src.fg6.integrated.rowtree.3,
     tree_src.fg6.integrated.rowtree.4, gene_src.fg6.integrated.rowtree.4,
     tree_src.fg6.integrated.rowtree.5, gene_src.fg6.integrated.rowtree.5,
     tree_src.fg6.integrated.rowtree.6, gene_src.fg6.integrated.rowtree.6,
     tf_gene_src.fg6.integrated.rowtree.6,
     cellorder_src.fg6.integrated.fg6,
     file = "figure.v08.07/organ_development_re_v240115/fg6_10x_parmeter.Rdata")

rm(markergene_src, tree_src.rowtree, gene_src.rowtree,
   src.filtergene, tree_src.rowtree.1, gene_src.rowtree.1,
   tree_src.coltree.2, tree_src.rowtree.2, gene_src.rowtree.2,
   tree_src.rowtree.3, gene_src.rowtree.3,
   tree_src.rowtree.4, gene_src.rowtree.4,
   tree_src.rowtree.5, gene_src.rowtree.5,
   tree_src.rowtree.6, gene_src.rowtree.6,
   tf_gene_src.rowtree.6,
   cellorder_src.fg6)
#================================================================================














