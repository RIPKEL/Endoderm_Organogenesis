#----------------
#  HG.1
#----------------
src.endoderm.hg1.ext_v1.1 = FindVariableFeatures(src.endoderm.hg1.ext_v1.1, nfeatures = 2000)
src.endoderm.hg1.ext_v1.1 = ScaleData(src.endoderm.hg1.ext_v1.1, split.by = "batch_phase",
                                      features = rownames(src.endoderm.hg1.ext_v1.1))
src.endoderm.hg1.ext_v1.1.filtergene = 
  Myfilter(as.matrix(src.endoderm.hg1.ext_v1.1@assays$RNA@data),
           gene = src.endoderm.hg1.ext_v1.1@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
#rm(src.endoderm.hg1.ext_v1.1.filtergene)
src.endoderm.hg1.ext.v1.1.selectgene = unique(c(
  select_gene_hg1, src.endoderm.hg1.ext_v1.1.filtergene))

#=============================
# Row
pdf("figure.v08.07/try.pdf")
#-----------------
hg1.ext.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.hg1.ext_v1.1@assays$RNA@data[
    src.endoderm.hg1.ext_v1.1.gene.fin,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.endoderm.hg1.ext_v1.1$Time,colors.time),
      MyName2Col(src.endoderm.hg1.ext_v1.1$cluster.extract.v1.1,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.hg1.ext_v1.1$cluster.v06.26.re,cluster.endoderm.color.v5)
    ),
    ColSideColorsSize = 1.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
#-------------------
dev.off()

pdf("figure.v08.07/try.tree.pdf",100,20)
plot(hg1.ext.re.rowtree)
dev.off()

hg1.ext.re.rowtree_1 = as.dendrogram(hg1.ext.re.rowtree)
src.endoderm.hg1.ext_v1.1.gene.fin = setdiff(
  src.endoderm.hg1.ext.v1.1.selectgene,
  c(labels(hg1.ext.re.rowtree_1[[2]][[1]]),
    labels(hg1.ext.re.rowtree_1[[2]][[2]][[1]][[2]][[1]])))

#--------------------------------------------------------------------------------
src.endoderm.hg1.ext_v1.1 = RunPCA(src.endoderm.hg1.ext_v1.1,
                                   features = src.endoderm.hg1.ext_v1.1.gene.fin)
src.endoderm.hg1.ext_v1.1 = 
  seurat_mnn(seurat_name = "src.endoderm.hg1.ext_v1.1", dims = 1:30, 
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.hg1.ext_v1.1.gene.fin)
src.endoderm.hg1.ext_v1.1 = 
  seurat_int(seurat_name = "src.endoderm.hg1.ext_v1.1", dims = 1:30,
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.hg1.ext_v1.1.gene.fin)

#---------------
src.endoderm.hg1.ext_v1.1 = 
  FindNeighbors(src.endoderm.hg1.ext_v1.1, dims=1:2, 
                reduction = "umap_mnn", assay = "RNA")
src.endoderm.hg1.ext_v1.1 = 
  FindClusters(src.endoderm.hg1.ext_v1.1, resolution = 0.8)

# src.endoderm.hg1.ext_v1.1@reductions$mnn@cell.embeddings = 
#   src.endoderm.hg1.ext_v1.1@reductions$mnn@cell.embeddings[colnames(src.endoderm.hg1.ext_v1.1),]
src.endoderm.hg1.ext_v1.1 = 
  FindNeighbors(src.endoderm.hg1.ext_v1.1, dims=1:30, 
                reduction = "mnn", assay = "integrated")
src.endoderm.hg1.ext_v1.1 = 
  FindClusters(src.endoderm.hg1.ext_v1.1, resolution = 2)
src.endoderm.hg1.ext_v1.1 = 
  FindClusters(src.endoderm.hg1.ext_v1.1, resolution = 4)
#---------------

pdf("figure.v08.07/try.pdf")
DimPlot(src.endoderm.hg1.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.0.8', reduction = "umap_mnn")
DimPlot(src.endoderm.hg1.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.2', reduction = "umap_mnn")
DimPlot(src.endoderm.hg1.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.4', reduction = "umap_mnn")

DimPlot(src.endoderm.hg1.ext_v1.1, group.by = 'Time',pt.size = 1.2,
        reduction = "umap_mnn", cols = colors.time)
DimPlot(src.endoderm.hg1.ext_v1.1, group.by = 'Time', pt.size = 1.2,
        reduction = "umap_integrated", cols = colors.time)
DimPlot(src.endoderm.hg1.ext_v1.1, group.by = 'batch', 
        reduction = "umap_integrated") #, cols = colors.time)

DimPlot(src.endoderm.hg1.ext_v1.1, group.by = 'cluster.extract.v1.1', 
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.hg1.ext_v1.1, group.by = 'cluster.extract.v1.1', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.hg1.ext_v1.1, group.by = 'cluster.v06.26.re', pt.size = 1.2,
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.hg1.ext_v1.1, group.by = 'cluster.v06.26.re', pt.size = 1.2,
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
dev.off()

src.endoderm.hg1.ext_v1.1@reductions$pca_integrated@cell.embeddings = 
  src.endoderm.hg1.ext_v1.1@reductions$pca_integrated@cell.embeddings[colnames(src.endoderm.hg1.ext_v1.1),]
#--------------------
src.endoderm.hg1.ext_v1.1$cluster.v06.26.re_raw = NA
src.endoderm.hg1.ext_v1.1@meta.data[
  intersect(colnames(src.endoderm.hg1.ext_v1.1),colnames(src.endoderm.hg1.ext.re)),]$cluster.v06.26.re_raw = 
  src.endoderm.hg1.ext.re@meta.data[
    intersect(colnames(src.endoderm.hg1.ext_v1.1),colnames(src.endoderm.hg1.ext.re)),]$cluster.v06.26.re
# src.endoderm.hg1.ext.re$cluster.v06.26.re = NA
# src.endoderm.hg1.ext.re@meta.data[colnames(src.endoderm.hg1.ext.re.sma.1_lar.1),]$cluster.v06.26.re =
#   src.endoderm.hg1.ext.re.sma.1_lar.1$cluster.v06.26.re
# src.endoderm.hg1.ext.re@meta.data[colnames(src.endoderm.hg1.ext.re.lar.2),]$cluster.v06.26.re =
#   src.endoderm.hg1.ext.re.lar.2$cluster.v06.26
#---------------------
src.endoderm.hg1.ext_v1.1$cluster.v06.26.re = NA
src.endoderm.hg1.ext_v1.1@meta.data[
  intersect(colnames(src.endoderm.hg1.ext_v1.1),colnames(src.endoderm.hg1.ext.re)),]$cluster.v06.26.re = 
  src.endoderm.hg1.ext.re@meta.data[
    intersect(colnames(src.endoderm.hg1.ext_v1.1),colnames(src.endoderm.hg1.ext.re)),]$cluster.v06.26.re


src.endoderm.hg1.ext_v1.1@meta.data[
  src.endoderm.hg1.ext_v1.1$RNA_snn_res.2%in%c(7,12),]$cluster.v06.26.re = "Large.intestine.2"
src.endoderm.hg1.ext_v1.1@meta.data[
  src.endoderm.hg1.ext_v1.1$RNA_snn_res.2%in%c(16,10,9,14,15),]$cluster.v06.26.re = "HG.1"
src.endoderm.hg1.ext_v1.1@meta.data[
  src.endoderm.hg1.ext_v1.1$RNA_snn_res.2%in%c(5,8,13,7),]$cluster.v06.26.re = "HG.1-Large.intestine.2"

src.endoderm.hg1.ext_v1.1@meta.data[
  src.endoderm.hg1.ext_v1.1$cluster.extract.v1.1%in%
    c("MG.2/3-HG.1"),]$cluster.v06.26.re = "Small.intestine.2"
src.endoderm.hg1.ext_v1.1@meta.data[
  src.endoderm.hg1.ext_v1.1$Time%in%c("ss24","ss27")&
    src.endoderm.hg1.ext_v1.1$cluster.revise.re.v1.30.re.fin%in%c("Small.intestine.2"),]$cluster.v06.26.re = "Small.intestine.2"
src.endoderm.hg1.ext_v1.1@meta.data[
  src.endoderm.hg1.ext_v1.1$Time%in%c("ss24","ss27")&
    src.endoderm.hg1.ext_v1.1$cluster.revise.re.v1.30.re.fin%in%c("Large.intestine.1"),]$cluster.v06.26.re = "Large.intestine.1"
src.endoderm.hg1.ext_v1.1@meta.data[
  src.endoderm.hg1.ext_v1.1$Time%in%c("ss24","ss27")&
    src.endoderm.hg1.ext_v1.1$cluster.revise.re.v1.30.re.fin%in%c("Large.intestine.2"),]$cluster.v06.26.re = "Large.intestine.2"

a = FNN::knn(
  src.endoderm.hg1.ext_v1.1@reductions$pca_integrated@cell.embeddings[
    rownames(src.endoderm.hg1.ext_v1.1@meta.data[!src.endoderm.hg1.ext_v1.1$cluster.v06.26.re%in%NA,]),],
  src.endoderm.hg1.ext_v1.1@reductions$pca_integrated@cell.embeddings[
    rownames(src.endoderm.hg1.ext_v1.1@meta.data[src.endoderm.hg1.ext_v1.1$cluster.v06.26.re%in%NA,]),],
  src.endoderm.hg1.ext_v1.1@meta.data[!src.endoderm.hg1.ext_v1.1$cluster.v06.26.re%in%NA,]$cluster.v06.26.re, k = 10)
src.endoderm.hg1.ext_v1.1@meta.data[
  src.endoderm.hg1.ext_v1.1$cluster.v06.26.re%in%NA,]$cluster.v06.26.re = as.character(a)

DimPlot(src.endoderm.hg1.ext_v1.1, group.by = 'cluster.v06.26.re', pt.size = 1.2,
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)

# FDL
src.endoderm.hg1.ext_v1.1 = 
  seurat_fdl("src.endoderm.hg1.ext_v1.1", "pca","RNA")
src.endoderm.hg1.ext_v1.1 = 
  seurat_fdl("src.endoderm.hg1.ext_v1.1", "mnn","RNA")
src.endoderm.hg1.ext_v1.1 = 
  seurat_fdl("src.endoderm.hg1.ext_v1.1", "pca_integrated","integrated")

# Graph
pdf("figure.v08.07/organ_development_re/hg1_summary_graph.re_snn_PI.pdf",12,7)
# SNN
for(red in c('pca',"mnn","pca_integrated")){
  
  assay = ifelse(red=="mnn","RNA","integrated")
  
  graph.hg1.ext_v1.1 = 
    seurat_GCN(seurat_name = "src.endoderm.hg1.ext_v1.1",
               src.endoderm.hg1.ext_v1.1.gene.fin,
               reduction = red, assay =  assay)
  graph.hg1_cluster = cluster_walktrap(graph.hg1.ext_v1.1, steps = 1)
  
  for(color in c("cluster.extract.v1.1","cluster.v06.26.re","cluster.v06.26.re_hc",
                 "Time")){
    for(edge.color in c(NA,"lightgray")){
      for(edge.color in c(NA,"lightgray")){
        set.seed(1)
        plot.igraph(graph.hg1.ext_v1.1,
                    layout = layout_with_fr(graph.hg1.ext_v1.1),
                    edge.color = edge.color,
                    vertex.size=2.5, vertex.label=NA,
                    vertex.label.color = "black", 
                    vertex.frame.color = NA, vertex.frame.width = 0.5,
                    #vertex.color =  colors.time[src.endoderm.hg1.ext_v1.1[["Time"]][graph.hg1_cluster$names,]],
                    #vertex.color =  cluster.endoderm.color.v5[src.endoderm.hg1.ext_v1.1[["cluster.v06.26.re"]][graph.hg1_cluster$names,]],
                    vertex.color =  colors_bg[src.endoderm.hg1.ext_v1.1[[color]][graph.hg1_cluster$names,]])
      }}}
  
}
dev.off()

pdf("figure.v08.07/organ_development_re/hg1_summary_graph.re_other.pdf",12,7)
# Reduction
for(red in c("umap_integrated","umap_mnn",
             "fdl_pca",'fdl_mnn',"fdl_pca_integrated")){
  p1 = DimPlot(src.endoderm.hg1.ext_v1.1,group.by = 'cluster.v06.26.re', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p2 = DimPlot(src.endoderm.hg1.ext_v1.1,group.by = 'cluster.extract.v1.1', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p3 = DimPlot(src.endoderm.hg1.ext_v1.1,group.by = 'Time', 
               reduction = red , cols = colors.time,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p4 = DimPlot(src.endoderm.hg1.ext_v1.1,group.by = 'cluster.v06.26.re_hc', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  print(p1); print(p2); print(p3);  print(p4)
}
dev.off()


pdf("figure.v08.07/organ_development_re/hg1_summary_proportion.re.pdf",12,7)
#-----------------------------
cell_type = c("HG.1","HG.1-Large.intestine.2",
              "Small.intestine.2","Large.intestine.1","Large.intestine.2")
meta_hg1 = src.endoderm.hg1.ext_v1.1@meta.data
meta_hg1$Time = factor(meta_hg1$Time, levels = names(colors.time))
meta_hg1$cluster.extract.v1.1 = factor(meta_hg1$cluster.extract.v1.1, 
                                       levels = c("HG.1","MG.2/3-HG.1"))
meta_hg1$cluster.v06.26.re = factor(meta_hg1$cluster.v06.26.re, levels = rev(cell_type))
meta_hg1$cluster.v06.26.re_hc = factor(meta_hg1$cluster.v06.26.re, levels = rev(cell_type))

meta_hg1_B0 = meta_hg1
meta_hg1_B1 = meta_hg1[meta_hg1$batch%in%1&!meta_hg1$Time%in%"ss9",]
meta_hg1_B2 = meta_hg1[meta_hg1$batch%in%2,]

for(data in c("meta_hg1_B0","meta_hg1_B1","meta_hg1_B2")){
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
  
  print(p + 
           geom_bar(data = get(data), 
                    mapping = aes(x = Time,
                                  group = cluster.v06.26.re_hc, fill = cluster.v06.26.re_hc),
                    stat = "count", position = 'fill') + 
           xlab("Time") + ylab("Cluster proportion of
  cell type")+ ggtitle(data))
  
}
#-----------------------------
dev.off()

save(src.endoderm.hg1.ext_v1.1,
     file = "figure.v08.07/organ_development_re/src.endoderm.hg1.ext_v1.1.re.Rdata")


# Marker
#------------------------
src.endoderm.hg1.ext_v1.1 = 
  SetIdent(src.endoderm.hg1.ext_v1.1,
           value = src.endoderm.hg1.ext_v1.1$cluster.v06.26.re)
markergene.hg1.ext_v1.1 = 
  FindAllMarkers(src.endoderm.hg1.ext_v1.1, assay = "RNA")
markergene.hg1.ext_v1.1 = markergene.hg1.ext_v1.1[
  markergene.hg1.ext_v1.1$avg_log2FC>0.2&
    markergene.hg1.ext_v1.1$p_val_adj<0.1, "gene"]

#------------------------
hg1.umap.embedding = src.endoderm.hg1.ext_v1.1[["fdl_pca_integrated"]]@cell.embeddings[,1:3]
hg1.umap.embedding = src.endoderm.hg1.ext_v1.1[["fdl_mnn"]]@cell.embeddings[,1:3]
hg1.umap.embedding = src.endoderm.hg1.ext_v1.1[["umap_mnn"]]@cell.embeddings[,1:3]
hg1.umap.embedding = src.endoderm.hg1.ext_v1.1[["umap_integrated"]]@cell.embeddings[,1:3]

# Cell order :: Princurve
#---------------------------
# HG.1
cellorder.hg1_hg1 = 
  rownames(src.endoderm.hg1.ext_v1.1@meta.data[
    src.endoderm.hg1.ext_v1.1$cluster.v06.26.re_hc%in%"HG.1",])
princurve.hg1_hg1 = princurve::principal_curve(
  x = hg1.umap.embedding[cellorder.hg1_hg1,],  smoother = "smooth.spline")
cellorder.hg1_hg1 =  names(
  princurve.hg1_hg1$lambda[cellorder.hg1_hg1][order(princurve.hg1_hg1$lambda[cellorder.hg1_hg1])])
# Small.intestine.2
cellorder.hg1_s2 = 
  rownames(src.endoderm.hg1.ext_v1.1@meta.data[
    src.endoderm.hg1.ext_v1.1$cluster.v06.26.re_hc%in%"Small.intestine.2",])
princurve.hg1_s2 = princurve::principal_curve(
  x = hg1.umap.embedding[cellorder.hg1_s2,],  smoother = "smooth.spline")
cellorder.hg1_s2 =  names(
  princurve.hg1_s2$lambda[cellorder.hg1_s2][order(princurve.hg1_s2$lambda[cellorder.hg1_s2])])
# Large.intestine.1
cellorder.hg1_l1 = 
  rownames(src.endoderm.hg1.ext_v1.1@meta.data[
    src.endoderm.hg1.ext_v1.1$cluster.v06.26.re_hc%in%"Large.intestine.1",])
princurve.hg1_l1 = princurve::principal_curve(
  x = hg1.umap.embedding[cellorder.hg1_l1,],  smoother = "smooth.spline")
cellorder.hg1_l1 =  names(
  princurve.hg1_l1$lambda[cellorder.hg1_l1][order(princurve.hg1_l1$lambda[cellorder.hg1_l1])])
# Large.intestine.2
cellorder.hg1_l2 = 
  rownames(src.endoderm.hg1.ext_v1.1@meta.data[
    src.endoderm.hg1.ext_v1.1$cluster.v06.26.re_hc%in%"Large.intestine.2",])
princurve.hg1_l2 = princurve::principal_curve(
  x = hg1.umap.embedding[cellorder.hg1_l2,],  smoother = "smooth.spline")
cellorder.hg1_l2 =  names(
  princurve.hg1_l2$lambda[cellorder.hg1_l2][order(princurve.hg1_l2$lambda[cellorder.hg1_l2])])
# HG.1-Large.intestine.2
cellorder.hg1_hg1l2 = 
  rownames(src.endoderm.hg1.ext_v1.1@meta.data[
    src.endoderm.hg1.ext_v1.1$cluster.v06.26.re_hc%in%"HG.1-Large.intestine.2",])
princurve.hg1_hg1l2 = princurve::principal_curve(
  x = hg1.umap.embedding[cellorder.hg1_hg1l2,],  smoother = "smooth.spline")
cellorder.hg1_hg1l2 =  names(
  princurve.hg1_hg1l2$lambda[cellorder.hg1_hg1l2][order(princurve.hg1_hg1l2$lambda[cellorder.hg1_hg1l2])])

#---------------------------

pdf("figure.v08.07/try.pdf")
#---------------------------
plot(c(1:length(cellorder.hg1_hg1)),
     c(1:length(cellorder.hg1_hg1)),
     col = colors.time[src.endoderm.hg1.ext_v1.1[,cellorder.hg1_hg1]$Time])
plot(c(1:length(cellorder.hg1_s2)),
     c(1:length(cellorder.hg1_s2)),
     col = colors.time[src.endoderm.hg1.ext_v1.1[,cellorder.hg1_s2]$Time])
plot(c(1:length(cellorder.hg1_l1)),
     c(1:length(cellorder.hg1_l1)),
     col = colors.time[src.endoderm.hg1.ext_v1.1[,cellorder.hg1_l1]$Time])
plot(c(1:length(cellorder.hg1_l2)),
     c(1:length(cellorder.hg1_l2)),
     col = colors.time[src.endoderm.hg1.ext_v1.1[,cellorder.hg1_l2]$Time])
plot(c(1:length(cellorder.hg1_hg1l2)),
     c(1:length(cellorder.hg1_hg1l2)),
     col = colors.time[src.endoderm.hg1.ext_v1.1[,cellorder.hg1_hg1l2]$Time])
#---------------------------
dev.off()

cellorder.hg1_hg1_FM = cellorder.hg1_hg1
cellorder.hg1_s2_FM = cellorder.hg1_s2
cellorder.hg1_l1_FM = cellorder.hg1_l1
cellorder.hg1_l2_FM = cellorder.hg1_l2
cellorder.hg1_hg1l2_UM = cellorder.hg1_hg1l2

#=================================================================================
#=================================================================================
pdf("figure.v08.07/organ_development_re/hg1_heatmap_marker.re_l2.pdf",6,9)
cellorder.hg1 = c(rev(cellorder.hg1_hg1_FM), 
                  cellorder.hg1_hg1l2_UM,
                  cellorder.hg1_l2_FM)

markergene.hg1.ext_v1.1 = 
  FindAllMarkers(src.endoderm.hg1.ext_v1.1[,cellorder.hg1], assay = "RNA")
markergene.hg1.ext_v1.1 = markergene.hg1.ext_v1.1[
  markergene.hg1.ext_v1.1$avg_log2FC>0.2&markergene.hg1.ext_v1.1$p_val_adj<0.1, "gene"]
selectgene = unique(c(markergene.hg1.ext_v1.1))

selectgene = src.endoderm.hg1.ext_v1.1.filtergene_fin
#-------------------
hg1.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.hg1.ext_v1.1@assays$RNA@data[selectgene, cellorder.hg1]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.hg1.ext_v1.1$Time[cellorder.hg1], colors.time),
              MyName2Col(src.endoderm.hg1.ext_v1.1$cluster.v06.26.re_hc[cellorder.hg1], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.hg1.ext_v1.1$cluster.extract.v1.1[cellorder.hg1], cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            # return.tree = "col",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.hg1.ext_v1.1@meta.data[cell_hg1_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

hg1.ext_v1.1.re.rowtree.1 = as.dendrogram(hg1.ext_v1.1.re.rowtree)
src.endoderm.hg1.ext_v1.1.filtergene_fin = c(
  setdiff(
    c(
      labels(hg1.ext_v1.1.re.rowtree.1[[1]])
    ),
    c()
  ),
  setdiff(
    c(
      labels(hg1.ext_v1.1.re.rowtree.1[[2]][[1]])
    ),
    c("Dkk1")
  ),
  setdiff(
   c(
      labels(hg1.ext_v1.1.re.rowtree.1[[2]][[2]][[2]][[2]][[1]]),
      labels(hg1.ext_v1.1.re.rowtree.1[[2]][[2]][[2]][[1]]),
      rev(labels(hg1.ext_v1.1.re.rowtree.1[[2]][[2]][[1]][[2]][[2]]))
      #labels(hg1.ext_v1.1.re.rowtree.1[[2]][[2]][[2]][[2]][[2]])
    ),
    c()
  )
)

names(src.endoderm.hg1.ext_v1.1.filtergene_fin) = c(
  rep(4,length(  
    setdiff(
      c(
        labels(hg1.ext_v1.1.re.rowtree.1[[1]])
      ),
      c()
    ))),
  rep(3,length(  
    setdiff(
      c(
        labels(hg1.ext_v1.1.re.rowtree.1[[2]][[1]])
      ),
      c("Dkk1")
    ))),
  rep(7,length(  
    setdiff(
      c(
        labels(hg1.ext_v1.1.re.rowtree.1[[2]][[2]][[1]][[2]][[2]]),
        labels(hg1.ext_v1.1.re.rowtree.1[[2]][[2]][[2]])
      ),
      c(
        labels(hg1.ext_v1.1.re.rowtree.1[[2]][[2]][[2]][[2]][[2]])
      )
    )))
)

src.endoderm.hg1.ext_v1.1.filtergene_l2  = 
  src.endoderm.hg1.ext_v1.1.filtergene_fin 



hg1.ext_v1.1.re.coltree.1 = as.dendrogram(hg1.ext_v1.1.re.rowtree)
src.endoderm.hg1.ext_v1.1$cluster.v06.26.re_hc = 
  src.endoderm.hg1.ext_v1.1$cluster.v06.26.re
src.endoderm.hg1.ext_v1.1@meta.data[
  labels(hg1.ext_v1.1.re.coltree.1[[2]][[1]]),]$cluster.v06.26.re_hc = "Large.intestine.2"
src.endoderm.hg1.ext_v1.1@meta.data[
  labels(hg1.ext_v1.1.re.coltree.1[[2]][[2]]),]$cluster.v06.26.re_hc = "HG.1-Large.intestine.2"
DimPlot(src.endoderm.hg1.ext_v1.1, group.by = 'cluster.v06.26.re_hc', pt.size = 1.2,
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)


#=================================================================================
#=================================================================================
pdf("figure.v08.07/organ_development_re/hg1_heatmap_marker.re_s2l1.pdf",6,9)
cellorder.hg1 = c(rev(cellorder.hg1_hg1_FM), 
                  cellorder.hg1_s2_FM,
                  cellorder.hg1_l1_FM)

markergene.hg1.ext_v1.1 = 
  FindAllMarkers(src.endoderm.hg1.ext_v1.1[,cellorder.hg1], assay = "RNA")
markergene.hg1.ext_v1.1 = markergene.hg1.ext_v1.1[
  markergene.hg1.ext_v1.1$avg_log2FC>0.2&markergene.hg1.ext_v1.1$p_val_adj<0.1, "gene"]
selectgene = unique(c(markergene.hg1.ext_v1.1))

selectgene = src.endoderm.hg1.ext_v1.1.filtergene_fin
#-------------------
hg1.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.hg1.ext_v1.1@assays$RNA@data[selectgene, cellorder.hg1]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.hg1.ext_v1.1$Time[cellorder.hg1], colors.time),
              MyName2Col(src.endoderm.hg1.ext_v1.1$cluster.v06.26.re_hc[cellorder.hg1], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.hg1.ext_v1.1$cluster.extract.v1.1[cellorder.hg1], cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            #return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.hg1.ext_v1.1@meta.data[cell_hg1_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

hg1.ext_v1.1.re.rowtree.2 = as.dendrogram(hg1.ext_v1.1.re.rowtree)
src.endoderm.hg1.ext_v1.1.filtergene_fin = c(
  setdiff(
    c(
      rev(labels(hg1.ext_v1.1.re.rowtree.2[[2]][[2]][[2]])),
      labels(hg1.ext_v1.1.re.rowtree.2[[2]][[2]][[1]])
    ),
    c()
  ),
  setdiff(
    c(
      labels(hg1.ext_v1.1.re.rowtree.2[[2]][[1]][[1]])
    ),
    c("Hist1h1b","Hist1h2ae")
  ),
  setdiff(
    c(
      labels(hg1.ext_v1.1.re.rowtree.2[[1]][[2]]),
      labels(hg1.ext_v1.1.re.rowtree.2[[1]][[1]])
    ),
    c()
  )
)

names(src.endoderm.hg1.ext_v1.1.filtergene_fin) = c(
  rep(4,length(  
    setdiff(
      c(
        labels(hg1.ext_v1.1.re.rowtree.2[[2]][[2]])
      ),
      c()
    ))),
  rep(3,length(  
    setdiff(
      c(
        labels(hg1.ext_v1.1.re.rowtree.2[[2]][[1]][[1]])
      ),
      c("Hist1h1b","Hist1h2ae")
    ))),
  rep(7,length(  
    setdiff(
      c(
        labels(hg1.ext_v1.1.re.rowtree.2[[1]][[2]]),
        labels(hg1.ext_v1.1.re.rowtree.2[[1]][[1]])
      ),
      c()
    )))
)


src.endoderm.hg1.ext_v1.1.filtergene_s2l1  = 
  src.endoderm.hg1.ext_v1.1.filtergene_fin 

#--------------------------------------------------------------------------------


save(
  # Gene 
  src.endoderm.hg1.ext.v1.1.selectgene,
  src.endoderm.hg1.ext_v1.1.gene.fin,
  #src.endoderm.hg1.ext_v1.1.filtergene_fin,
  src.endoderm.hg1.ext_v1.1.filtergene_l2,
  src.endoderm.hg1.ext_v1.1.filtergene_s2l1,
  # Cell type
  cellorder.hg1_hg1_FM, cellorder.hg1_hg1l2_UM, cellorder.hg1_l2_FM,
  cellorder.hg1_s2_FM, cellorder.hg1_l1_FM,
  file= "figure.v08.07/organ_development_re/hg1_heatmap_parameter.Rdata")




#---------------------
#-->> remake for HG.1
#---------------------
src.hg1.integrated.re = src.endoderm.hg1.ext_v1.1
src.hg1.integrated.re = NormalizeData(src.hg1.integrated.re, scale.factor = 10^5)
src.hg1.integrated.re = FindVariableFeatures(src.hg1.integrated.re, nfeatures = 2000)
src.hg1.integrated.re.filtergene = 
  Myfilter(as.matrix(src.hg1.integrated.re@assays$RNA@data),
           gene = src.hg1.integrated.re@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.hg1.integrated.re.filtergene.re = src.hg1.integrated.re.filtergene; rm(src.hg1.integrated.re.filtergene)
src.hg1.integrated.re.filtergene = src.hg1.integrated.re.filtergene.re; rm(src.hg1.integrated.re.filtergene.re)


src.hg1.integrated.re = SetIdent(src.hg1.integrated.re, value = src.hg1.integrated.re$cluster.v06.26.re_hc)
marker_src.hg1.integrated.re = FindAllMarkers(src.hg1.integrated.re)
marker_src.hg1.integrated.re$pct.ratio = marker_src.hg1.integrated.re$pct.1 / marker_src.hg1.integrated.re$pct.2
marker_src.hg1.integrated.re$rank = marker_src.hg1.integrated.re$pct.ratio * (-log(marker_src.hg1.integrated.re$p_val_adj))
marker_src.hg1.integrated.re = marker_src.hg1.integrated.re[order(marker_src.hg1.integrated.re$rank, decreasing = T),]
markergene_src.hg1.integrated.re = unique(marker_src.hg1.integrated.re$gene)
# markergene_src.hg1.integrated.re.raw = markergene_src.hg1.integrated.re

src.hg1.integrated.re = ScaleData(src.hg1.integrated.re, assay = "RNA",
                                  features = rownames(src.hg1.integrated.re@assays$integrated@data), 
                                  split.by = "batch_phase")
src.hg1.integrated.re = RunPCA(src.hg1.integrated.re, assay = "RNA",
                               features = rownames(src.hg1.integrated.re@assays$integrated@data))
src.hg1.integrated.re = RunUMAP(src.hg1.integrated.re, reduction = "mnn", 
                                dims = 1:30, n.neighbors = 100, # learning.rate = 0.001,
                                umap.method = "uwot-learn", seed.use = 30, # n.epochs = 500,
                                assay = "RNA", n.components = 3, min.dist = 0.01)

DimPlot(src.hg1.integrated.re, reduction = "umap", dims = c(1,3),
        group.by = "batch", label = T, label.size = 3)
DimPlot(src.hg1.integrated.re, reduction = "umap", dims = c(2,3), 
        group.by = "Time", cols = colors.time)
DimPlot(src.hg1.integrated.re, reduction = "umap", dims = c(2,3), 
        group.by = "cluster.v06.26.re_hc", cols = cluster.endoderm.color.v5)

src.hg1.integrated.re@reductions$umap_mnn_rot = src.hg1.integrated.re@reductions$umap
# src.hg1.integrated.re@reductions$umap_mnn_rot@cell.embeddings = 
#   t(view_src.hg1.integrated.re[1:3,1:3] %*% 
#       t(as.matrix(src.hg1.integrated.re@reductions$umap@cell.embeddings)))
src.hg1.integrated.re@reductions$umap_mnn_rot@cell.embeddings = cbind(
  src.hg1.integrated.re@reductions$umap@cell.embeddings[,2],
  cbind(src.hg1.integrated.re@reductions$umap@cell.embeddings[,3],
        src.hg1.integrated.re@reductions$umap@cell.embeddings[,1]))
colnames(src.hg1.integrated.re@reductions$umap_mnn_rot@cell.embeddings) = paste("UMAP_",c(1:3),sep="")
src.hg1.integrated.re@reductions$umap_mnn_rot@key = "UMAP_"

DimPlot(src.hg1.integrated.re, reduction = "umap_mnn_rot",
        group.by = "cluster.v06.26.re_hc", dims = c(1,2), 
        cols = cluster.endoderm.color.v5)
DimPlot(src.hg1.integrated.re, reduction = "umap_mnn_rot",
        group.by = "Time", dims = c(1,2), 
        cols = colors.time)

save(src.hg1.integrated.re, 
     file= "figure.v08.07/organ_development_re_v240115/src.hg1.integrated.re.Rdata")

src.hg1.integrated.merge.re = src.hg1.integrated.merge
src.hg1.integrated.merge.re@reductions$mnn_umap_fta = 
  CreateDimReducObject(embeddings = matrix(0, 
                                           nrow = length(colnames(src.hg1.integrated.merge.re)),
                                           ncol = 3),
                       # loadings = matrix(0), projected = matrix(0),
                       key = "UMAP_")
rownames(src.hg1.integrated.merge.re@reductions$mnn_umap_fta@cell.embeddings) = colnames(src.hg1.integrated.merge.re)
src.hg1.integrated.merge.re@reductions$mnn_umap_fta@cell.embeddings[colnames(src.hg1.integrated.re),] =
  src.hg1.integrated.re@reductions$umap_mnn_rot@cell.embeddings

# src.hg1.integrated.re = NormalizeData(src.hg1.integrated.re, scale.factor = 10^(4), assay = "RNA")
# src.hg1.integrated.re@assays$mnnRNA = 
#   CreateAssayObject(data = src.hg1.integrated.merge.re@assays$mnnRNA@data[,colnames(src.hg1.integrated.re)])
# src.hg1.integrated.merge.re@assays$mnnRNA[colnames(src.hg1.integrated.re)]

src.hg1.integrated = NormalizeData(src.hg1.integrated, scale.factor = 10^(4.7), assay = "RNA")
anchor.integrated = 
  FindTransferAnchors(reference = src.hg1.integrated, # src.hg1.integrated.merge.re[,colnames(src.hg1.integrated)],
                      query =  src.hg1.integrated.merge.re[,colnames(src.hg1.tracing)],
                      reference.assay = "RNA",
                      query.assay = "RNA", 
                      scale = F, 
                      features = src.endoderm.hg1.ext.v1.1.selectgene)
umap.transfer = TransferData(anchor.integrated, 
                             t(src.hg1.integrated.re@reductions$umap_mnn_rot@cell.embeddings))
cell_query = rownames(src.hg1.integrated.merge.re@meta.data[!src.hg1.integrated.merge.re$lineage%in%NA,])
src.hg1.integrated.merge.re@reductions$mnn_umap_fta@cell.embeddings[cell_query,] = as.matrix(t(umap.transfer@data))[cell_query,]

src.hg1.tracing@reductions$mnn_umap_fta = 
  CreateDimReducObject(embeddings =as.matrix(t(umap.transfer@data)), key = "UMAP_")

src.hg1.integrated.merge.re = 
  label_KNN_learn_Re(seurat = src.hg1.integrated.merge.re,
                     reduction = "mnn_umap_fta", 
                     label = "cluster.v06.26.re_hc", 
                     group = "Source_tech")

save(src.hg1.integrated.merge.re, file = "figure.v08.07/organ_development_re_v240115/src.hg1.integrated.merge.re.Rdata")
save(src.hg1.tracing, file = "figure.v08.07/organ_development_re_v240115/src.hg1.tracing.Rdata")


# == HG.1 Basic process 
#===============================================================================
pdf("figure.v08.07/organ_development_re_v240115/try.hg1_10x.pdf",9,7)
src.hg1.integrated.re.rowtree =
  MyHeatmap(as.matrix(src.hg1.integrated.re@assays$RNA@data[
    unique(c(markergene_src.hg1.integrated.re,
             src.hg1.integrated.re.filtergene,
             c()
    )),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.hg1.integrated.re$Time,
                 colors.time.2),
      MyName2Col(src.hg1.integrated.re$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.hg1.integrated.re$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.hg1.integrated.re$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()
pdf("figure.v08.07/organ_development_re_v240115/try.hg1_10x.tree.pdf",150,20)
plot(src.hg1.integrated.re.rowtree)
dev.off()
tree_src.hg1.integrated.re.rowtree = as.dendrogram(src.hg1.integrated.re.rowtree)
gene_src.hg1.integrated.re.rowtree = setdiff(
  unique(c(markergene_src.hg1.integrated.re, src.hg1.integrated.re.filtergene)),
  c(labels(tree_src.hg1.integrated.re.rowtree[[1]][[2]]),
    labels(tree_src.hg1.integrated.re.rowtree[[2]][[2]][[2]][[1]])))




