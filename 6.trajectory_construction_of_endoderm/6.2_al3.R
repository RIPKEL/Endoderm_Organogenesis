#----------------
#    AL.3
#----------------
src.endoderm.al3.ext_v1.1 = FindVariableFeatures(src.endoderm.al3.ext_v1.1, nfeatures = 2000)
src.endoderm.al3.ext_v1.1 = ScaleData(src.endoderm.al3.ext_v1.1, split.by = "batch_phase",
                                      features = rownames(src.endoderm.al3.ext_v1.1))
src.endoderm.al3.ext_v1.1.filtergene = 
  Myfilter(as.matrix(src.endoderm.al3.ext_v1.1@assays$RNA@data),
           gene = src.endoderm.al3.ext_v1.1@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
#rm(src.endoderm.al3.ext_v1.1.filtergene)
src.endoderm.al3.ext.v1.1.selectgene = unique(c(
  select_gene_al3, src.endoderm.al3.ext_v1.1.filtergene))

#=============================
# Row
pdf("figure.v08.07/try.pdf")
#-----------------
al3.ext.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.al3.ext_v1.1@assays$RNA@data[
    src.endoderm.al3.ext.v1.1.selectgene,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.endoderm.al3.ext_v1.1$Time,colors.time),
      MyName2Col(src.endoderm.al3.ext_v1.1$cluster.extract.v1.1,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.al3.ext_v1.1$cluster.v06.26.re,cluster.endoderm.color.v5)
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
plot(al3.ext.re.rowtree)
dev.off()

al3.ext.re.rowtree_1 = as.dendrogram(al3.ext.re.rowtree)
src.endoderm.al3.ext_v1.1.gene.fin = setdiff(
  src.endoderm.al3.ext.v1.1.selectgene,
  c(labels(al3.ext.re.rowtree_1[[2]][[1]]),
    labels(al3.ext.re.rowtree_1[[2]][[2]][[1]])))

#--------------------------------------------------------------------------------
src.endoderm.al3.ext_v1.1 = RunPCA(src.endoderm.al3.ext_v1.1,
                                   features = src.endoderm.al3.ext_v1.1.gene.fin)
src.endoderm.al3.ext_v1.1 = 
  seurat_mnn(seurat_name = "src.endoderm.al3.ext_v1.1", dims = 1:30, 
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.al3.ext_v1.1.gene.fin)
src.endoderm.al3.ext_v1.1 = 
  seurat_int(seurat_name = "src.endoderm.al3.ext_v1.1", dims = 1:30,
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.al3.ext_v1.1.gene.fin)

#---------------
src.endoderm.al3.ext_v1.1 = 
  FindNeighbors(src.endoderm.al3.ext_v1.1, dims=1:30, 
                reduction = "pca_integrated", assay = "integrated")
src.endoderm.al3.ext_v1.1 = 
  FindNeighbors(src.endoderm.al3.ext_v1.1, dims=1:30, 
                reduction = "pca", assay = "RNA")
src.endoderm.al3.ext_v1.1 = 
  FindClusters(src.endoderm.al3.ext_v1.1, resolution = 2)
src.endoderm.al3.ext_v1.1 = 
  FindClusters(src.endoderm.al3.ext_v1.1, resolution = 4)

src.endoderm.al3.ext_v1.1 = 
  FindNeighbors(src.endoderm.al3.ext_v1.1, dims=1:2, 
                reduction = "umap", assay = "RNA")
src.endoderm.al3.ext_v1.1 = 
  FindClusters(src.endoderm.al3.ext_v1.1, resolution = 1)
#---------------

pdf("figure.v08.07/try.pdf")
DimPlot(src.endoderm.al3.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.2.5', reduction = "umap_integrated")
DimPlot(src.endoderm.al3.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.1', reduction = "umap_integrated")

DimPlot(src.endoderm.al3.ext_v1.1, group.by = 'Time',
        reduction = "umap_mnn", cols = colors.time)
DimPlot(src.endoderm.al3.ext_v1.1, group.by = 'Time', 
        reduction = "umap_integrated", cols = colors.time)
DimPlot(src.endoderm.al3.ext_v1.1, group.by = 'batch', 
        reduction = "umap_integrated") #, cols = colors.time)

DimPlot(src.endoderm.al3.ext_v1.1, group.by = 'cluster.extract.v1.1', 
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.al3.ext_v1.1, group.by = 'cluster.extract.v1.1', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.al3.ext_v1.1, group.by = 'cluster.v06.26.re', 
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.al3.ext_v1.1, group.by = 'cluster.v06.26.re', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
dev.off()

src.endoderm.al3.ext_v1.1@reductions$pca_integrated@cell.embeddings = 
  src.endoderm.al3.ext_v1.1@reductions$pca_integrated@cell.embeddings[colnames(src.endoderm.al3.ext_v1.1),]

#--------------------
src.endoderm.al3.ext_v1.1$cluster.v06.26.re_raw = NA
src.endoderm.al3.ext_v1.1@meta.data[
  intersect(colnames(src.endoderm.al3.ext_v1.1),colnames(src.endoderm.al3.ext.re)),]$cluster.v06.26.re_raw = 
  src.endoderm.al3.ext.re@meta.data[
    intersect(colnames(src.endoderm.al3.ext_v1.1),colnames(src.endoderm.al3.ext.re)),]$cluster.v06.26
src.endoderm.al3.ext_v1.1@meta.data[
  src.endoderm.al3.ext_v1.1$cluster.extract.v1.1%in%c("FG.4-AL.1/2/3"),]$cluster.v06.26.re_raw = "Liver"

#---------------------
src.endoderm.al3.ext_v1.1$cluster.v06.26.re = NA
src.endoderm.al3.ext_v1.1@meta.data[
  intersect(colnames(src.endoderm.al3.ext_v1.1),colnames(src.endoderm.al3.ext.re)),]$cluster.v06.26.re = 
  src.endoderm.al3.ext.re@meta.data[
    intersect(colnames(src.endoderm.al3.ext_v1.1),colnames(src.endoderm.al3.ext.re)),]$cluster.v06.26
src.endoderm.al3.ext_v1.1@meta.data[
  src.endoderm.al3.ext_v1.1$cluster.extract.v1.1%in%c("FG.4-AL.1/2/3"),]$cluster.v06.26.re = "Liver"

src.endoderm.al3.ext_v1.1@meta.data[
  src.endoderm.al3.ext_v1.1$RNA_snn_res.2.5%in%c(4),]$cluster.v06.26.re = "AL.3-Small.intestine.1"
src.endoderm.al3.ext_v1.1@meta.data[
  src.endoderm.al3.ext_v1.1$RNA_snn_res.2.5%in%c(30,13,6,14),]$cluster.v06.26.re = "Small.intestine.1"
src.endoderm.al3.ext_v1.1@meta.data[
  src.endoderm.al3.ext_v1.1$RNA_snn_res.2.5%in%c(17),]$cluster.v06.26.re = "AL.3"
src.endoderm.al3.ext_v1.1@meta.data[
  src.endoderm.al3.ext_v1.1$RNA_snn_res.2.5%in%c(18)|
    src.endoderm.al3.ext_v1.1$RNA_snn_res.1%in%c(25,20),]$cluster.v06.26.re = "AL.3-EHBD/VP"
src.endoderm.al3.ext_v1.1@meta.data[
  src.endoderm.al3.ext_v1.1$RNA_snn_res.1%in%c(11,13,26,10,23)&
    src.endoderm.al3.ext_v1.1$cluster.extract.v1.1%in%c("AL.3"),]$cluster.v06.26.re = "AL.3-Liver"

a = FNN::knn(
  src.endoderm.al3.ext_v1.1@reductions$pca_integrated@cell.embeddings[
    rownames(src.endoderm.al3.ext_v1.1@meta.data[!src.endoderm.al3.ext_v1.1$cluster.v06.26.re%in%NA,]),],
  src.endoderm.al3.ext_v1.1@reductions$pca_integrated@cell.embeddings[
    rownames(src.endoderm.al3.ext_v1.1@meta.data[src.endoderm.al3.ext_v1.1$cluster.v06.26.re%in%NA,]),],
  src.endoderm.al3.ext_v1.1@meta.data[!src.endoderm.al3.ext_v1.1$cluster.v06.26.re%in%NA,]$cluster.v06.26.re, k = 10)
src.endoderm.al3.ext_v1.1@meta.data[
  src.endoderm.al3.ext_v1.1$cluster.v06.26.re%in%NA,]$cluster.v06.26.re = as.character(a)

DimPlot(src.endoderm.al3.ext_v1.1, group.by = 'cluster.v06.26.re', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)

# FDL
src.endoderm.al3.ext_v1.1 = 
  seurat_fdl("src.endoderm.al3.ext_v1.1", "pca","RNA")
src.endoderm.al3.ext_v1.1 = 
  seurat_fdl("src.endoderm.al3.ext_v1.1", "mnn","RNA")
src.endoderm.al3.ext_v1.1 = 
  seurat_fdl("src.endoderm.al3.ext_v1.1", "pca_integrated","integrated")


# Graph
pdf("figure.v08.07/organ_development_re/al3_summary_graph.re_snn_PI.pdf",12,7)
# SNN
for(red in c('pca',"mnn","pca_integrated")){
  
  assay = ifelse(red=="mnn","RNA","integrated")
  
  graph.al3.ext_v1.1 = 
    seurat_GCN(seurat_name = "src.endoderm.al3.ext_v1.1",
               src.endoderm.al3.ext_v1.1.gene.fin,
               reduction = red, assay =  assay)
  graph.al3_cluster = cluster_walktrap(graph.al3.ext_v1.1, steps = 1)
  
  for(color in c("cluster.extract.v1.1","cluster.v06.26.re",
                 "Time")){
    for(edge.color in c(NA,"lightgray")){
      for(edge.color in c(NA,"lightgray")){
        set.seed(1)
        plot.igraph(graph.al3.ext_v1.1,
                    layout = layout_with_fr(graph.al3.ext_v1.1),
                    edge.color = edge.color,
                    vertex.size=2.5, vertex.label=NA,
                    vertex.label.color = "black", 
                    vertex.frame.color = NA, vertex.frame.width = 0.5,
                    #vertex.color =  colors.time[src.endoderm.al3.ext_v1.1[["Time"]][graph.al3_cluster$names,]],
                    #vertex.color =  cluster.endoderm.color.v5[src.endoderm.al3.ext_v1.1[["cluster.v06.26.re"]][graph.al3_cluster$names,]],
                    vertex.color =  colors_bg[src.endoderm.al3.ext_v1.1[[color]][graph.al3_cluster$names,]])
      }}}
  
}
dev.off()

pdf("figure.v08.07/organ_development_re/al3_summary_graph.re_other.pdf",12,7)
# Reduction
for(red in c("umap_integrated","umap_mnn",
             "fdl_pca",'fdl_mnn',"fdl_pca_integrated")){
  p1 = DimPlot(src.endoderm.al3.ext_v1.1,group.by = 'cluster.v06.26.re', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p2 = DimPlot(src.endoderm.al3.ext_v1.1,group.by = 'cluster.extract.v1.1', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p3 = DimPlot(src.endoderm.al3.ext_v1.1,group.by = 'Time', 
               reduction = red , cols = colors.time,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  # p4 = DimPlot(src.endoderm.al3.ext_v1.1,group.by = 'cluster.v06.26.re_hc', 
  #              reduction = red, cols = cluster.endoderm.color.v5,
  #              pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  print(p1); print(p2); print(p3);# print(p4)
}
dev.off()


pdf("figure.v08.07/organ_development_re/al3_summary_proportion.re.pdf",12,7)
#-----------------------------
cell_type = c("AL.3","AL.3-Small.intestine.1","AL.3-EHBD/VP","AL.3-Liver",
              "Small.intestine.1","EHBD","VP","Liver")
meta_al3 = src.endoderm.al3.ext_v1.1@meta.data
meta_al3$Time = factor(meta_al3$Time, levels = names(colors.time))
meta_al3$cluster.extract.v1.1 = factor(meta_al3$cluster.extract.v1.1, 
                                       levels = c("AL.3","AL.3-MG.2","AL.3-MG.1/2/3","FG.4-AL.1/2/3"))
meta_al3$cluster.v06.26.re = factor(meta_al3$cluster.v06.26.re, levels = rev(cell_type))

meta_al3_B0 = meta_al3
meta_al3_B1 = meta_al3[meta_al3$batch%in%1&!meta_al3$Time%in%"ss9",]
meta_al3_B2 = meta_al3[meta_al3$batch%in%2,]

for(data in c("meta_al3_B0","meta_al3_B1","meta_al3_B2")){
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
  
  #  print(p + 
  #          geom_bar(data = get(data), 
  #                   mapping = aes(x = Time,
  #                                 group = cluster.v06.26.re_hc, fill = cluster.v06.26.re_hc),
  #                   stat = "count", position = 'fill') + 
  #          xlab("Time") + ylab("Cluster proportion of
  # cell type")+ ggtitle(data))
  
}
#-----------------------------
dev.off()

save(src.endoderm.al3.ext_v1.1,
     file = "figure.v08.07/organ_development_re/src.endoderm.al3.ext_v1.1.re.Rdata")



# Marker
#------------------------
src.endoderm.al3.ext_v1.1 = 
  SetIdent(src.endoderm.al3.ext_v1.1,
           value = src.endoderm.al3.ext_v1.1$cluster.v06.26.re)
markergene.al3.ext_v1.1 = 
  FindAllMarkers(src.endoderm.al3.ext_v1.1, assay = "RNA")
markergene.al3.ext_v1.1 = markergene.al3.ext_v1.1[
  markergene.al3.ext_v1.1$avg_log2FC>0.2&
    markergene.al3.ext_v1.1$p_val_adj<0.1, "gene"]

#------------------------
al3.umap.embedding = src.endoderm.al3.ext_v1.1[["fdl_pca_integrated"]]@cell.embeddings[,1:2]
al3.umap.embedding = src.endoderm.al3.ext_v1.1[["fdl_mnn"]]@cell.embeddings[,1:2]
al3.umap.embedding = src.endoderm.al3.ext_v1.1[["umap_mnn"]]@cell.embeddings[,1:3]
al3.umap.embedding = src.endoderm.al3.ext_v1.1[["umap_integrated"]]@cell.embeddings[,1:3]

# Cell order :: Princurve
#---------------------------
# AL.3
cellorder.al3_al3 = 
  rownames(src.endoderm.al3.ext_v1.1@meta.data[
    src.endoderm.al3.ext_v1.1$cluster.v06.26.re%in%"AL.3",])
princurve.al3_al3 = princurve::principal_curve(
  x = al3.umap.embedding[cellorder.al3_al3,],  smoother = "smooth.spline")
cellorder.al3_al3 =  names(
  princurve.al3_al3$lambda[cellorder.al3_al3][order(princurve.al3_al3$lambda[cellorder.al3_al3])])
# AL.3-Small.intestine.1
cellorder.al3_al3s1 = 
  rownames(src.endoderm.al3.ext_v1.1@meta.data[
    src.endoderm.al3.ext_v1.1$cluster.v06.26.re%in%"AL.3-Small.intestine.1",])
princurve.al3_al3s1 = princurve::principal_curve(
  x = al3.umap.embedding[cellorder.al3_al3s1,],  smoother = "smooth.spline")
cellorder.al3_al3s1 =  names(
  princurve.al3_al3s1$lambda[cellorder.al3_al3s1][order(princurve.al3_al3s1$lambda[cellorder.al3_al3s1])])
# AL.3-EHBD/VP
cellorder.al3_al3ev = 
  rownames(src.endoderm.al3.ext_v1.1@meta.data[
    src.endoderm.al3.ext_v1.1$cluster.v06.26.re%in%"AL.3-EHBD/VP",])
princurve.al3_al3ev = princurve::principal_curve(
  x = al3.umap.embedding[cellorder.al3_al3ev,],  smoother = "smooth.spline")
cellorder.al3_al3ev =  names(
  princurve.al3_al3ev$lambda[cellorder.al3_al3ev][order(princurve.al3_al3ev$lambda[cellorder.al3_al3ev])])
# AL.3-Liver
cellorder.al3_al3lv = 
  rownames(src.endoderm.al3.ext_v1.1@meta.data[
    src.endoderm.al3.ext_v1.1$cluster.v06.26.re%in%"AL.3-Liver",])
princurve.al3_al3lv = princurve::principal_curve(
  x = al3.umap.embedding[cellorder.al3_al3lv,],  smoother = "smooth.spline")
cellorder.al3_al3lv =  names(
  princurve.al3_al3lv$lambda[cellorder.al3_al3lv][order(princurve.al3_al3lv$lambda[cellorder.al3_al3lv])])
# Small.intestine.1
cellorder.al3_s1 = 
  rownames(src.endoderm.al3.ext_v1.1@meta.data[
    src.endoderm.al3.ext_v1.1$cluster.v06.26.re%in%"Small.intestine.1",])
princurve.al3_s1 = princurve::principal_curve(
  x = al3.umap.embedding[cellorder.al3_s1,],  smoother = "smooth.spline")
cellorder.al3_s1 =  names(
  princurve.al3_s1$lambda[cellorder.al3_s1][order(princurve.al3_s1$lambda[cellorder.al3_s1])])
# EHBD
cellorder.al3_eh = 
  rownames(src.endoderm.al3.ext_v1.1@meta.data[
    src.endoderm.al3.ext_v1.1$cluster.v06.26.re%in%"EHBD",])
princurve.al3_eh = princurve::principal_curve(
  x = al3.umap.embedding[cellorder.al3_eh,],  smoother = "smooth.spline")
cellorder.al3_eh =  names(
  princurve.al3_eh$lambda[cellorder.al3_eh][order(princurve.al3_eh$lambda[cellorder.al3_eh])])
# VP
cellorder.al3_vp = 
  rownames(src.endoderm.al3.ext_v1.1@meta.data[
    src.endoderm.al3.ext_v1.1$cluster.v06.26.re%in%"VP",])
princurve.al3_vp = princurve::principal_curve(
  x = al3.umap.embedding[cellorder.al3_vp,],  smoother = "smooth.spline")
cellorder.al3_vp =  names(
  princurve.al3_vp$lambda[cellorder.al3_vp][order(princurve.al3_vp$lambda[cellorder.al3_vp])])
# Liver
cellorder.al3_lv = 
  rownames(src.endoderm.al3.ext_v1.1@meta.data[
    src.endoderm.al3.ext_v1.1$cluster.v06.26.re%in%"Liver",])
princurve.al3_lv = princurve::principal_curve(
  x = al3.umap.embedding[cellorder.al3_lv,],  smoother = "smooth.spline")
cellorder.al3_lv =  names(
  princurve.al3_lv$lambda[cellorder.al3_lv][order(princurve.al3_lv$lambda[cellorder.al3_lv])])
#---------------------------

pdf("figure.v08.07/try.pdf")
#---------------------------
plot(c(1:length(cellorder.al3_al3)),
     c(1:length(cellorder.al3_al3)),
     col = colors.time[src.endoderm.al3.ext_v1.1[,cellorder.al3_al3]$Time])
plot(c(1:length(cellorder.al3_al3s1)),
     c(1:length(cellorder.al3_al3s1)),
     col = colors.time[src.endoderm.al3.ext_v1.1[,cellorder.al3_al3s1]$Time])
plot(c(1:length(cellorder.al3_al3ev)),
     c(1:length(cellorder.al3_al3ev)),
     col = colors.time[src.endoderm.al3.ext_v1.1[,cellorder.al3_al3ev]$Time])
plot(c(1:length(cellorder.al3_al3lv)),
     c(1:length(cellorder.al3_al3lv)),
     col = colors.time[src.endoderm.al3.ext_v1.1[,cellorder.al3_al3lv]$Time])
plot(c(1:length(cellorder.al3_s1)),
     c(1:length(cellorder.al3_s1)),
     col = colors.time[src.endoderm.al3.ext_v1.1[,cellorder.al3_s1]$Time])
plot(c(1:length(cellorder.al3_eh)),
     c(1:length(cellorder.al3_eh)),
     col = colors.time[src.endoderm.al3.ext_v1.1[,cellorder.al3_eh]$Time])
plot(c(1:length(cellorder.al3_vp)),
     c(1:length(cellorder.al3_vp)),
     col = colors.time[src.endoderm.al3.ext_v1.1[,cellorder.al3_vp]$Time])
plot(c(1:length(cellorder.al3_lv)),
     c(1:length(cellorder.al3_lv)),
     col = colors.time[src.endoderm.al3.ext_v1.1[,cellorder.al3_lv]$Time])
#---------------------------
dev.off()


cellorder.al3_al3_FPI = cellorder.al3_al3
cellorder.al3_al3s1_FPI = cellorder.al3_al3s1
cellorder.al3_al3ev_FPI = cellorder.al3_al3ev
cellorder.al3_al3lv_FPI = cellorder.al3_al3lv
cellorder.al3_s1_FPI = cellorder.al3_s1
#cellorder.al3_s1_FPI = rev(cellorder.al3_s1_FPI)
cellorder.al3_eh_FPI = cellorder.al3_eh
cellorder.al3_vp_FPI = cellorder.al3_vp
cellorder.al3_lv_FPI = cellorder.al3_lv

#=================================================================================
#=================================================================================

pdf("figure.v08.07/organ_development_re/al3_heatmap_marker.re_s1.pdf",6,9)
cellorder.al3 = c(cellorder.al3_al3_FPI, cellorder.al3_al3s1_FPI , cellorder.al3_s1_FPI)
markergene.al3.ext_v1.1 = 
  FindAllMarkers(src.endoderm.al3.ext_v1.1[,cellorder.al3], assay = "RNA")
markergene.al3.ext_v1.1 = markergene.al3.ext_v1.1[
  markergene.al3.ext_v1.1$avg_log2FC>0.2&markergene.al3.ext_v1.1$p_val_adj<0.1, "gene"]
# selectgene = unique(c(markergene.al3.ext_v1.1))
selectgene = src.endoderm.al3.ext_v1.1.filtergene_s1_fin
#-------------------
al3.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.al3.ext_v1.1@assays$RNA@data[selectgene, cellorder.al3]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.al3.ext_v1.1$Time[cellorder.al3], colors.time),
              MyName2Col(src.endoderm.al3.ext_v1.1$cluster.v06.26.re[cellorder.al3], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.al3.ext_v1.1$cluster.extract.v1.1[cellorder.al3], cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            #return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.al3.ext_v1.1@meta.data[cell_al3_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

al3.ext_v1.1.re.rowtree.1 = as.dendrogram(al3.ext_v1.1.re.rowtree)
src.endoderm.al3.ext_v1.1.filtergene_fin = c(
  labels(al3.ext_v1.1.re.rowtree.1[[2]][[2]][[2]]),
  labels(al3.ext_v1.1.re.rowtree.1[[2]][[2]][[1]][[2]][[2]][[2]]),
  labels(al3.ext_v1.1.re.rowtree.1[[2]][[2]][[1]][[2]][[1]]),
  labels(al3.ext_v1.1.re.rowtree.1[[2]][[2]][[1]][[1]]),
  
  labels(al3.ext_v1.1.re.rowtree.1[[2]][[1]][[1]]),
  labels(al3.ext_v1.1.re.rowtree.1[[2]][[1]][[2]][[1]]),
  labels(al3.ext_v1.1.re.rowtree.1[[2]][[1]][[2]][[2]]),
  
  rev(labels(al3.ext_v1.1.re.rowtree.1[[1]]))
  )

names(src.endoderm.al3.ext_v1.1.filtergene_fin) = c(
  rep(4,length(c(
    labels(al3.ext_v1.1.re.rowtree.1[[2]][[2]][[2]]),
    labels(al3.ext_v1.1.re.rowtree.1[[2]][[2]][[1]][[2]][[2]][[2]]),
    labels(al3.ext_v1.1.re.rowtree.1[[2]][[2]][[1]][[2]][[1]]),
    labels(al3.ext_v1.1.re.rowtree.1[[2]][[2]][[1]][[1]])
  ))),
  rep(3,length(c(
    labels(al3.ext_v1.1.re.rowtree.1[[2]][[1]][[1]]),
    labels(al3.ext_v1.1.re.rowtree.1[[2]][[1]][[2]][[1]]),
    labels(al3.ext_v1.1.re.rowtree.1[[2]][[1]][[2]][[2]])
  ))),
  rep(7,length(c(
    rev(labels(al3.ext_v1.1.re.rowtree.1[[1]]))
  )))
)
src.endoderm.al3.ext_v1.1.filtergene_s1_fin =
  src.endoderm.al3.ext_v1.1.filtergene_fin


pdf("figure.v08.07/organ_development_re/al3_heatmap_marker.re_ev.pdf",6,9)
cellorder.al3 = c(cellorder.al3_al3_FPI, cellorder.al3_al3ev_FPI , 
                  cellorder.al3_vp_FPI,cellorder.al3_eh_FPI)
markergene.al3.ext_v1.1 = 
  FindAllMarkers(src.endoderm.al3.ext_v1.1[,cellorder.al3], assay = "RNA")
markergene.al3.ext_v1.1 = markergene.al3.ext_v1.1[
  markergene.al3.ext_v1.1$avg_log2FC>0.2&markergene.al3.ext_v1.1$p_val_adj<0.1, "gene"]
# selectgene = unique(c(markergene.al3.ext_v1.1))
selectgene = src.endoderm.al3.ext_v1.1.filtergene_ev_fin
#-------------------
al3.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.al3.ext_v1.1@assays$RNA@data[selectgene, cellorder.al3]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.al3.ext_v1.1$Time[cellorder.al3], colors.time),
              MyName2Col(src.endoderm.al3.ext_v1.1$cluster.v06.26.re[cellorder.al3], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.al3.ext_v1.1$cluster.extract.v1.1[cellorder.al3], cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            #return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.al3.ext_v1.1@meta.data[cell_al3_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

al3.ext_v1.1.re.rowtree.2 = as.dendrogram(al3.ext_v1.1.re.rowtree)

src.endoderm.al3.ext_v1.1.filtergene_fin = c(
  setdiff(
    c(labels(al3.ext_v1.1.re.rowtree.2[[1]][[1]])),
    c("Dstn","Csrp2","Tpm1","Id1","Cdk2ap1")
  ),
  setdiff(
    c(
      labels(al3.ext_v1.1.re.rowtree.2[[1]][[2]][[1]]),
      labels(al3.ext_v1.1.re.rowtree.2[[1]][[2]][[2]][[1]][[1]]),
      labels(al3.ext_v1.1.re.rowtree.2[[1]][[2]][[2]][[1]][[2]][[1]]),
      labels(al3.ext_v1.1.re.rowtree.2[[1]][[2]][[2]][[2]][[1]][[2]]),
      labels(al3.ext_v1.1.re.rowtree.2[[1]][[2]][[2]][[2]][[2]])
    ),
    c()
  ),
  rev(setdiff(
    c(labels(al3.ext_v1.1.re.rowtree.2[[2]][[2]])),
    c("Parm1","Ptf1a","Cpa1")
  )),
  setdiff(
    c(labels(al3.ext_v1.1.re.rowtree.2[[2]][[1]])),
    c()
  )
)

names(src.endoderm.al3.ext_v1.1.filtergene_fin) = c(
  rep(4,length(  
    setdiff(
      c(labels(al3.ext_v1.1.re.rowtree.2[[1]][[1]])),
      c("Dstn","Csrp2","Tpm1","Id1","Cdk2ap1")
    ))),
  rep(3,length(  
    setdiff(
      c(
        labels(al3.ext_v1.1.re.rowtree.2[[1]][[2]][[1]]),
        labels(al3.ext_v1.1.re.rowtree.2[[1]][[2]][[2]][[1]][[1]]),
        labels(al3.ext_v1.1.re.rowtree.2[[1]][[2]][[2]][[1]][[2]][[1]]),
        labels(al3.ext_v1.1.re.rowtree.2[[1]][[2]][[2]][[2]][[1]][[2]]),
        labels(al3.ext_v1.1.re.rowtree.2[[1]][[2]][[2]][[2]][[2]])
      ),
      c()
    ))),
  rep(9,length(  
    setdiff(
      c(labels(al3.ext_v1.1.re.rowtree.2[[2]][[2]])),
      c("Parm1","Ptf1a","Cpa1")
    ))),
  rep(7,length(  
    setdiff(
      c(labels(al3.ext_v1.1.re.rowtree.2[[2]][[1]])),
      c()
    )))
)
src.endoderm.al3.ext_v1.1.filtergene_ev_fin =
  src.endoderm.al3.ext_v1.1.filtergene_fin


pdf("figure.v08.07/organ_development_re/al3_heatmap_marker.re_lv.pdf",6,9)
cellorder.al3 = c(cellorder.al3_al3_FPI, cellorder.al3_al3lv_FPI , 
                  cellorder.al3_lv_FPI)
markergene.al3.ext_v1.1 = 
  FindAllMarkers(src.endoderm.al3.ext_v1.1[,cellorder.al3], assay = "RNA")
markergene.al3.ext_v1.1 = markergene.al3.ext_v1.1[
  markergene.al3.ext_v1.1$avg_log2FC>0.2&markergene.al3.ext_v1.1$p_val_adj<0.1, "gene"]
# selectgene = unique(c(markergene.al3.ext_v1.1))
selectgene = src.endoderm.al3.ext_v1.1.filtergene_lv_fin
#-------------------
al3.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.al3.ext_v1.1@assays$RNA@data[selectgene, cellorder.al3]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.al3.ext_v1.1$Time[cellorder.al3], colors.time),
              MyName2Col(src.endoderm.al3.ext_v1.1$cluster.v06.26.re[cellorder.al3], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.al3.ext_v1.1$cluster.extract.v1.1[cellorder.al3], cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            #return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.al3.ext_v1.1@meta.data[cell_al3_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

al3.ext_v1.1.re.rowtree.3 = as.dendrogram(al3.ext_v1.1.re.rowtree)

src.endoderm.al3.ext_v1.1.filtergene_fin = c(
  setdiff(
    c(
      labels(al3.ext_v1.1.re.rowtree.3[[2]][[2]][[2]][[1]]),
      labels(al3.ext_v1.1.re.rowtree.3[[2]][[2]][[1]])
    ),
    c()
  ),
  setdiff(
    c( labels(al3.ext_v1.1.re.rowtree.3[[2]][[1]][[1]]),
      rev(labels(al3.ext_v1.1.re.rowtree.3[[2]][[1]][[2]][[1]])),
      labels(al3.ext_v1.1.re.rowtree.3[[2]][[1]][[2]][[2]][[1]]),
      labels(al3.ext_v1.1.re.rowtree.3[[2]][[1]][[2]][[2]][[2]][[1]])
    ),
    c()
  ),
  setdiff(
    c(
      labels(al3.ext_v1.1.re.rowtree.3[[1]][[1]][[1]]),
      labels(al3.ext_v1.1.re.rowtree.3[[1]][[1]]),
      rev(labels(al3.ext_v1.1.re.rowtree.3[[1]][[2]][[2]]))
    ),
    c()
  )
)

names(src.endoderm.al3.ext_v1.1.filtergene_fin) = c(
  rep(4,length(  
    setdiff(
      c(
        labels(al3.ext_v1.1.re.rowtree.3[[2]][[2]][[2]][[1]]),
        labels(al3.ext_v1.1.re.rowtree.3[[2]][[2]][[1]])
        #labels(al3.ext_v1.1.re.rowtree.3[[2]][[1]][[1]])
      ),
      c()
    ))),
  rep(3,length(  
    setdiff(
      c( labels(al3.ext_v1.1.re.rowtree.3[[2]][[1]][[1]]),
        rev(labels(al3.ext_v1.1.re.rowtree.3[[2]][[1]][[2]][[1]])),
        labels(al3.ext_v1.1.re.rowtree.3[[2]][[1]][[2]][[2]][[1]]),
        labels(al3.ext_v1.1.re.rowtree.3[[2]][[1]][[2]][[2]][[2]][[1]])
      ),
      c()
    ))),
  rep(7,length(  
    setdiff(
      c(
        labels(al3.ext_v1.1.re.rowtree.3[[1]][[1]][[1]]),
        labels(al3.ext_v1.1.re.rowtree.3[[1]][[1]]),
        rev(labels(al3.ext_v1.1.re.rowtree.3[[1]][[2]][[2]]))
      ),
      c()
    )))
)
src.endoderm.al3.ext_v1.1.filtergene_lv_fin =
  src.endoderm.al3.ext_v1.1.filtergene_fin


save(
  # Gene 
  src.endoderm.al3.ext.v1.1.selectgene,
  src.endoderm.al3.ext_v1.1.gene.fin,
  src.endoderm.al3.ext_v1.1.filtergene_s1_fin,
  src.endoderm.al3.ext_v1.1.filtergene_ev_fin,
  src.endoderm.al3.ext_v1.1.filtergene_lv_fin,
  # Cell type
  cellorder.al3_al3_FPI, cellorder.al3_al3s1_FPI , cellorder.al3_s1_FPI,
  cellorder.al3_al3ev_FPI, cellorder.al3_eh_FPI , cellorder.al3_vp_FPI,
  cellorder.al3_al3lv_FPI, cellorder.al3_lv_FPI , 
  file= "figure.v08.07/organ_development_re/al3_heatmap_parameter.Rdata")

