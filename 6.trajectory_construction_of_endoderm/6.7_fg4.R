#----------------------
#>>  FG.4, FG.4-sub
#----------------------

src.endoderm.fg4.ext_v1.1 = FindVariableFeatures(src.endoderm.fg4.ext_v1.1, nfeatures = 2000)
src.endoderm.fg4.ext_v1.1 = ScaleData(src.endoderm.fg4.ext_v1.1, split.by = "batch_phase",
                                      features = rownames(src.endoderm.fg4.ext_v1.1))
src.endoderm.fg4.ext_v1.1.filtergene = 
  Myfilter(as.matrix(src.endoderm.fg4.ext_v1.1@assays$RNA@data),
           gene = src.endoderm.fg4.ext_v1.1@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
#rm(src.endoderm.fg4.ext_v1.1.filtergene)
src.endoderm.fg4.ext.v1.1.selectgene = unique(c(
  select_gene_fg4,src.endoderm.fg4.ext_v1.1.filtergene))

#=============================
# Row
pdf("figure.v08.07/try.pdf")
#-----------------
fg4.ext.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.fg4.ext_v1.1@assays$RNA@data[
    src.endoderm.fg4.ext.v1.1.selectgene,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.endoderm.fg4.ext_v1.1$Time,colors.time),
      MyName2Col(src.endoderm.fg4.ext_v1.1$cluster.extract.v1.1,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.fg4.ext_v1.1$cluster.v06.26.re,cluster.endoderm.color.v5)
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
plot(fg4.ext.re.rowtree_1)
dev.off()

fg4.ext.re.rowtree_1 = as.dendrogram(fg4.ext.re.rowtree)
src.endoderm.fg4.ext_v1.1.gene.fin = setdiff(
  src.endoderm.fg4.ext.v1.1.selectgene,
  c(labels(fg4.ext.re.rowtree_1[[2]][[1]]),
    labels(fg4.ext.re.rowtree_1[[2]][[2]][[2]][[1]])))

#=============================
# Col
pdf("figure.v08.07/try.pdf")
#-----------------
fg4.ext.re.coltree =
  MyHeatmap(as.matrix(src.endoderm.fg4.ext_v1.1@assays$RNA@data[
    src.endoderm.fg4.ext_v1.1.gene.fin,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.endoderm.fg4.ext_v1.1$Time,colors.time),
      MyName2Col(src.endoderm.fg4.ext_v1.1$cluster.extract.v1.1,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.fg4.ext_v1.1$cluster.v06.26.re,cluster.endoderm.color.v5)
    ),
    ColSideColorsSize = 1.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "col",
    graph = T)
#-------------------
dev.off()

fg4.ext.re.coltree_1 = as.dendrogram(fg4.ext.re.coltree)

src.endoderm.fg4.ext_v1.1$cluster.v06.26.re_hc = NA
src.endoderm.fg4.ext_v1.1@meta.data[c(
  labels(fg4.ext.re.coltree_1[[1]][[2]][[2]][[1]])), 
  "cluster.v06.26.re_hc"] = "FG.3"

src.endoderm.fg4.ext_v1.1@meta.data[c(
  labels(fg4.ext.re.coltree_1[[1]][[2]][[2]][[2]][[2]])), 
  "cluster.v06.26.re_hc"] = "FG.3"

src.endoderm.fg4.ext_v1.1@meta.data[setdiff(
  labels(fg4.ext.re.coltree_1[[1]]),
  c(  
    labels(fg4.ext.re.coltree_1[[1]][[2]][[2]][[1]]),
    labels(fg4.ext.re.coltree_1[[1]][[2]][[2]][[2]][[2]])
  )), "cluster.v06.26.re_hc"] = "Pharynx.organ.4"

src.endoderm.fg4.ext_v1.1@meta.data[c(
  labels(fg4.ext.re.coltree_1[[2]][[2]][[2]][[1]])), 
  "cluster.v06.26.re_hc"] = "FG.3-Lung"

src.endoderm.fg4.ext_v1.1@meta.data[setdiff(
  labels(fg4.ext.re.coltree_1[[2]]),
  c(  
    labels(fg4.ext.re.coltree_1[[2]][[2]][[2]][[1]])
  )), "cluster.v06.26.re_hc"] = "Lung"


#--------------------------------------------------------------------------------
# src.endoderm.fg4.ext.re_v1.1.gene.fin =
#   rownames(src.endoderm.fg4.ext.re@reductions$pca@feature.loadings)

src.endoderm.fg4.ext_v1.1 = 
  seurat_mnn(seurat_name = "src.endoderm.fg4.ext_v1.1", dims = 1:30, 
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.fg4.ext_v1.1.gene.fin)
src.endoderm.fg4.ext_v1.1 = 
  seurat_int(seurat_name = "src.endoderm.fg4.ext_v1.1", dims = 1:30,
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.fg4.ext_v1.1.gene.fin)


# UMAP Integrated Rotation
#-------------------------------------------------
src.endoderm.fg4.ext_v1.1@reductions$umap_int_rotated = 
  src.endoderm.fg4.ext_v1.1@reductions$umap
src.endoderm.fg4.ext_v1.1@reductions$umap_int_rotated@cell.embeddings = 
  fg4_umap_integrated_rotatedEmbedding[colnames(src.endoderm.fg4.ext_v1.1),1:3]
src.endoderm.fg4.ext_v1.1@reductions$umap_int_rotated@feature.loadings.projected = 
  fg4_umap_integrated_rotation
colnames(src.endoderm.fg4.ext_v1.1@reductions$umap_int_rotated@cell.embeddings) = 
  c("UMAP_1","UMAP_2","UMAP_3")
#-------------------------------------------------


#---------------
src.endoderm.fg4.ext_v1.1 = 
  FindNeighbors(src.endoderm.fg4.ext_v1.1, dims=1:30, 
                reduction = "pca_integrated", assay = "integrated")
src.endoderm.fg4.ext_v1.1 = 
  FindNeighbors(src.endoderm.fg4.ext_v1.1, dims=1:30, 
                reduction = "pca", assay = "RNA")
src.endoderm.fg4.ext_v1.1 = 
  FindClusters(src.endoderm.fg4.ext_v1.1, 
               resolution = 2, graph.name = "RNA_snn",)

src.endoderm.fg4.ext_v1.1 = 
  FindNeighbors(src.endoderm.fg4.ext_v1.1, dims=1:3, 
                reduction = "umap_mnn", assay = "RNA")
src.endoderm.fg4.ext_v1.1 = 
  FindClusters(src.endoderm.fg4.ext_v1.1, 
               resolution = 1.5, graph.name = "RNA_snn",)
src.endoderm.fg4.ext_v1.1$cluster_snn_res.1.5 = 
  src.endoderm.fg4.ext_v1.1$RNA_snn_res.1.5
#---------------

DimPlot(src.endoderm.fg4.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.2', reduction = "umap_mnn")
DimPlot(src.endoderm.fg4.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'cluster_snn_res.1.5', reduction = "umap_mnn")
DimPlot(src.endoderm.fg4.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'cluster_snn_res.1', reduction = "umap_mnn")

DimPlot(src.endoderm.fg4.ext_v1.1, group.by = 'Time',
        reduction = "umap_mnn", cols = colors.time)
DimPlot(src.endoderm.fg4.ext_v1.1, group.by = 'Time', 
        reduction = "umap_integrated", cols = colors.time)
DimPlot(src.endoderm.fg4.ext_v1.1, group.by = 'batch', 
        reduction = "umap_integrated") #, cols = colors.time)

DimPlot(src.endoderm.fg4.ext_v1.1, group.by = 'cluster.extract.v1.1', 
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.fg4.ext_v1.1, group.by = 'cluster.extract.v1.1', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.fg4.ext_v1.1, group.by = 'cluster.v06.26.re', 
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.fg4.ext_v1.1, group.by = 'cluster.v06.26.re_v1.0', 
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.fg4.ext_v1.1, group.by = 'cluster.v06.26.re_hc',
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.fg4.ext_v1.1, group.by = 'cluster.v06.26.re_hc',
        reduction = "umap_int_rotated", cols = cluster.endoderm.color.v5)


meta = cbind(
  t(src.endoderm.fg4.ext_v1.1@assays$RNA["Hoxb1",]),
  src.endoderm.fg4.ext_v1.1@reductions$umap_mnn@cell.embeddings)
meta = as.data.frame(meta)
ggplot()+
  geom_point(data = meta,
             mapping = aes(x=UMAP_1,y=UMAP_2,color=Hoxb1))




src.endoderm.fg4.ext_v1.1$cluster.v06.26.re_v1.0 = NA
src.endoderm.fg4.ext_v1.1@meta.data[
  intersect(colnames(src.endoderm.fg4.ext_v1.1), colnames(src.endoderm.fg4.ext.re)),]$cluster.v06.26.re_v1.0  =
  src.endoderm.fg4.ext.re@meta.data[
    intersect(colnames(src.endoderm.fg4.ext_v1.1), colnames(src.endoderm.fg4.ext.re)),]$cluster.v06.26
a = FNN::knn(
  src.endoderm.fg4.ext_v1.1@reductions$umap_integrated@cell.embeddings[
    rownames(src.endoderm.fg4.ext_v1.1@meta.data[!src.endoderm.fg4.ext_v1.1$cluster.v06.26.re_v1.0%in%NA,]),],
  src.endoderm.fg4.ext_v1.1@reductions$umap_integrated@cell.embeddings[
    rownames(src.endoderm.fg4.ext_v1.1@meta.data[src.endoderm.fg4.ext_v1.1$cluster.v06.26.re_v1.0%in%NA,]),],
  src.endoderm.fg4.ext_v1.1@meta.data[!src.endoderm.fg4.ext_v1.1$cluster.v06.26.re_v1.0%in%NA,]$cluster.v06.26.re_v1.0, k = 10)
src.endoderm.fg4.ext_v1.1@meta.data[
  src.endoderm.fg4.ext_v1.1$cluster.v06.26.re_v1.0%in%NA,]$cluster.v06.26.re_v1.0 = as.character(a)

src.endoderm.fg4.ext_v1.1$cluster.v06.26.re = NA
src.endoderm.fg4.ext_v1.1@meta.data[
  src.endoderm.fg4.ext_v1.1$RNA_snn_res.2%in%c(12,4)|
    src.endoderm.fg4.ext_v1.1$cluster.extract.v1.1%in%c("FG.4-MG.1/3"),]$cluster.v06.26.re = "Stomach"
src.endoderm.fg4.ext_v1.1@meta.data[
  src.endoderm.fg4.ext_v1.1$RNA_snn_res.2%in%c(8,16,18)&
    src.endoderm.fg4.ext_v1.1$cluster.extract.v1.1%in%c("FG.4"),]$cluster.v06.26.re = "FG.4"
src.endoderm.fg4.ext_v1.1@meta.data[
  src.endoderm.fg4.ext_v1.1$RNA_snn_res.2%in%c(24,15,13,17,7,10,0)&
    src.endoderm.fg4.ext_v1.1$cluster.extract.v1.1%in%c("FG.4"),]$cluster.v06.26.re = "FG.4-Lung/Stomach"
src.endoderm.fg4.ext_v1.1@meta.data[
  src.endoderm.fg4.ext_v1.1$RNA_snn_res.2%in%c(7,10,17)&
    src.endoderm.fg4.ext_v1.1$cluster.extract.v1.1%in%c("FG.3/4"),]$cluster.v06.26.re = "Lung"
src.endoderm.fg4.ext_v1.1@meta.data[
  src.endoderm.fg4.ext_v1.1$RNA_snn_res.2%in%c(0,25,11,5)&
    src.endoderm.fg4.ext_v1.1$cluster.extract.v1.1%in%c("FG.3/4"),]$cluster.v06.26.re = "Lung" # "Trachea"
src.endoderm.fg4.ext_v1.1@meta.data[
  src.endoderm.fg4.ext_v1.1$cluster_snn_res.1%in%c(20)&
    !src.endoderm.fg4.ext_v1.1$RNA_snn_res.2%in%c(8)&
    src.endoderm.fg4.ext_v1.1$cluster.extract.v1.1%in%"FG.4",]$cluster.v06.26.re = "FG.4-Liver"
# src.endoderm.fg4.ext_v1.1@meta.data[
#   src.endoderm.fg4.ext_v1.1$cluster.v06.26.re_v1.0%in%"FG.4-Liver"&
#     src.endoderm.fg4.ext_v1.1$cluster.extract.v1.1%in%"FG.4",]$cluster.v06.26.re = "FG.4-Liver"
src.endoderm.fg4.ext_v1.1@meta.data[
    src.endoderm.fg4.ext_v1.1$cluster.extract.v1.1%in%c("FG.4-AL.1/2/3"),]$cluster.v06.26.re = "Liver"

a = FNN::knn(
  src.endoderm.fg4.ext_v1.1@reductions$umap_integrated@cell.embeddings[
    rownames(src.endoderm.fg4.ext_v1.1@meta.data[!src.endoderm.fg4.ext_v1.1$cluster.v06.26.re%in%NA,]),],
  src.endoderm.fg4.ext_v1.1@reductions$umap_integrated@cell.embeddings[
    rownames(src.endoderm.fg4.ext_v1.1@meta.data[src.endoderm.fg4.ext_v1.1$cluster.v06.26.re%in%NA,]),],
  src.endoderm.fg4.ext_v1.1@meta.data[!src.endoderm.fg4.ext_v1.1$cluster.v06.26.re%in%NA,]$cluster.v06.26.re, k = 10)
src.endoderm.fg4.ext_v1.1@meta.data[
  src.endoderm.fg4.ext_v1.1$cluster.v06.26.re%in%NA,]$cluster.v06.26.re = as.character(a)

src.endoderm.fg4.ext_v1.1@reductions$umap_mnn@cell.embeddings = 
  src.endoderm.fg4.ext_v1.1@reductions$umap_mnn@cell.embeddings[colnames(src.endoderm.fg4.ext_v1.1),]
src.endoderm.fg4.ext_v1.1@meta.data[
  src.endoderm.fg4.ext_v1.1@reductions$umap_mnn@cell.embeddings[,1]>(-1)&
    src.endoderm.fg4.ext_v1.1$cluster.v06.26.re%in%c("FG.4","FG.4-Lung/Stomach"),]$cluster.v06.26.re = "FG.4-Liver"

DimPlot(src.endoderm.fg4.ext_v1.1, group.by = 'cluster.v06.26.re', 
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)


# FDL
src.endoderm.fg4.ext_v1.1 = 
  seurat_fdl("src.endoderm.fg4.ext_v1.1", "pca","RNA")
src.endoderm.fg4.ext_v1.1 = 
  seurat_fdl("src.endoderm.fg4.ext_v1.1", "mnn","RNA")
src.endoderm.fg4.ext_v1.1 = 
  seurat_fdl("src.endoderm.fg4.ext_v1.1", "pca_integrated","integrated")


# Graph
pdf("figure.v08.07/organ_development_re/fg4_summary_graph.re_snn_PI.pdf",12,7)
# SNN
for(red in c("mnn","pca_integrated")){
  
  assay = ifelse(red=="mnn","RNA","integrated")
  
  graph.fg4.ext_v1.1 = 
    seurat_GCN(seurat_name = "src.endoderm.fg4.ext_v1.1",
               src.endoderm.fg4.ext_v1.1.gene.fin,
               reduction = red, assay =  assay)
  graph.fg4_cluster = cluster_walktrap(graph.fg4.ext_v1.1, steps = 1)
  
  for(color in c("cluster.extract.v1.1","cluster.v06.26.re",
                 "Time")){
    for(edge.color in c(NA,"lightgray")){
      for(edge.color in c(NA,"lightgray")){
        set.seed(1)
        plot.igraph(graph.fg4.ext_v1.1,
                    layout = layout_with_fr(graph.fg4.ext_v1.1),
                    edge.color = edge.color,
                    vertex.size=2.5, vertex.label=NA,
                    vertex.label.color = "black", 
                    vertex.frame.color = NA, vertex.frame.width = 0.5,
                    #vertex.color =  colors.time[src.endoderm.fg4.ext_v1.1[["Time"]][graph.fg4_cluster$names,]],
                    #vertex.color =  cluster.endoderm.color.v5[src.endoderm.fg4.ext_v1.1[["cluster.v06.26.re"]][graph.fg4_cluster$names,]],
                    vertex.color =  colors_bg[src.endoderm.fg4.ext_v1.1[[color]][graph.fg4_cluster$names,]])
      }}}

}
dev.off()

pdf("figure.v08.07/organ_development_re/fg4_summary_graph.re_other.pdf",12,7)
# Reduction
for(red in c("umap_integrated","umap_int_rotated","umap_mnn",
             "fdl_pca",'fdl_mnn',"fdl_pca_integrated")){
  p1 = DimPlot(src.endoderm.fg4.ext_v1.1,group.by = 'cluster.v06.26.re', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p2 = DimPlot(src.endoderm.fg4.ext_v1.1,group.by = 'cluster.extract.v1.1', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p3 = DimPlot(src.endoderm.fg4.ext_v1.1,group.by = 'Time', 
               reduction = red , cols = colors.time,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  # p4 = DimPlot(src.endoderm.fg4.ext_v1.1,group.by = 'cluster.v06.26.re_hc', 
  #              reduction = red, cols = cluster.endoderm.color.v5,
  #              pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  print(p1); print(p2); print(p3);# print(p4)
}
dev.off()


pdf("figure.v08.07/organ_development_re/fg4_summary_proportion.re.pdf",12,7)
#-----------------------------
cell_type = c("FG.4","FG.4-Lung/Stomach","FG.4-Liver","Lung","Stomach","Liver")
meta_fg4 = src.endoderm.fg4.ext_v1.1@meta.data
meta_fg4$Time = factor(meta_fg4$Time, levels = names(colors.time))
meta_fg4$cluster.extract.v1.1 = factor(meta_fg4$cluster.extract.v1.1, 
                                       levels = c("FG.4","FG.3/4","FG.4-MG.1/3",'FG.4-AL.1/2/3'))
meta_fg4$cluster.v06.26.re = factor(meta_fg4$cluster.v06.26.re, levels = rev(cell_type))
#meta_fg4$cluster.v06.26.re_hc = factor(meta_fg4$cluster.v06.26.re_hc, levels = rev(cell_type))

meta_fg4_B0 = meta_fg4
meta_fg4_B1 = meta_fg4[meta_fg4$batch%in%1&!meta_fg4$Time%in%"ss9",]
meta_fg4_B2 = meta_fg4[meta_fg4$batch%in%2,]

for(data in c("meta_fg4_B0","meta_fg4_B1","meta_fg4_B2")){
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

save(src.endoderm.fg4.ext_v1.1,
     file = "figure.v08.07/organ_development_re/src.endoderm.fg4.ext_v1.1.re.Rdata")


# Marker
#------------------------
src.endoderm.fg4.ext_v1.1 = 
  SetIdent(src.endoderm.fg4.ext_v1.1,
           value = src.endoderm.fg4.ext_v1.1$cluster.v06.26.re)
markergene.fg4.ext_v1.1 = 
  FindAllMarkers(src.endoderm.fg4.ext_v1.1, assay = "RNA")
markergene.fg4.ext_v1.1 = markergene.fg4.ext_v1.1[
  markergene.fg4.ext_v1.1$avg_log2FC>0.2&
    markergene.fg4.ext_v1.1$p_val_adj<0.1, "gene"]

markergene.fg4.ext_v1.1_nolv = 
  FindAllMarkers(
    src.endoderm.fg4.ext_v1.1[,src.endoderm.fg4.ext_v1.1$cluster.v06.26.re%in%c(
      "FG.4","FG.4-Lung/Stomach","Lung","Stomach")], assay = "RNA")
markergene.fg4.ext_v1.1_nolv = unique(
  markergene.fg4.ext_v1.1_nolv[,"gene"])

seurat = src.endoderm.fg4.ext_v1.1[,src.endoderm.fg4.ext_v1.1$cluster.v06.26.re%in%c(
  "FG.4","FG.4-Liver","Liver")]
seurat = SetIdent(seurat, value = seurat$cluster.v06.26.re)
markergene.fg4.ext_v1.1_lv = FindAllMarkers(seurat, assay = "RNA")
markergene.fg4.ext_v1.1_lv = markergene.fg4.ext_v1.1_lv[
  markergene.fg4.ext_v1.1_lv$avg_log2FC>0.2&
    markergene.fg4.ext_v1.1_lv$p_val_adj<0.1, "gene"]
rm(seurat)

seurat = src.endoderm.fg4.ext_v1.1[,src.endoderm.fg4.ext_v1.1$cluster.v06.26.re%in%c(
  "FG.4","FG.4-Lung/Stomach","Stomach")]
seurat = SetIdent(seurat, value = seurat$cluster.v06.26.re)
markergene.fg4.ext_v1.1_sto = FindAllMarkers(seurat, assay = "RNA")
markergene.fg4.ext_v1.1_sto = markergene.fg4.ext_v1.1_sto[
  markergene.fg4.ext_v1.1_sto$avg_log2FC>0.2&
    markergene.fg4.ext_v1.1_sto$p_val_adj<0.1, "gene"]
rm(seurat)


seurat = src.endoderm.fg4.ext_v1.1[,src.endoderm.fg4.ext_v1.1$cluster.v06.26.re%in%c(
  "FG.4","FG.4-Lung/Stomach","Lung")]
seurat = SetIdent(seurat, value = seurat$cluster.v06.26.re)
markergene.fg4.ext_v1.1_lu = FindAllMarkers(seurat, assay = "RNA")
markergene.fg4.ext_v1.1_lu = markergene.fg4.ext_v1.1_lu[
  markergene.fg4.ext_v1.1_lu$avg_log2FC>0.2&
    markergene.fg4.ext_v1.1_lu$p_val_adj<0.1, "gene"]
rm(seurat)
#------------------------


#------------------------
fg4.umap.embedding = src.endoderm.fg4.ext_v1.1[["fdl_pca_integrated"]]@cell.embeddings[,1:3]
fg4.umap.embedding = src.endoderm.fg4.ext_v1.1[["fdl_mnn"]]@cell.embeddings[,1:3]
fg4.umap.embedding = src.endoderm.fg4.ext_v1.1[["umap_mnn"]]@cell.embeddings[,1:3]


# Cell order :: Princurve
#---------------------------
#  FG.4
cellorder.fg4_fg4 = 
  rownames(src.endoderm.fg4.ext_v1.1@meta.data[
    src.endoderm.fg4.ext_v1.1$cluster.v06.26.re%in%"FG.4",])
princurve.fg4_fg4 = princurve::principal_curve(
  x = fg4.umap.embedding[cellorder.fg4_fg4,],  smoother = "smooth.spline")
cellorder.fg4_fg4 =  names(
  princurve.fg4_fg4$lambda[cellorder.fg4_fg4][order(princurve.fg4_fg4$lambda[cellorder.fg4_fg4])])

# FG.4-Lung/Stomach
cellorder.fg4_fg4lu = 
  rownames(src.endoderm.fg4.ext_v1.1@meta.data[
    src.endoderm.fg4.ext_v1.1$cluster.v06.26.re%in%"FG.4-Lung/Stomach",])
princurve.fg4_fg4lu = princurve::principal_curve(
  x = fg4.umap.embedding[cellorder.fg4_fg4lu,],  smoother = "smooth.spline")
cellorder.fg4_fg4lu =  names(
  princurve.fg4_fg4lu$lambda[cellorder.fg4_fg4lu][order(princurve.fg4_fg4lu$lambda[cellorder.fg4_fg4lu])])

# FG.4-Liver
cellorder.fg4_fg4lv = 
  rownames(src.endoderm.fg4.ext_v1.1@meta.data[
    src.endoderm.fg4.ext_v1.1$cluster.v06.26.re%in%"FG.4-Liver",])
princurve.fg4_fg4lv = princurve::principal_curve(
  x = fg4.umap.embedding[cellorder.fg4_fg4lv,],  smoother = "smooth.spline")
cellorder.fg4_fg4lv =  names(
  princurve.fg4_fg4lv$lambda[cellorder.fg4_fg4lv][order(princurve.fg4_fg4lv$lambda[cellorder.fg4_fg4lv])])


# Lung
cellorder.fg4_lung = 
  rownames(src.endoderm.fg4.ext_v1.1@meta.data[
    src.endoderm.fg4.ext_v1.1$cluster.v06.26.re%in%"Lung",])
princurve.fg4_lung = princurve::principal_curve(
  x = fg4.umap.embedding[cellorder.fg4_lung,],  smoother = "smooth.spline")
cellorder.fg4_lung =  names(
  princurve.fg4_lung$lambda[cellorder.fg4_lung][order(princurve.fg4_lung$lambda[cellorder.fg4_lung])])

# FG.4 - Liver
cellorder.fg4_lv = 
  rownames(src.endoderm.fg4.ext_v1.1@meta.data[
    src.endoderm.fg4.ext_v1.1$cluster.v06.26.re%in%"Liver",])
princurve.fg4_lv = princurve::principal_curve(
  x = fg4.umap.embedding[cellorder.fg4_lv,],  smoother = "smooth.spline")
cellorder.fg4_lv =  names(
  princurve.fg4_lv$lambda[cellorder.fg4_lv][order(princurve.fg4_lv$lambda[cellorder.fg4_lv])])

# FG.4 - Sto
cellorder.fg4_sto = 
  rownames(src.endoderm.fg4.ext_v1.1@meta.data[
    src.endoderm.fg4.ext_v1.1$cluster.v06.26.re%in%"Stomach",])
princurve.fg4_sto = princurve::principal_curve(
  x = fg4.umap.embedding[cellorder.fg4_sto,],  smoother = "smooth.spline")
cellorder.fg4_sto =  names(
  princurve.fg4_sto$lambda[cellorder.fg4_sto][order(princurve.fg4_sto$lambda[cellorder.fg4_sto])])
#---------------------------

pdf("figure.v08.07/try.pdf")
plot(c(1:length(cellorder.fg4_fg4)),
     c(1:length(cellorder.fg4_fg4)),
     col = colors.time[src.endoderm.fg4.ext_v1.1[,cellorder.fg4_fg4]$Time])
plot(c(1:length(cellorder.fg4_fg4lu)),
     c(1:length(cellorder.fg4_fg4lu)),
     col = colors.time[src.endoderm.fg4.ext_v1.1[,cellorder.fg4_fg4lu]$Time])
plot(c(1:length(cellorder.fg4_fg4lv)),
     c(1:length(cellorder.fg4_fg4lv)),
     col = colors.time[src.endoderm.fg4.ext_v1.1[,cellorder.fg4_fg4lv]$Time])
plot(c(1:length(cellorder.fg4_lung)),
     c(1:length(cellorder.fg4_lung)),
     col = colors.time[src.endoderm.fg4.ext_v1.1[,cellorder.fg4_lung]$Time])
plot(c(1:length(cellorder.fg4_sto)),
     c(1:length(cellorder.fg4_sto)),
     col = colors.time[src.endoderm.fg4.ext_v1.1[,cellorder.fg4_sto]$Time])
plot(c(1:length(cellorder.fg4_lv)),
     c(1:length(cellorder.fg4_lv)),
     col = colors.time[src.endoderm.fg4.ext_v1.1[,cellorder.fg4_lv]$Time])
dev.off()


# cellorder.fg4_fg4_FPI = cellorder.fg4_fg4# 1:3
# cellorder.fg4_lung_FPI = cellorder.fg4_lung # 1:3
# cellorder.fg4_lv_FPI = cellorder.fg4_lv # 1:3
# cellorder.fg4_sto_FPI = cellorder.fg4_sto # 1:3

cellorder.fg4_fg4lu_UM = cellorder.fg4_fg4lu # 1:3
cellorder.fg4_fg4lv_UM = cellorder.fg4_fg4lv # 1:3


# No-Lv-Filter
#---------------------------------------------------------------------------------
pdf("figure.v08.07/organ_development_re/fg4_heatmap_marker.re_nolv.pdf",6,9)
cellorder.fg4 = c(rev(cellorder.fg4_fg4_FPI), cellorder.fg4_fg4lu_UM,
                  cellorder.fg4_lung_FPI, cellorder.fg4_sto_FPI)
selectgene = src.endoderm.fg4.ext_v1.1.filtergene_fin
#-------------------
fg4.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.fg4.ext_v1.1@assays$RNA@data[selectgene, cellorder.fg4]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.fg4.ext_v1.1$Time[cellorder.fg4], colors.time),
              MyName2Col(src.endoderm.fg4.ext_v1.1$cluster.v06.26.re[cellorder.fg4], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.fg4.ext_v1.1$cluster.extract.v1.1[cellorder.fg4], cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            #return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.fg4.ext_v1.1@meta.data[cell_fg4_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

fg4.ext_v1.1.re.rowtree.1 = as.dendrogram(fg4.ext_v1.1.re.rowtree)
src.endoderm.fg4.ext_v1.1.filtergene_fin = c(
  labels(fg4.ext_v1.1.re.rowtree.1[[2]][[2]][[2]][[1]][[2]]), # FG.4
  labels(fg4.ext_v1.1.re.rowtree.1[[2]][[2]][[2]][[2]]), # FG.4
  labels(fg4.ext_v1.1.re.rowtree.1[[2]][[1]]), # FG.4-Lu/Sto
  labels(fg4.ext_v1.1.re.rowtree.1[[1]][[2]][[2]][[2]][[1]][[1]]),  # FG.4-Lu/Sto
  labels(fg4.ext_v1.1.re.rowtree.1[[1]][[2]][[2]][[1]][[2]][[2]]),  # FG.4-Lu/Sto
  labels(fg4.ext_v1.1.re.rowtree.1[[1]][[1]][[1]]), # FG.4-Lu/Sto
  # labels(fg4.ext_v1.1.re.rowtree.1[[1]][[1]][[2]]), # FG.4-Lu/Sto
  setdiff(
    rev(labels(fg4.ext_v1.1.re.rowtree.1[[2]][[2]][[1]])),
    c("Prtg")),  # Lung
  labels(fg4.ext_v1.1.re.rowtree.1[[1]][[2]][[2]][[2]][[2]]),  # Stomach
  labels(fg4.ext_v1.1.re.rowtree.1[[1]][[2]][[2]][[2]][[1]][[2]]),  # Stomach
  labels(fg4.ext_v1.1.re.rowtree.1[[1]][[2]][[2]][[1]][[2]][[1]]),  # Stomach
  labels(fg4.ext_v1.1.re.rowtree.1[[1]][[2]][[2]][[1]][[1]]) # Stomach
)

names(src.endoderm.fg4.ext_v1.1.filtergene_fin) = c(
  rep(4, length(c(
    labels(fg4.ext_v1.1.re.rowtree.1[[2]][[2]][[2]][[2]]), # FG.4
    labels(fg4.ext_v1.1.re.rowtree.1[[2]][[2]][[2]][[1]][[2]]) # FG.4
  ))),
  rep(3, length(c(
    labels(fg4.ext_v1.1.re.rowtree.1[[2]][[1]]), # FG.4-Lu/Sto
    labels(fg4.ext_v1.1.re.rowtree.1[[1]][[1]][[1]]), # FG.4-Lu/Sto
    # labels(fg4.ext_v1.1.re.rowtree.1[[1]][[1]][[2]]), # FG.4-Lu/Sto
    labels(fg4.ext_v1.1.re.rowtree.1[[1]][[2]][[2]][[2]][[1]][[1]]),  # FG.4-Lu/Sto
    labels(fg4.ext_v1.1.re.rowtree.1[[1]][[2]][[2]][[1]][[2]][[2]])  # FG.4-Lu/Sto
  ))),
  rep(9, length(c(
    setdiff(
      rev(labels(fg4.ext_v1.1.re.rowtree.1[[2]][[2]][[1]])),
      c("Prtg")) # Lung
  ))),
  rep(7, length(c(
    labels(fg4.ext_v1.1.re.rowtree.1[[1]][[2]][[2]][[2]][[2]]),  # Stomach
    labels(fg4.ext_v1.1.re.rowtree.1[[1]][[2]][[2]][[2]][[1]][[2]]),  # Stomach
    labels(fg4.ext_v1.1.re.rowtree.1[[1]][[2]][[2]][[1]][[2]][[1]]),  # Stomach
    labels(fg4.ext_v1.1.re.rowtree.1[[1]][[2]][[2]][[1]][[1]]) # Stomach
  )))
)

src.endoderm.fg4.ext_v1.1.filtergene_nolv = 
  src.endoderm.fg4.ext_v1.1.filtergene_fin
#---------------------------------------------------------------------------------


# Lv
#---------------------------------------------------------------------------------
pdf("figure.v08.07/organ_development_re/fg4_heatmap_marker.re_lv.pdf",6,9)
cellorder.fg4 = c(rev(cellorder.fg4_fg4_FPI), 
                  rev(cellorder.fg4_fg4lv_UM), cellorder.fg4_lv_FPI)
#selectgene = src.endoderm.fg4.ext_v1.1.filtergene_fin
selectgene = src.endoderm.fg4.ext_v1.1.filtergene_lv
#-------------------
fg4.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.fg4.ext_v1.1@assays$RNA@data[selectgene, cellorder.fg4]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.fg4.ext_v1.1$Time[cellorder.fg4], colors.time),
              MyName2Col(src.endoderm.fg4.ext_v1.1$cluster.v06.26.re[cellorder.fg4], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.fg4.ext_v1.1$cluster.extract.v1.1[cellorder.fg4], cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            #return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.fg4.ext_v1.1@meta.data[cell_fg4_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()


fg4.ext_v1.1.re.rowtree.2 = as.dendrogram(fg4.ext_v1.1.re.rowtree)
src.endoderm.fg4.ext_v1.1.filtergene_fin = c(
  c("Aard","Nkx2-3","Has2","Tgfb2","Fhl1","Rhou"),
  rev(labels(fg4.ext_v1.1.re.rowtree.2[[1]][[2]][[2]])),

  labels(fg4.ext_v1.1.re.rowtree.2[[1]][[1]][[2]][[2]][[2]]),
  labels(fg4.ext_v1.1.re.rowtree.2[[1]][[1]][[2]][[2]][[1]]),
  labels(fg4.ext_v1.1.re.rowtree.2[[1]][[1]][[2]][[1]]),
  rev(labels(fg4.ext_v1.1.re.rowtree.2[[1]][[2]][[1]])),
  setdiff(
    rev(labels(fg4.ext_v1.1.re.rowtree.2[[1]][[1]][[1]])),
    c("Aard","Nkx2-3","Has2","Tgfb2","Fhl1","Rhou","Krt18")
  ),
  
  labels(fg4.ext_v1.1.re.rowtree.2[[2]][[1]][[2]][[1]]),
  labels(fg4.ext_v1.1.re.rowtree.2[[2]][[1]][[1]]),
  labels(fg4.ext_v1.1.re.rowtree.2[[2]][[2]][[1]]),
  rev(labels(fg4.ext_v1.1.re.rowtree.2[[2]][[2]][[2]]))
)

names(src.endoderm.fg4.ext_v1.1.filtergene_fin) = c(
  rep(4,length(c(
    labels(fg4.ext_v1.1.re.rowtree.2[[1]][[2]][[2]]),
    c("Aard","Nkx2-3","Has2","Tgfb2","Fhl1","Rhou")
  ))),
  rep(3,length(
    c(
      labels(fg4.ext_v1.1.re.rowtree.2[[1]][[1]][[2]][[2]][[2]]),
      labels(fg4.ext_v1.1.re.rowtree.2[[1]][[1]][[2]][[2]][[1]]),
      labels(fg4.ext_v1.1.re.rowtree.2[[1]][[1]][[2]][[1]]),
      rev(labels(fg4.ext_v1.1.re.rowtree.2[[1]][[2]][[1]])),
      setdiff(
        rev(labels(fg4.ext_v1.1.re.rowtree.2[[1]][[1]][[1]])),
        c("Aard","Nkx2-3","Has2","Tgfb2","Fhl1","Rhou","Krt18")
      )
  ))),
  rep(7,length(c(
    labels(fg4.ext_v1.1.re.rowtree.2[[2]][[1]][[2]][[1]]),
    labels(fg4.ext_v1.1.re.rowtree.2[[2]][[1]][[1]]),
    labels(fg4.ext_v1.1.re.rowtree.2[[2]][[2]][[1]]),
    rev(labels(fg4.ext_v1.1.re.rowtree.2[[2]][[2]][[2]]))
  )))
)



# names(src.endoderm.fg4.ext_v1.1.filtergene_lv) = c(
#   rep(3,95),rep(4,70+11),rep(7,129))

src.endoderm.fg4.ext_v1.1.filtergene_lv = 
  src.endoderm.fg4.ext_v1.1.filtergene_fin
#---------------------------------------------------------------------------------


# Lu
#---------------------------------------------------------------------------------
pdf("figure.v08.07/organ_development_re/fg4_heatmap_marker.re_lu.pdf",6,9)
cellorder.fg4 = c(rev(cellorder.fg4_fg4_FPI), 
                  cellorder.fg4_fg4lu_UM, cellorder.fg4_lung_FPI)
selectgene = src.endoderm.fg4.ext_v1.1.filtergene_fin 
#-------------------
fg4.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.fg4.ext_v1.1@assays$RNA@data[selectgene, cellorder.fg4]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.fg4.ext_v1.1$Time[cellorder.fg4], colors.time),
              MyName2Col(src.endoderm.fg4.ext_v1.1$cluster.v06.26.re[cellorder.fg4], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.fg4.ext_v1.1$cluster.extract.v1.1[cellorder.fg4], cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            #return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.fg4.ext_v1.1@meta.data[cell_fg4_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

fg4.ext_v1.1.re.rowtree.3 = as.dendrogram(fg4.ext_v1.1.re.rowtree)
src.endoderm.fg4.ext_v1.1.filtergene_fin = c(
  labels(fg4.ext_v1.1.re.rowtree.3[[2]][[2]][[2]]),
  rev(labels(fg4.ext_v1.1.re.rowtree.3[[2]][[2]][[1]])),
  rev(labels(fg4.ext_v1.1.re.rowtree.3[[1]])),
  labels(fg4.ext_v1.1.re.rowtree.3[[2]][[1]])
)

names(src.endoderm.fg4.ext_v1.1.filtergene_fin) = c(
  rep(4, length(c(
    labels(fg4.ext_v1.1.re.rowtree.3[[2]][[2]][[2]]),
    rev(labels(fg4.ext_v1.1.re.rowtree.3[[2]][[2]][[1]]))
  ))),
  rep(3, length(c(
    rev(labels(fg4.ext_v1.1.re.rowtree.3[[1]]))
  ))),
  rep(7, length(c(
    labels(fg4.ext_v1.1.re.rowtree.3[[2]][[1]])
  )))
)

src.endoderm.fg4.ext_v1.1.filtergene_lu = 
  src.endoderm.fg4.ext_v1.1.filtergene_fin
#---------------------------------------------------------------------------------


# No-Lv-All
#---------------------------------------------------------------------------------
pdf("figure.v08.07/organ_development_re/fg4_heatmap_marker.re_nolv.pdf",6,9)
cellorder.fg4 = c(rev(cellorder.fg4_fg4_FPI), cellorder.fg4_fg4lu_UM,
                  cellorder.fg4_lung_FPI, cellorder.fg4_sto_FPI)
selectgene = markergene.fg4.ext_v1.1_nolv 
#-------------------
fg4.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.fg4.ext_v1.1@assays$RNA@data[selectgene, cellorder.fg4]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.fg4.ext_v1.1$Time[cellorder.fg4], colors.time),
              MyName2Col(src.endoderm.fg4.ext_v1.1$cluster.v06.26.re[cellorder.fg4], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.fg4.ext_v1.1$cluster.extract.v1.1[cellorder.fg4], cluster.endoderm.color.v5)
            ),
            # RowSideColors = t(cbind(
            #   MyName2Col(names(selectgene), colors.geneset)
            # )),
            ColSideColorsSize = 2.5,
            # Rowv = "none",
            # Colv = "none",
            return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.fg4.ext_v1.1@meta.data[cell_fg4_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

src.endoderm.fg4.ext_v1.1.filtergene_nolv_all = 
  src.endoderm.fg4.ext_v1.1.filtergene_fin
#---------------------------------------------------------------------------------


#----------------------------------------------
save(
  # Gene 
  src.endoderm.fg4.ext_v1.1.gene.fin,
  src.endoderm.fg4.ext_v1.1.selectgene,
  src.endoderm.fg4.ext_v1.1.filtergene_lv,
  src.endoderm.fg4.ext_v1.1.filtergene_nolv,
  src.endoderm.fg4.ext_v1.1.filtergene_lu,
  # Cell type
  cellorder.fg4_fg4_FPI,
  cellorder.fg4_fg4lu_UM,cellorder.fg4_lung_FPI,cellorder.fg4_sto_FPI,
  cellorder.fg4_fg4lv_UM,cellorder.fg4_lv_FPI,
  file= "figure.v08.07/organ_development_re/fg4_heatmap_parameter.Rdata")



#==============================================================================
#  src.endoderm.fg4.ext_v1.1.re: Extraction Pharynx.organ.4 from FG.3
#==============================================================================

src.endoderm.fg4.ext_v1.1$cluster.extract.v1.1.re = 
  src.endoderm.fg4.ext_v1.1$cluster.extract.v1.1
src.endoderm.fg4.ext_v1.1@meta.data[
  colnames(src.endoderm.fg3.ext_v1.1[,src.endoderm.fg3.ext_v1.1$cluster.extract.v1.1%in%"FG.3/4"]),]$cluster.extract.v1.1.re = 
  src.endoderm.fg3.ext_v1.1@meta.data[
    colnames(src.endoderm.fg3.ext_v1.1[,src.endoderm.fg3.ext_v1.1$cluster.extract.v1.1%in%"FG.3/4"]),]$cluster.extract.v1.1.re
DimPlot(src.endoderm.fg4.ext_v1.1, group.by = 'cluster.extract.v1.1.re', 
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)


src.endoderm.fg4.ext_v1.1.re = 
  src.endoderm.fg4.ext_v1.1[,!src.endoderm.fg4.ext_v1.1$cluster.extract.v1.1.re%in%"FG.3"]

# Correct
#---------------------------------------------

src.endoderm.fg4.ext_v1.1.re@meta.data[
  src.endoderm.fg4.ext_v1.1.re$Time%in%c("ss24","ss27")&
    src.endoderm.fg4.ext_v1.1.re$cluster.revise.re.v1.30.re.fin%in%c("Lung","Pharynx.organ.4"),]$cluster.v06.26.re = "Lung"
src.endoderm.fg4.ext_v1.1.re@meta.data[
  src.endoderm.fg4.ext_v1.1.re$Time%in%c("ss24","ss27")&
    src.endoderm.fg4.ext_v1.1.re$cluster.revise.re.v1.30.re.fin%in%c("Stomach"),]$cluster.v06.26.re = "Stomach"
src.endoderm.fg4.ext_v1.1.re@meta.data[
  src.endoderm.fg4.ext_v1.1.re$cluster.extract.v1.1%in%c("FG.4-AL.1/2/3"),]$cluster.v06.26.re = "Liver"
#---------------------------------------------


src.endoderm.fg4.ext_v1.1.re = 
  seurat_mnn(seurat_name = "src.endoderm.fg4.ext_v1.1.re", dims = 1:30, 
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.fg4.ext_v1.1.gene.fin)
src.endoderm.fg4.ext_v1.1.re = 
  seurat_int(seurat_name = "src.endoderm.fg4.ext_v1.1.re", dims = 1:30,
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.fg4.ext_v1.1.gene.fin)

#---------------
src.endoderm.fg4.ext_v1.1.re = 
  FindNeighbors(src.endoderm.fg4.ext_v1.1.re, dims=1:30, 
                reduction = "pca_integrated", assay = "integrated")
src.endoderm.fg4.ext_v1.1.re = 
  FindNeighbors(src.endoderm.fg4.ext_v1.1.re, dims=1:30, 
                reduction = "pca", assay = "RNA")
src.endoderm.fg4.ext_v1.1.re = 
  FindClusters(src.endoderm.fg4.ext_v1.1.re, 
               resolution = 2, graph.name = "RNA_snn",)
#---------------
src.endoderm.fg4.ext_v1.1.re@reductions$umap_mnn_raw = 
  src.endoderm.fg4.ext_v1.1[,colnames(src.endoderm.fg4.ext_v1.1.re)]@reductions$umap_mnn

DimPlot(src.endoderm.fg4.ext_v1.1.re, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.2', reduction = "umap_mnn_raw")
DimPlot(src.endoderm.fg4.ext_v1.1.re, label = T, label.size = 4, cols = colors.num,
        group.by = 'cluster_snn_res.1.5', reduction = "umap_mnn")
DimPlot(src.endoderm.fg4.ext_v1.1.re, label = T, label.size = 4, cols = colors.num,
        group.by = 'cluster_snn_res.1', reduction = "umap_mnn")

DimPlot(src.endoderm.fg4.ext_v1.1.re, group.by = 'Time',
        reduction = "umap_mnn", cols = colors.time)
DimPlot(src.endoderm.fg4.ext_v1.1.re, group.by = 'Time',
        reduction = "umap_mnn_raw", cols = colors.time)
DimPlot(src.endoderm.fg4.ext_v1.1.re, group.by = 'Time', 
        reduction = "umap_integrated", cols = colors.time)
DimPlot(src.endoderm.fg4.ext_v1.1.re, group.by = 'batch', 
        reduction = "umap_integrated") #, cols = colors.time)

DimPlot(src.endoderm.fg4.ext_v1.1.re, group.by = 'cluster.extract.v1.1', 
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.fg4.ext_v1.1.re, group.by = 'cluster.extract.v1.1', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.fg4.ext_v1.1.re, group.by = 'cluster.v06.26.re', 
        reduction = "umap_mnn_raw", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.fg4.ext_v1.1.re, group.by = 'cluster.v06.26.re_v1.0', 
        reduction = "umap_mnn_raw", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.fg4.ext_v1.1.re, group.by = 'cluster.v06.26.re_hc',
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.fg4.ext_v1.1.re, group.by = 'cluster.v06.26.re_hc',
        reduction = "umap_int_rotated", cols = cluster.endoderm.color.v5)

src.endoderm.fg4.ext_v1.1.re



# FDL
src.endoderm.fg4.ext_v1.1.re = 
  seurat_fdl("src.endoderm.fg4.ext_v1.1.re", "pca","RNA")
src.endoderm.fg4.ext_v1.1.re = 
  seurat_fdl("src.endoderm.fg4.ext_v1.1.re", "mnn","RNA")
src.endoderm.fg4.ext_v1.1.re = 
  seurat_fdl("src.endoderm.fg4.ext_v1.1.re", "pca_integrated","integrated")

# Graph
pdf("figure.v08.07/organ_development_re/fg4_summary_graph.re_snn_PI_ext.pdf",12,7)
# SNN
for(red in c("pca_integrated","mnn","pca")){
  
  assay = ifelse(red%in%c("pca_integrated"), "integrated", "RNA")
  graph.fg4.ext_v1.1.re = 
    seurat_GCN(seurat_name = "src.endoderm.fg4.ext_v1.1.re",
               src.endoderm.fg4.ext_v1.1.re.gene.fin,
               reduction = red, assay =  assay)
  graph.fg4_cluster = cluster_walktrap(graph.fg4.ext_v1.1.re, steps = 1)
  
  for(color in c("cluster.extract.v1.1","cluster.v06.26.re",
                 "Time")){
    for(edge.color in c(NA,"lightgray")){
      for(edge.color in c(NA,"lightgray")){
        set.seed(1)
        plot.igraph(graph.fg4.ext_v1.1.re,
                    layout = layout_with_fr(graph.fg4.ext_v1.1.re),
                    edge.color = edge.color,
                    vertex.size=2.5, vertex.label=NA,
                    vertex.label.color = "black", 
                    vertex.frame.color = NA, vertex.frame.width = 0.5,
                    #vertex.color =  colors.time[src.endoderm.fg4.ext_v1.1.re[["Time"]][graph.fg4_cluster$names,]],
                    #vertex.color =  cluster.endoderm.color.v5[src.endoderm.fg4.ext_v1.1.re[["cluster.v06.26.re"]][graph.fg4_cluster$names,]],
                    vertex.color =  colors_bg[src.endoderm.fg4.ext_v1.1.re[[color]][graph.fg4_cluster$names,]])
      }}}
  
}
dev.off()

pdf("figure.v08.07/organ_development_re/fg4_summary_graph.re_other_ext.pdf",12,7)
# Reduction
for(red in c("umap_integrated","umap_mnn","umap_mnn_raw",
             "fdl_pca",'fdl_mnn',"fdl_pca_integrated")){
  p1 = DimPlot(src.endoderm.fg4.ext_v1.1.re,group.by = 'cluster.v06.26.re', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p2 = DimPlot(src.endoderm.fg4.ext_v1.1.re,group.by = 'cluster.extract.v1.1', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p3 = DimPlot(src.endoderm.fg4.ext_v1.1.re,group.by = 'Time', 
               reduction = red , cols = colors.time,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  # p4 = DimPlot(src.endoderm.fg4.ext_v1.1.re,group.by = 'cluster.v06.26.re_hc', 
  #              reduction = red, cols = cluster.endoderm.color.v5,
  #              pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  print(p1); print(p2); print(p3);# print(p4)
}
dev.off()


pdf("figure.v08.07/organ_development_re/fg4_summary_proportion.re_ext.pdf",12,7)
#-----------------------------
cell_type = c("FG.4","FG.4-Lung/Stomach","FG.4-Liver","Lung","Stomach","Liver")
meta_fg4 = src.endoderm.fg4.ext_v1.1.re@meta.data
meta_fg4$Time = factor(meta_fg4$Time, levels = names(colors.time))
meta_fg4$cluster.extract.v1.1 = factor(meta_fg4$cluster.extract.v1.1, 
                                       levels = c("FG.4","FG.3/4","FG.4-MG.1/3",'FG.4-AL.1/2/3'))
meta_fg4$cluster.v06.26.re = factor(meta_fg4$cluster.v06.26.re, levels = rev(cell_type))
#meta_fg4$cluster.v06.26.re_hc = factor(meta_fg4$cluster.v06.26.re_hc, levels = rev(cell_type))

meta_fg4_B0 = meta_fg4
meta_fg4_B1 = meta_fg4[meta_fg4$batch%in%1&!meta_fg4$Time%in%"ss9",]
meta_fg4_B2 = meta_fg4[meta_fg4$batch%in%2,]

for(data in c("meta_fg4_B0","meta_fg4_B1","meta_fg4_B2")){
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

save(src.endoderm.fg4.ext_v1.1.re,
     file = "figure.v08.07/organ_development_re/src.endoderm.fg4.ext_v1.1.re.re.Rdata")

fg4.umap.embedding = src.endoderm.fg4.ext_v1.1.re[["fdl_pca_integrated"]]@cell.embeddings[,1:3]
fg4.umap.embedding = src.endoderm.fg4.ext_v1.1.re[["fdl_mnn"]]@cell.embeddings[,1:3]
fg4.umap.embedding = src.endoderm.fg4.ext_v1.1.re[["umap_mnn"]]@cell.embeddings[,1:3]

# Cell order :: Princurve
#---------------------------
#  FG.4
cellorder.fg4_fg4 = 
  rownames(src.endoderm.fg4.ext_v1.1.re@meta.data[
    src.endoderm.fg4.ext_v1.1.re$cluster.v06.26.re%in%"FG.4",])
princurve.fg4_fg4 = princurve::principal_curve(
  x = fg4.umap.embedding[cellorder.fg4_fg4,],  smoother = "smooth.spline")
cellorder.fg4_fg4 =  names(
  princurve.fg4_fg4$lambda[cellorder.fg4_fg4][order(princurve.fg4_fg4$lambda[cellorder.fg4_fg4])])

# FG.4-Lung/Stomach
cellorder.fg4_fg4lu = 
  rownames(src.endoderm.fg4.ext_v1.1.re@meta.data[
    src.endoderm.fg4.ext_v1.1.re$cluster.v06.26.re%in%"FG.4-Lung/Stomach",])
princurve.fg4_fg4lu = princurve::principal_curve(
  x = fg4.umap.embedding[cellorder.fg4_fg4lu,],  smoother = "smooth.spline")
cellorder.fg4_fg4lu =  names(
  princurve.fg4_fg4lu$lambda[cellorder.fg4_fg4lu][order(princurve.fg4_fg4lu$lambda[cellorder.fg4_fg4lu])])

# FG.4-Liver
cellorder.fg4_fg4lv = 
  rownames(src.endoderm.fg4.ext_v1.1.re@meta.data[
    src.endoderm.fg4.ext_v1.1.re$cluster.v06.26.re%in%"FG.4-Liver",])
princurve.fg4_fg4lv = princurve::principal_curve(
  x = fg4.umap.embedding[cellorder.fg4_fg4lv,],  smoother = "smooth.spline")
cellorder.fg4_fg4lv =  names(
  princurve.fg4_fg4lv$lambda[cellorder.fg4_fg4lv][order(princurve.fg4_fg4lv$lambda[cellorder.fg4_fg4lv])])


# Lung
cellorder.fg4_lung = 
  rownames(src.endoderm.fg4.ext_v1.1.re@meta.data[
    src.endoderm.fg4.ext_v1.1.re$cluster.v06.26.re%in%"Lung",])
princurve.fg4_lung = princurve::principal_curve(
  x = fg4.umap.embedding[cellorder.fg4_lung,],  smoother = "smooth.spline")
cellorder.fg4_lung =  names(
  princurve.fg4_lung$lambda[cellorder.fg4_lung][order(princurve.fg4_lung$lambda[cellorder.fg4_lung])])

# FG.4 - Liver
cellorder.fg4_lv = 
  rownames(src.endoderm.fg4.ext_v1.1.re@meta.data[
    src.endoderm.fg4.ext_v1.1.re$cluster.v06.26.re%in%"Liver",])
princurve.fg4_lv = princurve::principal_curve(
  x = fg4.umap.embedding[cellorder.fg4_lv,],  smoother = "smooth.spline")
cellorder.fg4_lv =  names(
  princurve.fg4_lv$lambda[cellorder.fg4_lv][order(princurve.fg4_lv$lambda[cellorder.fg4_lv])])

# FG.4 - Sto
cellorder.fg4_sto = 
  rownames(src.endoderm.fg4.ext_v1.1.re@meta.data[
    src.endoderm.fg4.ext_v1.1.re$cluster.v06.26.re%in%"Stomach",])
princurve.fg4_sto = princurve::principal_curve(
  x = fg4.umap.embedding[cellorder.fg4_sto,],  smoother = "smooth.spline")
cellorder.fg4_sto =  names(
  princurve.fg4_sto$lambda[cellorder.fg4_sto][order(princurve.fg4_sto$lambda[cellorder.fg4_sto])])
#---------------------------

pdf("figure.v08.07/try.pdf")
#---------------------------
plot(c(1:length(cellorder.fg4_fg4)),
     c(1:length(cellorder.fg4_fg4)),
     col = colors.time[src.endoderm.fg4.ext_v1.1.re[,cellorder.fg4_fg4]$Time])
plot(c(1:length(cellorder.fg4_fg4lu)),
     c(1:length(cellorder.fg4_fg4lu)),
     col = colors.time[src.endoderm.fg4.ext_v1.1.re[,cellorder.fg4_fg4lu]$Time])
plot(c(1:length(cellorder.fg4_fg4lv)),
     c(1:length(cellorder.fg4_fg4lv)),
     col = colors.time[src.endoderm.fg4.ext_v1.1.re[,cellorder.fg4_fg4lv]$Time])
plot(c(1:length(cellorder.fg4_lung)),
     c(1:length(cellorder.fg4_lung)),
     col = colors.time[src.endoderm.fg4.ext_v1.1.re[,cellorder.fg4_lung]$Time])
plot(c(1:length(cellorder.fg4_sto)),
     c(1:length(cellorder.fg4_sto)),
     col = colors.time[src.endoderm.fg4.ext_v1.1.re[,cellorder.fg4_sto]$Time])
plot(c(1:length(cellorder.fg4_lv)),
     c(1:length(cellorder.fg4_lv)),
     col = colors.time[src.endoderm.fg4.ext_v1.1.re[,cellorder.fg4_lv]$Time])
#---------------------------
dev.off()


cellorder.fg4_fg4_FPI = cellorder.fg4_fg4# 1:3
cellorder.fg4_fg4lu_FPI = cellorder.fg4_fg4lu # 1:3
cellorder.fg4_lung_FPI = cellorder.fg4_lung # 1:3
cellorder.fg4_sto_FPI = cellorder.fg4_sto # 1:3
cellorder.fg4_fg4lv_FPI = cellorder.fg4_fg4lv # 1:3
cellorder.fg4_lv_FPI = cellorder.fg4_lv # 1:3

src.endoderm.fg4.ext_v1.1.re = SetIdent(
  src.endoderm.fg4.ext_v1.1.re,
  value = src.endoderm.fg4.ext_v1.1.re$cluster.v06.26.re)


#=================================================================================
# Marker: Liver
#=================================================================================
pdf("figure.v08.07/organ_development_re/fg4_heatmap_marker.re_lv.pdf",6,9)
cellorder.fg4 = c(cellorder.fg4_fg4_FPI, 
                  cellorder.fg4_fg4lv_FPI, cellorder.fg4_lv_FPI)

markergene.fg4.ext_v1.1.re = 
  FindAllMarkers(src.endoderm.fg4.ext_v1.1.re[,cellorder.mg3], assay = "RNA")
markergene.fg4.ext_v1.1.re = markergene.fg4.ext_v1.1.re[
  markergene.fg4.ext_v1.1.re$avg_log2FC>0.2&markergene.fg4.ext_v1.1.re$p_val_adj<0.1, "gene"]
selectgene = unique(c(markergene.fg4.ext_v1.1.re))

selectgene = src.endoderm.fg4.ext_v1.1.filtergene_lv
#-------------------
fg4.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.fg4.ext_v1.1.re@assays$RNA@data[selectgene, cellorder.fg4]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.fg4.ext_v1.1.re$Time[cellorder.fg4], colors.time),
              MyName2Col(src.endoderm.fg4.ext_v1.1.re$cluster.v06.26.re[cellorder.fg4], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.fg4.ext_v1.1.re$cluster.extract.v1.1[cellorder.fg4], cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            #return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.fg4.ext_v1.1.re@meta.data[cell_fg4_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

src.endoderm.fg4.ext_v1.1.re.filtergene_lv = 
  src.endoderm.fg4.ext_v1.1.filtergene_lv

#=================================================================================
# Marker: Stomach
#=================================================================================
pdf("figure.v08.07/organ_development_re/fg4_heatmap_marker.re_sto.pdf",6,9)
cellorder.fg4 = c(cellorder.fg4_fg4_FPI, 
                  cellorder.fg4_fg4lu_FPI, cellorder.fg4_sto_FPI)

markergene.fg4.ext_v1.1.re = 
  FindAllMarkers(src.endoderm.fg4.ext_v1.1.re[,cellorder.fg4], assay = "RNA")
markergene.fg4.ext_v1.1.re = markergene.fg4.ext_v1.1.re[
  markergene.fg4.ext_v1.1.re$avg_log2FC>0.2&markergene.fg4.ext_v1.1.re$p_val_adj<0.1, "gene"]
selectgene = unique(c(markergene.fg4.ext_v1.1.re))

selectgene = src.endoderm.fg4.ext_v1.1.re.filtergene 
#-------------------
fg4.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.fg4.ext_v1.1.re@assays$RNA@data[selectgene, cellorder.fg4]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.fg4.ext_v1.1.re$Time[cellorder.fg4], colors.time),
              MyName2Col(src.endoderm.fg4.ext_v1.1.re$cluster.v06.26.re[cellorder.fg4], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.fg4.ext_v1.1.re$cluster.extract.v1.1[cellorder.fg4], cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            #return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.fg4.ext_v1.1.re@meta.data[cell_fg4_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

fg4.ext_v1.1.re.rowtree...1 = as.dendrogram(fg4.ext_v1.1.re.rowtree)
src.endoderm.fg4.ext_v1.1.re.filtergene = c(
  setdiff(
    c(
      labels(fg4.ext_v1.1.re.rowtree...1[[2]][[1]])
    ),
    c("Gm10260","Ddx21")
  ),
  setdiff(
    c(
      rev(labels(fg4.ext_v1.1.re.rowtree...1[[2]][[2]]))
    ),
    c()
  ),
  setdiff(
    c(
      labels(fg4.ext_v1.1.re.rowtree...1[[1]][[2]])
    ),
    c()
  ),
  setdiff(
    c(
      labels(fg4.ext_v1.1.re.rowtree...1[[1]][[1]])
    ),
    c()
  )
)

names(src.endoderm.fg4.ext_v1.1.re.filtergene) = c(
  rep(4,length(  
    setdiff(
      c(
        labels(fg4.ext_v1.1.re.rowtree...1[[2]][[1]])
      ),
      c("Gm10260","Ddx21")
    ))),
  rep(3,length(  
    setdiff(
      c(
        labels(fg4.ext_v1.1.re.rowtree...1[[2]][[2]])
      ),
      c()
    ))),
  rep(9,length(  
    setdiff(
      c(
        labels(fg4.ext_v1.1.re.rowtree...1[[1]][[2]])
      ),
      c()
    ))),
  rep(7,length(  
    setdiff(
      c(
        labels(fg4.ext_v1.1.re.rowtree...1[[1]][[1]])
      ),
      c()
    )))
)
src.endoderm.fg4.ext_v1.1.re.filtergene_sto = 
  src.endoderm.fg4.ext_v1.1.re.filtergene 

#=================================================================================
# Marker: Lung
#=================================================================================
pdf("figure.v08.07/organ_development_re/fg4_heatmap_marker.re_lu.pdf",6,9)
cellorder.fg4 = c(cellorder.fg4_fg4_FPI, 
                  cellorder.fg4_fg4lu_FPI, cellorder.fg4_lung_FPI)

markergene.fg4.ext_v1.1.re = 
  FindAllMarkers(src.endoderm.fg4.ext_v1.1.re[,cellorder.fg4], assay = "RNA")
markergene.fg4.ext_v1.1.re = markergene.fg4.ext_v1.1.re[
  markergene.fg4.ext_v1.1.re$avg_log2FC>0.1&markergene.fg4.ext_v1.1.re$p_val_adj<0.1, "gene"]
selectgene = unique(c(markergene.fg4.ext_v1.1.re))

selectgene = src.endoderm.fg4.ext_v1.1.re.filtergene 

#-------------------
fg4.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.fg4.ext_v1.1.re@assays$RNA@data[selectgene, cellorder.fg4]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.fg4.ext_v1.1.re$Time[cellorder.fg4], colors.time),
              MyName2Col(src.endoderm.fg4.ext_v1.1.re$cluster.v06.26.re[cellorder.fg4], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.fg4.ext_v1.1.re$cluster.extract.v1.1[cellorder.fg4], cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            # return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.fg4.ext_v1.1.re@meta.data[cell_fg4_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

fg4.ext_v1.1.re.rowtree...2 = as.dendrogram(fg4.ext_v1.1.re.rowtree)
src.endoderm.fg4.ext_v1.1.re.filtergene = c(
  setdiff(
    c(
      labels(fg4.ext_v1.1.re.rowtree...2[[1]][[2]]),
      labels(fg4.ext_v1.1.re.rowtree...2[[1]][[1]])
    ),
    c("Sfrp5")
  ),
  # setdiff(
  #   c(
  #     rev(labels(fg4.ext_v1.1.re.rowtree...2[[2]][[1]]))
  #   ),
  #   c()
  # ),
  setdiff(
    c(
      labels(fg4.ext_v1.1.re.rowtree...2[[2]][[2]])
    ),
    c()
  )
)

names(src.endoderm.fg4.ext_v1.1.re.filtergene) = c(
  rep(4,length(  
    setdiff(
      c(
        labels(fg4.ext_v1.1.re.rowtree...2[[1]])
      ),
      c("Sfrp5")
    ))),
  # rep(3,length(  
  #   setdiff(
  #     c(
  #       labels(fg4.ext_v1.1.re.rowtree...2[[2]][[1]])
  #     ),
  #     c()
  #   ))),
  rep(7,length(  
    setdiff(
      c(
        labels(fg4.ext_v1.1.re.rowtree...2[[2]][[2]])
      ),
      c()
    )))
)
src.endoderm.fg4.ext_v1.1.re.filtergene_lu = 
  src.endoderm.fg4.ext_v1.1.re.filtergene 


#---------------------------------------------------------------------------------
save(
  # Gene 
  src.endoderm.fg4.ext_v1.1.gene.fin,
  src.endoderm.fg4.ext_v1.1.selectgene,
  src.endoderm.fg4.ext_v1.1.filtergene_lv,
  src.endoderm.fg4.ext_v1.1.filtergene_nolv,
  src.endoderm.fg4.ext_v1.1.filtergene_lu,
  
  src.endoderm.fg4.ext_v1.1.re.filtergene_sto,
  src.endoderm.fg4.ext_v1.1.re.filtergene_lv,
  
  # Cell type
  cellorder.fg4_fg4_FPI,
  cellorder.fg4_fg4lu_FPI, cellorder.fg4_lung_FPI, cellorder.fg4_sto_FPI,
  cellorder.fg4_fg4lv_FPI, cellorder.fg4_lv_FPI,
  file= "figure.v08.07/organ_development_re/fg4_heatmap_parameter.Rdata")



