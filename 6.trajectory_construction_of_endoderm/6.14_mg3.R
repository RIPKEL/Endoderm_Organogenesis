#----------------
#     Remake
#----------------
src.endoderm.mg3.ext_v1.1 = FindVariableFeatures(src.endoderm.mg3.ext_v1.1, nfeatures = 2000)
src.endoderm.mg3.ext_v1.1 = ScaleData(src.endoderm.mg3.ext_v1.1, split.by = "batch_phase",
                                      features = rownames(src.endoderm.mg3.ext_v1.1))
src.endoderm.mg3.ext_v1.1.filtergene = 
  Myfilter(as.matrix(src.endoderm.mg3.ext_v1.1@assays$RNA@data),
           gene = src.endoderm.mg3.ext_v1.1@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
#rm(src.endoderm.mg3.ext_v1.1.filtergene)
src.endoderm.mg3.ext.v1.1.selectgene = unique(c(
  select_gene_mg3, src.endoderm.mg3.ext_v1.1.filtergene))

#=============================
# Row
pdf("figure.v08.07/try.pdf")
#-----------------
mg3.ext.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.mg3.ext_v1.1@assays$RNA@data[
    src.endoderm.mg3.ext_v1.1.gene.fin,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.endoderm.mg3.ext_v1.1$Time,colors.time),
      MyName2Col(src.endoderm.mg3.ext_v1.1$cluster.extract.v1.1,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.mg3.ext_v1.1$cluster.v06.26.re,cluster.endoderm.color.v5)
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
plot(mg3.ext.re.rowtree)
dev.off()

mg3.ext.re.rowtree_1 = as.dendrogram(mg3.ext.re.rowtree)
src.endoderm.mg3.ext_v1.1.gene.fin = setdiff(
  src.endoderm.mg3.ext.v1.1.selectgene,
  c(labels(mg3.ext.re.rowtree_1[[2]][[1]]),
    labels(mg3.ext.re.rowtree_1[[2]][[2]][[2]][[1]])))

#--------------------------------------------------------------------------------
src.endoderm.mg3.ext_v1.1 = RunPCA(src.endoderm.mg3.ext_v1.1,
                                   features = src.endoderm.mg3.ext_v1.1.gene.fin)
src.endoderm.mg3.ext_v1.1 = 
  seurat_mnn(seurat_name = "src.endoderm.mg3.ext_v1.1", dims = 1:25, 
             n.neighbors = 75,n.components = 2,
             seurat.selectgene = src.endoderm.mg3.ext_v1.1.gene.fin)
src.endoderm.mg3.ext_v1.1 = 
  seurat_int(seurat_name = "src.endoderm.mg3.ext_v1.1", dims = 1:25,
             n.neighbors = 75,n.components = 2,
             seurat.selectgene = src.endoderm.mg3.ext_v1.1.gene.fin)

#---------------
src.endoderm.mg3.ext_v1.1@reductions$mnn@cell.embeddings = 
  src.endoderm.mg3.ext_v1.1@reductions$mnn@cell.embeddings[colnames(src.endoderm.mg3.ext_v1.1),]
src.endoderm.mg3.ext_v1.1 = 
  FindNeighbors(src.endoderm.mg3.ext_v1.1, dims=1:30, 
                reduction = "mnn", assay = "RNA")
src.endoderm.mg3.ext_v1.1 = 
  FindClusters(src.endoderm.mg3.ext_v1.1, resolution = 2)
src.endoderm.mg3.ext_v1.1 = 
  FindClusters(src.endoderm.mg3.ext_v1.1, resolution = 4)
#---------------

pdf("figure.v08.07/try.pdf")
#---------------
DimPlot(src.endoderm.mg3.ext_v1.1, label = T, label.size = 4, cols = colors.num, pt.size = 1.2,
        group.by = 'RNA_snn_res.2', reduction = "umap_mnn")
DimPlot(src.endoderm.mg3.ext_v1.1, label = T, label.size = 4, cols = colors.num, pt.size = 1.2,
        group.by = 'RNA_snn_res.4', reduction = "umap_mnn")

DimPlot(src.endoderm.mg3.ext_v1.1, group.by = 'Time',pt.size = 1.2,
        reduction = "umap_mnn", cols = colors.time)
DimPlot(src.endoderm.mg3.ext_v1.1, group.by = 'Time', pt.size = 1.2,
        reduction = "umap_integrated", cols = colors.time)
DimPlot(src.endoderm.mg3.ext_v1.1, group.by = 'batch', 
        reduction = "umap_integrated") #, cols = colors.time)

DimPlot(src.endoderm.mg3.ext_v1.1, group.by = 'cluster.extract.v1.1', 
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.mg3.ext_v1.1, group.by = 'cluster.extract.v1.1', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.mg3.ext_v1.1, group.by = 'cluster.v06.26.re', pt.size = 1.2,
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.mg3.ext_v1.1, group.by = 'cluster.v06.26.re', pt.size = 1.2,
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
#---------------
dev.off()

src.endoderm.mg3.ext_v1.1@reductions$pca_integrated@cell.embeddings = 
  src.endoderm.mg3.ext_v1.1@reductions$pca_integrated@cell.embeddings[colnames(src.endoderm.mg3.ext_v1.1),]
src.endoderm.mg3.ext_v1.1@reductions$mnn@cell.embeddings = 
  src.endoderm.mg3.ext_v1.1@reductions$mnn@cell.embeddings[colnames(src.endoderm.mg3.ext_v1.1),]
#--------------------
src.endoderm.mg3.ext_v1.1$cluster.v06.26.re_raw = NA
src.endoderm.mg3.ext_v1.1@meta.data[
  intersect(colnames(src.endoderm.mg3.ext_v1.1),colnames(src.endoderm.mg3.ext.re)),]$cluster.v06.26.re_raw = 
  src.endoderm.mg3.ext.re@meta.data[
    intersect(colnames(src.endoderm.mg3.ext_v1.1),colnames(src.endoderm.mg3.ext.re)),]$cluster.v06.26

#---------------------
src.endoderm.mg3.ext_v1.1$cluster.v06.26.re = NA
src.endoderm.mg3.ext_v1.1@meta.data[
  intersect(colnames(src.endoderm.mg3.ext_v1.1),colnames(src.endoderm.mg3.ext.re)),]$cluster.v06.26.re = 
  src.endoderm.mg3.ext.re@meta.data[
    intersect(colnames(src.endoderm.mg3.ext_v1.1),colnames(src.endoderm.mg3.ext.re)),]$cluster.v06.26.re

src.endoderm.mg3.ext_v1.1@meta.data[
  src.endoderm.mg3.ext_v1.1$RNA_snn_res.2%in%c(16,18,8),]$cluster.v06.26.re = "MG.3"
src.endoderm.mg3.ext_v1.1@meta.data[
  src.endoderm.mg3.ext_v1.1$RNA_snn_res.2%in%c(11,9)|
    src.endoderm.mg3.ext_v1.1$RNA_snn_res.4%in%c(8,6),]$cluster.v06.26.re = "MG.3.A/M"
src.endoderm.mg3.ext_v1.1@meta.data[
  src.endoderm.mg3.ext_v1.1$RNA_snn_res.2%in%c(7)|
    src.endoderm.mg3.ext_v1.1$RNA_snn_res.4%in%c(32,23,20,19),]$cluster.v06.26.re = "MG.3.A/M"
src.endoderm.mg3.ext_v1.1@meta.data[
  src.endoderm.mg3.ext_v1.1$RNA_snn_res.2%in%c(3,4,12),]$cluster.v06.26.re = "MG.3.P"
src.endoderm.mg3.ext_v1.1@meta.data[
  src.endoderm.mg3.ext_v1.1$RNA_snn_res.4%in%c(12,24,26),]$cluster.v06.26.re = "DP"
src.endoderm.mg3.ext_v1.1@meta.data[
  src.endoderm.mg3.ext_v1.1$RNA_snn_res.4%in%c(2,3,7,21),]$cluster.v06.26.re = "Stomach"
src.endoderm.mg3.ext_v1.1@meta.data[
  src.endoderm.mg3.ext_v1.1$RNA_snn_res.2%in%c(6,13),]$cluster.v06.26.re = "Small.intestine.1"
src.endoderm.mg3.ext_v1.1@meta.data[
  src.endoderm.mg3.ext_v1.1$RNA_snn_res.2%in%c(1,0),]$cluster.v06.26.re = "Small.intestine.2"

src.endoderm.mg3.ext_v1.1@meta.data[src.endoderm.mg3.ext_v1.1$cluster.extract.v1.1%in%"MG.2/3",]$cluster.v06.26.re = 
  src.endoderm.mg2.ext_v1.1@meta.data[rownames(
    src.endoderm.mg3.ext_v1.1@meta.data[src.endoderm.mg3.ext_v1.1$cluster.extract.v1.1%in%"MG.2/3",]),]$cluster.v06.26.re
src.endoderm.mg3.ext_v1.1@meta.data[
  src.endoderm.mg3.ext_v1.1$cluster.extract.v1.1%in%c("FG.4-MG.1/3"),]$cluster.v06.26.re = "Stomach"
src.endoderm.mg3.ext_v1.1@meta.data[
  src.endoderm.mg3.ext_v1.1$cluster.extract.v1.1%in%c("MG.2/3-HG.1"),]$cluster.v06.26.re = "Small.intestine.2"
src.endoderm.mg3.ext_v1.1@meta.data[
  src.endoderm.mg3.ext_v1.1$cluster.extract.v1.1%in%c("MG.1/2/3","AL.3-MG.1/2/3"),]$cluster.v06.26.re = "Small.intestine.1"

src.endoderm.mg3.ext_v1.1@meta.data[
  src.endoderm.mg3.ext_v1.1$Time%in%c("ss24","ss27")&
    src.endoderm.mg3.ext_v1.1$cluster.revise.re.v1.30.re.fin%in%c("Small.intestine.2"),]$cluster.v06.26.re = "Small.intestine.1"
src.endoderm.mg3.ext_v1.1@meta.data[
  src.endoderm.mg3.ext_v1.1$Time%in%c("ss24","ss27")&
    src.endoderm.mg3.ext_v1.1$cluster.revise.re.v1.30.re.fin%in%c("Small.intestine.2"),]$cluster.v06.26.re = "Small.intestine.2"
src.endoderm.mg3.ext_v1.1@meta.data[
  src.endoderm.mg3.ext_v1.1$Time%in%c("ss24","ss27")&
    src.endoderm.mg3.ext_v1.1$cluster.revise.re.v1.30.re.fin%in%c("Stomach"),]$cluster.v06.26.re = "Stomach"
src.endoderm.mg3.ext_v1.1@meta.data[
  src.endoderm.mg3.ext_v1.1$Time%in%c("ss24","ss27")&
    src.endoderm.mg3.ext_v1.1$cluster.revise.re.v1.30.re.fin%in%c("DP","Pancreas"),]$cluster.v06.26.re = "DP"

src.endoderm.mg3.ext_v1.1@meta.data[
  src.endoderm.mg3.ext_v1.1$RNA_snn_res.4%in%c(27),]$cluster.v06.26.re = "EP.1"
src.endoderm.mg3.ext_v1.1@meta.data[
  src.endoderm.mg3.ext_v1.1$RNA_snn_res.4%in%c(33),]$cluster.v06.26.re = "EP.2"
src.endoderm.mg3.ext_v1.1@meta.data[
  intersect(colnames(src.endoderm.mg3.ext_v1.1),
            colnames(src.endoderm.mg3.ext.re[,src.endoderm.mg3.ext.re$cluster.v06.26.re%in%c("EP.1","EP.2","EP")])),]$cluster.v06.26.re = 
  src.endoderm.mg3.ext.re@meta.data[
    intersect(colnames(src.endoderm.mg3.ext_v1.1),
              colnames(src.endoderm.mg3.ext.re[,src.endoderm.mg3.ext.re$cluster.v06.26.re%in%c("EP.1","EP.2","EP")])),]$cluster.v06.26.re

src.endoderm.mg3.ext_v1.1@meta.data[
  src.endoderm.mg3.ext_v1.1$RNA_snn_res.4%in%c(16)|
    src.endoderm.mg3.ext_v1.1$RNA_snn_res.2%in%c(19),]$cluster.v06.26.re = NA
src.endoderm.mg3.ext_v1.1@meta.data[
  src.endoderm.mg3.ext_v1.1$cluster.v06.26.re%in%c("MG.2","MG.3.A/M-DP"),]$cluster.v06.26.re = NA

a = FNN::knn(
  src.endoderm.mg3.ext_v1.1@reductions$mnn@cell.embeddings[
    rownames(src.endoderm.mg3.ext_v1.1@meta.data[!src.endoderm.mg3.ext_v1.1$cluster.v06.26.re%in%NA,]),],
  src.endoderm.mg3.ext_v1.1@reductions$mnn@cell.embeddings[
    rownames(src.endoderm.mg3.ext_v1.1@meta.data[src.endoderm.mg3.ext_v1.1$cluster.v06.26.re%in%NA,]),],
  src.endoderm.mg3.ext_v1.1@meta.data[!src.endoderm.mg3.ext_v1.1$cluster.v06.26.re%in%NA,]$cluster.v06.26.re, k = 10)
src.endoderm.mg3.ext_v1.1@meta.data[
  src.endoderm.mg3.ext_v1.1$cluster.v06.26.re%in%NA,]$cluster.v06.26.re = as.character(a)

src.endoderm.mg3.ext_v1.1@meta.data[
  src.endoderm.mg3.ext_v1.1@reductions$umap_mnn@cell.embeddings[,1]>1&
    src.endoderm.mg3.ext_v1.1$cluster.v06.26.re%in%c("MG.3.P"),]$cluster.v06.26.re = "MG.3.A/M"

DimPlot(src.endoderm.mg3.ext_v1.1, group.by = 'cluster.v06.26.re', pt.size = 1.2,
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)

src.endoderm.mg3.ext_v1.1$cluster.v06.26.re_ep = src.endoderm.mg3.ext_v1.1$cluster.v06.26.re
src.endoderm.mg3.ext_v1.1@meta.data[
  src.endoderm.mg3.ext_v1.1$cluster.v06.26.re%in%c("EP.1","EP.2"),]$cluster.v06.26.re_ep = "EP"
#------------------------------------

# FDL
src.endoderm.mg3.ext_v1.1 = 
  seurat_fdl("src.endoderm.mg3.ext_v1.1", "pca","RNA")
src.endoderm.mg3.ext_v1.1 = 
  seurat_fdl("src.endoderm.mg3.ext_v1.1", "mnn","RNA")
src.endoderm.mg3.ext_v1.1 = 
  seurat_fdl("src.endoderm.mg3.ext_v1.1", "pca_integrated","integrated")

# Graph
pdf("figure.v08.07/organ_development_re/mg3_summary_graph.re_snn_PI.pdf",12,7)
# SNN
for(red in c('pca',"mnn","pca_integrated")){
  
  assay = ifelse(red=="mnn","RNA","integrated")
  
  graph.mg3.ext_v1.1 = 
    seurat_GCN(seurat_name = "src.endoderm.mg3.ext_v1.1",
               src.endoderm.mg3.ext_v1.1.gene.fin,
               reduction = red, assay =  assay)
  graph.mg3_cluster = cluster_walktrap(graph.mg3.ext_v1.1, steps = 1)
  
  for(color in c("cluster.extract.v1.1","cluster.v06.26.re","cluster.v06.26.re_ep",
                 "Time")){
    for(edge.color in c(NA,"lightgray")){
      for(edge.color in c(NA,"lightgray")){
        set.seed(1)
        plot.igraph(graph.mg3.ext_v1.1,
                    layout = layout_with_fr(graph.mg3.ext_v1.1),
                    edge.color = edge.color,
                    vertex.size=2.5, vertex.label=NA,
                    vertex.label.color = "black", 
                    vertex.frame.color = NA, vertex.frame.width = 0.5,
                    #vertex.color =  colors.time[src.endoderm.mg3.ext_v1.1[["Time"]][graph.mg3_cluster$names,]],
                    #vertex.color =  cluster.endoderm.color.v5[src.endoderm.mg3.ext_v1.1[["cluster.v06.26.re"]][graph.mg3_cluster$names,]],
                    vertex.color =  colors_bg[src.endoderm.mg3.ext_v1.1[[color]][graph.mg3_cluster$names,]])
      }}}
  
}
dev.off()

pdf("figure.v08.07/organ_development_re/mg3_summary_graph.re_other.pdf",12,7)
# Reduction
for(red in c("umap_integrated","umap_mnn",
             "fdl_pca",'fdl_mnn',"fdl_pca_integrated")){
  p1 = DimPlot(src.endoderm.mg3.ext_v1.1,group.by = 'cluster.v06.26.re', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p2 = DimPlot(src.endoderm.mg3.ext_v1.1,group.by = 'cluster.extract.v1.1', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p3 = DimPlot(src.endoderm.mg3.ext_v1.1,group.by = 'Time', 
               reduction = red , cols = colors.time,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p4 = DimPlot(src.endoderm.mg3.ext_v1.1,group.by = 'cluster.v06.26.re_ep', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  print(p1); print(p2); print(p3); print(p4)
}
dev.off()


pdf("figure.v08.07/organ_development_re/mg3_summary_proportion.re.pdf",12,7)
#-----------------------------
cell_type = c("MG.3","MG.3.A/M","MG.3.A/M","MG.3.P","Stomach","DP","EP","EP.1","EP.2","Small.intestine.1","Small.intestine.2")
meta_mg3 = src.endoderm.mg3.ext_v1.1@meta.data
meta_mg3$Time = factor(meta_mg3$Time, levels = names(colors.time))
meta_mg3$cluster.extract.v1.1 = factor(meta_mg3$cluster.extract.v1.1, 
                                       levels = c("MG.3","FG.4-MG.1/3","MG.2/3","MG.1/2/3",
                                                  "AL.3-MG.1/2/3","MG.2/3-HG.1"))
meta_mg3$cluster.v06.26.re = factor(meta_mg3$cluster.v06.26.re, levels = rev(cell_type))

meta_mg3_B0 = meta_mg3
meta_mg3_B1 = meta_mg3[meta_mg3$batch%in%1&!meta_mg3$Time%in%"ss9",]
meta_mg3_B2 = meta_mg3[meta_mg3$batch%in%2,]

for(data in c("meta_mg3_B0","meta_mg3_B1","meta_mg3_B2")){
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
                                 group = cluster.v06.26.re_ep, fill = cluster.v06.26.re_ep),
                   stat = "count", position = 'fill') + 
          xlab("Time") + ylab("Cluster proportion of
  cell type")+ ggtitle(data))
  
}
#-----------------------------
dev.off()

save(src.endoderm.mg3.ext_v1.1,
     file = "figure.v08.07/organ_development_re/src.endoderm.mg3.ext_v1.1.re.Rdata")


# Marker
#------------------------
src.endoderm.mg3.ext_v1.1 = 
  SetIdent(src.endoderm.mg3.ext_v1.1,
           value = src.endoderm.mg3.ext_v1.1$cluster.v06.26.re)
markergene.mg3.ext_v1.1 = 
  FindAllMarkers(src.endoderm.mg3.ext_v1.1, assay = "RNA")
markergene.mg3.ext_v1.1 = markergene.mg3.ext_v1.1[
  markergene.mg3.ext_v1.1$avg_log2FC>0.2&
    markergene.mg3.ext_v1.1$p_val_adj<0.1, "gene"]

#------------------------
mg3.umap.embedding = src.endoderm.mg3.ext_v1.1[["fdl_pca_integrated"]]@cell.embeddings[,1:3]
mg3.umap.embedding = src.endoderm.mg3.ext_v1.1[["fdl_mnn"]]@cell.embeddings[,1:3]
mg3.umap.embedding = src.endoderm.mg3.ext_v1.1[["umap_mnn"]]@cell.embeddings[,1:3]
mg3.umap.embedding = src.endoderm.mg3.ext_v1.1[["umap_integrated"]]@cell.embeddings[,1:3]


# Cell order :: Princurve
#---------------------------
# MG.3
cellorder.mg3_mg3 = 
  rownames(src.endoderm.mg3.ext_v1.1@meta.data[
    src.endoderm.mg3.ext_v1.1$cluster.v06.26.re%in%"MG.3",])
princurve.mg3_mg3 = princurve::principal_curve(
  x = mg3.umap.embedding[cellorder.mg3_mg3,],  smoother = "smooth.spline")
cellorder.mg3_mg3 =  names(
  princurve.mg3_mg3$lambda[cellorder.mg3_mg3][order(princurve.mg3_mg3$lambda[cellorder.mg3_mg3])])
# MG.3.A/M
cellorder.mg3_mg3a = 
  rownames(src.endoderm.mg3.ext_v1.1@meta.data[
    src.endoderm.mg3.ext_v1.1$cluster.v06.26.re%in%"MG.3.A/M",])
princurve.mg3_mg3a = princurve::principal_curve(
  x = mg3.umap.embedding[cellorder.mg3_mg3a,],  smoother = "smooth.spline")
cellorder.mg3_mg3a =  names(
  princurve.mg3_mg3a$lambda[cellorder.mg3_mg3a][order(princurve.mg3_mg3a$lambda[cellorder.mg3_mg3a])])
# MG.3.A/M
cellorder.mg3_mg3m = 
  rownames(src.endoderm.mg3.ext_v1.1@meta.data[
    src.endoderm.mg3.ext_v1.1$cluster.v06.26.re%in%"MG.3.A/M",])
princurve.mg3_mg3m = princurve::principal_curve(
  x = mg3.umap.embedding[cellorder.mg3_mg3m,],  smoother = "smooth.spline")
cellorder.mg3_mg3m =  names(
  princurve.mg3_mg3m$lambda[cellorder.mg3_mg3m][order(princurve.mg3_mg3m$lambda[cellorder.mg3_mg3m])])
# MG.3.P
cellorder.mg3_mg3p = 
  rownames(src.endoderm.mg3.ext_v1.1@meta.data[
    src.endoderm.mg3.ext_v1.1$cluster.v06.26.re%in%"MG.3.P",])
princurve.mg3_mg3p = princurve::principal_curve(
  x = mg3.umap.embedding[cellorder.mg3_mg3p,],  smoother = "smooth.spline")
cellorder.mg3_mg3p =  names(
  princurve.mg3_mg3p$lambda[cellorder.mg3_mg3p][order(princurve.mg3_mg3p$lambda[cellorder.mg3_mg3p])])
# Stomach
cellorder.mg3_sto = 
  rownames(src.endoderm.mg3.ext_v1.1@meta.data[
    src.endoderm.mg3.ext_v1.1$cluster.v06.26.re%in%"Stomach",])
princurve.mg3_sto = princurve::principal_curve(
  x = mg3.umap.embedding[cellorder.mg3_sto,],  smoother = "smooth.spline")
cellorder.mg3_sto =  names(
  princurve.mg3_sto$lambda[cellorder.mg3_sto][order(princurve.mg3_sto$lambda[cellorder.mg3_sto])])
# DP
cellorder.mg3_dp = 
  rownames(src.endoderm.mg3.ext_v1.1@meta.data[
    src.endoderm.mg3.ext_v1.1$cluster.v06.26.re%in%"DP",])
princurve.mg3_dp = princurve::principal_curve(
  x = mg3.umap.embedding[cellorder.mg3_dp,],  smoother = "smooth.spline")
cellorder.mg3_dp =  names(
  princurve.mg3_dp$lambda[cellorder.mg3_dp][order(princurve.mg3_dp$lambda[cellorder.mg3_dp])])
# EP.1
cellorder.mg3_ep1 = 
  rownames(src.endoderm.mg3.ext_v1.1@meta.data[
    src.endoderm.mg3.ext_v1.1$cluster.v06.26.re%in%"EP.1",])
princurve.mg3_ep1 = princurve::principal_curve(
  x = mg3.umap.embedding[cellorder.mg3_ep1,],  smoother = "smooth.spline")
cellorder.mg3_ep1 =  names(
  princurve.mg3_ep1$lambda[cellorder.mg3_ep1][order(princurve.mg3_ep1$lambda[cellorder.mg3_ep1])])
# EP.2
cellorder.mg3_ep2 = 
  rownames(src.endoderm.mg3.ext_v1.1@meta.data[
    src.endoderm.mg3.ext_v1.1$cluster.v06.26.re%in%"EP.2",])
princurve.mg3_ep2 = princurve::principal_curve(
  x = mg3.umap.embedding[cellorder.mg3_ep2,],  smoother = "smooth.spline")
cellorder.mg3_ep2 =  names(
  princurve.mg3_ep2$lambda[cellorder.mg3_ep2][order(princurve.mg3_ep2$lambda[cellorder.mg3_ep2])])
# Small.intestine.1
cellorder.mg3_s1 = 
  rownames(src.endoderm.mg3.ext_v1.1@meta.data[
    src.endoderm.mg3.ext_v1.1$cluster.v06.26.re%in%"Small.intestine.1",])
princurve.mg3_s1 = princurve::principal_curve(
  x = mg3.umap.embedding[cellorder.mg3_s1,],  smoother = "smooth.spline")
cellorder.mg3_s1 =  names(
  princurve.mg3_s1$lambda[cellorder.mg3_s1][order(princurve.mg3_s1$lambda[cellorder.mg3_s1])])
# Small.intestine.2
cellorder.mg3_s2 = 
  rownames(src.endoderm.mg3.ext_v1.1@meta.data[
    src.endoderm.mg3.ext_v1.1$cluster.v06.26.re%in%"Small.intestine.2",])
princurve.mg3_s2 = princurve::principal_curve(
  x = mg3.umap.embedding[cellorder.mg3_s2,],  smoother = "smooth.spline")
cellorder.mg3_s2 =  names(
  princurve.mg3_s2$lambda[cellorder.mg3_s2][order(princurve.mg3_s2$lambda[cellorder.mg3_s2])])
#---------------------------

pdf("figure.v08.07/try.pdf")
#----------------------------
plot(c(1:length(cellorder.mg3_mg3)),
     c(1:length(cellorder.mg3_mg3)),
     col = colors.time[src.endoderm.mg3.ext_v1.1[,cellorder.mg3_mg3]$Time])
plot(c(1:length(cellorder.mg3_mg3a)),
     c(1:length(cellorder.mg3_mg3a)),
     col = colors.time[src.endoderm.mg3.ext_v1.1[,cellorder.mg3_mg3a]$Time])
plot(c(1:length(cellorder.mg3_mg3m)),
     c(1:length(cellorder.mg3_mg3m)),
     col = colors.time[src.endoderm.mg3.ext_v1.1[,cellorder.mg3_mg3m]$Time])
plot(c(1:length(cellorder.mg3_mg3p)),
     c(1:length(cellorder.mg3_mg3p)),
     col = colors.time[src.endoderm.mg3.ext_v1.1[,cellorder.mg3_mg3p]$Time])
plot(c(1:length(cellorder.mg3_sto)),
     c(1:length(cellorder.mg3_sto)),
     col = colors.time[src.endoderm.mg3.ext_v1.1[,cellorder.mg3_sto]$Time])
plot(c(1:length(cellorder.mg3_dp)),
     c(1:length(cellorder.mg3_dp)),
     col = colors.time[src.endoderm.mg3.ext_v1.1[,cellorder.mg3_dp]$Time])
plot(c(1:length(cellorder.mg3_ep1)),
     c(1:length(cellorder.mg3_ep1)),
     col = colors.time[src.endoderm.mg3.ext_v1.1[,cellorder.mg3_ep1]$Time])
plot(c(1:length(cellorder.mg3_ep2)),
     c(1:length(cellorder.mg3_ep2)),
     col = colors.time[src.endoderm.mg3.ext_v1.1[,cellorder.mg3_ep2]$Time])
plot(c(1:length(cellorder.mg3_s1)),
     c(1:length(cellorder.mg3_s1)),
     col = colors.time[src.endoderm.mg3.ext_v1.1[,cellorder.mg3_s1]$Time])
plot(c(1:length(cellorder.mg3_s2)),
     c(1:length(cellorder.mg3_s2)),
     col = colors.time[src.endoderm.mg3.ext_v1.1[,cellorder.mg3_s2]$Time])
#---------------------------
dev.off()

cellorder.mg3_mg3_FM = cellorder.mg3_mg3
cellorder.mg3_mg3a_FM = cellorder.mg3_mg3a
cellorder.mg3_mg3m_FM = cellorder.mg3_mg3m
cellorder.mg3_mg3p_FM = cellorder.mg3_mg3p
cellorder.mg3_sto_FM = cellorder.mg3_sto
cellorder.mg3_dp_FM = cellorder.mg3_dp
cellorder.mg3_ep1_FM = cellorder.mg3_ep1
cellorder.mg3_ep2_FM = cellorder.mg3_ep2
cellorder.mg3_s1_FM = cellorder.mg3_s1
cellorder.mg3_s2_FM = cellorder.mg3_s2


src.endoderm.mg3.ext_v1.1 = 
  SetIdent(src.endoderm.mg3.ext_v1.1,
           value = src.endoderm.mg3.ext_v1.1$cluster.v06.26.re)
#=================================================================================
# Marker: MG.3.A/M
#=================================================================================
pdf("figure.v08.07/organ_development_re/mg3_heatmap_marker.re_mg3a.pdf",6,9)
cellorder.mg3 = c(cellorder.mg3_mg3_FM, cellorder.mg3_mg3a_FM, cellorder.mg3_sto_FM)

markergene.mg3.ext_v1.1 = 
  FindAllMarkers(src.endoderm.mg3.ext_v1.1[,cellorder.mg3], assay = "RNA")
markergene.mg3.ext_v1.1 = markergene.mg3.ext_v1.1[
  markergene.mg3.ext_v1.1$avg_log2FC>0.2&markergene.mg3.ext_v1.1$p_val_adj<0.1, "gene"]
# selectgene = unique(c(markergene.mg3.ext_v1.1))
selectgene = src.endoderm.mg3.ext_v1.1.filtergene_fin
#-------------------
mg3.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.mg3.ext_v1.1@assays$RNA@data[selectgene, cellorder.mg3]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.mg3.ext_v1.1$Time[cellorder.mg3], colors.time),
              MyName2Col(src.endoderm.mg3.ext_v1.1$cluster.v06.26.re[cellorder.mg3], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.mg3.ext_v1.1$cluster.extract.v1.1[cellorder.mg3], cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            # return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.mg3.ext_v1.1@meta.data[cell_mg3_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

mg3.ext_v1.1.re.rowtree.1 = as.dendrogram(mg3.ext_v1.1.re.rowtree)
src.endoderm.mg3.ext_v1.1.filtergene_fin = c(
  setdiff(
    c(
      rev(labels(mg3.ext_v1.1.re.rowtree.1[[2]][[1]][[1]][[1]])),
      labels(mg3.ext_v1.1.re.rowtree.1[[2]][[1]][[2]][[1]]),
      labels(mg3.ext_v1.1.re.rowtree.1[[2]][[1]][[2]][[2]]),
      labels(mg3.ext_v1.1.re.rowtree.1[[2]][[1]][[1]][[2]])
    ),
    c(
      labels(mg3.ext_v1.1.re.rowtree.1[[2]][[1]][[1]][[2]][[2]][[2]][[2]]),
      labels(mg3.ext_v1.1.re.rowtree.1[[2]][[1]][[1]][[2]][[2]][[2]][[1]][[2]]),
      c("Pdia3","Calr","Hspa5","Pdia6","Hsp90b1",'Emb',"Slc2a3","Ssr4","Erp29","Ssr2","Ppib")
    )
  ),
  setdiff(
    c(
      labels(mg3.ext_v1.1.re.rowtree.1[[2]][[2]])
    ),
    c(
      c("Myl6","Rbp1")
    )
  ),
  setdiff(
    c(
      labels(mg3.ext_v1.1.re.rowtree.1[[1]])
    ),
    c()
  )
)

names(src.endoderm.mg3.ext_v1.1.filtergene_fin) = c(
  rep(4,length(  
    setdiff(
      c(
        labels(mg3.ext_v1.1.re.rowtree.1[[2]][[1]][[1]][[1]]),
        labels(mg3.ext_v1.1.re.rowtree.1[[2]][[1]][[2]][[1]]),
        labels(mg3.ext_v1.1.re.rowtree.1[[2]][[1]][[2]][[2]]),
        labels(mg3.ext_v1.1.re.rowtree.1[[2]][[1]][[1]][[2]])
      ),
      c(
        labels(mg3.ext_v1.1.re.rowtree.1[[2]][[1]][[1]][[2]][[2]][[2]][[2]]),
        labels(mg3.ext_v1.1.re.rowtree.1[[2]][[1]][[1]][[2]][[2]][[2]][[1]][[2]]),
        c("Pdia3","Calr","Hspa5","Pdia6","Hsp90b1",'Emb',"Slc2a3","Ssr4","Erp29","Ssr2","Ppib")
      )
    ))),
  rep(3,length(  
    setdiff(
      c(
        labels(mg3.ext_v1.1.re.rowtree.1[[2]][[2]])
      ),
      c(
        c("Myl6","Rbp1")
      )
    ))),
  rep(7,length(  
    setdiff(
      c(
        labels(mg3.ext_v1.1.re.rowtree.1[[1]])
      ),
      c()
    )))
)

src.endoderm.mg3.ext_v1.1.filtergene_mg3a =
  src.endoderm.mg3.ext_v1.1.filtergene_fin

#=================================================================================
# Marker: MG.3.A/M: DP / EP
#=================================================================================
pdf("figure.v08.07/organ_development_re/mg3_heatmap_marker.re_mg3m.pdf",6,9)
cellorder.mg3 = c(cellorder.mg3_mg3_FM, cellorder.mg3_mg3m_FM, 
                  cellorder.mg3_dp_FM)
#cellorder.mg3_ep1_FM, cellorder.mg3_ep2_FM

markergene.mg3.ext_v1.1 = 
  FindAllMarkers(src.endoderm.mg3.ext_v1.1[,cellorder.mg3], assay = "RNA")
markergene.mg3.ext_v1.1 = markergene.mg3.ext_v1.1[
  markergene.mg3.ext_v1.1$avg_log2FC>0.2&markergene.mg3.ext_v1.1$p_val_adj<0.1, "gene"]
selectgene = unique(c(markergene.mg3.ext_v1.1))

selectgene = src.endoderm.mg3.ext_v1.1.filtergene_fin
#-------------------
mg3.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.mg3.ext_v1.1@assays$RNA@data[selectgene, cellorder.mg3]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.mg3.ext_v1.1$Time[cellorder.mg3], colors.time),
              MyName2Col(src.endoderm.mg3.ext_v1.1$cluster.v06.26.re[cellorder.mg3], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.mg3.ext_v1.1$cluster.extract.v1.1[cellorder.mg3], cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            # return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.mg3.ext_v1.1@meta.data[cell_mg3_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

mg3.ext_v1.1.re.rowtree.2 = as.dendrogram(mg3.ext_v1.1.re.rowtree)
src.endoderm.mg3.ext_v1.1.filtergene_fin = c(
  setdiff(
    c(
      rev(labels(mg3.ext_v1.1.re.rowtree.2[[1]]))
    ),
    c()
  ),
  setdiff(
    c(
      rev(labels(mg3.ext_v1.1.re.rowtree.2[[2]][[1]]))
    ),
    c()
  ),
  setdiff(
    c(
      labels(mg3.ext_v1.1.re.rowtree.2[[2]][[2]][[2]][[1]]),
      labels(mg3.ext_v1.1.re.rowtree.2[[2]][[2]][[2]][[2]]),
      labels(mg3.ext_v1.1.re.rowtree.2[[2]][[2]][[1]])
    ),
    c()
  )
)

names(src.endoderm.mg3.ext_v1.1.filtergene_fin) = c(
  rep(4,length(  
    setdiff(
      c(
        rev(labels(mg3.ext_v1.1.re.rowtree.2[[1]]))
      ),
      c()
    ))),
  rep(3,length(  
    setdiff(
      c(
        labels(mg3.ext_v1.1.re.rowtree.2[[2]][[1]])
      ),
      c()
    ))),
  rep(7,length(  
    setdiff(
      c(
        rev(labels(mg3.ext_v1.1.re.rowtree.2[[2]][[2]]))
      ),
      c()
    )))
)

src.endoderm.mg3.ext_v1.1.filtergene_mg3mdp =
  src.endoderm.mg3.ext_v1.1.filtergene_fin

#=================================================================================
# Marker: MG.3.A/M: DP / EP
#=================================================================================
pdf("figure.v08.07/organ_development_re/mg3_heatmap_marker.re_mg3m_ep.pdf",6,9)
cellorder.mg3 = c(cellorder.mg3_mg3_FM, cellorder.mg3_mg3m_FM, 
                  cellorder.mg3_ep1_FM,  cellorder.mg3_ep2_FM)
#cellorder.mg3_ep1_FM, cellorder.mg3_ep2_FM

markergene.mg3.ext_v1.1 = 
  FindAllMarkers(src.endoderm.mg3.ext_v1.1[,cellorder.mg3], assay = "RNA")
markergene.mg3.ext_v1.1 = markergene.mg3.ext_v1.1[
  markergene.mg3.ext_v1.1$avg_log2FC>0.2&markergene.mg3.ext_v1.1$p_val_adj<0.1, "gene"]
selectgene = unique(c(markergene.mg3.ext_v1.1))

selectgene = src.endoderm.mg3.ext_v1.1.filtergene_mg3mep
#-------------------
mg3.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.mg3.ext_v1.1@assays$RNA@data[selectgene, cellorder.mg3]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.mg3.ext_v1.1$Time[cellorder.mg3], colors.time),
              MyName2Col(src.endoderm.mg3.ext_v1.1$cluster.v06.26.re[cellorder.mg3], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.mg3.ext_v1.1$cluster.extract.v1.1[cellorder.mg3], cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            # return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.mg3.ext_v1.1@meta.data[cell_mg3_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

mg3.ext_v1.1.re.rowtree.3 = as.dendrogram(mg3.ext_v1.1.re.rowtree)
src.endoderm.mg3.ext_v1.1.filtergene_fin = c(
  setdiff(
    c(
      labels(mg3.ext_v1.1.re.rowtree.3[[1]])
    ),
    c()
  ),
  setdiff(
    c(
      labels(mg3.ext_v1.1.re.rowtree.3[[2]][[1]])
    ),
    c()
  ),
  setdiff(
    c(
      labels(mg3.ext_v1.1.re.rowtree.3[[2]][[2]][[2]][[2]]),
      labels(mg3.ext_v1.1.re.rowtree.3[[2]][[2]][[1]][[1]])
    ),
    c()
  ),
  setdiff(
    c(
      labels(mg3.ext_v1.1.re.rowtree.3[[2]][[2]][[1]][[2]])
    ),
    c()
  ),
  setdiff(
    c(
      labels(mg3.ext_v1.1.re.rowtree.3[[2]][[2]][[2]][[1]])
    ),
    c()
  )
)

names(src.endoderm.mg3.ext_v1.1.filtergene_fin) = c(
  rep(4,length(  
    setdiff(
      c(
        labels(mg3.ext_v1.1.re.rowtree.3[[1]])
      ),
      c()
    ))),
  rep(3,length(  
    setdiff(
      c(
        labels(mg3.ext_v1.1.re.rowtree.3[[2]][[1]])
      ),
      c()
    ))),
  rep(9,length(  
    setdiff(
      c(
        labels(mg3.ext_v1.1.re.rowtree.3[[2]][[2]][[2]][[2]]),
        labels(mg3.ext_v1.1.re.rowtree.3[[2]][[2]][[1]][[1]])
      ),
      c()
    ))),
  rep(2,length(  
    setdiff(
      c(
        labels(mg3.ext_v1.1.re.rowtree.3[[2]][[2]][[1]][[2]])
      ),
      c()
    ))),
  rep(7,length(  
    setdiff(
      c(
        labels(mg3.ext_v1.1.re.rowtree.3[[2]][[2]][[2]][[1]])
      ),
      c()
    )))
)

src.endoderm.mg3.ext_v1.1.filtergene_mg3mep =
  src.endoderm.mg3.ext_v1.1.filtergene_fin

#=================================================================================
# Marker: MG.3.p
#=================================================================================
pdf("figure.v08.07/organ_development_re/mg3_heatmap_marker.re_mg3p.pdf",6,9)
cellorder.mg3 = c(cellorder.mg3_mg3_FM, cellorder.mg3_mg3p_FM, 
                  rev(cellorder.mg3_s2_FM), cellorder.mg3_s1_FM)
#cellorder.mg3_ep1_FM, cellorder.mg3_ep2_FM

markergene.mg3.ext_v1.1 = 
  FindAllMarkers(src.endoderm.mg3.ext_v1.1[,cellorder.mg3], assay = "RNA")
markergene.mg3.ext_v1.1 = markergene.mg3.ext_v1.1[
  markergene.mg3.ext_v1.1$avg_log2FC>0.2&markergene.mg3.ext_v1.1$p_val_adj<0.1, "gene"]
selectgene = unique(c(markergene.mg3.ext_v1.1))

selectgene = src.endoderm.mg3.ext_v1.1.filtergene_fin
#-------------------
mg3.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.mg3.ext_v1.1@assays$RNA@data[selectgene, cellorder.mg3]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.mg3.ext_v1.1$Time[cellorder.mg3], colors.time),
              MyName2Col(src.endoderm.mg3.ext_v1.1$cluster.v06.26.re[cellorder.mg3], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.mg3.ext_v1.1$cluster.extract.v1.1[cellorder.mg3], cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            # return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.mg3.ext_v1.1@meta.data[cell_mg3_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

mg3.ext_v1.1.re.rowtree.4 = as.dendrogram(mg3.ext_v1.1.re.rowtree)
src.endoderm.mg3.ext_v1.1.filtergene_fin = c(
  setdiff(
    c(
      rev(labels(mg3.ext_v1.1.re.rowtree.4[[1]][[2]][[1]])),
      labels(mg3.ext_v1.1.re.rowtree.4[[1]][[2]][[2]])
    ),
    c(
      labels(mg3.ext_v1.1.re.rowtree.4[[1]][[2]][[1]][[2]]),
      labels(mg3.ext_v1.1.re.rowtree.4[[1]][[2]][[2]][[1]]),
      labels(mg3.ext_v1.1.re.rowtree.4[[1]][[2]][[1]][[1]][[1]])
    )
  ),
  setdiff(
    c(
      rev(labels(mg3.ext_v1.1.re.rowtree.4[[2]][[1]]))
    ),
    c()
  ),
  setdiff(
    c(
      labels(mg3.ext_v1.1.re.rowtree.4[[2]][[2]][[2]][[2]])
    ),
    c()
  ),
  setdiff(
    c(
      labels(mg3.ext_v1.1.re.rowtree.4[[2]][[2]][[2]][[1]]),
      labels(mg3.ext_v1.1.re.rowtree.4[[2]][[2]][[1]])
    ),
    c()
  )
)

names(src.endoderm.mg3.ext_v1.1.filtergene_fin) = c(
  rep(4,length(  
    setdiff(
      c(
        rev(labels(mg3.ext_v1.1.re.rowtree.4[[1]][[2]][[1]])),
        labels(mg3.ext_v1.1.re.rowtree.4[[1]][[2]][[2]])
      ),
      c(
        labels(mg3.ext_v1.1.re.rowtree.4[[1]][[2]][[1]][[2]]),
        labels(mg3.ext_v1.1.re.rowtree.4[[1]][[2]][[2]][[1]]),
        labels(mg3.ext_v1.1.re.rowtree.4[[1]][[2]][[1]][[1]][[1]])
      )
    ))),
  rep(3,length(  
    setdiff(
      c(
        rev(labels(mg3.ext_v1.1.re.rowtree.4[[2]][[1]]))
      ),
      c()
    ))),
  rep(9,length(  
    setdiff(
      c(
        labels(mg3.ext_v1.1.re.rowtree.4[[2]][[2]][[2]][[2]])
      ),
      c()
    ))),
  rep(7,length(  
    setdiff(
      c(
        labels(mg3.ext_v1.1.re.rowtree.4[[2]][[2]][[1]]),
        labels(mg3.ext_v1.1.re.rowtree.4[[2]][[2]][[2]][[1]])
      ),
      c()
    )))
)

src.endoderm.mg3.ext_v1.1.filtergene_mg3p =
  src.endoderm.mg3.ext_v1.1.filtergene_fin

#--------------------------------------------------------------------------------


save(
  # Gene 
  src.endoderm.mg3.ext.v1.1.selectgene,
  src.endoderm.mg3.ext_v1.1.gene.fin,
  src.endoderm.mg3.ext_v1.1.filtergene_fin,
  src.endoderm.mg3.ext_v1.1.filtergene_mg3a,
  src.endoderm.mg3.ext_v1.1.filtergene_mg3mdp,
  src.endoderm.mg3.ext_v1.1.filtergene_mg3mep,
  src.endoderm.mg3.ext_v1.1.filtergene_mg3p,
  # Cell type
  cellorder.mg3_mg3_FM, cellorder.mg3_mg3a_FM, cellorder.mg3_sto_FM,
  cellorder.mg3_mg3m_FM, cellorder.mg3_dp_FM, 
  cellorder.mg3_ep1_FM, cellorder.mg3_ep2_FM,
  cellorder.mg3_mg3p_FM, cellorder.mg3_s1_FM, cellorder.mg3_s2_FM,
  file= "figure.v08.07/organ_development_re/mg3_heatmap_parameter.Rdata")

