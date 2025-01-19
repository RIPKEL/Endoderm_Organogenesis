#----------------------------------------------------
#>>> AL.1/2
#----------------------------------------------------

src.endoderm.al12.ext_v1.1 = FindVariableFeatures(src.endoderm.al12.ext_v1.1, nfeatures = 2000)
src.endoderm.al12.ext_v1.1 = ScaleData(src.endoderm.al12.ext_v1.1, split.by = "batch_phase",
                                       features = rownames(src.endoderm.al12.ext_v1.1))
src.endoderm.al12.ext_v1.1.filtergene = 
  Myfilter(as.matrix(src.endoderm.al12.ext_v1.1@assays$RNA@data),
           gene = src.endoderm.al12.ext_v1.1@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
#rm(src.endoderm.al12.ext_v1.1.filtergene)
# src.endoderm.al12.ext.v1.1.selectgene = unique(c(select_gene_al12))

#> Row
pdf("figure.v08.07/try.pdf")
#-----------------
al12.ext.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.al12.ext_v1.1@assays$RNA@data[
    src.endoderm.al12.ext_v1.1.gene.fin ,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.endoderm.al12.ext_v1.1$Time,colors.time),
      MyName2Col(src.endoderm.al12.ext_v1.1$cluster.extract.v1.1,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.al12.ext_v1.1$cluster.v06.26.re,cluster.endoderm.color.v5)
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
plot(al12.ext.re.rowtree)
dev.off()

al12.ext.re.rowtree_1 = as.dendrogram(al12.ext.re.rowtree)
al12.ext.re.rowtree_2 = as.dendrogram(al12.ext.re.rowtree)
src.endoderm.al12.ext_v1.1.gene.fin = setdiff(
  src.endoderm.al12.ext.v1.1.selectgene,
  c(labels(al12.ext.re.rowtree_1[[2]][[1]]),
    labels(al12.ext.re.rowtree_1[[1]][[2]][[1]][[2]]),
    labels(al12.ext.re.rowtree_1[[2]][[2]][[2]][[1]]),
    labels(al12.ext.re.rowtree_2[[2]][[2]][[1]])))

#--------------------------------------------------------------------------------
src.endoderm.al12.ext_v1.1 = RunPCA(src.endoderm.al12.ext_v1.1,
                                    features = src.endoderm.al12.ext_v1.1.gene.fin)
src.endoderm.al12.ext_v1.1 = 
  seurat_mnn(seurat_name = "src.endoderm.al12.ext_v1.1", dims = 1:30, 
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.al12.ext_v1.1.gene.fin)
src.endoderm.al12.ext_v1.1 = 
  seurat_int(seurat_name = "src.endoderm.al12.ext_v1.1", dims = 1:30,
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.al12.ext_v1.1.gene.fin)

#---------------
src.endoderm.al12.ext_v1.1 = 
  FindNeighbors(src.endoderm.al12.ext_v1.1, dims=1:30, 
                reduction = "pca_integrated", assay = "integrated")
src.endoderm.al12.ext_v1.1 = 
  FindNeighbors(src.endoderm.al12.ext_v1.1, dims=1:30, 
                reduction = "pca", assay = "RNA")
src.endoderm.al12.ext_v1.1 = 
  FindClusters(src.endoderm.al12.ext_v1.1, resolution = 2)
#---------------

DimPlot(src.endoderm.al12.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.2', reduction = "umap_integrated")

DimPlot(src.endoderm.al12.ext_v1.1, group.by = 'Time',
        reduction = "umap_mnn", cols = colors.time)
DimPlot(src.endoderm.al12.ext_v1.1, group.by = 'Time', 
        reduction = "umap_integrated", cols = colors.time)
DimPlot(src.endoderm.al12.ext_v1.1, group.by = 'batch', 
        reduction = "umap_integrated") #, cols = colors.time)

DimPlot(src.endoderm.al12.ext_v1.1, group.by = 'cluster.extract.v1.1', 
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.al12.ext_v1.1, group.by = 'cluster.extract.v1.1', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.al12.ext_v1.1, group.by = 'cluster.v06.26.re', 
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.al12.ext_v1.1, group.by = 'cluster.v06.26.re', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
FeaturePlot(src.endoderm.al12.ext_v1.1, features = 'Xist', 
            reduction = "umap_integrated")


src.endoderm.al12.ext_v1.1@reductions$pca_integrated@cell.embeddings = 
  src.endoderm.al12.ext_v1.1@reductions$pca_integrated@cell.embeddings[colnames(src.endoderm.al12.ext_v1.1),]

src.endoderm.al12.ext_v1.1$cluster.v06.26.re = NA
src.endoderm.al12.ext_v1.1@meta.data[
  intersect(colnames(src.endoderm.al12.ext_v1.1),colnames(src.endoderm.al12.ext.re)),]$cluster.v06.26.re = 
  src.endoderm.al12.ext.re@meta.data[
    intersect(colnames(src.endoderm.al12.ext_v1.1),colnames(src.endoderm.al12.ext.re)),]$cluster.v06.26.re
a = FNN::knn(
  src.endoderm.al12.ext_v1.1@reductions$pca_integrated@cell.embeddings[
    rownames(src.endoderm.al12.ext_v1.1@meta.data[!src.endoderm.al12.ext_v1.1$cluster.v06.26.re%in%NA,]),],
  src.endoderm.al12.ext_v1.1@reductions$pca_integrated@cell.embeddings[
    rownames(src.endoderm.al12.ext_v1.1@meta.data[src.endoderm.al12.ext_v1.1$cluster.v06.26.re%in%NA,]),],
  src.endoderm.al12.ext_v1.1@meta.data[!src.endoderm.al12.ext_v1.1$cluster.v06.26.re%in%NA,]$cluster.v06.26.re, k = 10)
src.endoderm.al12.ext_v1.1@meta.data[
  src.endoderm.al12.ext_v1.1$cluster.v06.26.re%in%NA,]$cluster.v06.26.re = as.character(a)

src.endoderm.al12.ext_v1.1@meta.data[
  src.endoderm.al12.ext_v1.1$cluster.extract.v1.1%in%c("FG.4-AL.1/2/3"),]$cluster.v06.26.re = "Liver"
src.endoderm.al12.ext_v1.1@meta.data[
  (src.endoderm.al12.ext_v1.1$RNA_snn_res.2%in%c(2,4,10,1) #&
    #src.endoderm.al12.ext_v1.1$cluster.extract.v1.1%in%c("AL.1/2")
   )|
    src.endoderm.al12.ext_v1.1$cluster.v06.26.re%in%"AL.1/2",]$cluster.v06.26.re = "AL.1/2-Liver"
src.endoderm.al12.ext_v1.1@meta.data[
  src.endoderm.al12.ext_v1.1$RNA_snn_res.2%in%c(12,16),]$cluster.v06.26.re = "AL.2"
DimPlot(src.endoderm.al12.ext_v1.1, group.by = 'cluster.v06.26.re', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)


# FDL
src.endoderm.al12.ext_v1.1 = 
  seurat_fdl("src.endoderm.al12.ext_v1.1", "pca","RNA")
src.endoderm.al12.ext_v1.1 = 
  seurat_fdl("src.endoderm.al12.ext_v1.1", "mnn","RNA")
src.endoderm.al12.ext_v1.1 = 
  seurat_fdl("src.endoderm.al12.ext_v1.1", "pca_integrated","integrated")


# Graph
pdf("figure.v08.07/organ_development_re/al12_summary_graph.re_snn_PI.pdf",12,7)
# SNN
for(red in c("mnn","pca_integrated")){
  
  assay = ifelse(red=="mnn","RNA","integrated")
  
  graph.al12.ext_v1.1 = 
    seurat_GCN(seurat_name = "src.endoderm.al12.ext_v1.1",
               src.endoderm.al12.ext_v1.1.gene.fin,
               reduction = red, assay =  assay)
  graph.al12_cluster = cluster_walktrap(graph.al12.ext_v1.1, steps = 1)
  
  for(color in c("cluster.extract.v1.1","cluster.v06.26.re",
                 "cluster.v06.26.re_hc",
                 "Time")){
    for(edge.color in c(NA,"lightgray")){
      for(edge.color in c(NA,"lightgray")){
        set.seed(1)
        plot.igraph(graph.al12.ext_v1.1,
                    layout = layout_with_fr(graph.al12.ext_v1.1),
                    edge.color = edge.color,
                    vertex.size=2.5, vertex.label=NA,
                    vertex.label.color = "black", 
                    vertex.frame.color = NA, vertex.frame.width = 0.5,
                    #vertex.color =  colors.time[src.endoderm.al12.ext_v1.1[["Time"]][graph.al12_cluster$names,]],
                    #vertex.color =  cluster.endoderm.color.v5[src.endoderm.al12.ext_v1.1[["cluster.v06.26.re"]][graph.al12_cluster$names,]],
                    vertex.color =  colors_bg[src.endoderm.al12.ext_v1.1[[color]][graph.al12_cluster$names,]])
      }}}
  
}
dev.off()

pdf("figure.v08.07/organ_development_re/al12_summary_graph.re_other.pdf",12,7)
# Reduction
for(red in c("umap_integrated","umap_mnn",
             "fdl_pca",'fdl_mnn',"fdl_pca_integrated")){
  p1 = DimPlot(src.endoderm.al12.ext_v1.1,group.by = 'cluster.v06.26.re', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p2 = DimPlot(src.endoderm.al12.ext_v1.1,group.by = 'cluster.extract.v1.1', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p3 = DimPlot(src.endoderm.al12.ext_v1.1,group.by = 'Time', 
               reduction = red , cols = colors.time,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p4 = DimPlot(src.endoderm.al12.ext_v1.1,group.by = 'cluster.v06.26.re_hc', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  print(p1); print(p2); print(p3); print(p4)
}
dev.off()


pdf("figure.v08.07/organ_development_re/al12_summary_proportion.re.pdf",12,7)
#-----------------------------
cell_type = c("AL.1","AL.2","AL.1/2-Liver","Liver")
meta_al12 = src.endoderm.al12.ext_v1.1@meta.data
meta_al12$Time = factor(meta_al12$Time, levels = names(colors.time))
meta_al12$cluster.extract.v1.1 = factor(meta_al12$cluster.extract.v1.1, levels = c("AL.1","AL.2","AL.1/2","FG.4-AL.1/2/3"))
meta_al12$cluster.v06.26.re = factor(meta_al12$cluster.v06.26.re, levels = rev(cell_type))

meta_al12_B0 = meta_al12
meta_al12_B1 = meta_al12[meta_al12$batch%in%1&!meta_al12$Time%in%"ss9",]
meta_al12_B2 = meta_al12[meta_al12$batch%in%2,]

for(data in c("meta_al12_B0","meta_al12_B1","meta_al12_B2")){
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

save(src.endoderm.al12.ext_v1.1,
     file = "figure.v08.07/organ_development_re/src.endoderm.al12.ext_v1.1.re.Rdata")

# Correct 
#=================
pdf("figure.v08.07/organ_development_re/al12_heatmap_marker.re.pdf",6,9)
cellorder.al12 = c(cellorder.al12_al1_UPI, cellorder.al12_al2_UPI,
                   cellorder.al12_al12lv_UPI, cellorder.al12_lv_UPI)
selectgene = unique(c(markergene.al12.ext_v1.1,
                      src.endoderm.al12.ext_v1.1.gene.fin))
#-------------------
al12.ext_v1.1.re.coltree =
  MyHeatmap(as.matrix(src.endoderm.al12.ext_v1.1@assays$RNA@data[selectgene, cellorder.al12]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.al12.ext_v1.1$Time[cellorder.al12], colors.time),
              MyName2Col(src.endoderm.al12.ext_v1.1$cluster.v06.26.re[cellorder.al12], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.al12.ext_v1.1$cluster.extract.v1.1[cellorder.al12], cluster.endoderm.color.v5)
            ),
            # RowSideColors = t(cbind(
            #   MyName2Col(names(selectgene), colors.geneset)
            # )),
            ColSideColorsSize = 2.5,
            # Rowv = "none",
            # Colv = "none",
            return.tree = "col",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.al12.ext_v1.1@meta.data[cell_al12_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

al12.ext_v1.1.re.coltree = as.dendrogram(al12.ext_v1.1.re.coltree)
al12.ext_v1.1.re.colgene = selectgene
src.endoderm.al12.ext_v1.1$cluster.v06.26.re_hc = 
  src.endoderm.al12.ext_v1.1$cluster.v06.26.re
src.endoderm.al12.ext_v1.1@meta.data[labels(al12.ext_v1.1.re.coltree[[2]]),]$cluster.v06.26.re_hc = "Liver"
#=================


# Marker
#------------------------
src.endoderm.al12.ext_v1.1 = 
  SetIdent(src.endoderm.al12.ext_v1.1,
           value = src.endoderm.al12.ext_v1.1$cluster.v06.26.re_hc)
markergene.al12.ext_v1.1 = 
  FindAllMarkers(src.endoderm.al12.ext_v1.1, assay = "RNA")
markergene.al12.ext_v1.1 = markergene.al12.ext_v1.1[
  markergene.al12.ext_v1.1$avg_log2FC>0.2&
    markergene.al12.ext_v1.1$p_val_adj<0.1, "gene"]

#------------------------
al12.umap.embedding = src.endoderm.al12.ext_v1.1[["fdl_pca_integrated"]]@cell.embeddings[,1:2]
al12.umap.embedding = src.endoderm.al12.ext_v1.1[["fdl_mnn"]]@cell.embeddings[,1:2]
al12.umap.embedding = src.endoderm.al12.ext_v1.1[["umap_mnn"]]@cell.embeddings[,1:3]
al12.umap.embedding = src.endoderm.al12.ext_v1.1[["umap_integrated"]]@cell.embeddings[,1:3]

# Cell order :: Princurve
#---------------------------
# AL1
cellorder.al12_al1 = 
  rownames(src.endoderm.al12.ext_v1.1@meta.data[
    src.endoderm.al12.ext_v1.1$cluster.v06.26.re_hc%in%"AL.1",])
princurve.al12_al1 = princurve::principal_curve(
  x = al12.umap.embedding[cellorder.al12_al1,],  smoother = "smooth.spline")
cellorder.al12_al1 =  names(
  princurve.al12_al1$lambda[cellorder.al12_al1][order(princurve.al12_al1$lambda[cellorder.al12_al1])])

# AL2
cellorder.al12_al2 = 
  rownames(src.endoderm.al12.ext_v1.1@meta.data[
    src.endoderm.al12.ext_v1.1$cluster.v06.26.re_hc%in%"AL.2",])
princurve.al12_al2 = princurve::principal_curve(
  x = al12.umap.embedding[cellorder.al12_al2,],  smoother = "smooth.spline")
cellorder.al12_al2 =  names(
  princurve.al12_al2$lambda[cellorder.al12_al2][order(princurve.al12_al2$lambda[cellorder.al12_al2])])

# AL1/2
cellorder.al12_al12lv = 
  rownames(src.endoderm.al12.ext_v1.1@meta.data[
    src.endoderm.al12.ext_v1.1$cluster.v06.26.re_hc%in%"AL.1/2-Liver",])
princurve.al12_al12lv = princurve::principal_curve(
  x = al12.umap.embedding[cellorder.al12_al12lv,],  smoother = "smooth.spline")
cellorder.al12_al12lv =  names(
  princurve.al12_al12lv$lambda[cellorder.al12_al12lv][order(princurve.al12_al12lv$lambda[cellorder.al12_al12lv])])

# Lv
cellorder.al12_lv = 
  rownames(src.endoderm.al12.ext_v1.1@meta.data[
    src.endoderm.al12.ext_v1.1$cluster.v06.26.re_hc%in%"Liver",])
princurve.al12_lv = princurve::principal_curve(
  x = al12.umap.embedding[cellorder.al12_lv,],  smoother = "smooth.spline")
cellorder.al12_lv =  names(
  princurve.al12_lv$lambda[cellorder.al12_lv][order(princurve.al12_lv$lambda[cellorder.al12_lv])])
#---------------------------


pdf("figure.v08.07/try.pdf")
plot(c(1:length(cellorder.al12_al1)),
     c(1:length(cellorder.al12_al1)),
     col = colors.time[src.endoderm.al12.ext_v1.1[,cellorder.al12_al1]$Time])
plot(c(1:length(cellorder.al12_al2)),
     c(1:length(cellorder.al12_al2)),
     col = colors.time[src.endoderm.al12.ext_v1.1[,cellorder.al12_al2]$Time])
plot(c(1:length(cellorder.al12_al12lv)),
     c(1:length(cellorder.al12_al12lv)),
     col = colors.time[src.endoderm.al12.ext_v1.1[,cellorder.al12_al12lv]$Time])
plot(c(1:length(cellorder.al12_lv)),
     c(1:length(cellorder.al12_lv)),
     col = colors.time[src.endoderm.al12.ext_v1.1[,cellorder.al12_lv]$Time])
dev.off()


cellorder.al12_al1_UPI = cellorder.al12_al1
cellorder.al12_al2_UPI = cellorder.al12_al2
cellorder.al12_al12lv_UPI = cellorder.al12_al12lv
cellorder.al12_lv_UPI = cellorder.al12_lv


pdf("figure.v08.07/organ_development_re/al12_heatmap_marker.re.pdf",6,9)
cellorder.al12 = c(cellorder.al12_al1_UPI, cellorder.al12_al2_UPI,
                   cellorder.al12_al12lv_UPI, cellorder.al12_lv_UPI)
selectgene = src.endoderm.al12.ext_v1.1.filtergene_fin
#-------------------
al12.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.al12.ext_v1.1@assays$RNA@data[selectgene, cellorder.al12]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.al12.ext_v1.1$Time[cellorder.al12], colors.time),
              MyName2Col(src.endoderm.al12.ext_v1.1$cluster.v06.26.re[cellorder.al12], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.al12.ext_v1.1$cluster.extract.v1.1[cellorder.al12], cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            # return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.al12.ext_v1.1@meta.data[cell_al12_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()


al12.ext_v1.1.re.rowtree.1 = as.dendrogram(al12.ext_v1.1.re.rowtree)
src.endoderm.al12.ext_v1.1.filtergene_fin = c(
  labels(al12.ext_v1.1.re.rowtree.1[[2]][[2]][[2]][[2]]),
  labels(al12.ext_v1.1.re.rowtree.1[[2]][[2]][[1]]),
  labels(al12.ext_v1.1.re.rowtree.1[[2]][[2]][[2]][[1]]),

  labels(al12.ext_v1.1.re.rowtree.1[[2]][[1]]),
  
  labels(al12.ext_v1.1.re.rowtree.1[[1]][[2]][[2]][[2]]),
  labels(al12.ext_v1.1.re.rowtree.1[[1]][[1]]),
  rev(labels(al12.ext_v1.1.re.rowtree.1[[1]][[2]][[1]])),
  rev(labels(al12.ext_v1.1.re.rowtree.1[[1]][[2]][[2]][[1]]))
  )

names(src.endoderm.al12.ext_v1.1.filtergene_fin) = c(
  rep(4,length(c(
    labels(al12.ext_v1.1.re.rowtree.1[[2]][[2]][[2]][[2]]),
    labels(al12.ext_v1.1.re.rowtree.1[[2]][[2]][[1]]),
    labels(al12.ext_v1.1.re.rowtree.1[[2]][[2]][[2]][[1]])
  ))),
  rep(3,length(c(
    labels(al12.ext_v1.1.re.rowtree.1[[2]][[1]])
  ))),
  rep(7,length(c(
    labels(al12.ext_v1.1.re.rowtree.1[[1]][[2]][[2]][[2]]),
    labels(al12.ext_v1.1.re.rowtree.1[[1]][[1]]),
    rev(labels(al12.ext_v1.1.re.rowtree.1[[1]][[2]][[1]])),
    labels(al12.ext_v1.1.re.rowtree.1[[1]][[2]][[2]][[1]])
  )))
)



save(
  # Gene 
  src.endoderm.al12.ext.v1.1.selectgene,
  src.endoderm.al12.ext_v1.1.gene.fin,
  src.endoderm.al12.ext_v1.1.filtergene_fin,
  # Cell type
  cellorder.al12_al1_UPI, cellorder.al12_al2_UPI,
  cellorder.al12_al12lv_UPI, cellorder.al12_lv_UPI,
  file= "figure.v08.07/organ_development_re/al12_heatmap_parameter.Rdata")





#---------------------------------------------------------------------------
# --          AL.1/2  seprated
#---------------------------------------------------------------------------
src.al1.integrated.re = src.al1.integrated
src.al1.integrated.re$cluster.v06.26.re_correct_re = 
  src.endoderm.al12.ext_v1.1.re$cluster.v06.26.re_correct_re[colnames(src.al1.integrated.re)]
src.al1.integrated.re = src.al1.integrated.re[,!src.al1.integrated.re$cluster.v06.26.re_correct_re%in%"AL.2"]

src.al2.integrated.re = src.al2.integrated
src.al2.integrated.re$cluster.v06.26.re_correct_re = 
  src.endoderm.al12.ext_v1.1.re$cluster.v06.26.re_correct_re[colnames(src.al2.integrated.re)]
src.al2.integrated.re = src.al2.integrated.re[,!src.al2.integrated.re$cluster.v06.26.re_correct_re%in%"AL.1"]


for(i.seurat in c("src.al1.integrated.re",
                  "src.al2.integrated.re")){
  gene.used = src.endoderm.al12.ext_v1.1.gene.fin
  
  seurat.temp = get(i.seurat)
  seurat.temp = ScaleData(seurat.temp, features = rownames(seurat.temp))
  seurat.temp = RunPCA(seurat.temp,
                       features = seurat.temp.gene.fin)
  seurat.temp = 
    seurat_mnn(seurat_name = "seurat.temp", dims = 1:30, 
               n.neighbors = 100,n.components = 3,
               seurat.selectgene = seurat.temp.gene.fin)
  seurat.temp = 
    seurat_int(seurat_name = "seurat.temp", dims = 1:30,
               n.neighbors = 100,n.components = 3,
               seurat.selectgene = seurat.temp.gene.fin)
  
  assign(i.seurat, seurat.temp)
}

DimPlot(src.al1.integrated.re, reduction = "umap_integrated",
        group.by = "cluster.v06.26.re_correct_re")
DimPlot(src.al2.integrated.re, reduction = "umap_integrated",
        group.by = "cluster.v06.26.re_correct_re")

save(src.al1.integrated.re, file = "figure.v08.07/organ_development_re_v240115//src.al1.integrated.re.Rdata")
save(src.al2.integrated.re, file = "figure.v08.07/organ_development_re_v240115//src.al2.integrated.re.Rdata")


