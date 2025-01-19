#---------------------------------------------------------------------------
# --          AL.1 -- Re
#---------------------------------------------------------------------------
DimPlot(src.9ss.integrated.merge, group.by = "lineage", reduction = "umap_fta")
DimPlot(src.9ss.integrated.merge, group.by = "cluster.v06.26.re..merge_umap_fta", reduction = "umap_fta")


src.al1.tracing@meta.data[
  intersect(rownames(src.9ss.integrated.merge@meta.data[src.9ss.integrated.merge$cluster.v06.26.re..merge_umap_fta%in%"AL.1",]),
            rownames(src.al1.tracing@meta.data)),]$cluster.v06.26.re_mnn_umap_fta = "AL.1"
src.al2.tracing@meta.data[
  intersect(rownames(src.9ss.integrated.merge@meta.data[src.9ss.integrated.merge$cluster.v06.26.re..merge_umap_fta%in%"AL.2",]),
            rownames(src.al2.tracing@meta.data)),]$cluster.v06.26.re_mnn_umap_fta = "AL.2"

#------------------------
src.al1.tracing = FindVariableFeatures(src.al1.tracing, nfeatures = 2000)
src.al1.tracing.filtergene = 
  Myfilter(as.matrix(src.al1.tracing@assays$RNA@data),
           gene = src.al1.tracing@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.al1.tracing.filtergene.re = src.al1.tracing.filtergene; rm(src.al1.tracing.filtergene)
src.al1.tracing.filtergene = src.al1.tracing.filtergene.re; rm(src.al1.tracing.filtergene.re)


src.al1.tracing = SetIdent(src.al1.tracing, value = src.al1.tracing$cluster.v06.26.re_mnn_umap_fta)
marker_src.al1.tracing = FindAllMarkers(src.al1.tracing)
marker_src.al1.tracing$pct.ratio = marker_src.al1.tracing$pct.1 / marker_src.al1.tracing$pct.2
marker_src.al1.tracing$rank = marker_src.al1.tracing$pct.ratio * (-log(marker_src.al1.tracing$p_val_adj))
marker_src.al1.tracing = marker_src.al1.tracing[order(marker_src.al1.tracing$rank, decreasing = T),]
markergene_src.al1.tracing = unique(marker_src.al1.tracing$gene)
# markergene_src.al1.tracing.raw = markergene_src.al1.tracing

src.al1.tracing = ScaleData(src.al1.tracing, 
                            rownames(src.al1.tracing),split.by = "Phase")
src.al1.tracing = RunPCA(src.al1.tracing, 
                         features = src.al1.tracing.filtergene)
src.al1.tracing = RunUMAP(src.al1.tracing, reduction = "pca", 
                          dims = 1:30, n.neighbors = 150,
                          assay = "RNA", n.components = 3)

#====================================
# == AL.1 Basic process 
#===============================================================================
pdf("figure.v08.07/organ_development_re_v240115/try.al1.pdf",9,7)
src.al1.tracing.rowtree =
  MyHeatmap(as.matrix(src.al1.tracing@assays$RNA@data[
    unique(c(markergene_src.al1.tracing,
             src.al1.tracing.filtergene,
             c()
    )),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.al1.tracing$Time,
                 colors.time.2),
      MyName2Col(src.al1.tracing$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.al1.tracing$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.al1.tracing$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()
pdf("figure.v08.07/organ_development_re_v240115/try.al1.tree.pdf",150,20)
plot(src.al1.tracing.rowtree)
dev.off()

tree_src.al1.tracing.rowtree = as.dendrogram(src.al1.tracing.rowtree)
gene_src.al1.tracing.rowtree = setdiff(
  unique(c(markergene_src.al1.tracing, src.al1.tracing.filtergene)),
  c(labels(tree_src.al1.tracing.rowtree[[2]][[2]][[2]][[2]]),
    labels(tree_src.al1.tracing.rowtree[[2]][[2]][[1]])))


pdf("figure.v08.07/organ_development_re_v240115/try.al1.pdf",9,7)
src.al1.tracing.rowtree.1 =
  MyHeatmap(as.matrix(src.al1.tracing@assays$RNA@data[
    gene_src.al1.tracing.rowtree,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.al1.tracing$Time,
                 colors.time.2),
      MyName2Col(src.al1.tracing$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.al1.tracing$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.al1.tracing$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
src.al1.tracing.coltree.1 =
  MyHeatmap(as.matrix(src.al1.tracing@assays$RNA@data[
    gene_src.al1.tracing.rowtree,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.al1.tracing$Time,
                 colors.time.2),
      #MyName2Col(src.al1.tracing$tree.1,
      #           color.temp),
      MyName2Col(src.al1.tracing$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.al1.tracing$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.al1.tracing$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "col",
    graph = T)
dev.off()


tree_src.al1.tracing.rowtree.1 = as.dendrogram(src.al1.tracing.rowtree.1)
gene_src.al1.tracing.rowtree.1 = c(
  labels(tree_src.al1.tracing.rowtree.1[[1]][[1]]),
  setdiff(labels(tree_src.al1.tracing.rowtree.1[[1]][[2]]),
          c(labels(tree_src.al1.tracing.rowtree.1[[1]][[2]][[1]][[1]]),
            labels(tree_src.al1.tracing.rowtree.1[[1]][[2]][[1]][[2]][[1]]),
            labels(tree_src.al1.tracing.rowtree.1[[1]][[2]][[1]][[2]][[2]][[1]]))),
  labels(tree_src.al1.tracing.rowtree.1[[2]][[1]][[2]]),
  labels(tree_src.al1.tracing.rowtree.1[[2]][[2]][[1]][[1]]),
  rev(labels(tree_src.al1.tracing.rowtree.1[[2]][[2]][[2]])),
  labels(tree_src.al1.tracing.rowtree.1[[2]][[2]][[1]][[2]]))
names(gene_src.al1.tracing.rowtree.1) = c(
  rep(4, length(labels(tree_src.al1.tracing.rowtree.1[[1]][[1]]))),
  rep(3, length(setdiff(labels(tree_src.al1.tracing.rowtree.1[[1]][[2]]),
                        c(labels(tree_src.al1.tracing.rowtree.1[[1]][[2]][[1]][[1]]),
                          labels(tree_src.al1.tracing.rowtree.1[[1]][[2]][[1]][[2]][[1]]),
                          labels(tree_src.al1.tracing.rowtree.1[[1]][[2]][[1]][[2]][[2]][[1]]))))),
  rep(7,length(c(labels(tree_src.al1.tracing.rowtree.1[[2]][[1]][[2]]),
                 labels(tree_src.al1.tracing.rowtree.1[[2]][[2]][[1]][[1]]),
                 rev(labels(tree_src.al1.tracing.rowtree.1[[2]][[2]][[2]])),
                 labels(tree_src.al1.tracing.rowtree.1[[2]][[2]][[1]][[2]])))))
  
tree_src.al1.tracing.coltree.1 = as.dendrogram(src.al1.tracing.coltree.1)
src.al1.tracing$tree.1 = NA
src.al1.tracing@meta.data[labels(tree_src.al1.tracing.coltree.1[[1]][[1]][[1]]),]$tree.1 = 1
src.al1.tracing@meta.data[labels(tree_src.al1.tracing.coltree.1[[1]][[1]][[2]]),]$tree.1 = 2
src.al1.tracing@meta.data[labels(tree_src.al1.tracing.coltree.1[[1]][[2]][[1]]),]$tree.1 = 3
src.al1.tracing@meta.data[labels(tree_src.al1.tracing.coltree.1[[1]][[2]][[2]]),]$tree.1 = 4
src.al1.tracing@meta.data[labels(tree_src.al1.tracing.coltree.1[[2]][[1]][[1]]),]$tree.1 = 5
src.al1.tracing@meta.data[labels(tree_src.al1.tracing.coltree.1[[2]][[1]][[2]]),]$tree.1 = 6
src.al1.tracing@meta.data[labels(tree_src.al1.tracing.coltree.1[[2]][[2]][[1]]),]$tree.1 = 7
src.al1.tracing@meta.data[labels(tree_src.al1.tracing.coltree.1[[2]][[2]][[2]]),]$tree.1 = 8

src.al1.tracing$cluster.v06.26.re_correct = src.al1.tracing$cluster.v06.26.re_mnn_umap_fta
DimPlot(src.al1.tracing, reduction = "mnn_umap_fta", group.by = "tree.1", label = T) +
  DimPlot(src.al1.tracing, reduction = "mnn_umap_fta", group.by = "Time", label = T)
src.al1.tracing@meta.data[src.al1.tracing$tree.1%in%c(5,6),]$cluster.v06.26.re_correct = "AL.1"
src.al1.tracing@meta.data[src.al1.tracing$tree.1%in%c(7,8),]$cluster.v06.26.re_correct = "AL.1/2-Liver"
src.al1.tracing@meta.data[src.al1.tracing$tree.1%in%c(1,2,3,4),]$cluster.v06.26.re_correct = "Liver"
DimPlot(src.al1.tracing, reduction = "mnn_umap_fta", group.by = "cluster.v06.26.re_correct")

src.al1.tracing$cluster.v06.26.re_correct_re = NA
src.al1.tracing@meta.data[colnames(src.al1.tracing_re),]$cluster.v06.26.re_correct_re = 
  src.al1.tracing_re$cluster.v06.26.re_correct
DimPlot(src.al1.tracing, reduction = "mnn_umap_fta", group.by = "cluster.v06.26.re_correct_re")

#-- Pseudo time
cell_src.al1.tracing = rownames(src.al1.tracing@meta.data)
coord_src.al1.tracing = src.al1.tracing[["umap"]]@cell.embeddings[cell_src.al1.tracing, c(1,2)]
pcurve_src.al1.tracing = princurve::principal_curve(x = coord_src.al1.tracing, smoother = "smooth.spline")
src.al1.tracing$lambda = pcurve_src.al1.tracing$lambda[cell_src.al1.tracing]
src.al1.tracing@meta.data[cell_src.al1.tracing,]$lambda = 
  pcurve_src.al1.tracing$lambda[cell_src.al1.tracing]

cellorder_cell_src.al1.tracing = c(
  (intersect(names(src.al1.tracing$lambda[order(src.al1.tracing$lambda)]), 
             colnames(src.al1.tracing[,src.al1.tracing$cluster.v06.26.re_correct%in%c("AL.1"),]))),
  (intersect(names(src.al1.tracing$lambda[order(src.al1.tracing$lambda)]), 
             colnames(src.al1.tracing[,src.al1.tracing$cluster.v06.26.re_correct%in%c("AL.1/2-Liver"),]))),
  (intersect(names(src.al1.tracing$lambda[order(src.al1.tracing$lambda)]), 
             colnames(src.al1.tracing[,src.al1.tracing$cluster.v06.26.re_correct%in%c("Liver"),]))))

cellorder_cell_src.al1.tracing = intersect(
  cellorder_cell_src.al1.tracing,
  rownames(src.al1.tracing@meta.data[src.al1.tracing$lineage%in%c("Nkx2_3"),]))


#---------- Only Nkx2-3  --------
#===============================
src.al1.tracing_re = src.al1.tracing[,src.al1.tracing$lineage%in%c("Nkx2_3")]
src.al1.tracing_re = FindVariableFeatures(src.al1.tracing_re, nfeatures = 2000)
src.al1.tracing_re.filtergene = 
  Myfilter(as.matrix(src.al1.tracing_re@assays$RNA@data),
           gene = src.al1.tracing_re@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.al1.tracing_re.filtergene.re = src.al1.tracing_re.filtergene; rm(src.al1.tracing_re.filtergene)
src.al1.tracing_re.filtergene = src.al1.tracing_re.filtergene.re; rm(src.al1.tracing_re.filtergene.re)


src.al1.tracing_re = SetIdent(src.al1.tracing_re, value = src.al1.tracing_re$cluster.v06.26.re_mnn_umap_fta)
# src.al1.tracing_re = SetIdent(src.al1.tracing_re, value = src.al1.tracing_re$cluster.v06.26.re_correct_re)
marker_src.al1.tracing_re = FindAllMarkers(src.al1.tracing_re)
marker_src.al1.tracing_re$pct.ratio = marker_src.al1.tracing_re$pct.1 / marker_src.al1.tracing_re$pct.2
marker_src.al1.tracing_re$rank = marker_src.al1.tracing_re$pct.ratio * (-log(marker_src.al1.tracing_re$p_val_adj))
marker_src.al1.tracing_re = marker_src.al1.tracing_re[order(marker_src.al1.tracing_re$rank, decreasing = T),]
markergene_src.al1.tracing_re = unique(marker_src.al1.tracing_re$gene)
# markergene_src.al1.tracing_re.raw = markergene_src.al1.tracing_re

src.al1.tracing_re = ScaleData(src.al1.tracing_re, 
                               rownames(src.al1.tracing_re),split.by = "Phase")
src.al1.tracing_re = RunPCA(src.al1.tracing_re, 
                            features = src.al1.tracing_re.filtergene)
src.al1.tracing_re = RunUMAP(src.al1.tracing_re, reduction = "pca", 
                             dims = 1:30, n.neighbors = 20,
                             assay = "RNA", n.components = 3)
DimPlot(src.al1.tracing_re, group.by = "cluster.v06.26.re_correct")


pdf("figure.v08.07/organ_development_re_v240115/try.al1.re.pdf",9,7)
src.al1.tracing_re.rowtree.1 =
  MyHeatmap(as.matrix(src.al1.tracing_re@assays$RNA@data[
    unique(c(markergene_src.al1.tracing_re,
             src.al1.tracing_re.filtergene,
             c()
    )),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.al1.tracing_re$Time,
                 colors.time.2),
      MyName2Col(src.al1.tracing_re$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.al1.tracing_re$cluster.v06.26.re_correct_re,
                 cluster.endoderm.color.v5),
      MyName2Col(src.al1.tracing_re$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
src.al1.tracing_re.coltree.1 =
  MyHeatmap(as.matrix(src.al1.tracing_re@assays$RNA@data[
    unique(c(markergene_src.al1.tracing_re,
             src.al1.tracing_re.filtergene,
             c()
    )),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.al1.tracing_re$Time,
                 colors.time.2),
      MyName2Col(src.al1.tracing_re$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.al1.tracing_re$cluster.v06.26.re_correct_re,
                 cluster.endoderm.color.v5),
      MyName2Col(src.al1.tracing_re$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "col",
    graph = T)
dev.off()
pdf("figure.v08.07/organ_development_re_v240115/try.al1.re.tree.pdf",150,20)
plot(src.al1.tracing_re.rowtree.1)
dev.off()

tree_src.al1.tracing_re.rowtree.1 = as.dendrogram(src.al1.tracing_re.rowtree.1)
gene_src.al1.tracing_re.rowtree.1 = c(
  # labels(tree_src.al1.tracing_re.rowtree.1[[1]][[1]]),
  labels(tree_src.al1.tracing_re.rowtree.1[[1]][[2]][[1]]),
  labels(tree_src.al1.tracing_re.rowtree.1[[1]][[2]][[2]]),
  #labels(tree_src.al1.tracing_re.rowtree.1[[2]][[2]]),
  rev(c(
    labels(tree_src.al1.tracing_re.rowtree.1[[2]][[1]][[2]]),
    labels(tree_src.al1.tracing_re.rowtree.1[[2]][[1]][[1]])
  )))
names(gene_src.al1.tracing_re.rowtree.1) = c(
  rep(4, length(c(# labels(tree_src.al1.tracing_re.rowtree.1[[1]][[1]]),
                  labels(tree_src.al1.tracing_re.rowtree.1[[1]][[2]][[1]])))),
  rep(3, length(c(labels(tree_src.al1.tracing_re.rowtree.1[[1]][[2]][[2]])))),
  rep(7, length(c(# labels(tree_src.al1.tracing_re.rowtree.1[[2]][[2]]),
                  labels(tree_src.al1.tracing_re.rowtree.1[[2]][[1]])))))


tree_src.al1.tracing_re.coltree.1 = as.dendrogram(src.al1.tracing_re.coltree.1)
src.al1.tracing_re$tree.2 = NA
src.al1.tracing_re@meta.data[labels(tree_src.al1.tracing_re.coltree.1[[1]][[1]][[1]]),]$tree.2 = 1
src.al1.tracing_re@meta.data[labels(tree_src.al1.tracing_re.coltree.1[[1]][[1]][[2]]),]$tree.2 = 2
src.al1.tracing_re@meta.data[labels(tree_src.al1.tracing_re.coltree.1[[1]][[2]][[1]]),]$tree.2 = 3
src.al1.tracing_re@meta.data[labels(tree_src.al1.tracing_re.coltree.1[[1]][[2]][[2]]),]$tree.2 = 4
src.al1.tracing_re@meta.data[labels(tree_src.al1.tracing_re.coltree.1[[2]][[1]][[1]]),]$tree.2 = 5
src.al1.tracing_re@meta.data[labels(tree_src.al1.tracing_re.coltree.1[[2]][[1]][[2]]),]$tree.2 = 6
src.al1.tracing_re@meta.data[labels(tree_src.al1.tracing_re.coltree.1[[2]][[2]][[1]]),]$tree.2 = 7
src.al1.tracing_re@meta.data[labels(tree_src.al1.tracing_re.coltree.1[[2]][[2]][[2]]),]$tree.2 = 8
DimPlot(src.al1.tracing_re, reduction = "mnn_umap_fta", group.by = "tree.2", label = T, label.size = 5)

src.al1.tracing_re$cluster.v06.26.re_correct_re = src.al1.tracing_re$cluster.v06.26.re_mnn_umap_fta
src.al1.tracing_re@meta.data[src.al1.tracing_re$tree.2%in%c(7,8),]$cluster.v06.26.re_correct_re = "AL.1"
src.al1.tracing_re@meta.data[src.al1.tracing_re$tree.2%in%c(5,6),]$cluster.v06.26.re_correct_re = "AL.1/2-Liver"
src.al1.tracing_re@meta.data[src.al1.tracing_re$tree.2%in%c(1,2,3,4),]$cluster.v06.26.re_correct_re = "Liver"
DimPlot(src.al1.tracing_re, reduction = "mnn_umap_fta", group.by = "cluster.v06.26.re_correct_re", label = T, label.size = 5)

#-- pseudo time
cell_src.al1.tracing_re = rownames(src.al1.tracing_re@meta.data)
coord_src.al1.tracing_re = src.al1.tracing_re[["mnn_umap_fta"]]@cell.embeddings[cell_src.al1.tracing_re, c(1,2)]
pcurve_src.al1.tracing_re = princurve::principal_curve(x = coord_src.al1.tracing_re, smoother = "smooth.spline")
src.al1.tracing_re$lambda = pcurve_src.al1.tracing_re$lambda[cell_src.al1.tracing_re]
src.al1.tracing_re@meta.data[cell_src.al1.tracing_re,]$lambda = 
  pcurve_src.al1.tracing_re$lambda[cell_src.al1.tracing_re]
src.al1.tracing_re@meta.data[!src.al1.tracing_re$lambda%in%NA,]$lambda = 
  norm_range(src.al1.tracing_re@meta.data[!src.al1.tracing_re$lambda%in%NA,]$lambda)

cellorder_cell_src.al1.tracing_re = c(
  rev(intersect(names(src.al1.tracing_re$lambda[order(src.al1.tracing_re$lambda)]), 
                colnames(src.al1.tracing_re[,src.al1.tracing_re$cluster.v06.26.re_correct_re%in%c("AL.1"),]))),
  rev(intersect(names(src.al1.tracing_re$lambda[order(src.al1.tracing_re$lambda)]), 
                colnames(src.al1.tracing_re[,src.al1.tracing_re$cluster.v06.26.re_correct_re%in%c("AL.1/2-Liver"),]))),
  rev(intersect(names(src.al1.tracing_re$lambda[order(src.al1.tracing_re$lambda)]), 
                colnames(src.al1.tracing_re[,src.al1.tracing_re$cluster.v06.26.re_correct_re%in%c("Liver"),]))))

pdf("figure.v08.07/organ_development_re_v240115/Heatmap_for_trajectory/try.al1.re.pdf", 10, 10)
src.al1.tracing_re.coltree.1 =
  MyHeatmap(as.matrix(src.al1.tracing_re@assays$RNA@data[
    gene_src.al1.tracing_re.rowtree.1,
    cellorder_cell_src.al1.tracing_re]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.al1.tracing_re$Time[cellorder_cell_src.al1.tracing_re],
                 colors.time.2),
      MyName2Col(src.al1.tracing_re$lineage[cellorder_cell_src.al1.tracing_re],
                 color.lineage),
      MyName2Col(src.al1.tracing_re$cluster.v06.26.re_correct_re[cellorder_cell_src.al1.tracing_re],
                 cluster.endoderm.color.v5)),
    RowSideColors = t(cbind(
      MyName2Col(names(gene_src.al1.tracing_re.rowtree.1),
                 colors.geneset))),
    ColSideColorsSize = 4.8,
    RowSideColorsSize = 1.2,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    Rowv = "none", Colv = "none",
    margins = c(10,10),
    #return.tree = "col",
    graph = T)
dev.off()

save(src.al1.tracing, file = 'figure.v08.07/organ_development_re_v240115/src.al1.tracing.Rdata')
save(src.al1.tracing_re,
     gene_src.al1.tracing_re.rowtree.1, file = 'figure.v08.07/organ_development_re_v240115/src.al1.tracing_re.Rdata')
save(src.al1.tracing.filtergene,
     marker_src.al1.tracing,
     markergene_src.al1.tracing,
     tree_src.al1.tracing.rowtree,
     tree_src.al1.tracing.rowtree.1,
     tree_src.al1.tracing.coltree.1,
     gene_src.al1.tracing.rowtree,
     gene_src.al1.tracing.rowtree.1,
     
     src.al1.tracing_re.filtergene,
     markergene_src.al1.tracing_re,
     tree_src.al1.tracing_re.coltree.1,
     tree_src.al1.tracing_re.rowtree.1,
     gene_src.al1.tracing_re.rowtree.1, 
     
     cellorder_cell_src.al1.tracing,
     cellorder_cell_src.al1.tracing_re,
     file = 'figure.v08.07/organ_development_re_v240115/src.al1.tracing.parameter.Rdata')

#-- Correct to src.al1.tracing ------------------

pdf("figure.v08.07/organ_development_re_v240115/Heatmap_for_trajectory/try.al1.pdf",10,10)
src.al1.tracing.coltree.1 =
  MyHeatmap(as.matrix(src.al1.tracing_re@assays$RNA@data[
    gene_src.al1.tracing.rowtree.1,
    cellorder_cell_src.al1.tracing_re]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.al1.tracing_re$Time[cellorder_cell_src.al1.tracing_re],
                 colors.time.2),
      MyName2Col(src.al1.tracing_re$lineage[cellorder_cell_src.al1.tracing_re],
                 color.lineage),
      MyName2Col(src.al1.tracing_re$cluster.v06.26.re_correct_re[cellorder_cell_src.al1.tracing_re],
                 cluster.endoderm.color.v5)),
    RowSideColors = t(cbind(
      MyName2Col(names(gene_src.al1.tracing.rowtree.1),
                 colors.geneset))),
    ColSideColorsSize = 4.8,
    RowSideColorsSize = 1.2,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    Rowv = "none", Colv = "none",
    margins = c(10,10), 
    #return.tree = "col",
    graph = T)
dev.off()


#===============================================================================


#---------------------------------------------------------------------------
# --          AL.2 -- Re
#---------------------------------------------------------------------------
DimPlot(src.9ss.integrated.merge, group.by = "lineage", reduction = "umap_fta")
DimPlot(src.9ss.integrated.merge, group.by = "cluster.predict.umap_int.ext.v1.1", reduction = "umap_fta")

src.al2.tracing@meta.data[
  intersect(rownames(src.9ss.integrated.merge@meta.data[src.9ss.integrated.merge$cluster.predict.umap_int.ext.v1.1%in%"AL.2",]),
            rownames(src.al2.tracing@meta.data)),]$cluster.v06.26.re_mnn_umap_fta = "AL.2"
#------------------------
src.al2.tracing = FindVariableFeatures(src.al2.tracing, nfeatures = 2000)
src.al2.tracing.filtergene = 
  Myfilter(as.matrix(src.al2.tracing@assays$RNA@data),
           gene = src.al2.tracing@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.al2.tracing.filtergene.re = src.al2.tracing.filtergene; rm(src.al2.tracing.filtergene)
src.al2.tracing.filtergene = src.al2.tracing.filtergene.re; rm(src.al2.tracing.filtergene.re)


src.al2.tracing = SetIdent(src.al2.tracing, value = src.al2.tracing$cluster.v06.26.re_mnn_umap_fta)
marker_src.al2.tracing = FindAllMarkers(src.al2.tracing)
marker_src.al2.tracing$pct.ratio = marker_src.al2.tracing$pct.1 / marker_src.al2.tracing$pct.2
marker_src.al2.tracing$rank = marker_src.al2.tracing$pct.ratio * (-log(marker_src.al2.tracing$p_val_adj))
marker_src.al2.tracing = marker_src.al2.tracing[order(marker_src.al2.tracing$rank, decreasing = T),]
markergene_src.al2.tracing = unique(marker_src.al2.tracing$gene)
# markergene_src.al2.tracing.raw = markergene_src.al2.tracing

src.al2.tracing = ScaleData(src.al2.tracing, 
                            rownames(src.al2.tracing),split.by = "Phase")
src.al2.tracing = RunPCA(src.al2.tracing, 
                         features = src.al2.tracing.filtergene)
src.al2.tracing = RunUMAP(src.al2.tracing, reduction = "pca", 
                          dims = 1:30, n.neighbors = 150,
                          assay = "RNA", n.components = 3)
# == AL.2 Basic process 
#===============================================================================
pdf("figure.v08.07/organ_development_re_v240115/try.al2.pdf",9,7)
src.al2.tracing.rowtree =
  MyHeatmap(as.matrix(src.al2.tracing@assays$RNA@data[
    unique(c(markergene_src.al2.tracing,
             src.al2.tracing.filtergene,
             c()
    )),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.al2.tracing$Time,
                 colors.time.2),
      MyName2Col(src.al2.tracing$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.al2.tracing$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.al2.tracing$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()
pdf("figure.v08.07/organ_development_re_v240115/try.al2.tree.pdf",150,20)
plot(src.al2.tracing.rowtree)
dev.off()

tree_src.al2.tracing.rowtree = as.dendrogram(src.al2.tracing.rowtree)
gene_src.al2.tracing.rowtree = setdiff(
  unique(c(markergene_src.al2.tracing, src.al2.tracing.filtergene)),
  c(labels(tree_src.al2.tracing.rowtree[[2]][[2]][[2]][[1]][[1]]),
    labels(tree_src.al2.tracing.rowtree[[2]][[2]][[1]])))


pdf("figure.v08.07/organ_development_re_v240115/try.al2.pdf",9,7)
src.al2.tracing.rowtree.1 =
  MyHeatmap(as.matrix(src.al2.tracing@assays$RNA@data[
    gene_src.al2.tracing.rowtree,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.al2.tracing$Time,
                 colors.time.2),
      MyName2Col(src.al2.tracing$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.al2.tracing$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.al2.tracing$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
src.al2.tracing.coltree.1 =
  MyHeatmap(as.matrix(src.al2.tracing@assays$RNA@data[
    gene_src.al2.tracing.rowtree,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.al2.tracing$Time,
                 colors.time.2),
      # MyName2Col(src.al2.tracing$tree.1,
      #            color.temp),
      MyName2Col(src.al2.tracing$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.al2.tracing$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.al2.tracing$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "col",
    graph = T)
dev.off()


tree_src.al2.tracing.rowtree.1 = as.dendrogram(src.al2.tracing.rowtree.1)
gene_src.al2.tracing.rowtree.1 = c(
  labels(tree_src.al2.tracing.rowtree.1[[1]][[1]][[2]][[2]]),
  labels(tree_src.al2.tracing.rowtree.1[[1]][[1]][[1]]),
  labels(tree_src.al2.tracing.rowtree.1[[1]][[1]][[2]][[1]]),
  
  labels(tree_src.al2.tracing.rowtree.1[[2]][[2]][[2]][[1]]),
  labels(tree_src.al2.tracing.rowtree.1[[2]][[2]][[2]][[2]][[1]]),
  labels(tree_src.al2.tracing.rowtree.1[[2]][[2]][[1]][[2]]),
  
  labels(tree_src.al2.tracing.rowtree.1[[2]][[2]][[1]][[1]]),
  labels(tree_src.al2.tracing.rowtree.1[[2]][[2]][[2]][[2]][[2]]))

names(gene_src.al2.tracing.rowtree.1) = c(
  rep(4, length(c(labels(tree_src.al2.tracing.rowtree.1[[1]][[1]])))),
  rep(3,length(c(labels(tree_src.al2.tracing.rowtree.1[[2]][[2]][[2]][[1]]),
                 labels(tree_src.al2.tracing.rowtree.1[[2]][[2]][[2]][[2]][[1]]),
                 labels(tree_src.al2.tracing.rowtree.1[[2]][[2]][[1]][[1]])))),
  rep(7,length(c(labels(tree_src.al2.tracing.rowtree.1[[2]][[2]][[1]][[2]]),
                 labels(tree_src.al2.tracing.rowtree.1[[2]][[2]][[2]][[2]][[2]])))))


tree_src.al2.tracing.coltree.1 = as.dendrogram(src.al2.tracing.coltree.1)
src.al2.tracing$tree.1 = NA
src.al2.tracing@meta.data[labels(tree_src.al2.tracing.coltree.1[[1]][[1]][[1]]),]$tree.1 = 1
src.al2.tracing@meta.data[labels(tree_src.al2.tracing.coltree.1[[1]][[1]][[2]]),]$tree.1 = 2
src.al2.tracing@meta.data[labels(tree_src.al2.tracing.coltree.1[[1]][[2]][[1]]),]$tree.1 = 3
src.al2.tracing@meta.data[labels(tree_src.al2.tracing.coltree.1[[1]][[2]][[2]]),]$tree.1 = 4
src.al2.tracing@meta.data[labels(tree_src.al2.tracing.coltree.1[[2]][[1]][[1]]),]$tree.1 = 5
src.al2.tracing@meta.data[labels(tree_src.al2.tracing.coltree.1[[2]][[1]][[2]]),]$tree.1 = 6
src.al2.tracing@meta.data[labels(tree_src.al2.tracing.coltree.1[[2]][[2]][[1]]),]$tree.1 = 7
src.al2.tracing@meta.data[labels(tree_src.al2.tracing.coltree.1[[2]][[2]][[2]]),]$tree.1 = 8
DimPlot(src.al2.tracing, reduction = "mnn_umap_fta", group.by = "tree.1", label = T, label.size = 5) + 
  DimPlot(src.al2.tracing, reduction = "mnn_umap_fta", group.by = "Time")

src.al2.tracing$cluster.v06.26.re_correct = src.al2.tracing$cluster.v06.26.re_mnn_umap_fta
src.al2.tracing@meta.data[src.al2.tracing$tree.1%in%c(5,6),]$cluster.v06.26.re_correct = "AL.2"
src.al2.tracing@meta.data[src.al2.tracing$tree.1%in%c(7,8),]$cluster.v06.26.re_correct = "AL.1/2-Liver"
src.al2.tracing@meta.data[src.al2.tracing$tree.1%in%c(1,2,3,4),]$cluster.v06.26.re_correct = "Liver"


#-- Pseudo time
cell_src.al2.tracing = rownames(src.al2.tracing@meta.data)
coord_src.al2.tracing = src.al2.tracing[["umap"]]@cell.embeddings[cell_src.al2.tracing, c(1,2)]
pcurve_src.al2.tracing = princurve::principal_curve(x = coord_src.al2.tracing, smoother = "smooth.spline")
src.al2.tracing$lambda = pcurve_src.al2.tracing$lambda[cell_src.al2.tracing]
src.al2.tracing@meta.data[cell_src.al2.tracing,]$lambda = 
  pcurve_src.al2.tracing$lambda[cell_src.al2.tracing]
src.al2.tracing@meta.data[!src.al2.tracing$lambda%in%NA,]$lambda = 
  norm_range(src.al2.tracing@meta.data[!src.al2.tracing$lambda%in%NA,]$lambda)

cellorder_cell_src.al2.tracing = c(
  (intersect(names(src.al2.tracing$lambda[order(src.al2.tracing$lambda)]), 
                colnames(src.al2.tracing[,src.al2.tracing$cluster.v06.26.re_correct%in%c("AL.2"),]))),
  (intersect(names(src.al2.tracing$lambda[order(src.al2.tracing$lambda)]), 
                colnames(src.al2.tracing[,src.al2.tracing$cluster.v06.26.re_correct%in%c("AL.1/2-Liver"),]))),
  (intersect(names(src.al2.tracing$lambda[order(src.al2.tracing$lambda)]), 
                colnames(src.al2.tracing[,src.al2.tracing$cluster.v06.26.re_correct%in%c("Liver"),]))))


pdf("figure.v08.07/organ_development_re_v240115/Heatmap_for_trajectory/try.al2.pdf",10,10)
src.al2.tracing.rowtree.2 =
  MyHeatmap(as.matrix(src.al2.tracing@assays$RNA@data[
    gene_src.al2.tracing.rowtree.1,
    cellorder_cell_src.al2.tracing]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.al2.tracing$Time[cellorder_cell_src.al2.tracing],
                 colors.time.2),
      MyName2Col(src.al2.tracing$lineage[cellorder_cell_src.al2.tracing],
                 color.lineage),
      MyName2Col(src.al2.tracing$cluster.v06.26.re_correct[cellorder_cell_src.al2.tracing],
                 cluster.endoderm.color.v5)),
    RowSideColors = t(cbind(
      MyName2Col(names(gene_src.al2.tracing.rowtree.1),
                 colors.geneset))),
    ColSideColorsSize = 4.8,
    RowSideColorsSize = 1.2,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    Rowv = "none", Colv = "none",
    margins = c(10,10), 
    # return.tree = "row",
    graph = T)
dev.off()

save(src.al2.tracing, file = 'figure.v08.07/organ_development_re_v240115/src.al2.tracing.Rdata')

save(src.al2.tracing.filtergene,
     marker_src.al2.tracing,
     markergene_src.al1.tracing,
     tree_src.al2.tracing.rowtree,
     tree_src.al2.tracing.rowtree.1,
     tree_src.al2.tracing.coltree.1,
     gene_src.al2.tracing.rowtree,
     gene_src.al2.tracing.rowtree.1,
     cellorder_cell_src.al2.tracing,
     file = 'figure.v08.07/organ_development_re_v240115/src.al2.tracing.parameter.Rdata')
#===============================================================================




#---------------------------------------------------------------------------
# --          AL.1/2 Check !!!
#---------------------------------------------------------------------------
src.endoderm.al12.ext_v1.1.re = src.endoderm.al12.ext_v1.1
src.endoderm.al12.ext_v1.1.re = NormalizeData(src.endoderm.al12.ext_v1.1.re, scale.factor = 10^5)
src.endoderm.al12.ext_v1.1.re = ScaleData(src.endoderm.al12.ext_v1.1.re, 
                                          features = rownames(src.endoderm.al12.ext_v1.1.re),
                                          split.by = "batch_phase")
#---------------
src.endoderm.al12.ext_v1.1.re = FindVariableFeatures(src.endoderm.al12.ext_v1.1.re, assay = "RNA", nfeatures = 2000)
src.endoderm.al12.ext_v1.1.re_filtergene = 
  Myfilter(as.matrix(src.endoderm.al12.ext_v1.1.re@assays$RNA@data),
           gene = src.endoderm.al12.ext_v1.1.re@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.endoderm.al12.ext_v1.1.re_filtergene_re = src.endoderm.al12.ext_v1.1.re_filtergene
rm(src.endoderm.al12.ext_v1.1.re_filtergene)
src.endoderm.al12.ext_v1.1.re_filtergene = src.endoderm.al12.ext_v1.1.re_filtergene_re

src.endoderm.al12.ext_v1.1.re = RunPCA(src.endoderm.al12.ext_v1.1.re, 
                                       features = src.endoderm.al12.ext_v1.1.re_filtergene)
src.endoderm.al12.ext_v1.1.re = RunUMAP(src.endoderm.al12.ext_v1.1.re, reduction = "pca",
                                        n.neighbors = 150, dims = c(1:25), n.components = 2)
#---------------

cellorder_cell_src.endoderm.al12.ext_v1.1.re = 
  c(cellorder.al12_al1_UPI, cellorder.al12_al2_UPI,
    cellorder.al12_al12lv_UPI, cellorder.al12_lv_UPI)
  
pdf("figure.v08.07/organ_development_re_v240115/try.al1_10x.pdf",9,7)
src.endoderm.al12.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.al12.ext_v1.1.re@assays$RNA@data[
    unique(c(gene_src.al2.tracing.rowtree.1,
             gene_src.al2.tracing.rowtree.1)),
    cellorder_cell_src.endoderm.al12.ext_v1.1.re]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.endoderm.al12.ext_v1.1.re$Time[cellorder_cell_src.endoderm.al12.ext_v1.1.re],
                 colors.time),
      # MyName2Col(src.endoderm.al12.ext_v1.1.re$tree[cellorder_cell_src.endoderm.al12.ext_v1.1.re],
      #            color.temp),
      MyName2Col(src.endoderm.al12.ext_v1.1.re$cluster.extract.v1.1[cellorder_cell_src.endoderm.al12.ext_v1.1.re], 
                 cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.al12.ext_v1.1.re$cluster.v06.26.re_hc[cellorder_cell_src.endoderm.al12.ext_v1.1.re],
                 cluster.endoderm.color.v5)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
src.endoderm.al12.ext_v1.1.re.coltree =
  MyHeatmap(as.matrix(src.endoderm.al12.ext_v1.1.re@assays$RNA@data[
    unique(c(gene_src.al2.tracing.rowtree.1,
             gene_src.al2.tracing.rowtree.1)),
    cellorder_cell_src.endoderm.al12.ext_v1.1.re]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.endoderm.al12.ext_v1.1.re$Time[cellorder_cell_src.endoderm.al12.ext_v1.1.re],
                 colors.time.2),
      MyName2Col(src.endoderm.al12.ext_v1.1.re$cluster.extract.v1.1[cellorder_cell_src.endoderm.al12.ext_v1.1.re], 
                 cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.al12.ext_v1.1.re$cluster.v06.26.re_hc[cellorder_cell_src.endoderm.al12.ext_v1.1.re],
                 cluster.endoderm.color.v5)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "col",
    graph = T)
dev.off()


tree_src.endoderm.al12.ext_v1.1.re.coltree = as.dendrogram(src.endoderm.al12.ext_v1.1.re.coltree)
src.endoderm.al12.ext_v1.1.re$tree = NA
src.endoderm.al12.ext_v1.1.re@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re.coltree[[1]][[1]][[1]]),]$tree = 1
src.endoderm.al12.ext_v1.1.re@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re.coltree[[1]][[1]][[2]]),]$tree = 2
src.endoderm.al12.ext_v1.1.re@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re.coltree[[1]][[2]][[1]]),]$tree = 3
src.endoderm.al12.ext_v1.1.re@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re.coltree[[1]][[2]][[2]]),]$tree = 4
src.endoderm.al12.ext_v1.1.re@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re.coltree[[2]][[1]][[1]]),]$tree = 5
src.endoderm.al12.ext_v1.1.re@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re.coltree[[2]][[1]][[2]]),]$tree = 6
src.endoderm.al12.ext_v1.1.re@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re.coltree[[2]][[2]][[1]][[1]]),]$tree = 7
src.endoderm.al12.ext_v1.1.re@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re.coltree[[2]][[2]][[1]][[2]]),]$tree = 10
src.endoderm.al12.ext_v1.1.re@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re.coltree[[2]][[2]][[2]][[1]]),]$tree = 8
src.endoderm.al12.ext_v1.1.re@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re.coltree[[2]][[2]][[2]][[2]]),]$tree = 9
DimPlot(src.endoderm.al12.ext_v1.1.re, reduction = 'umap_integrated', group.by = 'tree', label = T, label.size = 5) + 
  DimPlot(src.endoderm.al12.ext_v1.1.re, reduction = 'umap_integrated', group.by = 'Time') # , label = T, label.size = 5)

src.endoderm.al12.ext_v1.1.re$cluster.v06.26.re_correct = src.endoderm.al12.ext_v1.1.re$cluster.v06.26.re_hc
src.endoderm.al12.ext_v1.1.re@meta.data[src.endoderm.al12.ext_v1.1.re$tree%in%c(5,6),]$cluster.v06.26.re_correct = "AL.2"
src.endoderm.al12.ext_v1.1.re@meta.data[src.endoderm.al12.ext_v1.1.re$tree%in%c(7,9,10),]$cluster.v06.26.re_correct = "AL.1/2-Liver"
src.endoderm.al12.ext_v1.1.re@meta.data[src.endoderm.al12.ext_v1.1.re$tree%in%c(1,2,3,4,8),]$cluster.v06.26.re_correct = "Liver"
DimPlot(src.endoderm.al12.ext_v1.1.re, reduction = "umap_integrated", group.by = "cluster.v06.26.re_correct")

save(src.endoderm.al12.ext_v1.1.re, file = "figure.v08.07/organ_development_re_v240115/src.endoderm.al12.ext_v1.1.re.Rdata")
save(src.endoderm.al12.ext_v1.1.re_filtergene,
     tree_src.endoderm.al12.ext_v1.1.re.coltree,
     cellorder_cell_src.endoderm.al12.ext_v1.1.re,
     file = "figure.v08.07/organ_development_re_v240115/src.endoderm.al12.ext_v1.1.re.parameter.Rdata")



#-----  Set AL.1/2
#------------------------------
cell_extract = rownames(src.endoderm.al12.ext_v1.1.re@meta.data[
  src.endoderm.al12.ext_v1.1.re$cluster.v06.26.re_correct%in%"AL.2",])

src.endoderm.al12.ext_v1.1.re_al12 = src.endoderm.al12.ext_v1.1.re[, cell_extract]
src.endoderm.al12.ext_v1.1.re_al12 = NormalizeData(src.endoderm.al12.ext_v1.1.re_al12, scale.factor = 10^5)
src.endoderm.al12.ext_v1.1.re_al12 = ScaleData(src.endoderm.al12.ext_v1.1.re_al12, 
                                               features = rownames(src.endoderm.al12.ext_v1.1.re_al12),
                                               split.by = "batch_phase")

src.endoderm.al12.ext_v1.1.re_al12 = FindVariableFeatures(src.endoderm.al12.ext_v1.1.re_al12, assay = "RNA", nfeatures = 2000)
src.endoderm.al12.ext_v1.1.re_al12_filtergene = 
  Myfilter(as.matrix(src.endoderm.al12.ext_v1.1.re_al12@assays$RNA@data),
           gene = src.endoderm.al12.ext_v1.1.re_al12@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.endoderm.al12.ext_v1.1.re_al12_filtergene_re = src.endoderm.al12.ext_v1.1.re_al12_filtergene
rm(src.endoderm.al12.ext_v1.1.re_al12_filtergene)
src.endoderm.al12.ext_v1.1.re_al12_filtergene = src.endoderm.al12.ext_v1.1.re_al12_filtergene_re
cellorder_cell_src.endoderm.al12.ext_v1.1.re_al12 = colnames(src.endoderm.al12.ext_v1.1.re_al12)

cell_extract_al1 = c(
  paste("ss9_", rownames(src.9ss.integrated@meta.data[src.9ss.integrated$cluster.v06.26.re..merge%in%"AL.1",]), sep=''),
  paste('ss12_', rownames(src.12ss.integrated@meta.data[src.12ss.integrated$cluster.extract.v1.1_re%in%"AL.1",]), sep=''))
src.endoderm.al12.ext_v1.1.re_al12$check = "AL.2"
src.endoderm.al12.ext_v1.1.re_al12@meta.data[cell_extract_al1 ,]$check = "AL.1"

src.endoderm.al12.ext_v1.1.re_al12 = RunPCA(src.endoderm.al12.ext_v1.1.re_al12, 
                                            features = src.endoderm.al12.ext_v1.1.re_al12_filtergene)
src.endoderm.al12.ext_v1.1.re_al12 = RunUMAP(src.endoderm.al12.ext_v1.1.re_al12, reduction = "pca", dims = 1:20)
src.endoderm.al12.ext_v1.1.re_al12@reductions$umap_filter = src.endoderm.al12.ext_v1.1.re_al12@reductions$umap

DimPlot(src.endoderm.al12.ext_v1.1.re_al12, reduction = "pca", group.by = "Time", c(1,3)) +
  DimPlot(src.endoderm.al12.ext_v1.1.re_al12, reduction = "pca", group.by = "cluster.v06.26.re_hc") +
  DimPlot(src.endoderm.al12.ext_v1.1.re_al12, reduction = "pca", group.by = "Phase") +
  DimPlot(src.endoderm.al12.ext_v1.1.re_al12, reduction = "pca", group.by = "batch") +
  DimPlot(src.endoderm.al12.ext_v1.1.re_al12, reduction = "pca", group.by = "tree", 
          label = T, label.size = 5, dims = c(1,3)) +
  DimPlot(src.endoderm.al12.ext_v1.1.re_al12, reduction = "pca", group.by = "check", dims = c(4,3)) 


library(dplyr)
pc_src.endoderm.al12.ext_v1.1.re_al12 = src.endoderm.al12.ext_v1.1.re_al12[['pca']]@feature.loadings
pc_src.endoderm.al12.ext_v1.1.re_al12 = as.data.frame(pc_src.endoderm.al12.ext_v1.1.re_al12)
pc_src.endoderm.al12.ext_v1.1.re_al12$PC_3 = as.numeric(pc_src.endoderm.al12.ext_v1.1.re_al12$PC_3)
pc_src.endoderm.al12.ext_v1.1.re_al12 = pc_src.endoderm.al12.ext_v1.1.re_al12[
  order(pc_src.endoderm.al12.ext_v1.1.re_al12$PC_3),]
gene_src.endoderm.al12.ext_v1.1.re_al12_pc3 = c(
  rownames(pc_src.endoderm.al12.ext_v1.1.re_al12)[1:100],
  rownames(pc_src.endoderm.al12.ext_v1.1.re_al12)[(nrow(pc_src.endoderm.al12.ext_v1.1.re_al12)-100):nrow(pc_src.endoderm.al12.ext_v1.1.re_al12)])

src.endoderm.al12.ext_v1.1.re_al12 = RunPCA(src.endoderm.al12.ext_v1.1.re_al12, 
                                            features = gene_src.endoderm.al12.ext_v1.1.re_al12_pc3)
src.endoderm.al12.ext_v1.1.re_al12 = RunUMAP(src.endoderm.al12.ext_v1.1.re_al12, reduction = "pca", dims = 1:10)

pdf("figure.v08.07/organ_development_re_v240115/try.al1_10x.pdf",9,7)
src.endoderm.al12.ext_v1.1.re_al12.rowtree =
  MyHeatmap(as.matrix(src.endoderm.al12.ext_v1.1.re_al12@assays$RNA@data[
    src.endoderm.al12.ext_v1.1.re_al12_filtergene,
    cellorder_cell_src.endoderm.al12.ext_v1.1.re_al12]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.endoderm.al12.ext_v1.1.re_al12$Time[cellorder_cell_src.endoderm.al12.ext_v1.1.re_al12],
                 colors.time),
      # MyName2Col(src.endoderm.al12.ext_v1.1.re_al12$tree[cellorder_cell_src.endoderm.al12.ext_v1.1.re_al12],
      #            color.temp),
      MyName2Col(src.endoderm.al12.ext_v1.1.re_al12$cluster.extract.v1.1[cellorder_cell_src.endoderm.al12.ext_v1.1.re_al12], 
                 cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.al12.ext_v1.1.re_al12$cluster.v06.26.re_hc[cellorder_cell_src.endoderm.al12.ext_v1.1.re_al12],
                 cluster.endoderm.color.v5)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
src.endoderm.al12.ext_v1.1.re_al12.coltree =
  MyHeatmap(as.matrix(src.endoderm.al12.ext_v1.1.re_al12@assays$RNA@data[
    src.endoderm.al12.ext_v1.1.re_al12_filtergene,
    cellorder_cell_src.endoderm.al12.ext_v1.1.re_al12]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.endoderm.al12.ext_v1.1.re_al12$Time[cellorder_cell_src.endoderm.al12.ext_v1.1.re_al12],
                 colors.time.2),
      MyName2Col(src.endoderm.al12.ext_v1.1.re_al12$cluster.extract.v1.1[cellorder_cell_src.endoderm.al12.ext_v1.1.re_al12], 
                 cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.al12.ext_v1.1.re_al12$cluster.v06.26.re_hc[cellorder_cell_src.endoderm.al12.ext_v1.1.re_al12],
                 cluster.endoderm.color.v5)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "col",
    graph = T)
dev.off()

tree_src.endoderm.al12.ext_v1.1.re_al12.rowtree = as.dendrogram(src.endoderm.al12.ext_v1.1.re_al12.rowtree)
tree_src.endoderm.al12.ext_v1.1.re_al12.coltree = as.dendrogram(src.endoderm.al12.ext_v1.1.re_al12.coltree)

src.endoderm.al12.ext_v1.1.re_al12$tree = NA
src.endoderm.al12.ext_v1.1.re_al12@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re_al12.coltree[[1]][[1]][[1]]),]$tree = 1
src.endoderm.al12.ext_v1.1.re_al12@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re_al12.coltree[[1]][[1]][[2]]),]$tree = 2
src.endoderm.al12.ext_v1.1.re_al12@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re_al12.coltree[[1]][[2]][[1]]),]$tree = 3
src.endoderm.al12.ext_v1.1.re_al12@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re_al12.coltree[[1]][[2]][[2]]),]$tree = 4
src.endoderm.al12.ext_v1.1.re_al12@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re_al12.coltree[[2]][[1]][[1]]),]$tree = 5
src.endoderm.al12.ext_v1.1.re_al12@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re_al12.coltree[[2]][[1]][[2]]),]$tree = 6
src.endoderm.al12.ext_v1.1.re_al12@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re_al12.coltree[[2]][[2]][[1]]),]$tree = 7
src.endoderm.al12.ext_v1.1.re_al12@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re_al12.coltree[[2]][[2]][[2]]),]$tree = 8


pdf("figure.v08.07/organ_development_re_v240115/try.al1_10x.pdf",9,7)
src.endoderm.al12.ext_v1.1.re_al12.rowtree.1 =
  MyHeatmap(as.matrix(src.endoderm.al12.ext_v1.1.re_al12@assays$RNA@data[
    gene_src.endoderm.al12.ext_v1.1.re_al12_pc3,
    cellorder_cell_src.endoderm.al12.ext_v1.1.re_al12]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.endoderm.al12.ext_v1.1.re_al12$Time[cellorder_cell_src.endoderm.al12.ext_v1.1.re_al12],
                 colors.time),
      MyName2Col(src.endoderm.al12.ext_v1.1.re_al12$check[cellorder_cell_src.endoderm.al12.ext_v1.1.re_al12],
                 cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.al12.ext_v1.1.re_al12$tree.1[cellorder_cell_src.endoderm.al12.ext_v1.1.re_al12],
                 color.temp),
      MyName2Col(src.endoderm.al12.ext_v1.1.re_al12$cluster.extract.v1.1[cellorder_cell_src.endoderm.al12.ext_v1.1.re_al12], 
                 cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.al12.ext_v1.1.re_al12$cluster.v06.26.re_hc[cellorder_cell_src.endoderm.al12.ext_v1.1.re_al12],
                 cluster.endoderm.color.v5)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)

src.endoderm.al12.ext_v1.1.re_al12.coltree.1 =
  MyHeatmap(as.matrix(src.endoderm.al12.ext_v1.1.re_al12@assays$RNA@data[
    gene_src.endoderm.al12.ext_v1.1.re_al12_pc3,
    cellorder_cell_src.endoderm.al12.ext_v1.1.re_al12]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.endoderm.al12.ext_v1.1.re_al12$Time[cellorder_cell_src.endoderm.al12.ext_v1.1.re_al12],
                 colors.time.2),
      MyName2Col(src.endoderm.al12.ext_v1.1.re_al12$cluster.extract.v1.1[cellorder_cell_src.endoderm.al12.ext_v1.1.re_al12], 
                 cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.al12.ext_v1.1.re_al12$cluster.v06.26.re_hc[cellorder_cell_src.endoderm.al12.ext_v1.1.re_al12],
                 cluster.endoderm.color.v5)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "col",
    graph = T)
dev.off()

tree_src.endoderm.al12.ext_v1.1.re_al12.rowtree.1 = as.dendrogram(src.endoderm.al12.ext_v1.1.re_al12.rowtree.1)
tree_src.endoderm.al12.ext_v1.1.re_al12.coltree.1 = as.dendrogram(src.endoderm.al12.ext_v1.1.re_al12.coltree.1)  

src.endoderm.al12.ext_v1.1.re_al12$tree.1 = NA
src.endoderm.al12.ext_v1.1.re_al12@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re_al12.coltree.1[[1]][[1]][[1]]),]$tree.1 = 1
src.endoderm.al12.ext_v1.1.re_al12@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re_al12.coltree.1[[1]][[1]][[2]]),]$tree.1 = 2
src.endoderm.al12.ext_v1.1.re_al12@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re_al12.coltree.1[[1]][[2]][[1]]),]$tree.1 = 3
src.endoderm.al12.ext_v1.1.re_al12@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re_al12.coltree.1[[1]][[2]][[2]]),]$tree.1 = 4
src.endoderm.al12.ext_v1.1.re_al12@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re_al12.coltree.1[[2]][[1]][[1]]),]$tree.1 = 5
src.endoderm.al12.ext_v1.1.re_al12@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re_al12.coltree.1[[2]][[1]][[2]]),]$tree.1 = 6
src.endoderm.al12.ext_v1.1.re_al12@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re_al12.coltree.1[[2]][[2]][[1]][[1]]),]$tree.1 = 7
src.endoderm.al12.ext_v1.1.re_al12@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re_al12.coltree.1[[2]][[2]][[1]][[2]]),]$tree.1 = 8
src.endoderm.al12.ext_v1.1.re_al12@meta.data[labels(tree_src.endoderm.al12.ext_v1.1.re_al12.coltree.1[[2]][[2]][[2]]),]$tree.1 = 9

DimPlot(src.endoderm.al12.ext_v1.1.re_al12, group.by = "Time", label = T, label.size = 5) +
  DimPlot(src.endoderm.al12.ext_v1.1.re_al12, group.by = "check", label = T, label.size = 5) +
  DimPlot(src.endoderm.al12.ext_v1.1.re_al12, group.by = "tree.1", label = T, label.size = 5)

src.endoderm.al12.ext_v1.1.re_al12$cluster.v06.26.re_correct_re = src.endoderm.al12.ext_v1.1.re_al12$cluster.v06.26.re_correct
src.endoderm.al12.ext_v1.1.re_al12@meta.data[src.endoderm.al12.ext_v1.1.re_al12$tree.1%in%c(5,6,8,1),]$cluster.v06.26.re_correct_re = "AL.1"
src.endoderm.al12.ext_v1.1.re_al12@meta.data[src.endoderm.al12.ext_v1.1.re_al12$tree.1%in%c(9,7,2,3,4),]$cluster.v06.26.re_correct_re = "AL.2"
src.endoderm.al12.ext_v1.1.re_al12@meta.data[
  paste("ss9_", rownames(src.9ss.integrated@meta.data[src.9ss.integrated$cluster.v06.26.re..merge%in%"AL.1",]), sep=''),]$cluster.v06.26.re_correct_re = "AL.1"


src.endoderm.al12.ext_v1.1.re$cluster.v06.26.re_correct_re = src.endoderm.al12.ext_v1.1.re$cluster.v06.26.re_correct
src.endoderm.al12.ext_v1.1.re@meta.data[colnames(src.endoderm.al12.ext_v1.1.re_al12),]$cluster.v06.26.re_correct_re =
  src.endoderm.al12.ext_v1.1.re_al12$cluster.v06.26.re_correct_re
DimPlot(src.endoderm.al12.ext_v1.1.re, group.by = "cluster.v06.26.re_correct_re", reduction = "umap_integrated")


save(src.endoderm.al12.ext_v1.1.re, file = 'figure.v08.07/organ_development_re_v240115/src.endoderm.al12.ext_v1.1.re.Rdata')
save(src.endoderm.al12.ext_v1.1.re_al12, file = 'figure.v08.07/organ_development_re_v240115/src.endoderm.al12.ext_v1.1.re_al12.Rdata')

save(tree_src.endoderm.al12.ext_v1.1.re.coltree,
     src.endoderm.al12.ext_v1.1.re_al12_filtergene,
     gene_src.endoderm.al12.ext_v1.1.re_al12_pc3,
     tree_src.endoderm.al12.ext_v1.1.re_al12.rowtree,
     tree_src.endoderm.al12.ext_v1.1.re_al12.rowtree.1,
     tree_src.endoderm.al12.ext_v1.1.re_al12.coltree,
     tree_src.endoderm.al12.ext_v1.1.re_al12.coltree.1,
     file = 'figure.v08.07/organ_development_re_v240115/src.endoderm.al12.ext_v1.1.re.parameter.Rdata')


#-- 12SS: cluster.extract.v1.1_re for AL.1 / AL.2 10X
src.12ss.integrated = FindNeighbors(src.12ss.integrated, reduction = "pca", dims = 1:30)
src.12ss.integrated = FindClusters(src.12ss.integrated, resolution = 15)
FeaturePlot(src.12ss.integrated, features = "Nkx2-3")
DimPlot(src.12ss.integrated, group.by = "cluster.extract.v1.1")
DimPlot(src.12ss.integrated, group.by = "RNA_snn_res.10", label = T, label.size = 4)
DimPlot(src.12ss.integrated, group.by = "RNA_snn_res.15", label = T, label.size = 4)

src.12ss.integrated$cluster.extract.v1.1_re = src.12ss.integrated$cluster.extract.v1.1
src.12ss.integrated@meta.data[
  (src.12ss.integrated$RNA_snn_res.10%in%c(3,16,80,113)|
     src.12ss.integrated$RNA_snn_res.15%in%c(56,59,78)) & 
    src.12ss.integrated$cluster.extract.v1.1%in%"AL.1/2",]$cluster.extract.v1.1_re = "AL.1"
src.12ss.integrated@meta.data[
  (!(src.12ss.integrated$RNA_snn_res.10%in%c(3,16,80,113)|
     src.12ss.integrated$RNA_snn_res.15%in%c(56,59,78))) & 
    src.12ss.integrated$cluster.extract.v1.1%in%"AL.1/2",]$cluster.extract.v1.1_re = "AL.2"
DimPlot(src.12ss.integrated, group.by = "cluster.extract.v1.1_re", 
        cols = cluster.endoderm.color.v5)
save(src.12ss.integrated, file = "figure.v08.07/organ_development_re_v240115//src.12ss.integrated.Rdata")



#-----  Correct for AL.1 / AL.2 10X
src.al1.integrated.re = src.al1.integrated
src.al1.integrated.re$cluster.v06.26.re_correct_re = 
  src.endoderm.al12.ext_v1.1.re$cluster.v06.26.re_correct_re[colnames(src.al1.integrated.re)]
src.al1.integrated.re = src.al1.integrated.re[,!src.al1.integrated.re$cluster.v06.26.re_correct_re%in%"AL.2"]
DimPlot(src.al1.integrated.re, reduction = "umap_integrated",
        group.by = "cluster.v06.26.re_correct_re")

src.al2.integrated.re = src.al2.integrated
src.al2.integrated.re$cluster.v06.26.re_correct_re = 
  src.endoderm.al12.ext_v1.1.re$cluster.v06.26.re_correct_re[colnames(src.al2.integrated.re)]
src.al2.integrated.re = src.al2.integrated.re[,!src.al2.integrated.re$cluster.v06.26.re_correct_re%in%"AL.1"]
DimPlot(src.al2.integrated.re, reduction = "umap_integrated",
        group.by = "cluster.v06.26.re_correct_re")

save(src.al1.integrated.re, file = "figure.v08.07/organ_development_re_v240115//src.al1.integrated.re.Rdata")
save(src.al2.integrated.re, file = "figure.v08.07/organ_development_re_v240115//src.al2.integrated.re.Rdata")


#-----  Correct for AL.3 10X
src.al3.integrated$cluster.v06.26.re_correct = src.al3.integrated$cluster.v06.26.re
src.al3.integrated@meta.data[
  intersect(colnames(src.al3.integrated),
            colnames(src.endoderm.al12.ext_v1.1.re)),]$cluster.v06.26.re_correct =
  src.endoderm.al12.ext_v1.1.re@meta.data[
    intersect(colnames(src.al3.integrated),
              colnames(src.endoderm.al12.ext_v1.1.re)),]$cluster.v06.26.re_correct 

src.al3.integrated.re = src.al3.integrated[, !src.al3.integrated$cluster.v06.26.re_correct%in%'AL.2']
DimPlot(src.al3.integrated.re, reduction = "umap_integrated", group.by = "cluster.v06.26.re_correct")
save(src.al3.integrated.re, file = "figure.v08.07/organ_development_re_v240115/src.al3.integrated.re.Rdata")

#-----  Correct for Liver 10X
src.liver.integrated = src.endoderm.liver.ext_v1.1.re

src.liver.integrated$cluster.v06.26.re_correct = src.liver.integrated$cluster.v06.26.re..merge
src.liver.integrated@meta.data[
  intersect(colnames(src.liver.integrated),
            colnames(src.endoderm.al12.ext_v1.1.re)),]$cluster.v06.26.re_correct =
  src.endoderm.al12.ext_v1.1.re@meta.data[
    intersect(colnames(src.liver.integrated),
              colnames(src.endoderm.al12.ext_v1.1.re)),]$cluster.v06.26.re_correct 

src.liver.integrated@meta.data[
  intersect(colnames(src.liver.integrated),
            colnames(src.fg4.integrated_re_fg4)),]$cluster.v06.26.re_correct =
  src.fg4.integrated_re_fg4@meta.data[
    intersect(colnames(src.liver.integrated),
              colnames(src.fg4.integrated_re_fg4)),]$cluster.v06.26.re_correct_re.tree3

src.liver.integrated = src.liver.integrated[,!src.liver.integrated$cluster.v06.26.re_correct%in%c("FG.4-Lung/Stomach")]
DimPlot(src.liver.integrated, reduction = "umap_mnn", group.by = "cluster.v06.26.re_correct")
save(src.liver.integrated, file = 'figure.v08.07/organ_development_re_v240115/src.liver.integrated.Rdata')


#---------------------------------
# -- AL.3 tracing
#=================================================================================
src.al3.integrated.merge.re = src.al3.integrated.merge[, c(colnames(src.al3.integrated.re),
                                                           colnames((src.al3.tracing)))]

DimPlot(src.al3.integrated.merge.re, reduction = "mnn_umap_fta", group.by = 'Time')
DimPlot(src.al3.integrated.merge.re, reduction = "mnn_umap_fta", 
        group.by = 'lineage', cols = color.lineage, na.value = "#eeeeee")


# Set for function
#----------------------------
src.al3.integrated.selectgene = src.endoderm.al3.ext_v1.1.gene.fin

seurat = src.al3.integrated.merge.re
# seurat.filtergene = src.endoderm.al3.ext_v1.1.gene.fin
# seurat = ScaleData(seurat, features = seurat.filtergene, split.by = 'Source_tech')
# seurat = RunPCA(seurat, seurat.filtergene, assay = 'RNA')
# seurat = RunUMAP(seurat, reduction = "pca", dims = 1:30, n.components = 3, assay = "RNA")
cell_refer1 = rownames(seurat@meta.data[seurat$Source_tech%in%"refer" & seurat$batch%in%1, ])
cell_refer2 = rownames(seurat@meta.data[seurat$Source_tech%in%"refer" & seurat$batch%in%2, ])
cell_query = rownames(seurat@meta.data[seurat$Source_tech%in%"query",])
seurat.selectgene = src.al3.integrated.selectgene
src.al3.integrated = NormalizeData(src.al3.integrated, assay = "RNA", scale.factor = 10^5)


for(assay in c("RNA")){
  
  anchor.integrated.merge = 
    FindIntegrationAnchors(c(src.al3.integrated[,cell_refer1],
                             src.al3.integrated[,cell_refer2],
                             seurat[,cell_query]
    ),
    assay = c(assay, assay, assay), 
    reduction = "rpca", 
    anchor.features = seurat.selectgene,
    dims = 1:30)
  
  seurat.re = IntegrateData(anchorset = anchor.integrated.merge, dims = 1:30)
  DefaultAssay(seurat.re) = "integrated"
  seurat.re = ScaleData(seurat.re, features = rownames(seurat.re))
  seurat.re = RunPCA(seurat.re, features = seurat.selectgene)
  seurat.re = RunUMAP(seurat.re, dims = 1:30, n.neighbors = 75,
                      n.components=3)
  
  # seurat.re@reductions$umap
  DimPlot(seurat.re, reduction = 'umap',
          # group.by = "Treatment",
          # group.by = "Time",
          # group.by = "lineage",
          group.by = "cluster.v06.26.re_correct",
          # group.by = "cluster.v06.26.re_mnn_umap_fta",
          pt.size = 1.5, dims = c(1,2),
          na.value = "#eeeeee",
          cols = c(# colors.time,
            colors.time.2, color.lineage, cluster.endoderm.color.v5,colors.type))
  
  
  assay_name_re = gsub('RNA',"", gsub("mnn_","mnn",assay))
  seurat.re[[paste(assay_name_re,"pca_integrated",sep="")]] = seurat.re[["pca"]]
  seurat.re[[paste(assay_name_re,"umap_integrated",sep="")]] = seurat.re[["umap"]]
  
  seurat = RunUMAP(seurat.re, dims = 1:30, n.neighbors = 75, n.components=3)
  
  seurat[[paste(assay_name_re,"pca_integrated",sep="")]] = seurat[["pca"]]
  seurat[[paste(assay_name_re,"pca_integrated",sep="")]]@cell.embeddings[colnames(seurat),] = 
    seurat.re[[paste(assay_name_re,"pca_integrated",sep="")]]@cell.embeddings[colnames(seurat),] 
  
  seurat[[paste(assay_name_re,"umap_integrated",sep="")]] = seurat[["umap"]]
  seurat[[paste(assay_name_re,"umap_integrated",sep="")]]@cell.embeddings[colnames(seurat),] = 
    seurat.re[[paste(assay_name_re,"umap_integrated",sep="")]]@cell.embeddings[colnames(seurat),] 
}

label_learning = FNN::knn(seurat[["umap_integrated"]]@cell.embeddings[c(cell_refer1,cell_refer2),], 
                          seurat[["umap_integrated"]]@cell.embeddings[cell_query,],
                          seurat$cluster.v06.26.re_correct[c(cell_refer1,cell_refer2)], k = 10)
seurat$cluster.v06.26.re_correct_mnn_umap_fta = NA
seurat$cluster.v06.26.re_correct_mnn_umap_fta[c(cell_refer1,cell_refer2)] =
  seurat$cluster.v06.26.re_correct[c(cell_refer1,cell_refer2)]
seurat$cluster.v06.26.re_correct_mnn_umap_fta[cell_query]=
  as.character(label_learning)

seurat@reductions$mnn_umap_fta = seurat@reductions$umap
seurat@reductions$mnn_umap_fta@cell.embeddings[,c(1,2)] = 
  src.al3.integrated.merge.re@reductions$mnn_umap_fta@cell.embeddings[colnames(seurat),c(1,2)]

DimPlot(seurat[,cell_query], reduction = "umap_integrated", 
        group.by = 'Time', 
        cols = colors.time.2, na.value = "#eeeeee")
DimPlot(seurat[,cell_query], reduction = "mnn_umap_fta", 
        group.by = 'cluster.v06.26.re_correct_mnn_umap_fta', 
        cols = cluster.endoderm.color.v5, na.value = "#eeeeee")

# set for src.al3.integrated.merge.re
src.al3.integrated.merge.re@reductions$umap_integrated = src.al3.integrated.merge.re@reductions$umap
src.al3.integrated.merge.re@reductions$umap_integrated@cell.embeddings = 
  seurat@reductions$umap_integrated@cell.embeddings
colnames(src.al3.integrated.merge.re@reductions$umap_integrated@cell.embeddings) = paste("UMAP_",c(1:3),sep="")
src.al3.integrated.merge.re@reductions$umap_integrated@key = "UMAP_"

src.al3.integrated.merge.re$cluster.v06.26.re_correct_mnn_umap_fta =
  seurat$cluster.v06.26.re_correct_mnn_umap_fta

DimPlot(src.al3.integrated.merge.re, reduction = "mnn_umap_fta",
        # group.by = "cluster.v06.26.re_correct",
        group.by = "cluster.v06.26.re_mnn_umap_fta",
        pt.size = 1.5, dims = c(1,2),
        na.value = "#eeeeee",
        cols = c(# colors.time,
          colors.time.2, color.lineage, cluster.endoderm.color.v5,colors.type))
save(src.al3.integrated.merge.re, file = "figure.v08.07/organ_development_re_v240115/src.al3.integrated.merge.re.Rdata")

#----------------------------

src.al3.tracing$cluster.v06.26.re_correct_mnn_umap_fta = 
  src.al3.integrated.merge.re$cluster.v06.26.re_correct_mnn_umap_fta[colnames(src.al3.tracing)]
DimPlot(src.al3.tracing, reduction = "mnn_umap_fta", 
        group.by = 'Time', cols = colors.time.2, na.value = "#eeeeee")
DimPlot(src.al3.tracing, reduction = "mnn_umap_fta", 
        group.by = 'lineage', cols = color.lineage, na.value = "#eeeeee")
DimPlot(src.al3.tracing, reduction = "mnn_umap_fta", 
        group.by = 'cluster.v06.26.re_mnn_umap_fta', 
        cols = cluster.endoderm.color.v5, na.value = "#eeeeee")


src.al3.tracing_al3 = src.al3.tracing[,src.al3.tracing$cluster.v06.26.re_correct_mnn_umap_fta%in%c(
  "AL.3","AL.3-Small.intestine.1","AL.3-EHBD/VP","AL.3-Liver")]
src.al3.tracing_al3 = ScaleData(src.al3.tracing_al3, 
                                features = rownames(src.al3.tracing_al3), split.by = "Phase")
src.al3.tracing_al3 = RunPCA(src.al3.tracing_al3, 
                             features = gene_src.al3.integrated_re_al3.rowtree.2)
src.al3.tracing_al3 = RunUMAP(src.al3.tracing_al3, reduction = 'pca', dims = 1:30)
DimPlot(src.al3.tracing_al3, reduction = "pca",
        group.by = 'Time')
DimPlot(src.al3.tracing_al3, reduction = "pca",
        group.by = 'cluster.v06.26.re_correct_mnn_umap_fta')
DimPlot(src.al3.tracing_al3, reduction = "umap",
        group.by = 'cluster.v06.26.re_correct_mnn_umap_fta')



src.al3.tracing_al3 = FindVariableFeatures(src.al3.tracing_al3)
src.al3.tracing_al3_filtergene = 
  Myfilter(as.matrix(src.al3.tracing_al3@assays$RNA@data),
           gene = unique(c(
             src.al3.tracing_al3@assays$RNA@var.features)),
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.al3.tracing_al3_filtergene_re = src.al3.tracing_al3_filtergene
rm(src.al3.tracing_al3_filtergene)
src.al3.tracing_al3_filtergene = src.al3.tracing_al3_filtergene_re


pdf("figure.v08.07/organ_development_re_v240115/try.al3_al3_sm3_tree.pdf",9,7)
src.al3.tracing_al3.rowtree =
  MyHeatmap(as.matrix(src.al3.tracing_al3@assays$RNA@data[
    src.al3.tracing_al3_filtergene, ]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.al3.tracing_al3$Time,
                 colors.time.2),
      MyName2Col(src.al3.tracing_al3$cluster.v06.26.re_correct_mnn_umap_fta,
                 cluster.endoderm.color.v5)),
    #RowSideColors = t(cbind(
    #  MyName2Col(names(src.al3.tracing_al3_filtergene_re), colors.geneset))),
    ColSideColorsSize = 3,
    RowSideColorsSize = 1.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    return.tree = "row",
    # Colv="none", Rowv = "none",
    graph = T)
src.al3.tracing_al3.coltree =
  MyHeatmap(as.matrix(src.al3.tracing_al3@assays$RNA@data[
    src.al3.tracing_al3_filtergene, ]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.al3.tracing_al3$Time,
                 colors.time.2),
      MyName2Col(src.al3.tracing_al3$cluster.v06.26.re_correct_mnn_umap_fta,
                 cluster.endoderm.color.v5)),
    #RowSideColors = t(cbind(
    #  MyName2Col(names(src.al3.tracing_al3_filtergene_re), colors.geneset))),
    ColSideColorsSize = 3,
    RowSideColorsSize = 1.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    return.tree = "col",
    # Colv="none", Rowv = "none",
    graph = T)
dev.off()

tree_src.al3.tracing_al3.rowtree = as.dendrogram(src.al3.tracing_al3.rowtree)
tree_src.al3.tracing_al3.coltree = as.dendrogram(src.al3.tracing_al3.coltree)

gene_src.al3.tracing_al3.rowtree = setdiff(
  src.al3.tracing_al3_filtergene,
  c(labels(tree_src.al3.tracing_al3.rowtree[[2]][[2]][[2]][[2]][[1]]),
    labels(tree_src.al3.tracing_al3.rowtree[[2]][[2]][[1]][[2]][[2]]),
    labels(tree_src.al3.tracing_al3.rowtree[[2]][[2]][[1]][[2]][[1]][[2]])))

pdf("figure.v08.07/organ_development_re_v240115/try.al3_al3_sm3.pdf",9,7)
src.al3.tracing_al3.rowtree.1 =
  MyHeatmap(as.matrix(src.al3.tracing_al3@assays$RNA@data[
    gene_src.al3.tracing_al3.rowtree, ]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.al3.tracing_al3$Time,
                 colors.time.2),
      MyName2Col(src.al3.tracing_al3$cluster.v06.26.re_correct_mnn_umap_fta,
                 cluster.endoderm.color.v5)),
    #RowSideColors = t(cbind(
    #  MyName2Col(names(src.al3.tracing_al3_filtergene_re), colors.geneset))),
    ColSideColorsSize = 3,
    RowSideColorsSize = 1.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    return.tree = "row",
    # Colv="none", Rowv = "none",
    graph = T)
src.al3.tracing_al3.coltree.1 =
  MyHeatmap(as.matrix(src.al3.tracing_al3@assays$RNA@data[
    gene_src.al3.tracing_al3.rowtree, ]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.al3.tracing_al3$Time,
                 colors.time.2),
      MyName2Col(src.al3.tracing_al3$cluster.v06.26.re_correct_mnn_umap_fta,
                 cluster.endoderm.color.v5)),
    #RowSideColors = t(cbind(
    #  MyName2Col(names(src.al3.tracing_al3_filtergene_re), colors.geneset))),
    ColSideColorsSize = 3,
    RowSideColorsSize = 1.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    return.tree = "col",
    # Colv="none", Rowv = "none",
    graph = T)
dev.off()

pdf("figure.v08.07/organ_development_re_v240115/try.al3_al3_sm3_tree.pdf",100,20)
plot(src.al3.tracing_al3.rowtree.1)
dev.off()

tree_src.al3.tracing_al3.rowtree.1 = as.dendrogram(src.al3.tracing_al3.rowtree.1)
tree_src.al3.tracing_al3.coltree.1 = as.dendrogram(src.al3.tracing_al3.coltree.1)
gene_src.al3.tracing_al3.rowtree.1 = setdiff(
  gene_src.al3.tracing_al3.rowtree,
  c(labels(tree_src.al3.tracing_al3.rowtree.1[[2]][[1]])))

pdf("figure.v08.07/organ_development_re_v240115/try.al3_al3_sm3.pdf",9,7)
src.al3.tracing_al3.rowtree.2 =
  MyHeatmap(as.matrix(src.al3.tracing_al3@assays$RNA@data[
    gene_src.al3.tracing_al3.rowtree.1, ]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.al3.tracing_al3$Time,
                 colors.time.2),
      MyName2Col(src.al3.tracing_al3$cluster.v06.26.re_correct_mnn_umap_fta,
                 cluster.endoderm.color.v5)),
    #RowSideColors = t(cbind(
    #  MyName2Col(names(src.al3.tracing_al3_filtergene_re), colors.geneset))),
    ColSideColorsSize = 3,
    RowSideColorsSize = 1.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    return.tree = "row",
    # Colv="none", Rowv = "none",
    graph = T)
src.al3.tracing_al3.coltree.2 =
  MyHeatmap(as.matrix(src.al3.tracing_al3@assays$RNA@data[
    gene_src.al3.tracing_al3.rowtree.1, ]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.al3.tracing_al3$Time,
                 colors.time.2),
      MyName2Col(src.al3.tracing_al3$cluster.v06.26.re_correct_mnn_umap_fta,
                 cluster.endoderm.color.v5)),
    #RowSideColors = t(cbind(
    #  MyName2Col(names(src.al3.tracing_al3_filtergene_re), colors.geneset))),
    ColSideColorsSize = 3,
    RowSideColorsSize = 1.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    return.tree = "col",
    # Colv="none", Rowv = "none",
    graph = T)
dev.off()













