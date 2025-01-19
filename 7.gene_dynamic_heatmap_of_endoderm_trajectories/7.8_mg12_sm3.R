#---------------------------------------------------------------------------
# --          MG.2 -- Re
#---------------------------------------------------------------------------
DimPlot(src.9ss.integrated.merge, group.by = "lineage", reduction = "umap_fta")
DimPlot(src.9ss.integrated.merge, group.by = "cluster.v06.26.re..merge_umap_fta", reduction = "umap_fta")

src.mg2.tracing@meta.data[
  intersect(rownames(src.9ss.integrated.merge@meta.data[src.9ss.integrated.merge$cluster.predict.umap_int.ext.v1.1%in%"MG.2",]),
            rownames(src.mg2.tracing@meta.data)),]$cluster.v06.26.re_mnn_umap_fta = "MG.2"

#------------------------
src.mg2.tracing = FindVariableFeatures(src.mg2.tracing, nfeatures = 2000)
src.mg2.tracing.filtergene = 
  Myfilter(as.matrix(src.mg2.tracing@assays$RNA@data),
           gene = src.mg2.tracing@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.mg2.tracing.filtergene.re = src.mg2.tracing.filtergene; rm(src.mg2.tracing.filtergene)
src.mg2.tracing.filtergene = src.mg2.tracing.filtergene.re; rm(src.mg2.tracing.filtergene.re)


src.mg2.tracing = SetIdent(src.mg2.tracing, value = src.mg2.tracing$cluster.v06.26.re_mnn_umap_fta)
marker_src.mg2.tracing = FindAllMarkers(src.mg2.tracing)
marker_src.mg2.tracing$pct.ratio = marker_src.mg2.tracing$pct.1 / marker_src.mg2.tracing$pct.2
marker_src.mg2.tracing$rank = marker_src.mg2.tracing$pct.ratio * (-log(marker_src.mg2.tracing$p_val_adj))
marker_src.mg2.tracing = marker_src.mg2.tracing[order(marker_src.mg2.tracing$rank, decreasing = T),]
markergene_src.mg2.tracing = unique(marker_src.mg2.tracing$gene)
# markergene_src.mg2.tracing.raw = markergene_src.mg2.tracing

src.mg2.tracing = ScaleData(src.mg2.tracing, 
                            rownames(src.mg2.tracing),split.by = "Phase")
src.mg2.tracing = RunPCA(src.mg2.tracing, 
                         features = src.mg2.tracing.filtergene)
src.mg2.tracing = RunUMAP(src.mg2.tracing, reduction = "pca", 
                          dims = 1:30, n.neighbors = 150,
                          assay = "RNA", n.components = 3)


# == MG.2 Basic process 
#===============================================================================
pdf("figure.v08.07/organ_development_re_v240115/try.mg2.pdf",9,7)
src.mg2.tracing.rowtree =
  MyHeatmap(as.matrix(src.mg2.tracing@assays$RNA@data[
    unique(c(markergene_src.mg2.tracing,
             src.mg2.tracing.filtergene,
             c()
    )),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.mg2.tracing$Time,
                 colors.time.2),
      MyName2Col(src.mg2.tracing$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.mg2.tracing$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.mg2.tracing$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()
pdf("figure.v08.07/organ_development_re_v240115/try.mg2.tree.pdf",150,20)
plot(src.mg2.tracing.rowtree)
dev.off()
tree_src.mg2.tracing.rowtree = as.dendrogram(src.mg2.tracing.rowtree)
gene_src.mg2.tracing.rowtree = setdiff(
  unique(c(markergene_src.mg2.tracing, src.mg2.tracing.filtergene)),
  c(labels(tree_src.mg2.tracing.rowtree[[1]][[2]]),
    labels(tree_src.mg2.tracing.rowtree[[2]][[2]][[2]][[1]])))


pdf("figure.v08.07/organ_development_re_v240115/try.mg2.pdf",9,7)
src.mg2.tracing.rowtree.1 =
  MyHeatmap(as.matrix(src.mg2.tracing@assays$RNA@data[
    gene_src.mg2.tracing.rowtree,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.mg2.tracing$Time,
                 colors.time.2),
      MyName2Col(src.mg2.tracing$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.mg2.tracing$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.mg2.tracing$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
src.mg2.tracing.coltree.1 =
  MyHeatmap(as.matrix(src.mg2.tracing@assays$RNA@data[
    gene_src.mg2.tracing.rowtree,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.mg2.tracing$Time,
                 colors.time.2),
      # MyName2Col(src.mg2.tracing$tree.1,
      #            color.temp),
      MyName2Col(src.mg2.tracing$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.mg2.tracing$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.mg2.tracing$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "col",
    graph = T)
dev.off()
pdf("figure.v08.07/organ_development_re_v240115/try.mg2.tree.pdf",150,20)
plot(src.mg2.tracing.rowtree.1)
dev.off()

tree_src.mg2.tracing.rowtree.1 = as.dendrogram(src.mg2.tracing.rowtree.1)
gene_src.mg2.tracing.rowtree.1 = c(
  labels(tree_src.mg2.tracing.rowtree.1[[2]][[1]][[2]]),
  labels(tree_src.mg2.tracing.rowtree.1[[2]][[1]][[1]][[2]]),
  
  labels(tree_src.mg2.tracing.rowtree.1[[2]][[2]][[1]][[1]]),
  labels(tree_src.mg2.tracing.rowtree.1[[2]][[2]][[1]][[2]][[1]]),
  labels(tree_src.mg2.tracing.rowtree.1[[2]][[2]][[2]][[2]][[1]]),

  #labels(tree_src.mg2.tracing.rowtree.1[[2]][[2]][[2]][[1]]),
  rev(c(
    labels(tree_src.mg2.tracing.rowtree.1[[1]][[2]][[2]][[2]]),
    labels(tree_src.mg2.tracing.rowtree.1[[1]][[1]]),
    labels(tree_src.mg2.tracing.rowtree.1[[1]][[2]][[1]][[1]])
  )))

names(gene_src.mg2.tracing.rowtree.1) = c(
  rep(4, length(c(
    labels(tree_src.mg2.tracing.rowtree.1[[2]][[1]][[2]]),
    labels(tree_src.mg2.tracing.rowtree.1[[2]][[1]][[1]][[2]])
  ))),
  rep(3,length(c(
    labels(tree_src.mg2.tracing.rowtree.1[[2]][[2]][[2]][[2]][[1]]),
    labels(tree_src.mg2.tracing.rowtree.1[[2]][[2]][[1]][[1]]),
    labels(tree_src.mg2.tracing.rowtree.1[[2]][[2]][[1]][[2]][[1]])
  ))),
  rep(9, length(c(
    labels(tree_src.mg2.tracing.rowtree.1[[1]][[1]]),
    labels(tree_src.mg2.tracing.rowtree.1[[1]][[2]][[1]][[1]])
  ))),
  rep(7, length(c(
    labels(tree_src.mg2.tracing.rowtree.1[[1]][[2]][[2]][[2]])
  ))))

tree_src.mg2.tracing.coltree.1 = as.dendrogram(src.mg2.tracing.coltree.1)
src.mg2.tracing$tree.1 = NA
src.mg2.tracing@meta.data[labels(tree_src.mg2.tracing.coltree.1[[1]][[1]][[1]]),]$tree.1 = 1
src.mg2.tracing@meta.data[labels(tree_src.mg2.tracing.coltree.1[[1]][[1]][[2]]),]$tree.1 = 2
src.mg2.tracing@meta.data[labels(tree_src.mg2.tracing.coltree.1[[1]][[2]][[1]]),]$tree.1 = 3
src.mg2.tracing@meta.data[labels(tree_src.mg2.tracing.coltree.1[[1]][[2]][[2]]),]$tree.1 = 4
src.mg2.tracing@meta.data[labels(tree_src.mg2.tracing.coltree.1[[2]][[1]][[1]]),]$tree.1 = 5
src.mg2.tracing@meta.data[labels(tree_src.mg2.tracing.coltree.1[[2]][[1]][[2]]),]$tree.1 = 6
src.mg2.tracing@meta.data[labels(tree_src.mg2.tracing.coltree.1[[2]][[2]][[1]]),]$tree.1 = 7
src.mg2.tracing@meta.data[labels(tree_src.mg2.tracing.coltree.1[[2]][[2]][[2]]),]$tree.1 = 8

DimPlot(src.mg2.tracing, reduction = "mnn_umap_fta", group.by = "tree.1", label = T, label.size = 5)

src.mg2.tracing$cluster.v06.26.re_correct = src.mg2.tracing$cluster.v06.26.re_mnn_umap_fta
src.mg2.tracing@meta.data[src.mg2.tracing$tree.1%in%c(1,2,3,4),]$cluster.v06.26.re_correct = "MG.2"
src.mg2.tracing@meta.data[src.mg2.tracing$tree.1%in%c(7,8),]$cluster.v06.26.re_correct = "Small.intestine.1"
src.mg2.tracing@meta.data[src.mg2.tracing$tree.1%in%c(5,6),]$cluster.v06.26.re_correct = "Small.intestine.2"
DimPlot(src.mg2.tracing, reduction = "mnn_umap_fta", group.by = "Time")

#-- Pseudo time
cell_src.mg2.tracing_sm1 = rownames(src.mg2.tracing@meta.data[src.mg2.tracing$cluster.v06.26.re_correct%in%c("MG.2","Small.intestine.1"),])
coord_src.mg2.tracing_sm1 = src.mg2.tracing[["umap"]]@cell.embeddings[cell_src.mg2.tracing_sm1, c(1,2)]
pcurve_src.mg2.tracing_sm1 = princurve::principal_curve(x = coord_src.mg2.tracing_sm1, smoother = "smooth.spline")
src.mg2.tracing$lambda_sm1 = NA
src.mg2.tracing@meta.data[cell_src.mg2.tracing_sm1,]$lambda_sm1 = pcurve_src.mg2.tracing_sm1$lambda[cell_src.mg2.tracing_sm1]
src.mg2.tracing@meta.data[!src.mg2.tracing$lambda_sm1%in%NA,]$lambda_sm1 = 
  norm_range(src.mg2.tracing@meta.data[!src.mg2.tracing$lambda_sm1%in%NA,]$lambda_sm1)

cell_src.mg2.tracing_sm2 = rownames(src.mg2.tracing@meta.data[src.mg2.tracing$cluster.v06.26.re_correct%in%c("MG.2","Small.intestine.2"),])
coord_src.mg2.tracing_sm2 = src.mg2.tracing[["umap"]]@cell.embeddings[cell_src.mg2.tracing_sm2, c(1,2)]
pcurve_src.mg2.tracing_sm2 = princurve::principal_curve(x = coord_src.mg2.tracing_sm2, smoother = "smooth.spline")
src.mg2.tracing$lambda_sm2 = NA
src.mg2.tracing@meta.data[cell_src.mg2.tracing_sm2,]$lambda_sm2 = pcurve_src.mg2.tracing_sm2$lambda[cell_src.mg2.tracing_sm2]
src.mg2.tracing@meta.data[!src.mg2.tracing$lambda_sm2%in%NA,]$lambda_sm2 = 
  norm_range(src.mg2.tracing@meta.data[!src.mg2.tracing$lambda_sm2%in%NA,]$lambda_sm2)


cell_src.mg2.tracing_mg2 = rownames(src.mg2.tracing@meta.data[src.mg2.tracing$cluster.v06.26.re_correct%in%c("MG.2"),])
w1 = length(cell_src.mg2.tracing_sm1) / (length(cell_src.mg2.tracing_sm1) + length(cell_src.mg2.tracing_sm2))
w2 = length(cell_src.mg2.tracing_sm1) / (length(cell_src.mg2.tracing_sm1) + length(cell_src.mg2.tracing_sm2))
src.mg2.tracing$lambda_mg2 = NA
src.mg2.tracing@meta.data[cell_src.mg2.tracing_mg2, ]$lambda_mg2 = 
  w1*src.mg2.tracing@meta.data[cell_src.mg2.tracing_mg2, ]$lambda_sm1 + 
  w2*src.mg2.tracing@meta.data[cell_src.mg2.tracing_mg2, ]$lambda_sm2
src.mg2.tracing@meta.data[!src.mg2.tracing$lambda_mg2%in%NA,]$lambda_sm2 = 
  norm_range(src.mg2.tracing@meta.data[!src.mg2.tracing$lambda_mg2%in%NA,]$lambda_sm2)


cellorder_cell_src.mg2.tracing = c(
  (intersect(names(src.mg2.tracing$lambda_mg2[order(src.mg2.tracing$lambda_mg2)]), 
                colnames(src.mg2.tracing[,src.mg2.tracing$cluster.v06.26.re_correct%in%c("MG.2"),]))),
  (intersect(names(src.mg2.tracing$lambda_sm1[order(src.mg2.tracing$lambda_sm1)]), 
                colnames(src.mg2.tracing[,src.mg2.tracing$cluster.v06.26.re_correct%in%c("Small.intestine.1"),]))),
  rev(intersect(names(src.mg2.tracing$lambda_sm2[order(src.mg2.tracing$lambda_sm2)]), 
                colnames(src.mg2.tracing[,src.mg2.tracing$cluster.v06.26.re_correct%in%c("Small.intestine.2"),]))))

gene_src.mg2.tracing.rowtree.2 = c(
  gene_src.mg2.tracing.rowtree.1[names(gene_src.mg2.tracing.rowtree.1)%in%c(4)],
  gene_src.mg2.tracing.rowtree.1[names(gene_src.mg2.tracing.rowtree.1)%in%c(3)],
  rev(gene_src.mg2.tracing.rowtree.1[names(gene_src.mg2.tracing.rowtree.1)%in%c(7)]))

pdf("figure.v08.07/organ_development_re_v240115/Heatmap_for_trajectory/try.mg2.pdf", 10, 10)
src.mg2.tracing.coltree.2 =
  MyHeatmap(as.matrix(src.mg2.tracing@assays$RNA@data[
    # gene_src.mg2.tracing.rowtree.1,
    gene_src.mg2.tracing.rowtree.2,
    cellorder_cell_src.mg2.tracing]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.mg2.tracing$Time[cellorder_cell_src.mg2.tracing],
                 colors.time.2),
      MyName2Col(src.mg2.tracing$lineage[cellorder_cell_src.mg2.tracing],
                 color.lineage),
      MyName2Col(src.mg2.tracing$cluster.v06.26.re_correct[cellorder_cell_src.mg2.tracing],
                 cluster.endoderm.color.v5)
    ),
    RowSideColors = t(cbind(
      MyName2Col(names(
        gene_src.mg2.tracing.rowtree.2),
        # gene_src.mg2.tracing.rowtree.1),
                 colors.geneset)
    )),
    ColSideColorsSize = 4.8,
    RowSideColorsSize = 1.2,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    Rowv = "none", Colv = "none",
    margins = c(10, 10),
    #return.tree = "col",
    graph = T)


cellorder_cell_src.mg2.tracing = 
  intersect(cellorder_cell_src.mg2.tracing,
            rownames(src.mg2.tracing@meta.data[src.mg2.tracing$lineage%in%"Nepn",]))
src.mg2.tracing.coltree.2 =
  MyHeatmap(as.matrix(src.mg2.tracing@assays$RNA@data[
    gene_src.mg2.tracing.rowtree.2,
    # gene_src.mg2.tracing.rowtree.1,
    cellorder_cell_src.mg2.tracing]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.mg2.tracing$Time[cellorder_cell_src.mg2.tracing],
                 colors.time.2),
      MyName2Col(src.mg2.tracing$lineage[cellorder_cell_src.mg2.tracing],
                 color.lineage),
      MyName2Col(src.mg2.tracing$cluster.v06.26.re_correct[cellorder_cell_src.mg2.tracing],
                 cluster.endoderm.color.v5)
    ),
    RowSideColors = t(cbind(
      MyName2Col(names(
        gene_src.mg2.tracing.rowtree.2),
        # gene_src.mg2.tracing.rowtree.1),
        colors.geneset)
    )),
    ColSideColorsSize = 4.8,
    RowSideColorsSize = 1.2,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    Rowv = "none", Colv = "none",
    margins = c(10,10),
    #return.tree = "col",
    graph = T)
dev.off()
#=======================

save(src.mg2.tracing, file = "figure.v08.07/organ_development_re_v240115/src.mg2.tracing.Rdata")
save(src.mg2.tracing.filtergene,
     marker_src.mg2.tracing,
     markergene_src.mg2.tracing,
     tree_src.mg2.tracing.rowtree,
     tree_src.mg2.tracing.rowtree.1,
     tree_src.mg2.tracing.coltree.1,
     
     gene_src.mg2.tracing.rowtree,
     gene_src.mg2.tracing.rowtree.1,
     cellorder_cell_src.mg2.tracing,
     
     file = "figure.v08.07/organ_development_re_v240115/src.mg2.tracing.parameter.Rdata")
     



#---------------------------------------------------------------------------
# --          MG.1 -- Re
#---------------------------------------------------------------------------
DimPlot(src.9ss.integrated.merge, group.by = "lineage", reduction = "umap_fta")
DimPlot(src.9ss.integrated.merge, group.by = "cluster.v06.26.re..merge_umap_fta", reduction = "umap_fta")

src.mg1.tracing@meta.data[
  intersect(rownames(src.9ss.integrated.merge@meta.data[src.9ss.integrated.merge$cluster.predict.umap_int.ext.v1.1%in%"MG.1",]),
            rownames(src.mg1.tracing@meta.data)),]$cluster.v06.26.re_mnn_umap_fta = "MG.1"

#------------------------
src.mg1.tracing = FindVariableFeatures(src.mg1.tracing, nfeatures = 2000)
src.mg1.tracing.filtergene = 
  Myfilter(as.matrix(src.mg1.tracing@assays$RNA@data),
           gene = src.mg1.tracing@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.mg1.tracing.filtergene.re = src.mg1.tracing.filtergene; rm(src.mg1.tracing.filtergene)
src.mg1.tracing.filtergene = src.mg1.tracing.filtergene.re; rm(src.mg1.tracing.filtergene.re)


src.mg1.tracing = SetIdent(src.mg1.tracing, value = src.mg1.tracing$cluster.v06.26.re_mnn_umap_fta)
marker_src.mg1.tracing = FindAllMarkers(src.mg1.tracing)
marker_src.mg1.tracing$pct.ratio = marker_src.mg1.tracing$pct.1 / marker_src.mg1.tracing$pct.2
marker_src.mg1.tracing$rank = marker_src.mg1.tracing$pct.ratio * (-log(marker_src.mg1.tracing$p_val_adj))
marker_src.mg1.tracing = marker_src.mg1.tracing[order(marker_src.mg1.tracing$rank, decreasing = T),]
markergene_src.mg1.tracing = unique(marker_src.mg1.tracing$gene)
# markergene_src.mg1.tracing.raw = markergene_src.mg1.tracing

src.mg1.tracing = ScaleData(src.mg1.tracing, 
                            rownames(src.mg1.tracing),split.by = "Phase")
src.mg1.tracing = RunPCA(src.mg1.tracing, 
                         features = src.mg1.tracing.filtergene)
src.mg1.tracing = RunUMAP(src.mg1.tracing, reduction = "pca", 
                          dims = 1:30, n.neighbors = 150,
                          assay = "RNA", n.components = 3)
DimPlot(src.mg1.tracing, reduction = "umap", group.by = "Time")


# == MG.1 Basic process 
#===============================================================================
pdf("figure.v08.07/organ_development_re_v240115/try.mg1.pdf",9,7)
src.mg1.tracing.rowtree =
  MyHeatmap(as.matrix(src.mg1.tracing@assays$RNA@data[
    unique(c(markergene_src.mg1.tracing)
             #src.mg1.tracing.filtergene,
           #c(src.endoderm.mg2.ext_v1.1.filtergene_fin)
    ),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.mg1.tracing$Time,
                 colors.time.2),
      MyName2Col(src.mg1.tracing$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.mg1.tracing$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.mg1.tracing$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none", Colv = "none",
    return.tree = "row",
    graph = T)

src.mg1.tracing.coltree =
  MyHeatmap(as.matrix(src.mg1.tracing@assays$RNA@data[
    unique(c(markergene_src.mg1.tracing)
           #src.mg1.tracing.filtergene,
           #c(src.endoderm.mg2.ext_v1.1.filtergene_fin)
    ),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.mg1.tracing$Time,
                 colors.time.2),
      MyName2Col(src.mg1.tracing$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.mg1.tracing$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.mg1.tracing$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none", Colv = "none",
    return.tree = "col",
    graph = T)
dev.off()

pdf("figure.v08.07/organ_development_re_v240115/try.mg1.tree.pdf",150,20)
plot(src.mg1.tracing.rowtree)
dev.off()
tree_src.mg1.tracing.rowtree = as.dendrogram(src.mg1.tracing.rowtree)
gene_src.mg1.tracing.rowtree = c(
  labels(tree_src.mg1.tracing.rowtree[[2]][[2]][[2]][[2]][[2]][[1]]),
  labels(tree_src.mg1.tracing.rowtree[[2]][[2]][[2]][[2]][[2]][[2]]),
  labels(tree_src.mg1.tracing.rowtree[[2]][[2]][[1]]),
  labels(tree_src.mg1.tracing.rowtree[[2]][[2]][[2]][[1]]))

tree_src.mg1.tracing.coltree = as.dendrogram(src.mg1.tracing.coltree)
src.mg1.tracing$tree = NA
src.mg1.tracing@meta.data[labels(tree_src.mg1.tracing.coltree[[1]][[1]][[1]]),]$tree = 1
src.mg1.tracing@meta.data[labels(tree_src.mg1.tracing.coltree[[1]][[1]][[2]]),]$tree = 2
src.mg1.tracing@meta.data[labels(tree_src.mg1.tracing.coltree[[1]][[2]][[1]]),]$tree = 3
src.mg1.tracing@meta.data[labels(tree_src.mg1.tracing.coltree[[1]][[2]][[2]]),]$tree = 4
src.mg1.tracing@meta.data[labels(tree_src.mg1.tracing.coltree[[2]][[1]][[1]]),]$tree = 5
src.mg1.tracing@meta.data[labels(tree_src.mg1.tracing.coltree[[2]][[1]][[2]]),]$tree = 6
src.mg1.tracing@meta.data[labels(tree_src.mg1.tracing.coltree[[2]][[2]][[1]][[1]]),]$tree = 7.1
src.mg1.tracing@meta.data[labels(tree_src.mg1.tracing.coltree[[2]][[2]][[1]][[2]]),]$tree = 7.2
src.mg1.tracing@meta.data[labels(tree_src.mg1.tracing.coltree[[2]][[2]][[2]][[1]]),]$tree = 8.1
src.mg1.tracing@meta.data[labels(tree_src.mg1.tracing.coltree[[2]][[2]][[2]][[2]]),]$tree = 8.2
DimPlot(src.mg1.integrated, reduction = "umap_mnn", group.by = "Time", label = T, label.size = 5)

src.mg1.integrated = NormalizeData(src.mg1.integrated, scale.factor = 10*5)
anchor.integrated = 
  FindTransferAnchors(reference = src.mg1.integrated,
                      query = src.mg1.tracing,
                      reference.assay =  "mnnRNA",
                      query.assay =  'mnnRNA', 
                      scale = T, dims = 1:30, 
                      features = src.endoderm.mg1.ext_v1.1.gene.fin)
umap_embedding = src.mg1.integrated[["umap_mnn"]]@cell.embeddings
umap.transfer = TransferData(anchor.integrated, t(umap_embedding), k.weight = 5)
src.mg1.tracing@reductions$mnn_umap_fta@cell.embeddings = as.matrix(t(umap.transfer@data))
src.mg1.tracing@reductions$mnn_umap_fta@key = "UMAP_"
colnames(src.mg1.tracing@reductions$mnn_umap_fta@cell.embeddings) = paste("UMAP_",c(1:3),sep='')

DimPlot(src.mg1.tracing, reduction = "mnn_umap_fta", group.by = "Time", label = T, label.size = 5) +
  DimPlot(src.mg1.tracing, reduction = "mnn_umap_fta", group.by = "tree", label = T, label.size = 5) +
  DimPlot(src.mg1.tracing, reduction = "mnn_umap_fta", group.by = "cluster.v06.26.re_mnn_umap_fta", label = T, label.size = 5)

src.mg1.tracing$cluster.v06.26.re_correct = src.mg1.tracing$cluster.v06.26.re_mnn_umap_fta
src.mg1.tracing@meta.data[src.mg1.tracing$tree%in%c(3,4),]$cluster.v06.26.re_correct = "MG.1"
src.mg1.tracing@meta.data[src.mg1.tracing$tree%in%c(1,2,8.1,8.2),]$cluster.v06.26.re_correct = "Stomach"
src.mg1.tracing@meta.data[src.mg1.tracing$tree%in%c(5,6,7.1,7.2),]$cluster.v06.26.re_correct = "Small.intestine.1"
src.mg1.tracing@meta.data[src.mg1.tracing$Time%in%c("9ss","12ss","15ss"),]$cluster.v06.26.re_correct = "MG.1"
src.mg1.tracing@meta.data[src.mg1.tracing$cluster.v06.26.re_mnn_umap_fta%in%c("Small.intestine.1") &
                            src.mg1.tracing$cluster.v06.26.re_correct%in%c("Stomach"),]$cluster.v06.26.re_correct = "Small.intestine.1"
DimPlot(src.mg1.tracing, reduction = "mnn_umap_fta", group.by = "cluster.v06.26.re_correct", label = T, label.size = 5) +
  DimPlot(src.mg1.tracing, reduction = "mnn_umap_fta", group.by = "cluster.v06.26.re_mnn_umap_fta")


src.mg1.tracing = SetIdent(src.mg1.tracing, 
                           value = src.mg1.tracing$cluster.v06.26.re_correct)
marker_src.mg1.tracing = FindAllMarkers(src.mg1.tracing)
marker_src.mg1.tracing$pct.ratio = marker_src.mg1.tracing$pct.1 / marker_src.mg1.tracing$pct.2
marker_src.mg1.tracing$rank = marker_src.mg1.tracing$pct.ratio * (-log(marker_src.mg1.tracing$p_val_adj))
marker_src.mg1.tracing = marker_src.mg1.tracing[order(marker_src.mg1.tracing$rank, decreasing = T),]
markergene_src.mg1.tracing.1 = unique(marker_src.mg1.tracing$gene)
# markergene_src.mg1.tracing.raw = markergene_src.mg1.tracing


pdf("figure.v08.07/organ_development_re_v240115/try.mg1.pdf",9,7)
src.mg1.tracing.rowtree.1 =
  MyHeatmap(as.matrix(src.mg1.tracing@assays$RNA@data[
    unique(c(markergene_src.mg1.tracing.1,
             src.mg1.tracing.filtergene)
           #src.mg1.tracing.filtergene,
           #c(src.endoderm.mg2.ext_v1.1.filtergene_fin)
    ),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.mg1.tracing$Time,
                 colors.time.2),
      MyName2Col(src.mg1.tracing$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.mg1.tracing$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.mg1.tracing$cluster.v06.26.re_correct,
                 cluster.endoderm.color.v5),
      MyName2Col(src.mg1.tracing$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none", Colv = "none",
    return.tree = "row",
    graph = T)
src.mg1.tracing.coltree.1 =
  MyHeatmap(as.matrix(src.mg1.tracing@assays$RNA@data[
    unique(c(markergene_src.mg1.tracing.1,
             src.mg1.tracing.filtergene)
           #src.mg1.tracing.filtergene,
           #c(src.endoderm.mg2.ext_v1.1.filtergene_fin)
    ),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.mg1.tracing$Time,
                 colors.time.2),
      # MyName2Col(src.mg1.tracing$tree,
      #            color.temp),
      MyName2Col(src.mg1.tracing$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.mg1.tracing$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.mg1.tracing$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none", Colv = "none",
    return.tree = "col",
    graph = T)
dev.off()

pdf("figure.v08.07/organ_development_re_v240115/try.mg1.tree.pdf",150,20)
plot(src.mg1.tracing.rowtree.1)
dev.off()

tree_src.mg1.tracing.rowtree.1 = as.dendrogram(src.mg1.tracing.rowtree.1)
gene_src.mg1.tracing.rowtree.1 = c(
  labels(tree_src.mg1.tracing.rowtree.1[[1]][[2]][[1]][[2]]), # MG.1
  
  labels(tree_src.mg1.tracing.rowtree.1[[1]][[1]][[2]][[2]]),
  labels(tree_src.mg1.tracing.rowtree.1[[1]][[1]][[1]]),
  labels(tree_src.mg1.tracing.rowtree.1[[1]][[1]][[2]][[1]]), # Stomach,
  
  labels(tree_src.mg1.tracing.rowtree.1[[2]][[1]]), # Stomach & Sm1
  
  rev(c(
    labels(tree_src.mg1.tracing.rowtree.1[[2]][[2]][[1]][[2]]), # Sm1,
    rev(labels(tree_src.mg1.tracing.rowtree.1[[2]][[2]][[2]][[2]][[2]][[2]])) # Sm1
  ))
)

names(gene_src.mg1.tracing.rowtree.1) = c(
  rep(4, length(labels(tree_src.mg1.tracing.rowtree.1[[1]][[2]][[1]][[2]]))),
  
  rep(3, length(c(
    labels(tree_src.mg1.tracing.rowtree.1[[1]][[1]][[2]][[2]]),
    labels(tree_src.mg1.tracing.rowtree.1[[1]][[1]][[1]]),
    labels(tree_src.mg1.tracing.rowtree.1[[1]][[1]][[2]][[1]])
  ))),
  
  rep(5, length(labels(tree_src.mg1.tracing.rowtree.1[[2]][[1]]))),
  
  rep(7, length(c(
    labels(tree_src.mg1.tracing.rowtree.1[[2]][[2]][[1]][[2]]), # Sm1,
    rev(labels(tree_src.mg1.tracing.rowtree.1[[2]][[2]][[2]][[2]][[2]][[2]])) 
  )))
)

DimPlot(src.mg1.tracing)
DimPlot(src.mg1.tracing, group.by = "Time", reduction = 'mnn_umap_fta')

# -- Pseudo time
#-------------------------------------------------------------------------------
cell_src.mg1.tracing_sto = rownames(src.mg1.tracing@meta.data[src.mg1.tracing$cluster.v06.26.re_correct%in%c("MG.1","Stomach"),])
coord_src.mg1.tracing_sto = src.mg1.tracing[["mnn_umap_fta"]]@cell.embeddings[cell_src.mg1.tracing_sto, c(2,2)]
pcurve_src.mg1.tracing_sto = princurve::principal_curve(x = coord_src.mg1.tracing_sto, smoother = "smooth.spline")
src.mg1.tracing$lambda_sto = NA
src.mg1.tracing@meta.data[cell_src.mg1.tracing_sto,]$lambda_sto = pcurve_src.mg1.tracing_sto$lambda[cell_src.mg1.tracing_sto]
src.mg1.tracing@meta.data[!src.mg1.tracing$lambda_sto%in%NA,]$lambda_sto = 
  norm_range(src.mg1.tracing@meta.data[!src.mg1.tracing$lambda_sto%in%NA,]$lambda_sto)

cell_src.mg1.tracing_sm1 = rownames(src.mg1.tracing@meta.data[src.mg1.tracing$cluster.v06.26.re_correct%in%c("MG.1","Small.intestine.1"),])
coord_src.mg1.tracing_sm1 = src.mg1.tracing[["mnn_umap_fta"]]@cell.embeddings[cell_src.mg1.tracing_sm1, c(2,1)]
pcurve_src.mg1.tracing_sm1 = princurve::principal_curve(x = coord_src.mg1.tracing_sm1, smoother = "smooth.spline")
src.mg1.tracing$lambda_sm1 = NA
src.mg1.tracing@meta.data[cell_src.mg1.tracing_sm1,]$lambda_sm1 = pcurve_src.mg1.tracing_sm1$lambda[cell_src.mg1.tracing_sm1]
src.mg1.tracing@meta.data[!src.mg1.tracing$lambda_sm1%in%NA,]$lambda_sm1 = 
  norm_range(src.mg1.tracing@meta.data[!src.mg1.tracing$lambda_sm1%in%NA,]$lambda_sm1)

cell_src.mg1.tracing_mg1 = rownames(src.mg1.tracing@meta.data[src.mg1.tracing$cluster.v06.26.re_correct%in%c("MG.1","Stomach"),])
coord_src.mg1.tracing_mg1 = src.mg1.tracing[["mnn_umap_fta"]]@cell.embeddings[cell_src.mg1.tracing_mg1, c(1,2)]
pcurve_src.mg1.tracing_mg1 = princurve::principal_curve(x = coord_src.mg1.tracing_mg1, smoother = "smooth.spline")
src.mg1.tracing$lambda_mg1 = NA
src.mg1.tracing@meta.data[cell_src.mg1.tracing_mg1,]$lambda_mg1 = pcurve_src.mg1.tracing_mg1$lambda[cell_src.mg1.tracing_mg1]
src.mg1.tracing@meta.data[!src.mg1.tracing$lambda_mg1%in%NA,]$lambda_mg1 = 
  norm_range(src.mg1.tracing@meta.data[!src.mg1.tracing$lambda_mg1%in%NA,]$lambda_mg1)
#-------------------------------------------------------------------------------

cellorder_cell_src.mg1.tracing = c(
  (intersect(names(src.mg1.tracing$lambda_mg1[order(src.mg1.tracing$lambda_mg1)]), 
             colnames(src.mg1.tracing[,src.mg1.tracing$cluster.v06.26.re_correct%in%c("MG.1"),]))),
  (intersect(names(src.mg1.tracing$lambda_sm1[order(src.mg1.tracing$lambda_sm1)]), 
             colnames(src.mg1.tracing[,src.mg1.tracing$cluster.v06.26.re_correct%in%c("Small.intestine.1"),]))),
  (intersect(names(src.mg1.tracing$lambda_sto[order(src.mg1.tracing$lambda_sto)]), 
             colnames(src.mg1.tracing[,src.mg1.tracing$cluster.v06.26.re_correct%in%c("Stomach"),]))))


pdf("figure.v08.07/organ_development_re_v240115/try.mg1.pdf",10,10)
src.mg1.tracing.rowtree.2.fin =
  MyHeatmap(as.matrix(src.mg1.tracing@assays$RNA@data[
    gene_src.mg1.tracing.rowtree.1,
    # gene_src.mg1.tracing.rowtree.1,
    cellorder_cell_src.mg1.tracing]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.mg1.tracing$Time[cellorder_cell_src.mg1.tracing],
                 colors.time.2),
      MyName2Col(src.mg1.tracing$lineage[cellorder_cell_src.mg1.tracing],
                 color.lineage),
      MyName2Col(src.mg1.tracing$cluster.v06.26.re_correct[cellorder_cell_src.mg1.tracing],
                 cluster.endoderm.color.v5)
    ),
    RowSideColors = t(cbind(MyName2Col(
      names(gene_src.mg1.tracing.rowtree.1), colors.geneset))),
    ColSideColorsSize = 4.8,
    RowSideColorsSize = 1.2,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    Rowv = "none", Colv = "none",
    margins = c(10,10),
    # return.tree = "row",
    graph = T)
dev.off()

# save(src.mg1.tracing, )
DimPlot(src.mg1.tracing, group.by = "cluster.v06.26.re_correct", reduction = "mnn_umap_fta")

save(src.mg1.tracing, file = 'figure.v08.07/organ_development_re_v240115/src.mg1.tracing.Rdata')
save(tree_src.mg1.tracing.rowtree,
     tree_src.mg1.tracing.rowtree.1,
     tree_src.mg1.tracing.coltree,
     gene_src.mg1.tracing.rowtree,
     gene_src.mg1.tracing.rowtree.1,
     src.mg1.tracing.filtergene, 
     markergene_src.mg1.tracing,
     markergene_src.mg1.tracing.1,
     cellorder_cell_src.mg1.tracing,
     file = 'figure.v08.07/organ_development_re_v240115/src.mg1.tracing.parameter.Rdata')


src.mg1.tracing.re = src.mg1.tracing[,src.mg1.tracing$lineage%in%"Sox2"]
save(src.mg1.tracing.re, file = 'figure.v08.07/organ_development_re_v240115/src.mg1.tracing.re.Rdata')
cellorder_cell_src.mg1.tracing.re = intersect(cellorder_cell_src.mg1.tracing, colnames(src.mg1.tracing.re))
DimPlot(src.mg1.tracing.re, reduction = "mnn_umap_fta")

pdf("figure.v08.07/organ_development_re_v240115/Heatmap_for_trajectory/try.mg1.re.pdf",10,10)
src.mg1.tracing.re.rowtree =
  MyHeatmap(as.matrix(src.mg1.tracing.re@assays$RNA@data[
    gene_src.mg1.tracing.rowtree.1,
    # gene_src.mg1.tracing.re.rowtree.1,
    cellorder_cell_src.mg1.tracing.re]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.mg1.tracing.re$Time[cellorder_cell_src.mg1.tracing.re],
                 colors.time.2),
      MyName2Col(src.mg1.tracing.re$lineage[cellorder_cell_src.mg1.tracing.re],
                 color.lineage),
      MyName2Col(src.mg1.tracing.re$cluster.v06.26.re_correct[cellorder_cell_src.mg1.tracing.re],
                 cluster.endoderm.color.v5)
    ),
    RowSideColors = t(cbind(MyName2Col(
      names(gene_src.mg1.tracing.rowtree.1), colors.geneset))),
    ColSideColorsSize = 4,8,
    RowSideColorsSize = 1,2,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    Rowv = "none", Colv = "none",
    # return.tree = "row",
    margins = c(10,10),
    graph = T)
dev.off()

save(gene_src.mg1.tracing.rowtree.1,
     cellorder_cell_src.mg1.tracing.re,
     cellorder_cell_src.mg1.tracing,
     src.mg1.tracing, src.mg1.tracing.re,
     file = "figure.v08.07/organ_development_re_v240115/src.mg1.tracing.re.parameter.v0508.Rdata")




