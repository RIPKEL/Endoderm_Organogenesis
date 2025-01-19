#-----------
# HG.2-Re
#-----------
src.hg2.tracing = src.hg2.integrated.merge[,!src.hg2.integrated.merge$lineage%in%NA]
src.hg2.tracing$cluster.v06.26.re_mnn_umap_fta = 
  src.hg2.integrated.merge$cluster.v06.26.re_mnn_umap_fta[colnames(src.hg2.tracing)]
src.hg2.tracing@meta.data[
  intersect(rownames(src.9ss.integrated.merge@meta.data[
    src.9ss.integrated.merge$cluster.predict.umap_int.ext.v1.1%in%"HG.2",]),
    rownames(src.hg2.tracing@meta.data)),]$cluster.v06.26.re_mnn_umap_fta = "HG.2"

src.hg2.tracing = FindVariableFeatures(src.hg2.tracing, nfeatures = 2000)
src.hg2.tracing.filtergene = 
  Myfilter(as.matrix(src.hg2.tracing@assays$RNA@data),
           gene = src.hg2.tracing@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.hg2.tracing.filtergene.re = src.hg2.tracing.filtergene; rm(src.hg2.tracing.filtergene)
src.hg2.tracing.filtergene = src.hg2.tracing.filtergene.re; rm(src.hg2.tracing.filtergene.re)

src.hg2.tracing = SetIdent(src.hg2.tracing, value = src.hg2.tracing$cluster.v06.26.re_mnn_umap_fta)
marker_src.hg2.tracing = FindAllMarkers(src.hg2.tracing)
marker_src.hg2.tracing$pct.ratio = marker_src.hg2.tracing$pct.1 / marker_src.hg2.tracing$pct.2
marker_src.hg2.tracing$rank = marker_src.hg2.tracing$pct.ratio * (-log(marker_src.hg2.tracing$p_val_adj))
marker_src.hg2.tracing = marker_src.hg2.tracing[order(marker_src.hg2.tracing$rank, decreasing = T),]
markergene_src.hg2.tracing = unique(marker_src.hg2.tracing$gene)

src.hg2.tracing = RunPCA(src.hg2.tracing, features = src.hg2.tracing.filtergene)
src.hg2.tracing = RunUMAP(src.hg2.tracing, dims = 1:30, reduction = 'pca', 
                          n.neighbors = 100, n.components = 2)

src.hg2.tracing[["mnn_umap_fta"]] = src.hg2.tracing[["umap"]]
src.hg2.tracing@reductions$mnn_umap_fta@key = "Coord_"
src.hg2.tracing@reductions$mnn_umap_fta@cell.embeddings = 
  src.hg2.integrated.merge@reductions$mnn_umap_fta@cell.embeddings[colnames(src.hg2.tracing), c(1:2)]
colnames(src.hg2.tracing@reductions$mnn_umap_fta@cell.embeddings) = c("Coord_1","Coord_2")


pdf("figure.v08.07/organ_development_re_v240115/try.hg2.re.pdf",10,10)
#--------------------
src.hg2.tracing.rowtree  =  MyHeatmap(as.matrix(
    src.hg2.tracing@assays$RNA@data[
      unique(c(markergene_src.hg2.tracing,
               src.hg2.tracing.filtergene))
      ,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.hg2.tracing$Time, #[cellorder_cell_src.hg2.tracing],
                 colors.time.2),
      MyName2Col(src.hg2.tracing$lineage, #[cellorder_cell_src.hg2.tracing],
                 color.lineage),
      MyName2Col(src.hg2.tracing$cluster.v06.26.re_mnn_umap_fta, #[cellorder_cell_src.hg2.tracing],
                 cluster.endoderm.color.v5)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    margins = c(10,10),
    graph = T)

src.hg2.tracing.coltree  =  MyHeatmap(as.matrix(
  src.hg2.tracing@assays$RNA@data[
    unique(c(markergene_src.hg2.tracing,
             src.hg2.tracing.filtergene))
    ,]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.hg2.tracing$Time, #[cellorder_cell_src.hg2.tracing],
               colors.time.2),
    MyName2Col(src.hg2.tracing$lineage, #[cellorder_cell_src.hg2.tracing],
               color.lineage),
    MyName2Col(src.hg2.tracing$cluster.v06.26.re_mnn_umap_fta, #[cellorder_cell_src.hg2.tracing],
               cluster.endoderm.color.v5)
  ),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  #Rowv = "none",
  return.tree = "col",
  margins = c(10,10),
  graph = T)
#--------------------
dev.off()

tree_src.hg2.tracing.rowtree = as.dendrogram(src.hg2.tracing.rowtree)
gene_src.hg2.tracing.rowtree = c(
  labels(tree_src.hg2.tracing.rowtree[[1]][[2]][[1]][[1]]),
  labels(tree_src.hg2.tracing.rowtree[[1]][[2]][[1]][[2]][[2]]),
  labels(tree_src.hg2.tracing.rowtree[[1]][[2]][[2]][[1]]),
  labels(tree_src.hg2.tracing.rowtree[[1]][[2]][[2]][[2]]),
  
  rev(c(
    labels(tree_src.hg2.tracing.rowtree[[1]][[1]])
  )),
  
  rev(c(
    labels(tree_src.hg2.tracing.rowtree[[2]][[1]][[1]]),
    labels(tree_src.hg2.tracing.rowtree[[2]][[1]][[2]][[1]][[1]]),
    labels(tree_src.hg2.tracing.rowtree[[2]][[1]][[2]][[2]])
  ))
)
names(gene_src.hg2.tracing.rowtree) = c(
  rep(4, length(c(
    labels(tree_src.hg2.tracing.rowtree[[1]][[2]][[1]][[1]]),
    labels(tree_src.hg2.tracing.rowtree[[1]][[2]][[1]][[2]][[2]]),
    labels(tree_src.hg2.tracing.rowtree[[1]][[2]][[2]][[1]]),
    labels(tree_src.hg2.tracing.rowtree[[1]][[2]][[2]][[2]])
  ))),
  rep(3, length(c(
    labels(tree_src.hg2.tracing.rowtree[[1]][[1]])
  ))),
  rep(7, length(c(
    rev(labels(tree_src.hg2.tracing.rowtree[[2]][[1]][[1]])),
    labels(tree_src.hg2.tracing.rowtree[[2]][[1]][[2]][[1]][[1]]),
    labels(tree_src.hg2.tracing.rowtree[[2]][[1]][[2]][[2]])
  )))
)



pdf("figure.v08.07/organ_development_re_v240115/try.hg2.re.pdf",10,10)
src.hg2.tracing.coltree1  =  MyHeatmap(as.matrix(
  src.hg2.tracing@assays$RNA@data[
    gene_src.hg2.tracing.rowtree,
    cellorder_src.hg2.tracing]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.hg2.tracing$Time[cellorder_src.hg2.tracing],
               colors.time.2),
    MyName2Col(src.hg2.tracing$lineage[cellorder_src.hg2.tracing],
               color.lineage),
    MyName2Col(src.hg2.tracing$cluster.v06.26.re_correct[cellorder_src.hg2.tracing],
               cluster.endoderm.color.v5)
  ),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  # Rowv = "none", Colv = "none",
  return.tree = "col",
  margins = c(10,10),
  graph = T)
dev.off()


tree_src.hg2.tracing.coltree1 = as.dendrogram(src.hg2.tracing.coltree1)
src.hg2.tracing$tree = NA
for(i1 in c(2)){for(i2 in c(1,2)){
  for(i3 in c(1,2)){for(i4 in c(1,2)){for(i5 in c(1,2)){
    src.hg2.tracing@meta.data[
      labels(tree_src.hg2.tracing.coltree1[[i1]][[i2]][[i3]][[i4]][[i5]]),]$tree = paste(i1,i2,i3,i4,i5,sep=".")
  }}}}
}

DimPlot(src.hg2.tracing, group.by = 'tree', reduction = 'umap', label = T) +
  DimPlot(src.hg2.tracing, group.by = "cluster.v06.26.re_mnn_umap_fta", reduction = "umap", label = T) +
  DimPlot(src.hg2.tracing, group.by = "Time", reduction = "umap")

src.hg2.tracing$cluster.v06.26.re_correct = src.hg2.tracing$cluster.v06.26.re_mnn_umap_fta
src.hg2.tracing@meta.data[
  src.hg2.tracing$tree%in%c("2.1.1.1.1", '2.1.1.1.2', "2.1.1.2.1", "2.1.1.2.2",
                            "2.1.2.2.1", "2.1.2.2.2",
                            '2.2.1.1.1', '2.2.1.1.2', "2.2.1.2.1", '2.2.1.2.2'),]$cluster.v06.26.re_correct = "HG.2"
src.hg2.tracing@meta.data[src.hg2.tracing$tree%in%c(NA,"2.1.2.1.1", "2.1.2.1.2"),]$cluster.v06.26.re_correct = "Large.intestine.1"
src.hg2.tracing@meta.data[
  src.hg2.tracing$tree%in%c("2.2.2.1.1", '2.2.2.1.2',"2.2.2.2.1", "2.2.2.2.2"),]$cluster.v06.26.re_correct = "Large.intestine.3"
#-- correct for 9SS and 24/27SS
src.hg2.tracing@meta.data[
  src.hg2.tracing$cluster.v06.26.re_correct%in%"HG.2" &
    src.hg2.tracing$Time%in%c("24ss", "27ss"),]$cluster.v06.26.re_correct = "Large.intestine.3"
src.hg2.tracing@meta.data[
  src.hg2.tracing$cluster.v06.26.re_correct%in%c("Large.intestine.1", "Large.intestine.3") &
    src.hg2.tracing$Time%in%c("9ss"),]$cluster.v06.26.re_correct = "HG.2"
DimPlot(src.hg2.tracing, group.by = "cluster.v06.26.re_correct", reduction = "umap") 
table(src.hg2.tracing$cluster.v06.26.re_correct,
      src.hg2.tracing$Time)

#-- Pseudo time
# type_define = "cluster.v06.26.re_mnn_umap_fta"
type_define = "cluster.v06.26.re_correct"
#---------------------
#-- lar1 --
cell_src.hg2.tracing_lar1 = 
  rownames(src.hg2.tracing@meta.data[unlist(src.hg2.tracing[[type_define]]) %in%c('HG.2', "Large.intestine.1"),])
coord_src.hg2.tracing_lar1 = src.hg2.tracing[["umap"]]@cell.embeddings[cell_src.hg2.tracing_lar1,c(2,1)]
pcurve_src.hg2.tracing_lar1 = princurve::principal_curve(x = coord_src.hg2.tracing_lar1, smoother = "smooth.spline")
src.hg2.tracing$lambda_lar1 = pcurve_src.hg2.tracing_lar1$lambda[cell_src.hg2.tracing_lar1]
src.hg2.tracing@meta.data[cell_src.hg2.tracing_lar1,]$lambda_lar1 = 
  pcurve_src.hg2.tracing_lar1$lambda[cell_src.hg2.tracing_lar1]
src.hg2.tracing@meta.data[!src.hg2.tracing$lambda_lar1%in%NA,]$lambda_lar1 = 
  norm_range(src.hg2.tracing@meta.data[!src.hg2.tracing$lambda_lar1%in%NA,]$lambda_lar1)

#-- lar3 --
cell_src.hg2.tracing_lar3 = 
  rownames(src.hg2.tracing@meta.data[unlist(src.hg2.tracing[[type_define]]) %in%c('HG.2', "Large.intestine.3"),])
coord_src.hg2.tracing_lar3 = src.hg2.tracing[["umap"]]@cell.embeddings[cell_src.hg2.tracing_lar3, c(2,1)]
pcurve_src.hg2.tracing_lar3 = princurve::principal_curve(x = coord_src.hg2.tracing_lar3, smoother = "smooth.spline")
src.hg2.tracing$lambda_lar3 = pcurve_src.hg2.tracing_lar3$lambda[cell_src.hg2.tracing_lar3]
src.hg2.tracing@meta.data[cell_src.hg2.tracing_lar3,]$lambda_lar3 = 
  pcurve_src.hg2.tracing_lar3$lambda[cell_src.hg2.tracing_lar3]
src.hg2.tracing@meta.data[!src.hg2.tracing$lambda_lar3%in%NA,]$lambda_lar3 = 
  norm_range(src.hg2.tracing@meta.data[!src.hg2.tracing$lambda_lar3%in%NA,]$lambda_lar3)

cell_src.hg2.tracing_hg2 = rownames(src.hg2.tracing@reductions$umap@cell.embeddings[
  order(src.hg2.tracing@reductions$umap@cell.embeddings[,2]),])

cellorder_src.hg2.tracing = c(
  (intersect(cell_src.hg2.tracing_hg2,
             colnames(src.hg2.tracing[, unlist(src.hg2.tracing[[type_define]]) %in%"HG.2"]))),
  (intersect(cell_src.hg2.tracing_hg2,
             colnames(src.hg2.tracing[, unlist(src.hg2.tracing[[type_define]]) %in%"Large.intestine.3"]))),
  (intersect(names(src.hg2.tracing$lambda_lar1[order(src.hg2.tracing$lambda_lar1)]), 
             colnames(src.hg2.tracing[, unlist(src.hg2.tracing[[type_define]]) %in%"Large.intestine.1"])))
  )
#---------------------


pdf("figure.v08.07/organ_development_re_v240115/Heatmap_for_trajectory/try.hg2.re.pdf",10,10)
src.hg2.tracing.rowtree1  =  MyHeatmap(as.matrix(
  src.hg2.tracing@assays$RNA@data[
    gene_src.hg2.tracing.rowtree,
    cellorder_src.hg2.tracing]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.hg2.tracing$Time[cellorder_src.hg2.tracing],
               colors.time.2),
    MyName2Col(src.hg2.tracing$lineage[cellorder_src.hg2.tracing],
               color.lineage),
    MyName2Col(src.hg2.tracing$cluster.v06.26.re_correct[cellorder_src.hg2.tracing],
               cluster.endoderm.color.v5)
  ),
  RowSideColors = t(cbind(MyName2Col(
    names(gene_src.hg2.tracing.rowtree), colors.geneset))),
  ColSideColorsSize = 4,8,
  RowSideColorsSize = 1,2,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  Rowv = "none", Colv = "none",
  # return.tree = "row",
  margins = c(10,10),
  graph = T)
dev.off()

save(src.hg2.tracing,
     gene_src.hg2.tracing.rowtree,
     cellorder_src.hg2.tracing, 
     file = "figure.v08.07/organ_development_re_v240115/src.hg2.tracing.parameter.Rdata")

















