#---------------------------------------------------------------------------
# --          HG.1 -- Re
#---------------------------------------------------------------------------

#-----------
# HG.1-Re
#-----------
src.hg1.tracing = src.hg1.integrated.merge[,!src.hg1.integrated.merge$lineage%in%NA]
src.hg1.tracing$cluster.v06.26.re_mnn_umap_fta = 
  src.hg1.integrated.merge$cluster.v06.26.re_mnn_umap_fta[colnames(src.hg1.tracing)]
src.hg1.tracing@meta.data[
  intersect(rownames(src.9ss.integrated.merge@meta.data[
    src.9ss.integrated.merge$cluster.predict.umap_int.ext.v1.1%in%"HG.1",]),
    rownames(src.hg1.tracing@meta.data)),]$cluster.v06.26.re_mnn_umap_fta = "HG.1"

src.hg1.tracing = FindVariableFeatures(src.hg1.tracing, nfeatures = 2000)
src.hg1.tracing.filtergene = 
  Myfilter(as.matrix(src.hg1.tracing@assays$RNA@data),
           gene = src.hg1.tracing@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.hg1.tracing.filtergene.re = src.hg1.tracing.filtergene; rm(src.hg1.tracing.filtergene)
src.hg1.tracing.filtergene = src.hg1.tracing.filtergene.re; rm(src.hg1.tracing.filtergene.re)

src.hg1.tracing = SetIdent(src.hg1.tracing, value = src.hg1.tracing$cluster.v06.26.re_mnn_umap_fta)
marker_src.hg1.tracing = FindAllMarkers(src.hg1.tracing)
marker_src.hg1.tracing$pct.ratio = marker_src.hg1.tracing$pct.1 / marker_src.hg1.tracing$pct.2
marker_src.hg1.tracing$rank = marker_src.hg1.tracing$pct.ratio * (-log(marker_src.hg1.tracing$p_val_adj))
marker_src.hg1.tracing = marker_src.hg1.tracing[order(marker_src.hg1.tracing$rank, decreasing = T),]
markergene_src.hg1.tracing = unique(marker_src.hg1.tracing$gene)

src.hg1.tracing = RunPCA(src.hg1.tracing, features = src.hg1.tracing.filtergene)
src.hg1.tracing = RunUMAP(src.hg1.tracing, dims = 1:30, reduction = 'pca', 
                          n.neighbors = 100, n.components = 2)

src.hg1.tracing[["mnn_umap_fta"]] = src.hg1.tracing[["umap"]]
src.hg1.tracing@reductions$mnn_umap_fta@key = "Coord_"
src.hg1.tracing@reductions$mnn_umap_fta@cell.embeddings = 
  src.hg1.integrated.merge@reductions$mnn_umap_fta@cell.embeddings[colnames(src.hg1.tracing), c(1:2)]
colnames(src.hg1.tracing@reductions$mnn_umap_fta@cell.embeddings) = c("Coord_1","Coord_2")

DimPlot(src.hg1.tracing, reduction = 'mnn_umap_fta', 
        group.by = "cluster.v06.26.re_hc_mnn_umap_fta")


#--- hg1 - lar2
src.hg1.tracing.lar2 = src.hg1.tracing[, src.hg1.tracing$cluster.v06.26.re_hc_mnn_umap_fta%in%c(
  "HG.1", "HG.1-Large.intestine.2", "Large.intestine.2")]
src.hg1.tracing.lar2 = FindVariableFeatures(src.hg1.tracing.lar2, nfeatures = 2000)
src.hg1.tracing.lar2.filtergene = 
  Myfilter(as.matrix(src.hg1.tracing.lar2@assays$RNA@data),
           gene = src.hg1.tracing.lar2@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.hg1.tracing.lar2.filtergene.re = src.hg1.tracing.lar2.filtergene; rm(src.hg1.tracing.lar2.filtergene)
src.hg1.tracing.lar2.filtergene = src.hg1.tracing.lar2.filtergene.re; rm(src.hg1.tracing.lar2.filtergene.re)

src.hg1.tracing.lar2 = SetIdent(src.hg1.tracing.lar2, value = src.hg1.tracing.lar2$cluster.v06.26.re_mnn_umap_fta)
marker_src.hg1.tracing.lar2 = FindAllMarkers(src.hg1.tracing.lar2)
marker_src.hg1.tracing.lar2$pct.ratio = marker_src.hg1.tracing.lar2$pct.1 / marker_src.hg1.tracing.lar2$pct.2
marker_src.hg1.tracing.lar2$rank = marker_src.hg1.tracing.lar2$pct.ratio * (-log(marker_src.hg1.tracing.lar2$p_val_adj))
marker_src.hg1.tracing.lar2 = marker_src.hg1.tracing.lar2[order(marker_src.hg1.tracing.lar2$rank, decreasing = T),]
markergene_src.hg1.tracing.lar2 = unique(marker_src.hg1.tracing.lar2$gene)


pdf("figure.v08.07/organ_development_re_v240115/try.hg1.pdf",10,10)
#--------------------
src.hg1.tracing.lar2.rowtree  =  MyHeatmap(as.matrix(
  src.hg1.tracing.lar2@assays$RNA@data[
    unique(c(markergene_src.hg1.tracing.lar2,
             src.hg1.tracing.lar2.filtergene))
    ,]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.hg1.tracing.lar2$Time, #[cellorder_cell_src.hg1.tracing.lar2],
               colors.time.2),
    MyName2Col(src.hg1.tracing.lar2$lineage, #[cellorder_cell_src.hg1.tracing.lar2],
               color.lineage),
    MyName2Col(src.hg1.tracing.lar2$cluster.v06.26.re_hc_mnn_umap_fta, #[cellorder_cell_src.hg1.tracing.lar2],
               cluster.endoderm.color.v5)
  ),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  #Rowv = "none",
  return.tree = "row",
  margins = c(10,10),
  graph = T)

src.hg1.tracing.lar2.coltree  =  MyHeatmap(as.matrix(
  src.hg1.tracing.lar2@assays$RNA@data[
    unique(c(markergene_src.hg1.tracing.lar2,
             src.hg1.tracing.lar2.filtergene))
    ,]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.hg1.tracing.lar2$Time, #[cellorder_cell_src.hg1.tracing.lar2],
               colors.time.2),
    MyName2Col(src.hg1.tracing.lar2$lineage, #[cellorder_cell_src.hg1.tracing.lar2],
               color.lineage),
    MyName2Col(src.hg1.tracing.lar2$cluster.v06.26.re_hc_mnn_umap_fta, #[cellorder_cell_src.hg1.tracing.lar2],
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

tree_src.hg1.tracing.lar2.rowtree = as.dendrogram(src.hg1.tracing.lar2.rowtree)
gene_src.hg1.tracing.lar2.rowtree = c(
  labels(tree_src.hg1.tracing.lar2.rowtree[[2]][[1]][[1]][[1]]),
  labels(tree_src.hg1.tracing.lar2.rowtree[[2]][[1]][[1]][[2]][[1]]),
  labels(tree_src.hg1.tracing.lar2.rowtree[[2]][[1]][[2]]),
  
  labels(tree_src.hg1.tracing.lar2.rowtree[[1]][[2]][[2]][[2]][[1]]),
  
  labels(tree_src.hg1.tracing.lar2.rowtree[[2]][[2]])
)


pdf("figure.v08.07/organ_development_re_v240115/try.hg1.pdf",10,10)
#--------------------
src.hg1.tracing.lar2.rowtree1  =  MyHeatmap(as.matrix(
  src.hg1.tracing.lar2@assays$RNA@data[
    gene_src.hg1.tracing.lar2.rowtree
    ,]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.hg1.tracing.lar2$Time, #[cellorder_cell_src.hg1.tracing.lar2],
               colors.time.2),
    MyName2Col(src.hg1.tracing.lar2$lineage, #[cellorder_cell_src.hg1.tracing.lar2],
               color.lineage),
    MyName2Col(src.hg1.tracing.lar2$cluster.v06.26.re_hc_mnn_umap_fta, #[cellorder_cell_src.hg1.tracing.lar2],
               cluster.endoderm.color.v5)
  ),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  #Rowv = "none",
  return.tree = "row",
  margins = c(10,10),
  graph = T)

src.hg1.tracing.lar2.coltree1  =  MyHeatmap(as.matrix(
  src.hg1.tracing.lar2@assays$RNA@data[
    gene_src.hg1.tracing.lar2.rowtree
    ,]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.hg1.tracing.lar2$Time, #[cellorder_cell_src.hg1.tracing.lar2],
               colors.time.2),
    MyName2Col(src.hg1.tracing.lar2$lineage, #[cellorder_cell_src.hg1.tracing.lar2],
               color.lineage),
    MyName2Col(src.hg1.tracing.lar2$cluster.v06.26.re_hc_mnn_umap_fta, #[cellorder_cell_src.hg1.tracing.lar2],
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
tree_src.hg1.tracing.lar2.rowtree1 = as.dendrogram(src.hg1.tracing.lar2.rowtree1)
gene_src.hg1.tracing.lar2.rowtree1 = c(
  rev(labels(tree_src.hg1.tracing.lar2.rowtree1[[2]][[1]])),
  
  labels(tree_src.hg1.tracing.lar2.rowtree1[[2]][[2]][[2]][[2]]),
  labels(tree_src.hg1.tracing.lar2.rowtree1[[2]][[2]][[2]][[1]][[1]][[1]]),
  # labels(tree_src.hg1.tracing.lar2.rowtree1[[2]][[2]][[2]][[1]][[2]][[1]]),
  labels(tree_src.hg1.tracing.lar2.rowtree1[[2]][[2]][[2]][[1]][[2]][[2]][[1]]),
  
  rev(labels(tree_src.hg1.tracing.lar2.rowtree1[[1]][[1]])),
  
  rev(c(
    # labels(tree_src.hg1.tracing.lar2.rowtree1[[1]][[1]]),
    labels(tree_src.hg1.tracing.lar2.rowtree1[[1]][[2]][[1]]),
    labels(tree_src.hg1.tracing.lar2.rowtree1[[1]][[2]][[2]][[1]]),
    labels(tree_src.hg1.tracing.lar2.rowtree1[[1]][[2]][[2]][[2]][[1]])
  ))
)
names(gene_src.hg1.tracing.lar2.rowtree1) = c(
  rep(4, length(c(
    labels(tree_src.hg1.tracing.lar2.rowtree1[[2]][[1]])
  ))),
  rep(3, length(c(
    labels(tree_src.hg1.tracing.lar2.rowtree1[[2]][[2]][[2]][[2]]),
    labels(tree_src.hg1.tracing.lar2.rowtree1[[2]][[2]][[2]][[1]][[1]][[1]]),
    # labels(tree_src.hg1.tracing.lar2.rowtree1[[2]][[2]][[2]][[1]][[2]][[1]]),
    labels(tree_src.hg1.tracing.lar2.rowtree1[[2]][[2]][[2]][[1]][[2]][[2]][[1]])
  ))),
  rep(5, length(c(
    labels(tree_src.hg1.tracing.lar2.rowtree1[[1]][[1]])
  ))),
  rep(7, length(c(
    labels(tree_src.hg1.tracing.lar2.rowtree1[[1]][[2]][[1]]),
    # labels(tree_src.hg1.tracing.lar2.rowtree1[[1]][[1]]),
    labels(tree_src.hg1.tracing.lar2.rowtree1[[1]][[2]][[2]][[1]]),
    labels(tree_src.hg1.tracing.lar2.rowtree1[[1]][[2]][[2]][[2]][[1]])
  )))
)


gene_src.hg1.tracing.lar2.rowtree1.4 = gene_src.hg1.tracing.lar2.rowtree1[
  names(gene_src.hg1.tracing.lar2.rowtree1)%in%4]
data_temp = src.hg1.tracing.lar2@assays$RNA@data[
  gene_src.hg1.tracing.lar2.rowtree1.4,
  src.hg1.tracing.lar2$cluster.v06.26.re_correct%in%"HG.1"]
data_temp = rowSums(data_temp)
gene_src.hg1.tracing.lar2.rowtree1.4 = unlist(rev(names(data_temp)[order(data_temp)]))

gene_src.hg1.tracing.lar2.rowtree1.7 = gene_src.hg1.tracing.lar2.rowtree1[
  names(gene_src.hg1.tracing.lar2.rowtree1)%in%7]
data_temp = src.hg1.tracing.lar2@assays$RNA@data[
  gene_src.hg1.tracing.lar2.rowtree1.7,
  src.hg1.tracing.lar2$cluster.v06.26.re_correct%in%"Large.intestine.2"]
data_temp = rowSums(data_temp)
gene_src.hg1.tracing.lar2.rowtree1.7 = unlist(rev(names(data_temp)[order(data_temp)]))


gene_src.hg1.tracing.lar2.rowtree2 = c(
  gene_src.hg1.tracing.lar2.rowtree1.4,
  gene_src.hg1.tracing.lar2.rowtree1[names(
    gene_src.hg1.tracing.lar2.rowtree1)%in%c(3,5)],
  gene_src.hg1.tracing.lar2.rowtree1.7)
names(gene_src.hg1.tracing.lar2.rowtree2) = c(
  rep(4, length(gene_src.hg1.tracing.lar2.rowtree1.4)),
  names(gene_src.hg1.tracing.lar2.rowtree1[names(
    gene_src.hg1.tracing.lar2.rowtree1)%in%c(3,5)]),
  rep(7, length(gene_src.hg1.tracing.lar2.rowtree1.7)))

src.hg1.tracing.lar2$cluster.v06.26.re_correct = 
  src.hg1.tracing.lar2$cluster.v06.26.re_hc_mnn_umap_fta
src.hg1.tracing.lar2@meta.data[
  src.hg1.tracing.lar2$cluster.v06.26.re_hc_mnn_umap_fta%in%c("HG.1-Large.intestine.2") &
    src.hg1.tracing.lar2$Time%in%c("9ss"),]$cluster.v06.26.re_correct = "HG.1"

src.hg1.tracing$cluster.v06.26.re_correct = 
  src.hg1.tracing$cluster.v06.26.re_hc_mnn_umap_fta
src.hg1.tracing@meta.data[colnames(src.hg1.tracing.lar2),]$cluster.v06.26.re_correct =
  src.hg1.tracing.lar2$cluster.v06.26.re_correct
src.hg1.tracing@meta.data[
  src.hg1.tracing$cluster.v06.26.re_correct%in%"HG.1" &
    src.hg1.tracing$Time%in%c("24ss","27ss"),]$cluster.v06.26.re_correct = "Large.intestine.1"


gene_order = function(gene, type){
  data_temp = src.hg1.tracing.lar2@assays$RNA@data[
    gene,
    src.hg1.tracing.lar2$cluster.v06.26.re_correct%in%type]
  data_temp = rowSums(data_temp)
  gene = unlist(rev(names(data_temp)[order(data_temp)]))
  return(gene)
}

gene_src.hg1.tracing.lar2.rowtree2.fin = c(
  gene_order(gene_src.hg1.tracing.lar2.rowtree2[names(gene_src.hg1.tracing.lar2.rowtree2)==4], type = "HG.1-Large.intestine.2"),
  # gene_order(gene_src.hg1.tracing.lar2.rowtree2[names(gene_src.hg1.tracing.lar2.rowtree2)==3], type = "Large.intestine.2"),
  rev(gene_order(gene_src.hg1.tracing.lar2.rowtree2[names(gene_src.hg1.tracing.lar2.rowtree2)==5], type = "Large.intestine.2")),
  (gene_order(gene_src.hg1.tracing.lar2.rowtree2[names(gene_src.hg1.tracing.lar2.rowtree2)==7], type = "Large.intestine.2"))
)
names(gene_src.hg1.tracing.lar2.rowtree2.fin) = c(
  rep(4, length(gene_src.hg1.tracing.lar2.rowtree2[names(gene_src.hg1.tracing.lar2.rowtree2)==4])),
  # rep(3, length(gene_src.hg1.tracing.lar2.rowtree2[names(gene_src.hg1.tracing.lar2.rowtree2)==3])),
  rep(3, length(gene_src.hg1.tracing.lar2.rowtree2[names(gene_src.hg1.tracing.lar2.rowtree2)==5])),
  rep(7, length(gene_src.hg1.tracing.lar2.rowtree2[names(gene_src.hg1.tracing.lar2.rowtree2)==7]))
)


# type_define = "cluster.v06.26.re_hc_mnn_umap_fta"
type_define = 'cluster.v06.26.re_correct'
#---------------------
cell_src.hg1.tracing_hg1 = 
  rownames(src.hg1.tracing@meta.data[unlist(src.hg1.tracing[[type_define]]) %in%c('HG.1'),])
coord_src.hg1.tracing_hg1 = src.hg1.tracing[["mnn_umap_fta"]]@cell.embeddings[cell_src.hg1.tracing_hg1,c(2,1)]
pcurve_src.hg1.tracing_hg1 = princurve::principal_curve(x = coord_src.hg1.tracing_hg1, smoother = "smooth.spline")
src.hg1.tracing$lambda_hg1 = pcurve_src.hg1.tracing_hg1$lambda[cell_src.hg1.tracing_hg1]
src.hg1.tracing@meta.data[cell_src.hg1.tracing_hg1,]$lambda_hg1 = 
  pcurve_src.hg1.tracing_hg1$lambda[cell_src.hg1.tracing_hg1]
src.hg1.tracing@meta.data[!src.hg1.tracing$lambda_hg1%in%NA,]$lambda_hg1 = 
  norm_range(src.hg1.tracing@meta.data[!src.hg1.tracing$lambda_hg1%in%NA,]$lambda_hg1)

#-- pre-lar2 --
cell_src.hg1.tracing_prelar2 = 
  rownames(src.hg1.tracing@meta.data[unlist(src.hg1.tracing[[type_define]]) %in%c("HG.1-Large.intestine.2"),])
coord_src.hg1.tracing_prelar2 = src.hg1.tracing[["mnn_umap_fta"]]@cell.embeddings[cell_src.hg1.tracing_prelar2,c(2,1)]
pcurve_src.hg1.tracing_prelar2 = princurve::principal_curve(x = coord_src.hg1.tracing_prelar2, smoother = "smooth.spline")
src.hg1.tracing$lambda_prelar2 = pcurve_src.hg1.tracing_prelar2$lambda[cell_src.hg1.tracing_prelar2]
src.hg1.tracing@meta.data[cell_src.hg1.tracing_prelar2,]$lambda_prelar2 = 
  pcurve_src.hg1.tracing_prelar2$lambda[cell_src.hg1.tracing_prelar2]
src.hg1.tracing@meta.data[!src.hg1.tracing$lambda_prelar2%in%NA,]$lambda_prelar2 = 
  norm_range(src.hg1.tracing@meta.data[!src.hg1.tracing$lambda_prelar2%in%NA,]$lambda_prelar2)

#-- lar2 --
cell_src.hg1.tracing_lar2 = 
  rownames(src.hg1.tracing@meta.data[unlist(src.hg1.tracing[[type_define]]) %in%c("Large.intestine.2"),])
coord_src.hg1.tracing_lar2 = src.hg1.tracing[["mnn_umap_fta"]]@cell.embeddings[cell_src.hg1.tracing_lar2,c(2,1)]
pcurve_src.hg1.tracing_lar2 = princurve::principal_curve(x = coord_src.hg1.tracing_lar2, smoother = "smooth.spline")
src.hg1.tracing$lambda_lar2 = pcurve_src.hg1.tracing_lar2$lambda[cell_src.hg1.tracing_lar2]
src.hg1.tracing@meta.data[cell_src.hg1.tracing_lar2,]$lambda_lar2 = 
  pcurve_src.hg1.tracing_lar2$lambda[cell_src.hg1.tracing_lar2]
src.hg1.tracing@meta.data[!src.hg1.tracing$lambda_lar2%in%NA,]$lambda_lar2 = 
  norm_range(src.hg1.tracing@meta.data[!src.hg1.tracing$lambda_lar2%in%NA,]$lambda_lar2)

cellorder_src.hg1.tracing.lar2 = c(
  rev(intersect(names(src.hg1.tracing$lambda_hg1[order(src.hg1.tracing$lambda_hg1)]), 
                colnames(src.hg1.tracing[, unlist(src.hg1.tracing[[type_define]]) %in%"HG.1"]))),
  
  rev(intersect(names(src.hg1.tracing$lambda_prelar2[order(src.hg1.tracing$lambda_prelar2)]), 
                colnames(src.hg1.tracing[, unlist(src.hg1.tracing[[type_define]]) %in%"HG.1-Large.intestine.2"]))),
  
  rev(intersect(names(src.hg1.tracing$lambda_lar2[order(src.hg1.tracing$lambda_lar2)]), 
                colnames(src.hg1.tracing[, unlist(src.hg1.tracing[[type_define]]) %in%"Large.intestine.2"])))
)
#---------------------

pdf("figure.v08.07/organ_development_re_v240115/try.hg1.pdf",10,10)
gene = gene_src.hg1.tracing.lar2.rowtree2.fin
#--------------------
src.hg1.tracing.lar2.rowtree2  =  MyHeatmap(as.matrix(
  src.hg1.tracing.lar2@assays$RNA@data[
    gene,
    cellorder_src.hg1.tracing.lar2]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.hg1.tracing.lar2$Time[cellorder_src.hg1.tracing.lar2],
               colors.time.2),
    MyName2Col(src.hg1.tracing.lar2$lineage[cellorder_src.hg1.tracing.lar2],
               color.lineage),
    #MyName2Col(src.hg1.tracing.lar2$cluster.v06.26.re_hc_mnn_umap_fta[cellorder_src.hg1.tracing.lar2],
    #           cluster.endoderm.color.v5),
    MyName2Col(src.hg1.tracing.lar2$cluster.v06.26.re_correct[cellorder_src.hg1.tracing.lar2],
               cluster.endoderm.color.v5)
  ),
  RowSideColors = t(cbind(MyName2Col(
    names(gene), colors.geneset))),
  ColSideColorsSize = 4,8,
  RowSideColorsSize = 1,2,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  Rowv = "none", Colv = "none",
  # return.tree = "row",
  margins = c(10,10),
  graph = T)
#--------------------
dev.off()


save(src.hg1.tracing,
     src.hg1.tracing.filtergene,
     markergene_src.hg1.tracing.lar2,
     gene_src.hg1.tracing.lar2.rowtree1,
     gene_src.hg1.tracing.lar2.rowtree2,
     cellorder_src.hg1.tracing.lar2, 
     file = "figure.v08.07/organ_development_re_v240115/src.hg1.tracing.lar2.parameter.Rdata")


#--- hg1 - sm2lar1
src.hg1.tracing.sm2lar1= src.hg1.tracing[, src.hg1.tracing$cluster.v06.26.re_correct%in%c(
  "HG.1", "Large.intestine.1", "Small.intestine.2")]
src.hg1.tracing.sm2lar1= FindVariableFeatures(src.hg1.tracing.sm2lar1, nfeatures = 2000)
src.hg1.tracing.sm2lar1.filtergene = 
  Myfilter(as.matrix(src.hg1.tracing.sm2lar1@assays$RNA@data),
           gene = src.hg1.tracing.sm2lar1@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.hg1.tracing.sm2lar1.filtergene.re = src.hg1.tracing.sm2lar1.filtergene; rm(src.hg1.tracing.sm2lar1.filtergene)
src.hg1.tracing.sm2lar1.filtergene = src.hg1.tracing.sm2lar1.filtergene.re; rm(src.hg1.tracing.sm2lar1.filtergene.re)

src.hg1.tracing.sm2lar1= SetIdent(src.hg1.tracing.sm2lar1, value = src.hg1.tracing.sm2lar1$cluster.v06.26.re_correct)
marker_src.hg1.tracing.sm2lar1= FindAllMarkers(src.hg1.tracing.sm2lar1)
marker_src.hg1.tracing.sm2lar1$pct.ratio = marker_src.hg1.tracing.sm2lar1$pct.1 / marker_src.hg1.tracing.sm2lar1$pct.2
marker_src.hg1.tracing.sm2lar1$rank = marker_src.hg1.tracing.sm2lar1$pct.ratio * (-log(marker_src.hg1.tracing.sm2lar1$p_val_adj))
marker_src.hg1.tracing.sm2lar1= marker_src.hg1.tracing.sm2lar1[order(marker_src.hg1.tracing.sm2lar1$rank, decreasing = T),]
markergene_src.hg1.tracing.sm2lar1= unique(marker_src.hg1.tracing.sm2lar1$gene)

pdf("figure.v08.07/organ_development_re_v240115/try.hg1.pdf",10,10)
#--------------------
src.hg1.tracing.sm2lar1.rowtree  =  MyHeatmap(as.matrix(
  src.hg1.tracing.sm2lar1@assays$RNA@data[
    unique(c(markergene_src.hg1.tracing.sm2lar1,
             src.hg1.tracing.sm2lar1.filtergene))
    ,]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.hg1.tracing.sm2lar1$Time, #[cellorder_cell_src.hg1.tracing.sm2lar1],
               colors.time.2),
    MyName2Col(src.hg1.tracing.sm2lar1$lineage, #[cellorder_cell_src.hg1.tracing.sm2lar1],
               color.lineage),
    MyName2Col(src.hg1.tracing.sm2lar1$cluster.v06.26.re_hc_mnn_umap_fta, #[cellorder_cell_src.hg1.tracing.sm2lar1],
               cluster.endoderm.color.v5)
  ),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  #Rowv = "none",
  return.tree = "row",
  margins = c(10,10),
  graph = T)

src.hg1.tracing.sm2lar1.coltree  =  MyHeatmap(as.matrix(
  src.hg1.tracing.sm2lar1@assays$RNA@data[
    unique(c(markergene_src.hg1.tracing.sm2lar1,
             src.hg1.tracing.sm2lar1.filtergene))
    ,]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.hg1.tracing.sm2lar1$Time, #[cellorder_cell_src.hg1.tracing.sm2lar1],
               colors.time.2),
    MyName2Col(src.hg1.tracing.sm2lar1$lineage, #[cellorder_cell_src.hg1.tracing.sm2lar1],
               color.lineage),
    MyName2Col(src.hg1.tracing.sm2lar1$cluster.v06.26.re_hc_mnn_umap_fta, #[cellorder_cell_src.hg1.tracing.sm2lar1],
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


gene_order = function(gene){
  data_temp = src.hg1.tracing.sm2lar1@assays$RNA@data[
    gene,
    src.hg1.tracing.sm2lar1$cluster.v06.26.re_correct%in%"HG.1"]
  data_temp = rowSums(data_temp)
  gene = unlist(rev(names(data_temp)[order(data_temp)]))
  return(gene)
}

tree_src.hg1.tracing.sm2lar1.rowtree = as.dendrogram(src.hg1.tracing.sm2lar1.rowtree)
gene_src.hg1.tracing.sm2lar1.rowtree = c(
  gene_order(c(
    labels(tree_src.hg1.tracing.sm2lar1.rowtree[[2]][[1]][[1]][[2]][[2]]),
    labels(tree_src.hg1.tracing.sm2lar1.rowtree[[2]][[1]][[2]][[2]][[2]])
  )))
names(gene_src.hg1.tracing.sm2lar1.rowtree) = c(
  rep(4, length(c(
    labels(tree_src.hg1.tracing.sm2lar1.rowtree[[2]][[1]][[1]][[2]][[2]])
  ))),
  rep(4, length(c(
    labels(tree_src.hg1.tracing.sm2lar1.rowtree[[2]][[1]][[2]][[2]][[2]])
  )))
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gene_src.hg1.tracing.sm2lar1.rowtree1 = setdiff(
  unique(marker_src.hg1.tracing.sm2lar1[
    marker_src.hg1.tracing.sm2lar1$cluster%in%c("Large.intestine.1", "Small.intestine.2"), "gene"]),
  c(labels(tree_src.hg1.tracing.sm2lar1.rowtree[[2]][[1]][[1]]),
    labels(tree_src.hg1.tracing.sm2lar1.rowtree[[2]][[1]][[2]][[2]][[2]])))
cellist_src.hg1.tracing.sm2lar1 = rownames(src.hg1.tracing.sm2lar1@meta.data[
  src.hg1.tracing.sm2lar1$cluster.v06.26.re_correct%in%c("Large.intestine.1", "Small.intestine.2"),])

pdf("figure.v08.07/organ_development_re_v240115/try.hg1.pdf",10,10)
cell_list = cellist_src.hg1.tracing.sm2lar1 
#--------------------
src.hg1.tracing.sm2lar1.rowtree2  =  MyHeatmap(as.matrix(
  src.hg1.tracing.sm2lar1@assays$RNA@data[
    gene_src.hg1.tracing.sm2lar1.rowtree1,
    cellist_src.hg1.tracing.sm2lar1]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.hg1.tracing.sm2lar1$Time[cell_list],
               colors.time.2),
    MyName2Col(src.hg1.tracing.sm2lar1$lineage[cell_list],
               color.lineage),
    MyName2Col(src.hg1.tracing.sm2lar1$cluster.v06.26.re_hc_mnn_umap_fta[cell_list],
               cluster.endoderm.color.v5)
  ),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  #Rowv = "none",
  return.tree = "row",
  margins = c(10,10),
  graph = T)

src.hg1.tracing.sm2lar1.coltree2 = MyHeatmap(as.matrix(
  src.hg1.tracing.sm2lar1@assays$RNA@data[
    gene_src.hg1.tracing.sm2lar1.rowtree1,
    cellist_src.hg1.tracing.sm2lar1]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.hg1.tracing.sm2lar1$Time[cell_list],
               colors.time.2),
    MyName2Col(src.hg1.tracing.sm2lar1$lineage[cell_list],
               color.lineage),
    MyName2Col(src.hg1.tracing.sm2lar1$cluster.v06.26.re_hc_mnn_umap_fta[cell_list],
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

tree_src.hg1.tracing.sm2lar1.coltree2 = as.dendrogram(src.hg1.tracing.sm2lar1.coltree2)
src.hg1.tracing.sm2lar1$tree = NA
for(i1 in c(1:2)){ for(i2 in c(1:2)){ for(i3 in c(1:2)){
  src.hg1.tracing.sm2lar1@meta.data[
    labels(tree_src.hg1.tracing.sm2lar1.coltree2[[i1]][[i2]][[i3]]),]$tree = paste(i1,i2,i3,sep=".")
}}
}
DimPlot(src.hg1.tracing.sm2lar1, group.by = "tree",
        reduction = "mnn_umap_fta", cols = tree.hg1.sm2lar1) +
  DimPlot(src.hg1.tracing.sm2lar1, group.by = "Time", reduction = "mnn_umap_fta")

src.hg1.tracing.sm2lar1$cluster.v06.26.re_correct.refine = 
  src.hg1.tracing.sm2lar1$cluster.v06.26.re_correct 
src.hg1.tracing.sm2lar1@meta.data[
  src.hg1.tracing.sm2lar1$tree%in%c("1.1.1","1.1.2","1.2.1","1.2.2","2.2.1"),]$cluster.v06.26.re_correct.refine = 'Large.intestine.1'
src.hg1.tracing.sm2lar1@meta.data[
  src.hg1.tracing.sm2lar1$tree%in%c("2.2.2","2.1.1","2.1.2"),]$cluster.v06.26.re_correct.refine = 'Small.intestine.intestine.2'

src.hg1.tracing$cluster.v06.26.re_correct.refine = src.hg1.tracing$cluster.v06.26.re_correct
src.hg1.tracing@meta.data[colnames(src.hg1.tracing.sm2lar1),]$cluster.v06.26.re_correct.refine = 
  src.hg1.tracing.sm2lar1$cluster.v06.26.re_correct.refine
DimPlot(src.hg1.tracing, reduction = "mnn_umap_fta", group.by = "cluster.v06.26.re_correct.refine")


gene_order = function(gene){
  data_temp = src.hg1.tracing.sm2lar1@assays$RNA@data[
    gene,
    src.hg1.tracing.sm2lar1$cluster.v06.26.re_correct.refine%in%"HG.1"]
  data_temp = rowSums(data_temp)
  gene = unlist(rev(names(data_temp)[order(data_temp)]))
  return(gene)
}

tree_src.hg1.tracing.sm2lar1.rowtree2 = as.dendrogram(src.hg1.tracing.sm2lar1.rowtree2)
gene_src.hg1.tracing.sm2lar1.rowtree2 = c(
  gene_order(
    labels(tree_src.hg1.tracing.sm2lar1.rowtree2[[2]][[1]][[1]])
  ))
names(gene_src.hg1.tracing.sm2lar1.rowtree2) = c(
  rep(7, length(c(
    labels(tree_src.hg1.tracing.sm2lar1.rowtree2[[2]][[1]][[1]])
  )))
)

gene_src.hg1.tracing.sm2lar1.rowtree2.fin = c(
  rev(gene_src.hg1.tracing.sm2lar1.rowtree),
  gene_src.hg1.tracing.sm2lar1.rowtree2
)

save(src.hg1.tracing.sm2lar1,
     src.hg1.tracing,
     markergene_src.hg1.tracing.sm2lar1,
     gene_src.hg1.tracing.sm2lar1.rowtree,
     gene_src.hg1.tracing.sm2lar1.rowtree1,
     gene_src.hg1.tracing.sm2lar1.rowtree2,
     gene_src.hg1.tracing.sm2lar1.rowtree2.fin,
     file = "figure.v08.07/organ_development_re_v240115/src.hg1.tracing.sm2lar1.parameter.Rdata")


# type_define = "cluster.v06.26.re_hc_mnn_umap_fta"
type_define = 'cluster.v06.26.re_correct.refine'
#---------------------
cell_src.hg1.tracing_hg1 = 
  rownames(src.hg1.tracing@meta.data[unlist(src.hg1.tracing[[type_define]]) %in%c('HG.1'),])
coord_src.hg1.tracing_hg1 = src.hg1.tracing[["mnn_umap_fta"]]@cell.embeddings[cell_src.hg1.tracing_hg1,c(2,1)]
pcurve_src.hg1.tracing_hg1 = princurve::principal_curve(x = coord_src.hg1.tracing_hg1, smoother = "smooth.spline")
src.hg1.tracing$lambda_hg1 = pcurve_src.hg1.tracing_hg1$lambda[cell_src.hg1.tracing_hg1]
src.hg1.tracing@meta.data[cell_src.hg1.tracing_hg1,]$lambda_hg1 = 
  pcurve_src.hg1.tracing_hg1$lambda[cell_src.hg1.tracing_hg1]
src.hg1.tracing@meta.data[!src.hg1.tracing$lambda_hg1%in%NA,]$lambda_hg1 = 
  norm_range(src.hg1.tracing@meta.data[!src.hg1.tracing$lambda_hg1%in%NA,]$lambda_hg1)

#-- pre-sm2 --
cell_src.hg1.tracing_lar1 = 
  rownames(src.hg1.tracing@meta.data[unlist(src.hg1.tracing[[type_define]]) %in%c("Large.intestine.1"),])
coord_src.hg1.tracing_lar1 = src.hg1.tracing[["mnn_umap_fta"]]@cell.embeddings[cell_src.hg1.tracing_lar1,c(2,1)]
pcurve_src.hg1.tracing_lar1 = princurve::principal_curve(x = coord_src.hg1.tracing_lar1, smoother = "smooth.spline")
src.hg1.tracing$lambda_lar1 = pcurve_src.hg1.tracing_lar1$lambda[cell_src.hg1.tracing_lar1]
src.hg1.tracing@meta.data[cell_src.hg1.tracing_lar1,]$lambda_lar1 = 
  pcurve_src.hg1.tracing_lar1$lambda[cell_src.hg1.tracing_lar1]
src.hg1.tracing@meta.data[!src.hg1.tracing$lambda_lar1%in%NA,]$lambda_lar1 = 
  norm_range(src.hg1.tracing@meta.data[!src.hg1.tracing$lambda_lar1%in%NA,]$lambda_lar1)

#-- sm2 --
cell_src.hg1.tracing_sm2 = 
  rownames(src.hg1.tracing@meta.data[unlist(src.hg1.tracing[[type_define]]) %in%c("Small.intestine.2"),])
coord_src.hg1.tracing_sm2 = src.hg1.tracing[["mnn_umap_fta"]]@cell.embeddings[cell_src.hg1.tracing_sm2,c(2,1)]
pcurve_src.hg1.tracing_sm2 = princurve::principal_curve(x = coord_src.hg1.tracing_sm2, smoother = "smooth.spline")
src.hg1.tracing$lambda_sm2 = pcurve_src.hg1.tracing_sm2$lambda[cell_src.hg1.tracing_sm2]
src.hg1.tracing@meta.data[cell_src.hg1.tracing_sm2,]$lambda_sm2 = 
  pcurve_src.hg1.tracing_sm2$lambda[cell_src.hg1.tracing_sm2]
src.hg1.tracing@meta.data[!src.hg1.tracing$lambda_sm2%in%NA,]$lambda_sm2 = 
  norm_range(src.hg1.tracing@meta.data[!src.hg1.tracing$lambda_sm2%in%NA,]$lambda_sm2)

cellorder_src.hg1.tracing.sm2 = c(
  rev(intersect(names(src.hg1.tracing$lambda_hg1[order(src.hg1.tracing$lambda_hg1)]), 
                colnames(src.hg1.tracing[, unlist(src.hg1.tracing[[type_define]]) %in%"HG.1"]))),
  
  (intersect(names(src.hg1.tracing$lambda_lar1[order(src.hg1.tracing$lambda_lar1)]), 
             colnames(src.hg1.tracing[, unlist(src.hg1.tracing[[type_define]]) %in%"Large.intestine.1"]))),
  
  (intersect(names(src.hg1.tracing$lambda_sm2[order(src.hg1.tracing$lambda_sm2)]), 
             colnames(src.hg1.tracing[, unlist(src.hg1.tracing[[type_define]]) %in%"Small.intestine.2"])))
)
#---------------------

pdf("figure.v08.07/organ_development_re_v240115/try.hg1.pdf",10,10)
cell_list = cellorder_src.hg1.tracing.sm2 
gene_list = gene_src.hg1.tracing.sm2lar1.rowtree2.fin
#--------------------
src.hg1.tracing.sm2lar1.coltree2.re = MyHeatmap(as.matrix(
  src.hg1.tracing.sm2lar1@assays$RNA@data[
    gene_list,
    cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.hg1.tracing.sm2lar1$Time[cell_list],
               colors.time.2),
    MyName2Col(src.hg1.tracing.sm2lar1$lineage[cell_list],
               color.lineage),
    MyName2Col(src.hg1.tracing.sm2lar1$cluster.v06.26.re_correct.refine[cell_list],
               cluster.endoderm.color.v5)
  ),
  RowSideColors = t(cbind(MyName2Col(
    names(gene_list), colors.geneset))),
  ColSideColorsSize = 4,8,
  RowSideColorsSize = 1,2,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  Rowv = "none", Colv = "none",
  # return.tree = "row",
  margins = c(10,10),
  graph = T)
#--------------------
dev.off()

