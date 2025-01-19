#-------------
#  MG.3 - re
#-------------
src.mg3.tracing = src.mg3.integrated.merge[,!src.mg3.integrated.merge$lineage%in%NA]
src.mg3.tracing$cluster.v06.26.re_mnn_umap_fta = 
  src.mg3.integrated.merge$cluster.v06.26.re_mnn_umap_fta[colnames(src.mg3.tracing)]

src.mg3.tracing@meta.data[
  intersect(rownames(src.9ss.integrated.merge@meta.data[
    src.9ss.integrated.merge$cluster.predict.umap_int.ext.v1.1%in%"MG.3",]),
    rownames(src.mg3.tracing@meta.data)),]$cluster.v06.26.re_mnn_umap_fta = "MG.3"

src.mg3.tracing = FindVariableFeatures(src.mg3.tracing, nfeatures = 2000)
src.mg3.tracing.filtergene = 
  Myfilter(as.matrix(src.mg3.tracing@assays$RNA@data),
           gene = src.mg3.tracing@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.mg3.tracing.filtergene.re = src.mg3.tracing.filtergene; rm(src.mg3.tracing.filtergene)
src.mg3.tracing.filtergene = src.mg3.tracing.filtergene.re; rm(src.mg3.tracing.filtergene.re)

src.mg3.tracing = SetIdent(src.mg3.tracing, value = src.mg3.tracing$cluster.v06.26.re_mnn_umap_fta)
marker_src.mg3.tracing = FindAllMarkers(src.mg3.tracing)
marker_src.mg3.tracing$pct.ratio = marker_src.mg3.tracing$pct.1 / marker_src.mg3.tracing$pct.2
marker_src.mg3.tracing$rank = marker_src.mg3.tracing$pct.ratio * (-log(marker_src.mg3.tracing$p_val_adj))
marker_src.mg3.tracing = marker_src.mg3.tracing[order(marker_src.mg3.tracing$rank, decreasing = T),]
markergene_src.mg3.tracing = unique(marker_src.mg3.tracing$gene)

src.mg3.tracing = RunPCA(src.mg3.tracing, features = src.mg3.tracing.filtergene)
src.mg3.tracing = RunUMAP(src.mg3.tracing, dims = 1:30, reduction = 'pca', 
                          n.neighbors = 100, n.components = 2)

src.mg3.tracing[["mnn_umap_fta"]] = src.mg3.tracing[["umap"]]
src.mg3.tracing@reductions$mnn_umap_fta@key = "Coord_"
src.mg3.tracing@reductions$mnn_umap_fta@cell.embeddings = 
  src.mg3.integrated.merge@reductions$mnn_umap_fta@cell.embeddings[colnames(src.mg3.tracing), c(1:2)]
colnames(src.mg3.tracing@reductions$mnn_umap_fta@cell.embeddings) = c("Coord_1","Coord_2")
DimPlot(src.mg3.tracing, reduction = "mnn_umap_fta")

src.mg3.tracing$cluster.v06.26.re_correct = 
  src.mg3.tracing$cluster.v06.26.re_mnn_umap_fta


src.mg3.tracing.sto = src.mg3.tracing[,src.mg3.tracing$cluster.v06.26.re_mnn_umap_fta%in%c("MG.3","MG.3.A/M","Stomach")]
src.mg3.tracing.pan = src.mg3.tracing[,src.mg3.tracing$cluster.v06.26.re_mnn_umap_fta%in%c("MG.3","MG.3.A/M","DP","EP.1","EP.2")]
src.mg3.tracing.int = src.mg3.tracing[,src.mg3.tracing$cluster.v06.26.re_mnn_umap_fta%in%c("MG.3","MG.3.P",'Small.intestine.1',"Small.intestine.2")]

for(i.seurat in c("src.mg3.tracing.sto", "src.mg3.tracing.pan", "src.mg3.tracing.int")){
  seurat = get(i.seurat)
  seurat = FindVariableFeatures(seurat, nfeatures = 2000)
  seurat.filtergene = 
    Myfilter(as.matrix(seurat@assays$RNA@data),
             gene = seurat@assays$RNA@var.features,
             pearson.threshold = 0.15,partner.threshold = 5,
             bottom.dispersion.interval = 0.1)
  seurat.filtergene.re = seurat.filtergene; rm(seurat.filtergene)
  seurat.filtergene = seurat.filtergene.re; rm(seurat.filtergene.re)
  
  seurat = SetIdent(seurat, value = seurat$cluster.v06.26.re_mnn_umap_fta)
  marker_seurat = FindAllMarkers(seurat)
  marker_seurat$pct.ratio = marker_seurat$pct.1 / marker_seurat$pct.2
  marker_seurat$rank = marker_seurat$pct.ratio * (-log(marker_seurat$p_val_adj))
  marker_seurat = marker_seurat[order(marker_seurat$rank, decreasing = T),]
  markergene_seurat = unique(marker_seurat$gene)
  
  assign(i.seurat, seurat)
  assign(paste(i.seurat,".filtergene", sep=""), seurat.filtergene)
  assign(paste("marker_",i.seurat, sep=""), marker_seurat)
  assign(paste("markergene_",i.seurat, sep=""), markergene_seurat)
}


# MG.3 Sto
#-----------------
src.mg3.tracing.sto@meta.data[src.mg3.tracing.sto$cluster.v06.26.re_correct%in%"MG.3.A/M",]$cluster.v06.26.re_correct = "MG.3.A/M"

type_define = "cluster.v06.26.re_correct"
#---------------------
#-- mg3 --
cell_src.mg3.tracing.sto_mg3 = 
  rownames(src.mg3.tracing.sto@meta.data[unlist(src.mg3.tracing.sto[[type_define]]) %in%c("MG.3"),])
coord_src.mg3.tracing.sto_mg3 = src.mg3.tracing.sto[["mnn_umap_fta"]]@cell.embeddings[cell_src.mg3.tracing.sto_mg3,c(2,1)]
pcurve_src.mg3.tracing.sto_mg3 = princurve::principal_curve(x = coord_src.mg3.tracing.sto_mg3, smoother = "smooth.spline")
src.mg3.tracing.sto$lambda_mg3 = pcurve_src.mg3.tracing.sto_mg3$lambda[cell_src.mg3.tracing.sto_mg3]
src.mg3.tracing.sto@meta.data[cell_src.mg3.tracing.sto_mg3,]$lambda_mg3 = 
  pcurve_src.mg3.tracing.sto_mg3$lambda[cell_src.mg3.tracing.sto_mg3]
src.mg3.tracing.sto@meta.data[!src.mg3.tracing.sto$lambda_mg3%in%NA,]$lambda_mg3 = 
  norm_range(src.mg3.tracing.sto@meta.data[!src.mg3.tracing.sto$lambda_mg3%in%NA,]$lambda_mg3)

#-- mg3m --
cell_src.mg3.tracing.sto_mg3m = 
  rownames(src.mg3.tracing.sto@meta.data[unlist(src.mg3.tracing.sto[[type_define]]) %in%c("MG.3.A/M","MG.3"),])
coord_src.mg3.tracing.sto_mg3m = src.mg3.tracing.sto[["mnn_umap_fta"]]@cell.embeddings[cell_src.mg3.tracing.sto_mg3m,c(2,1)]
pcurve_src.mg3.tracing.sto_mg3m = princurve::principal_curve(x = coord_src.mg3.tracing.sto_mg3m, smoother = "smooth.spline")
src.mg3.tracing.sto$lambda_mg3m = pcurve_src.mg3.tracing.sto_mg3m$lambda[cell_src.mg3.tracing.sto_mg3m]
src.mg3.tracing.sto@meta.data[cell_src.mg3.tracing.sto_mg3m,]$lambda_mg3m = 
  pcurve_src.mg3.tracing.sto_mg3m$lambda[cell_src.mg3.tracing.sto_mg3m]
src.mg3.tracing.sto@meta.data[!src.mg3.tracing.sto$lambda_mg3m%in%NA,]$lambda_mg3m = 
  norm_range(src.mg3.tracing.sto@meta.data[!src.mg3.tracing.sto$lambda_mg3m%in%NA,]$lambda_mg3m)

#-- sto --
cell_src.mg3.tracing.sto_sto = 
  rownames(src.mg3.tracing.sto@meta.data[unlist(src.mg3.tracing.sto[[type_define]]) %in%c("Stomach"),])
coord_src.mg3.tracing.sto_sto = src.mg3.tracing.sto[["mnn_umap_fta"]]@cell.embeddings[cell_src.mg3.tracing.sto_sto, c(2,1)]
pcurve_src.mg3.tracing.sto_sto = princurve::principal_curve(x = coord_src.mg3.tracing.sto_sto, smoother = "smooth.spline")
src.mg3.tracing.sto$lambda_sto = pcurve_src.mg3.tracing.sto_sto$lambda[cell_src.mg3.tracing.sto_sto]
src.mg3.tracing.sto@meta.data[cell_src.mg3.tracing.sto_sto,]$lambda_sto = 
  pcurve_src.mg3.tracing.sto_sto$lambda[cell_src.mg3.tracing.sto_sto]
src.mg3.tracing.sto@meta.data[!src.mg3.tracing.sto$lambda_sto%in%NA,]$lambda_sto = 
  norm_range(src.mg3.tracing.sto@meta.data[!src.mg3.tracing.sto$lambda_sto%in%NA,]$lambda_sto)

cellorder_src.mg3.tracing.sto = c(
  rev(intersect(names(src.mg3.tracing.sto$lambda_mg3m[order(src.mg3.tracing.sto$lambda_mg3m)]), 
                colnames(src.mg3.tracing.sto[, unlist(src.mg3.tracing.sto[[type_define]]) %in%"MG.3"]))),
  rev(intersect(names(src.mg3.tracing.sto$lambda_mg3m[order(src.mg3.tracing.sto$lambda_mg3m)]), 
                colnames(src.mg3.tracing.sto[, unlist(src.mg3.tracing.sto[[type_define]]) %in%"MG.3.A/M"]))),
  rev(intersect(names(src.mg3.tracing.sto$lambda_sto[order(src.mg3.tracing.sto$lambda_sto)]), 
                colnames(src.mg3.tracing.sto[, unlist(src.mg3.tracing.sto[[type_define]]) %in%"Stomach"])))
)
#---------------------


pdf("figure.v08.07/organ_development_re_v240115/try.mg3.sto.pdf",10,10)
gene_list = unique(c(markergene_src.mg3.tracing.sto,
                     src.mg3.tracing.sto.filtergene))
cell_list = colnames(src.mg3.tracing.sto)
#--------------------
src.mg3.tracing.sto.rowtree  =  MyHeatmap(as.matrix(
  src.mg3.tracing.sto@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.mg3.tracing.sto$Time[cell_list], #[cellorder_cell_src.mg3.tracing.sto],
               colors.time.2),
    MyName2Col(src.mg3.tracing.sto$lineage[cell_list], #[cellorder_cell_src.mg3.tracing.sto],
               color.lineage),
    MyName2Col(src.mg3.tracing.sto$cluster.v06.26.re_correct[cell_list], #[cellorder_cell_src.mg3.tracing.sto],
               cluster.endoderm.color.v5)
  ),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  #Rowv = "none",
  return.tree = "row",
  margins = c(10,10),
  graph = T)

src.mg3.tracing.sto.coltree  =  MyHeatmap(as.matrix(
  src.mg3.tracing.sto@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.mg3.tracing.sto$Time[cell_list], #[cellorder_cell_src.mg3.tracing.sto],
               colors.time.2),
    MyName2Col(src.mg3.tracing.sto$lineage[cell_list], #[cellorder_cell_src.mg3.tracing.sto],
               color.lineage),
    MyName2Col(src.mg3.tracing.sto$cluster.v06.26.re_correct[cell_list], #[cellorder_cell_src.mg3.tracing.sto],
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
tree_src.mg3.tracing.sto.rowtree = as.dendrogram(src.mg3.tracing.sto.rowtree)
tree_src.mg3.tracing.sto.coltree = as.dendrogram(src.mg3.tracing.sto.coltree)

gene_src.mg3.tracing.sto.rowtree = c(
  labels(tree_src.mg3.tracing.sto.rowtree[[1]][[2]][[1]][[1]]),
  labels(tree_src.mg3.tracing.sto.rowtree[[1]][[2]][[1]][[2]][[1]]),
  
  labels(tree_src.mg3.tracing.sto.rowtree[[1]][[2]][[2]][[2]]),
  labels(tree_src.mg3.tracing.sto.rowtree[[1]][[1]]),
  
  labels(tree_src.mg3.tracing.sto.rowtree[[2]][[2]][[2]][[2]][[1]]),
  
  labels(tree_src.mg3.tracing.sto.rowtree[[2]][[1]])
)


pdf("figure.v08.07/organ_development_re_v240115/try.mg3.sto.pdf",10,10)
gene_list = gene_src.mg3.tracing.sto.rowtree
cell_list = cellorder_src.mg3.tracing.sto
#--------------------
src.mg3.tracing.sto.rowtree1  =  MyHeatmap(as.matrix(
  src.mg3.tracing.sto@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.mg3.tracing.sto$Time[cell_list], #[cellorder_cell_src.mg3.tracing.sto],
               colors.time.2),
    MyName2Col(src.mg3.tracing.sto$lineage[cell_list], #[cellorder_cell_src.mg3.tracing.sto],
               color.lineage),
    MyName2Col(src.mg3.tracing.sto$cluster.v06.26.re_correct[cell_list], #[cellorder_cell_src.mg3.tracing.sto],
               cluster.endoderm.color.v5)
  ),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  #Rowv = "none",
  return.tree = "row",
  margins = c(10,10),
  graph = T)

src.mg3.tracing.sto.coltree1  =  MyHeatmap(as.matrix(
  src.mg3.tracing.sto@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.mg3.tracing.sto$Time[cell_list], #[cellorder_cell_src.mg3.tracing.sto],
               colors.time.2),
    MyName2Col(src.mg3.tracing.sto$lineage[cell_list], #[cellorder_cell_src.mg3.tracing.sto],
               color.lineage),
    MyName2Col(src.mg3.tracing.sto$cluster.v06.26.re_correct[cell_list], #[cellorder_cell_src.mg3.tracing.sto],
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
tree_src.mg3.tracing.sto.rowtree1 = as.dendrogram(src.mg3.tracing.sto.rowtree1)
tree_src.mg3.tracing.sto.coltree1 = as.dendrogram(src.mg3.tracing.sto.coltree1)

gene_src.mg3.tracing.sto.rowtree1 = c(
  (labels(tree_src.mg3.tracing.sto.rowtree1[[2]][[1]][[2]])),
  
  labels(tree_src.mg3.tracing.sto.rowtree1[[2]][[1]][[1]]),
  labels(tree_src.mg3.tracing.sto.rowtree1[[2]][[2]][[2]]),
  
  rev(labels(tree_src.mg3.tracing.sto.rowtree1[[2]][[2]][[1]])),
  
  labels(tree_src.mg3.tracing.sto.rowtree1[[1]][[1]][[2]]),
  labels(tree_src.mg3.tracing.sto.rowtree1[[1]][[1]][[1]]),
  labels(tree_src.mg3.tracing.sto.rowtree1[[1]][[2]][[2]])
)
names(gene_src.mg3.tracing.sto.rowtree1) = c(
  rep(4, length(c(
    rev(labels(tree_src.mg3.tracing.sto.rowtree1[[2]][[1]][[2]]))
  ))),
  rep(3, length(c(
    labels(tree_src.mg3.tracing.sto.rowtree1[[2]][[1]][[1]]),
    labels(tree_src.mg3.tracing.sto.rowtree1[[2]][[2]][[2]])
  ))),
  rep(5, length(c(
    labels(tree_src.mg3.tracing.sto.rowtree1[[2]][[2]][[1]])
  ))),
  rep(7, length(c(
    labels(tree_src.mg3.tracing.sto.rowtree1[[1]][[1]][[2]]),
    labels(tree_src.mg3.tracing.sto.rowtree1[[1]][[1]][[1]]),
    labels(tree_src.mg3.tracing.sto.rowtree1[[1]][[2]][[2]])
  )))
)

pdf("figure.v08.07/organ_development_re_v240115/try.mg3.sto.pdf",10,10)
gene_list = gene_src.mg3.tracing.sto.rowtree1
cell_list = cellorder_src.mg3.tracing.sto
#--------------------
src.mg3.tracing.sto.rowtree.fin = MyHeatmap(as.matrix(
  src.mg3.tracing.sto@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.mg3.tracing.sto$Time[cell_list], #[cellorder_cell_src.mg3.tracing.sto],
               colors.time.2),
    MyName2Col(src.mg3.tracing.sto$lineage[cell_list], #[cellorder_cell_src.mg3.tracing.sto],
               color.lineage),
    MyName2Col(src.mg3.tracing.sto$cluster.v06.26.re_correct[cell_list], #[cellorder_cell_src.mg3.tracing.sto],
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

save(src.mg3.tracing.sto,
     gene_src.mg3.tracing.sto.rowtree1,
     cellorder_src.mg3.tracing.sto,
     file = "figure.v08.07/organ_development_re_v240115/src.mg3.tracing.sto.parameter.Rdata")


# MG.3 Pan
#-----------------
src.mg3.tracing.pan@meta.data[src.mg3.tracing.pan$cluster.v06.26.re_correct%in%"MG.3.A/M",]$cluster.v06.26.re_correct = "MG.3.A/M"

type_define = "cluster.v06.26.re_correct"
#---------------------
#-- mg3 --
cell_src.mg3.tracing.pan_mg3 = 
  rownames(src.mg3.tracing.pan@meta.data[unlist(src.mg3.tracing.pan[[type_define]]) %in%c("MG.3"),])
coord_src.mg3.tracing.pan_mg3 = src.mg3.tracing.pan[["mnn_umap_fta"]]@cell.embeddings[cell_src.mg3.tracing.pan_mg3,c(2,1)]
pcurve_src.mg3.tracing.pan_mg3 = princurve::principal_curve(x = coord_src.mg3.tracing.pan_mg3, smoother = "smooth.spline")
src.mg3.tracing.pan$lambda_mg3 = pcurve_src.mg3.tracing.pan_mg3$lambda[cell_src.mg3.tracing.pan_mg3]
src.mg3.tracing.pan@meta.data[cell_src.mg3.tracing.pan_mg3,]$lambda_mg3 = 
  pcurve_src.mg3.tracing.pan_mg3$lambda[cell_src.mg3.tracing.pan_mg3]
src.mg3.tracing.pan@meta.data[!src.mg3.tracing.pan$lambda_mg3%in%NA,]$lambda_mg3 = 
  norm_range(src.mg3.tracing.pan@meta.data[!src.mg3.tracing.pan$lambda_mg3%in%NA,]$lambda_mg3)

#-- mg3m --
cell_src.mg3.tracing.pan_mg3m = 
  rownames(src.mg3.tracing.pan@meta.data[unlist(src.mg3.tracing.pan[[type_define]]) %in%c("MG.3.A/M","MG.3"),])
coord_src.mg3.tracing.pan_mg3m = src.mg3.tracing.pan[["mnn_umap_fta"]]@cell.embeddings[cell_src.mg3.tracing.pan_mg3m,c(2,1)]
pcurve_src.mg3.tracing.pan_mg3m = princurve::principal_curve(x = coord_src.mg3.tracing.pan_mg3m, smoother = "smooth.spline")
src.mg3.tracing.pan$lambda_mg3m = pcurve_src.mg3.tracing.pan_mg3m$lambda[cell_src.mg3.tracing.pan_mg3m]
src.mg3.tracing.pan@meta.data[cell_src.mg3.tracing.pan_mg3m,]$lambda_mg3m = 
  pcurve_src.mg3.tracing.pan_mg3m$lambda[cell_src.mg3.tracing.pan_mg3m]
src.mg3.tracing.pan@meta.data[!src.mg3.tracing.pan$lambda_mg3m%in%NA,]$lambda_mg3m = 
  norm_range(src.mg3.tracing.pan@meta.data[!src.mg3.tracing.pan$lambda_mg3m%in%NA,]$lambda_mg3m)

#-- dp--
cell_src.mg3.tracing.pan_dp = 
  rownames(src.mg3.tracing.pan@meta.data[unlist(src.mg3.tracing.pan[[type_define]]) %in%c("MG.3.A/M","DP"),])
coord_src.mg3.tracing.pan_dp = src.mg3.tracing.pan[["mnn_umap_fta"]]@cell.embeddings[cell_src.mg3.tracing.pan_dp, c(2,1)]
pcurve_src.mg3.tracing.pan_dp = princurve::principal_curve(x = coord_src.mg3.tracing.pan_dp, smoother = "smooth.spline")
src.mg3.tracing.pan$lambda_dp = pcurve_src.mg3.tracing.pan_dp$lambda[cell_src.mg3.tracing.pan_dp]
src.mg3.tracing.pan@meta.data[cell_src.mg3.tracing.pan_dp,]$lambda_dp = 
  pcurve_src.mg3.tracing.pan_dp$lambda[cell_src.mg3.tracing.pan_dp]
src.mg3.tracing.pan@meta.data[!src.mg3.tracing.pan$lambda_dp%in%NA,]$lambda_dp = 
  norm_range(src.mg3.tracing.pan@meta.data[!src.mg3.tracing.pan$lambda_dp%in%NA,]$lambda_dp)

#-- ep--
cell_src.mg3.tracing.pan_ep = 
  rownames(src.mg3.tracing.pan@meta.data[unlist(src.mg3.tracing.pan[[type_define]]) %in%c("DP"),])
coord_src.mg3.tracing.pan_ep = src.mg3.tracing.pan[["mnn_umap_fta"]]@cell.embeddings[cell_src.mg3.tracing.pan_ep, c(2,1)]
pcurve_src.mg3.tracing.pan_ep = princurve::principal_curve(x = coord_src.mg3.tracing.pan_ep, smoother = "smooth.spline")
src.mg3.tracing.pan$lambda_ep = pcurve_src.mg3.tracing.pan_ep$lambda[cell_src.mg3.tracing.pan_ep]
src.mg3.tracing.pan@meta.data[cell_src.mg3.tracing.pan_ep,]$lambda_ep = 
  pcurve_src.mg3.tracing.pan_ep$lambda[cell_src.mg3.tracing.pan_ep]
src.mg3.tracing.pan@meta.data[!src.mg3.tracing.pan$lambda_ep%in%NA,]$lambda_ep = 
  norm_range(src.mg3.tracing.pan@meta.data[!src.mg3.tracing.pan$lambda_ep%in%NA,]$lambda_ep)


cellorder_src.mg3.tracing.pan = c(
  rev(intersect(names(src.mg3.tracing.pan$lambda_mg3m[order(src.mg3.tracing.pan$lambda_mg3m)]), 
                colnames(src.mg3.tracing.pan[, unlist(src.mg3.tracing.pan[[type_define]]) %in%"MG.3"]))),
  rev(intersect(names(src.mg3.tracing.pan$lambda_mg3m[order(src.mg3.tracing.pan$lambda_mg3m)]), 
                colnames(src.mg3.tracing.pan[, unlist(src.mg3.tracing.pan[[type_define]]) %in%"MG.3.A/M"]))),
  rev(intersect(names(src.mg3.tracing.pan$lambda_dp[order(src.mg3.tracing.pan$lambda_dp)]), 
                colnames(src.mg3.tracing.pan[, unlist(src.mg3.tracing.pan[[type_define]]) %in%"DP"]))),
  (intersect(names(src.mg3.tracing.pan$lambda_ep[order(src.mg3.tracing.pan$lambda_ep)]), 
             colnames(src.mg3.tracing.pan[, unlist(src.mg3.tracing.pan[[type_define]]) %in%"EP.1"]))),
  (intersect(names(src.mg3.tracing.pan$lambda_ep[order(src.mg3.tracing.pan$lambda_ep)]), 
             colnames(src.mg3.tracing.pan[, unlist(src.mg3.tracing.pan[[type_define]]) %in%"EP.2"])))
)
#---------------------


pdf("figure.v08.07/organ_development_re_v240115/try.mg3.pan.pdf",10,10)
gene_list = unique(marker_src.mg3.tracing.pan[
  marker_src.mg3.tracing.pan$cluster%in%c("MG.3","MG.3.A/M","DP"), "gene"])
cell_list = colnames(src.mg3.tracing.pan)
#--------------------
src.mg3.tracing.pan.rowtree  =  MyHeatmap(as.matrix(
  src.mg3.tracing.pan@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.mg3.tracing.pan$Time[cell_list], #[cellorder_cell_src.mg3.tracing.pan],
               colors.time.2),
    MyName2Col(src.mg3.tracing.pan$lineage[cell_list], #[cellorder_cell_src.mg3.tracing.pan],
               color.lineage),
    MyName2Col(src.mg3.tracing.pan$cluster.v06.26.re_correct[cell_list], #[cellorder_cell_src.mg3.tracing.pan],
               cluster.endoderm.color.v5)
  ),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  #Rowv = "none",
  return.tree = "row",
  margins = c(10,10),
  graph = T)

src.mg3.tracing.pan.coltree  =  MyHeatmap(as.matrix(
  src.mg3.tracing.pan@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.mg3.tracing.pan$Time[cell_list], #[cellorder_cell_src.mg3.tracing.pan],
               colors.time.2),
    MyName2Col(src.mg3.tracing.pan$lineage[cell_list], #[cellorder_cell_src.mg3.tracing.pan],
               color.lineage),
    MyName2Col(src.mg3.tracing.pan$cluster.v06.26.re_correct[cell_list], #[cellorder_cell_src.mg3.tracing.pan],
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
tree_src.mg3.tracing.pan.rowtree = as.dendrogram(src.mg3.tracing.pan.rowtree)
tree_src.mg3.tracing.pan.coltree = as.dendrogram(src.mg3.tracing.pan.coltree)
gene_src.mg3.tracing.pan.rowtree = c(
  labels(tree_src.mg3.tracing.pan.rowtree[[1]][[1]][[2]][[2]][[2]]),
  
  labels(tree_src.mg3.tracing.pan.rowtree[[1]][[1]][[2]][[2]][[1]]),
  labels(tree_src.mg3.tracing.pan.rowtree[[1]][[1]][[2]][[1]]),
  
  rev(c(
    labels(tree_src.mg3.tracing.pan.rowtree[[2]][[2]][[1]]),
    labels(tree_src.mg3.tracing.pan.rowtree[[2]][[2]][[2]][[2]]),
    labels(tree_src.mg3.tracing.pan.rowtree[[2]][[2]][[2]][[1]])
  )),
  
  labels(tree_src.mg3.tracing.pan.rowtree[[2]][[1]])
)
names(gene_src.mg3.tracing.pan.rowtree) = c(
  rep(4, length(c(
    labels(tree_src.mg3.tracing.pan.rowtree[[1]][[1]][[2]][[2]][[2]])
  ))),
  rep(3, length(c(
    labels(tree_src.mg3.tracing.pan.rowtree[[1]][[1]][[2]][[2]][[1]]),
    labels(tree_src.mg3.tracing.pan.rowtree[[1]][[1]][[2]][[1]])
  ))),
  rep(5, length(c(
    labels(tree_src.mg3.tracing.pan.rowtree[[2]][[2]][[1]]),
    labels(tree_src.mg3.tracing.pan.rowtree[[2]][[2]][[2]][[2]]),
    labels(tree_src.mg3.tracing.pan.rowtree[[2]][[2]][[2]][[1]])
  ))),
  rep(2, length(c(
    labels(tree_src.mg3.tracing.pan.rowtree[[2]][[1]])
  )))
)

pdf("figure.v08.07/organ_development_re_v240115/try.mg3.pan.pdf",10,10)
gene_list = setdiff(
  unique(marker_src.mg3.tracing.pan[
    marker_src.mg3.tracing.pan$cluster%in%c("EP.1","EP.2"), "gene"]),
  labels(tree_src.mg3.tracing.pan.rowtree)
)
cell_list = cellorder_src.mg3.tracing.pan
#--------------------
src.mg3.tracing.pan.rowtree1  =  MyHeatmap(as.matrix(
  src.mg3.tracing.pan@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.mg3.tracing.pan$Time[cell_list], #[cellorder_cell_src.mg3.tracing.pan],
               colors.time.2),
    MyName2Col(src.mg3.tracing.pan$lineage[cell_list], #[cellorder_cell_src.mg3.tracing.pan],
               color.lineage),
    MyName2Col(src.mg3.tracing.pan$cluster.v06.26.re_correct[cell_list], #[cellorder_cell_src.mg3.tracing.pan],
               cluster.endoderm.color.v5)
  ),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  #Rowv = "none",
  return.tree = "row",
  margins = c(10,10),
  graph = T)

src.mg3.tracing.pan.coltree1  =  MyHeatmap(as.matrix(
  src.mg3.tracing.pan@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.mg3.tracing.pan$Time[cell_list], #[cellorder_cell_src.mg3.tracing.pan],
               colors.time.2),
    MyName2Col(src.mg3.tracing.pan$lineage[cell_list], #[cellorder_cell_src.mg3.tracing.pan],
               color.lineage),
    MyName2Col(src.mg3.tracing.pan$cluster.v06.26.re_correct[cell_list], #[cellorder_cell_src.mg3.tracing.pan],
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
tree_src.mg3.tracing.pan.rowtree1 = as.dendrogram(src.mg3.tracing.pan.rowtree1)
tree_src.mg3.tracing.pan.coltree1 = as.dendrogram(src.mg3.tracing.pan.coltree1)
gene_src.mg3.tracing.pan.rowtree1 = c(
  labels(tree_src.mg3.tracing.pan.rowtree1[[1]][[1]][[1]]),
  labels(tree_src.mg3.tracing.pan.rowtree1[[1]][[2]][[1]]),
  labels(tree_src.mg3.tracing.pan.rowtree1[[1]][[2]][[2]][[1]])
)
names(gene_src.mg3.tracing.pan.rowtree1) = c(
  rep(8, length(c(
    labels(tree_src.mg3.tracing.pan.rowtree1[[1]][[1]][[1]]),
    labels(tree_src.mg3.tracing.pan.rowtree1[[1]][[2]][[1]])
  ))),
  rep(7, length(c(
    labels(tree_src.mg3.tracing.pan.rowtree1[[1]][[2]][[2]][[1]])
  )))
)

gene_src.mg3.tracing.pan.rowtree.fin = c(
  gene_src.mg3.tracing.pan.rowtree,
  gene_src.mg3.tracing.pan.rowtree1
)

pdf("figure.v08.07/organ_development_re_v240115/try.mg3.pan.pdf",10,10)
gene_list = gene_src.mg3.tracing.pan.rowtree.fin
cell_list = cellorder_src.mg3.tracing.pan
#--------------------
src.mg3.tracing.pan.rowtree1  =  MyHeatmap(as.matrix(
  src.mg3.tracing.pan@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.mg3.tracing.pan$Time[cell_list], #[cellorder_cell_src.mg3.tracing.pan],
               colors.time.2),
    MyName2Col(src.mg3.tracing.pan$lineage[cell_list], #[cellorder_cell_src.mg3.tracing.pan],
               color.lineage),
    MyName2Col(src.mg3.tracing.pan$cluster.v06.26.re_correct[cell_list], #[cellorder_cell_src.mg3.tracing.pan],
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

save(src.mg3.tracing.pan,
     gene_src.mg3.tracing.pan.rowtree.fin,
     gene_src.mg3.tracing.pan.rowtree,
     gene_src.mg3.tracing.pan.rowtree1,
     cellorder_src.mg3.tracing.pan,
     file = "figure.v08.07/organ_development_re_v240115/src.mg3.tracing.pan.parameter.Rdata")



# MG.3 Intestine
#-----------------

type_define = "cluster.v06.26.re_correct"
#---------------------
#-- mg3 --
cell_src.mg3.tracing.int_mg3 = 
  rownames(src.mg3.tracing.int@meta.data[unlist(src.mg3.tracing.int[[type_define]]) %in%c("MG.3"),])
coord_src.mg3.tracing.int_mg3 = src.mg3.tracing.int[["mnn_umap_fta"]]@cell.embeddings[cell_src.mg3.tracing.int_mg3,c(2,1)]
pcurve_src.mg3.tracing.int_mg3 = princurve::principal_curve(x = coord_src.mg3.tracing.int_mg3, smoother = "smooth.spline")
src.mg3.tracing.int$lambda_mg3 = pcurve_src.mg3.tracing.int_mg3$lambda[cell_src.mg3.tracing.int_mg3]
src.mg3.tracing.int@meta.data[cell_src.mg3.tracing.int_mg3,]$lambda_mg3 = 
  pcurve_src.mg3.tracing.int_mg3$lambda[cell_src.mg3.tracing.int_mg3]
src.mg3.tracing.int@meta.data[!src.mg3.tracing.int$lambda_mg3%in%NA,]$lambda_mg3 = 
  norm_range(src.mg3.tracing.int@meta.data[!src.mg3.tracing.int$lambda_mg3%in%NA,]$lambda_mg3)

#-- mg3m --
cell_src.mg3.tracing.int_mg3m = 
  rownames(src.mg3.tracing.int@meta.data[unlist(src.mg3.tracing.int[[type_define]]) %in%c("MG.3.P","MG.3"),])
coord_src.mg3.tracing.int_mg3m = src.mg3.tracing.int[["mnn_umap_fta"]]@cell.embeddings[cell_src.mg3.tracing.int_mg3m,c(2,1)]
pcurve_src.mg3.tracing.int_mg3m = princurve::principal_curve(x = coord_src.mg3.tracing.int_mg3m, smoother = "smooth.spline")
src.mg3.tracing.int$lambda_mg3m = pcurve_src.mg3.tracing.int_mg3m$lambda[cell_src.mg3.tracing.int_mg3m]
src.mg3.tracing.int@meta.data[cell_src.mg3.tracing.int_mg3m,]$lambda_mg3m = 
  pcurve_src.mg3.tracing.int_mg3m$lambda[cell_src.mg3.tracing.int_mg3m]
src.mg3.tracing.int@meta.data[!src.mg3.tracing.int$lambda_mg3m%in%NA,]$lambda_mg3m = 
  norm_range(src.mg3.tracing.int@meta.data[!src.mg3.tracing.int$lambda_mg3m%in%NA,]$lambda_mg3m)

#-- dp--
cell_src.mg3.tracing.int_dp = 
  rownames(src.mg3.tracing.int@meta.data[unlist(src.mg3.tracing.int[[type_define]]) %in%c("MG.3.P","Small.intestine.1"),])
coord_src.mg3.tracing.int_dp = src.mg3.tracing.int[["mnn_umap_fta"]]@cell.embeddings[cell_src.mg3.tracing.int_dp, c(2,1)]
pcurve_src.mg3.tracing.int_dp = princurve::principal_curve(x = coord_src.mg3.tracing.int_dp, smoother = "smooth.spline")
src.mg3.tracing.int$lambda_dp = pcurve_src.mg3.tracing.int_dp$lambda[cell_src.mg3.tracing.int_dp]
src.mg3.tracing.int@meta.data[cell_src.mg3.tracing.int_dp,]$lambda_dp = 
  pcurve_src.mg3.tracing.int_dp$lambda[cell_src.mg3.tracing.int_dp]
src.mg3.tracing.int@meta.data[!src.mg3.tracing.int$lambda_dp%in%NA,]$lambda_dp = 
  norm_range(src.mg3.tracing.int@meta.data[!src.mg3.tracing.int$lambda_dp%in%NA,]$lambda_dp)

#-- ep--
cell_src.mg3.tracing.int_ep = 
  rownames(src.mg3.tracing.int@meta.data[unlist(src.mg3.tracing.int[[type_define]]) %in%c("MG.3.P","Small.intestine.2"),])
coord_src.mg3.tracing.int_ep = src.mg3.tracing.int[["mnn_umap_fta"]]@cell.embeddings[cell_src.mg3.tracing.int_ep, c(2,1)]
pcurve_src.mg3.tracing.int_ep = princurve::principal_curve(x = coord_src.mg3.tracing.int_ep, smoother = "smooth.spline")
src.mg3.tracing.int$lambda_ep = pcurve_src.mg3.tracing.int_ep$lambda[cell_src.mg3.tracing.int_ep]
src.mg3.tracing.int@meta.data[cell_src.mg3.tracing.int_ep,]$lambda_ep = 
  pcurve_src.mg3.tracing.int_ep$lambda[cell_src.mg3.tracing.int_ep]
src.mg3.tracing.int@meta.data[!src.mg3.tracing.int$lambda_ep%in%NA,]$lambda_ep = 
  norm_range(src.mg3.tracing.int@meta.data[!src.mg3.tracing.int$lambda_ep%in%NA,]$lambda_ep)


cellorder_src.mg3.tracing.int = c(
  (intersect(names(src.mg3.tracing.int$lambda_mg3m[order(src.mg3.tracing.int$lambda_mg3m)]), 
             colnames(src.mg3.tracing.int[, unlist(src.mg3.tracing.int[[type_define]]) %in%"MG.3"]))),
  (intersect(names(src.mg3.tracing.int$lambda_mg3m[order(src.mg3.tracing.int$lambda_mg3m)]), 
             colnames(src.mg3.tracing.int[, unlist(src.mg3.tracing.int[[type_define]]) %in%"MG.3.P"]))),
  rev(intersect(names(src.mg3.tracing.int$lambda_ep[order(src.mg3.tracing.int$lambda_ep)]), 
                colnames(src.mg3.tracing.int[, unlist(src.mg3.tracing.int[[type_define]]) %in%"Small.intestine.2"]))),
  rev(intersect(names(src.mg3.tracing.int$lambda_dp[order(src.mg3.tracing.int$lambda_dp)]), 
                colnames(src.mg3.tracing.int[, unlist(src.mg3.tracing.int[[type_define]]) %in%"Small.intestine.1"])))
)
#---------------------


pdf("figure.v08.07/organ_development_re_v240115/try.mg3.int.pdf",10,10)
gene_list = unique(marker_src.mg3.tracing.pan[
  marker_src.mg3.tracing.pan$cluster%in%c("MG.3","MG.3.P"), "gene"])
cell_list = colnames(src.mg3.tracing.int)
#--------------------
src.mg3.tracing.int.rowtree  =  MyHeatmap(as.matrix(
  src.mg3.tracing.int@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.mg3.tracing.int$Time[cell_list], #[cellorder_cell_src.mg3.tracing.int],
               colors.time.2),
    MyName2Col(src.mg3.tracing.int$lineage[cell_list], #[cellorder_cell_src.mg3.tracing.int],
               color.lineage),
    MyName2Col(src.mg3.tracing.int$cluster.v06.26.re_correct[cell_list], #[cellorder_cell_src.mg3.tracing.int],
               cluster.endoderm.color.v5)
  ),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  #Rowv = "none",
  return.tree = "row",
  margins = c(10,10),
  graph = T)

src.mg3.tracing.int.coltree  =  MyHeatmap(as.matrix(
  src.mg3.tracing.int@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.mg3.tracing.int$Time[cell_list], #[cellorder_cell_src.mg3.tracing.int],
               colors.time.2),
    MyName2Col(src.mg3.tracing.int$lineage[cell_list], #[cellorder_cell_src.mg3.tracing.int],
               color.lineage),
    MyName2Col(src.mg3.tracing.int$cluster.v06.26.re_correct[cell_list], #[cellorder_cell_src.mg3.tracing.int],
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
tree_src.mg3.tracing.int.rowtree = as.dendrogram(src.mg3.tracing.int.rowtree)
tree_src.mg3.tracing.int.coltree = as.dendrogram(src.mg3.tracing.int.coltree)
gene_src.mg3.tracing.int.rowtree = c(
  labels(tree_src.mg3.tracing.int.rowtree[[2]][[2]][[2]][[1]]),
  labels(tree_src.mg3.tracing.int.rowtree[[2]][[2]][[2]][[2]][[1]]),
  labels(tree_src.mg3.tracing.int.rowtree[[2]][[2]][[2]][[2]][[2]][[1]]),
  labels(tree_src.mg3.tracing.int.rowtree[[1]][[2]][[2]]),
  labels(tree_src.mg3.tracing.int.rowtree[[1]][[2]][[1]][[1]])
)


pdf("figure.v08.07/organ_development_re_v240115/try.mg3.int.pdf",10,10)
gene_list = setdiff(
  unique(marker_src.mg3.tracing.int[
    marker_src.mg3.tracing.int$cluster%in%c("Small.intestine.1","Small.intestine.2"), "gene"]),
  labels(tree_src.mg3.tracing.int.rowtree))
cell_list = cellorder_src.mg3.tracing.int
#--------------------
src.mg3.tracing.int.rowtree1  =  MyHeatmap(as.matrix(
  src.mg3.tracing.int@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.mg3.tracing.int$Time[cell_list], #[cellorder_cell_src.mg3.tracing.int],
               colors.time.2),
    MyName2Col(src.mg3.tracing.int$lineage[cell_list], #[cellorder_cell_src.mg3.tracing.int],
               color.lineage),
    MyName2Col(src.mg3.tracing.int$cluster.v06.26.re_correct[cell_list], #[cellorder_cell_src.mg3.tracing.int],
               cluster.endoderm.color.v5)
  ),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  #Rowv = "none",
  return.tree = "row",
  margins = c(10,10),
  graph = T)

src.mg3.tracing.int.coltree1  =  MyHeatmap(as.matrix(
  src.mg3.tracing.int@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.mg3.tracing.int$Time[cell_list], #[cellorder_cell_src.mg3.tracing.int],
               colors.time.2),
    MyName2Col(src.mg3.tracing.int$lineage[cell_list], #[cellorder_cell_src.mg3.tracing.int],
               color.lineage),
    MyName2Col(src.mg3.tracing.int$cluster.v06.26.re_correct[cell_list], #[cellorder_cell_src.mg3.tracing.int],
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
tree_src.mg3.tracing.int.rowtree1 = as.dendrogram(src.mg3.tracing.int.rowtree1)
tree_src.mg3.tracing.int.coltree1 = as.dendrogram(src.mg3.tracing.int.coltree1)
gene_src.mg3.tracing.int.rowtree1 = c(
  labels(tree_src.mg3.tracing.int.rowtree1[[1]]),
  labels(tree_src.mg3.tracing.int.rowtree1[[2]][[2]][[2]][[2]][[1]]),
  labels(tree_src.mg3.tracing.int.rowtree1[[2]][[1]][[2]]),
  labels(tree_src.mg3.tracing.int.rowtree1[[2]][[2]][[1]]),
  labels(tree_src.mg3.tracing.int.rowtree1[[2]][[2]][[2]][[2]][[2]])
)


pdf("figure.v08.07/organ_development_re_v240115/try.mg3.int.pdf",10,10)
gene_list = c(gene_src.mg3.tracing.int.rowtree,
              gene_src.mg3.tracing.int.rowtree1)
cell_list = cellorder_src.mg3.tracing.int
#--------------------
src.mg3.tracing.int.rowtree2  =  MyHeatmap(as.matrix(
  src.mg3.tracing.int@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.mg3.tracing.int$Time[cell_list], #[cellorder_cell_src.mg3.tracing.int],
               colors.time.2),
    MyName2Col(src.mg3.tracing.int$lineage[cell_list], #[cellorder_cell_src.mg3.tracing.int],
               color.lineage),
    MyName2Col(src.mg3.tracing.int$cluster.v06.26.re_correct[cell_list], #[cellorder_cell_src.mg3.tracing.int],
               cluster.endoderm.color.v5)
  ),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  #Rowv = "none",
  return.tree = "row",
  margins = c(10,10),
  graph = T)

src.mg3.tracing.int.coltree2  =  MyHeatmap(as.matrix(
  src.mg3.tracing.int@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.mg3.tracing.int$Time[cell_list], #[cellorder_cell_src.mg3.tracing.int],
               colors.time.2),
    MyName2Col(src.mg3.tracing.int$lineage[cell_list], #[cellorder_cell_src.mg3.tracing.int],
               color.lineage),
    MyName2Col(src.mg3.tracing.int$cluster.v06.26.re_correct[cell_list], #[cellorder_cell_src.mg3.tracing.int],
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
tree_src.mg3.tracing.int.rowtree2 = as.dendrogram(src.mg3.tracing.int.rowtree2)
tree_src.mg3.tracing.int.coltree2 = as.dendrogram(src.mg3.tracing.int.coltree2)


gene_order = function(gene, type = "MG.3"){
  data_temp = src.mg3.tracing.int@assays$RNA@data[
    gene,
    src.mg3.tracing.int$cluster.v06.26.re_correct%in%type]
  data_temp = rowSums(data_temp)
  gene = unlist(rev(names(data_temp)[order(data_temp)]))
  return(gene)
}

gene_src.mg3.tracing.int.rowtree2 = c(
  gene_order(
    gene = c(labels(tree_src.mg3.tracing.int.rowtree2[[1]][[1]]),
             rev(labels(tree_src.mg3.tracing.int.rowtree2[[1]][[2]]))),
    type = "MG.3"),
  
  labels(tree_src.mg3.tracing.int.rowtree2[[1]][[2]][[2]]),
  
  labels(tree_src.mg3.tracing.int.rowtree2[[2]][[2]][[1]]),
  
  rev(labels(tree_src.mg3.tracing.int.rowtree2[[2]][[1]])),
  
  labels(tree_src.mg3.tracing.int.rowtree2[[2]][[2]][[2]])
)
names(gene_src.mg3.tracing.int.rowtree2) = c(
  rep(4, length(c(
    labels(tree_src.mg3.tracing.int.rowtree2[[1]][[1]]),
    rev(labels(tree_src.mg3.tracing.int.rowtree2[[1]][[2]]))
  ))),
  rep(3, length(c(
    labels(tree_src.mg3.tracing.int.rowtree2[[1]][[2]][[2]])
  ))),
  rep(5, length(c(
    labels(tree_src.mg3.tracing.int.rowtree2[[2]][[2]][[1]])
  ))),
  rep(2, length(c(
    labels(tree_src.mg3.tracing.int.rowtree2[[2]][[1]])
  ))),
  rep(7, length(c(
    labels(tree_src.mg3.tracing.int.rowtree2[[2]][[2]][[2]])
  )))
)


pdf("figure.v08.07/organ_development_re_v240115/try.mg3.int.pdf",10,10)
gene_list = gene_src.mg3.tracing.int.rowtree2
cell_list = cellorder_src.mg3.tracing.int
#--------------------
src.mg3.tracing.int.rowtree1  =  MyHeatmap(as.matrix(
  src.mg3.tracing.int@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.mg3.tracing.int$Time[cell_list], #[cellorder_cell_src.mg3.tracing.int],
               colors.time.2),
    MyName2Col(src.mg3.tracing.int$lineage[cell_list], #[cellorder_cell_src.mg3.tracing.int],
               color.lineage),
    MyName2Col(src.mg3.tracing.int$cluster.v06.26.re_correct[cell_list], #[cellorder_cell_src.mg3.tracing.int],
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

save(src.mg3.tracing.int,
     gene_src.mg3.tracing.int.rowtree,
     gene_src.mg3.tracing.int.rowtree1,
     gene_src.mg3.tracing.int.rowtree2,
     cellorder_src.mg3.tracing.int,
     file = "figure.v08.07/organ_development_re_v240115/src.mg3.tracing.int.parameter.Rdata")