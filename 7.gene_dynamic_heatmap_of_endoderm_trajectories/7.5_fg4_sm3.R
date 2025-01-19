#-------------
#  FG.4 - re
#-------------
# src.fg4.tracing.re = src.mg3.integrated.merge[,!src.mg3.integrated.merge$lineage%in%NA]
src.fg4.tracing.re$cluster.v06.26.re_mnn_umap_fta = 
  src.fg4.integrated.merge.re$cluster.v06.26.re_correct_refine_mnn_umap_fta[colnames(src.fg4.tracing.re)]

src.fg4.tracing.re@meta.data[
  intersect(rownames(src.9ss.integrated.merge@meta.data[
    src.9ss.integrated.merge$cluster.predict.umap_int.ext.v1.1%in%"FG.4",]),
    rownames(src.fg4.tracing.re@meta.data)),]$cluster.v06.26.re_mnn_umap_fta = "FG.4"

src.fg4.tracing.re = FindVariableFeatures(src.fg4.tracing.re, nfeatures = 2000)
src.fg4.tracing.re.filtergene = 
  Myfilter(as.matrix(src.fg4.tracing.re@assays$RNA@data),
           gene = src.fg4.tracing.re@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.fg4.tracing.re.filtergene.re = src.fg4.tracing.re.filtergene; rm(src.fg4.tracing.re.filtergene)
src.fg4.tracing.re.filtergene = src.fg4.tracing.re.filtergene.re; rm(src.fg4.tracing.re.filtergene.re)

src.fg4.tracing.re = SetIdent(src.fg4.tracing.re, value = src.fg4.tracing.re$cluster.v06.26.re_mnn_umap_fta)
marker_src.fg4.tracing.re = FindAllMarkers(src.fg4.tracing.re)
marker_src.fg4.tracing.re$pct.ratio = marker_src.fg4.tracing.re$pct.1 / marker_src.fg4.tracing.re$pct.2
marker_src.fg4.tracing.re$rank = marker_src.fg4.tracing.re$pct.ratio * (-log(marker_src.fg4.tracing.re$p_val_adj))
marker_src.fg4.tracing.re = marker_src.fg4.tracing.re[order(marker_src.fg4.tracing.re$rank, decreasing = T),]
markergene_src.fg4.tracing.re = unique(marker_src.fg4.tracing.re$gene)

src.fg4.tracing.re = RunPCA(src.fg4.tracing.re, features = src.fg4.tracing.re.filtergene)
src.fg4.tracing.re = RunUMAP(src.fg4.tracing.re, dims = 1:30, reduction = 'pca', 
                             n.neighbors = 100, n.components = 2)

src.fg4.tracing.re[["mnn_umap_fta"]] = src.fg4.tracing.re[["umap"]]
src.fg4.tracing.re@reductions$mnn_umap_fta@key = "Coord_"
src.fg4.tracing.re@reductions$mnn_umap_fta@cell.embeddings = 
  src.fg4.integrated.merge.re@reductions$mnn_umap_fta@cell.embeddings[colnames(src.fg4.tracing.re), c(1:2)]
colnames(src.fg4.tracing.re@reductions$mnn_umap_fta@cell.embeddings) = c("Coord_1","Coord_2")
DimPlot(src.fg4.tracing.re, reduction = "mnn_umap_fta")


src.fg4.tracing.re.liv = src.fg4.tracing.re[,src.fg4.tracing.re$cluster.v06.26.re_mnn_umap_fta%in%c("FG.4","FG.4-Liver",'AL.1/2-Liver',"Liver")]
src.fg4.tracing.re.lun = src.fg4.tracing.re[,src.fg4.tracing.re$cluster.v06.26.re_mnn_umap_fta%in%c("FG.4","FG.4-Lung/Stomach",'Lung',"Stomach")]
src.fg4.tracing.re.pha = src.fg4.tracing.re[,src.fg4.tracing.re$cluster.v06.26.re_mnn_umap_fta%in%c("FG.4","Pharynx.organ.5")]

for(i.seurat in c("src.fg4.tracing.re.liv", "src.fg4.tracing.re.lun", "src.fg4.tracing.re.pha")){
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


# FG.4 pha
#-----------------

type_define = "cluster.v06.26.re_mnn_umap_fta"
#---------------------
#-- mg3 --
cell_src.fg4.tracing.re.pha_mg3 = 
  rownames(src.fg4.tracing.re.pha@meta.data[unlist(src.fg4.tracing.re.pha[[type_define]]) %in%c("FG.4"),])
coord_src.fg4.tracing.re.pha_mg3 = src.fg4.tracing.re.pha[["mnn_umap_fta"]]@cell.embeddings[cell_src.fg4.tracing.re.pha_mg3,c(1,2)]
pcurve_src.fg4.tracing.re.pha_mg3 = princurve::principal_curve(x = coord_src.fg4.tracing.re.pha_mg3, smoother = "smooth.spline")
src.fg4.tracing.re.pha$lambda_mg3 = pcurve_src.fg4.tracing.re.pha_mg3$lambda[cell_src.fg4.tracing.re.pha_mg3]
src.fg4.tracing.re.pha@meta.data[cell_src.fg4.tracing.re.pha_mg3,]$lambda_mg3 = 
  pcurve_src.fg4.tracing.re.pha_mg3$lambda[cell_src.fg4.tracing.re.pha_mg3]
src.fg4.tracing.re.pha@meta.data[!src.fg4.tracing.re.pha$lambda_mg3%in%NA,]$lambda_mg3 = 
  norm_range(src.fg4.tracing.re.pha@meta.data[!src.fg4.tracing.re.pha$lambda_mg3%in%NA,]$lambda_mg3)

#-- sto --
cell_src.fg4.tracing.re.pha_sto = 
  rownames(src.fg4.tracing.re.pha@meta.data[unlist(src.fg4.tracing.re.pha[[type_define]]) %in%c("FG.4","Pharynx.organ.5"),])
coord_src.fg4.tracing.re.pha_sto = src.fg4.tracing.re.pha[["mnn_umap_fta"]]@cell.embeddings[cell_src.fg4.tracing.re.pha_sto, c(1,2)]
pcurve_src.fg4.tracing.re.pha_sto = princurve::principal_curve(x = coord_src.fg4.tracing.re.pha_sto, smoother = "smooth.spline")
src.fg4.tracing.re.pha$lambda_sto = pcurve_src.fg4.tracing.re.pha_sto$lambda[cell_src.fg4.tracing.re.pha_sto]
src.fg4.tracing.re.pha@meta.data[cell_src.fg4.tracing.re.pha_sto,]$lambda_sto = 
  pcurve_src.fg4.tracing.re.pha_sto$lambda[cell_src.fg4.tracing.re.pha_sto]
src.fg4.tracing.re.pha@meta.data[!src.fg4.tracing.re.pha$lambda_sto%in%NA,]$lambda_sto = 
  norm_range(src.fg4.tracing.re.pha@meta.data[!src.fg4.tracing.re.pha$lambda_sto%in%NA,]$lambda_sto)

cellorder_src.fg4.tracing.re.pha = c(
  (intersect(names(src.fg4.tracing.re.pha$lambda_mg3[order(src.fg4.tracing.re.pha$lambda_mg3)]), 
                colnames(src.fg4.tracing.re.pha[, unlist(src.fg4.tracing.re.pha[[type_define]]) %in%"FG.4"]))),
  rev(intersect(names(src.fg4.tracing.re.pha$lambda_sto[order(src.fg4.tracing.re.pha$lambda_sto)]), 
                colnames(src.fg4.tracing.re.pha[, unlist(src.fg4.tracing.re.pha[[type_define]]) %in%"Pharynx.organ.5"])))
)
#---------------------


pdf("figure.v08.07/organ_development_re_v240115/try.fg4.pha.pdf",10,10)
gene_list = unique(c(markergene_src.fg4.tracing.re.pha,
                     src.fg4.tracing.re.pha.filtergene))
cell_list = colnames(src.fg4.tracing.re.pha)
#--------------------
src.fg4.tracing.re.pha.rowtree  =  MyHeatmap(as.matrix(
  src.fg4.tracing.re.pha@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.fg4.tracing.re.pha$Time[cell_list], #[cellorder_cell_src.fg4.tracing.re.pha],
               colors.time.2),
    MyName2Col(src.fg4.tracing.re.pha$lineage[cell_list], #[cellorder_cell_src.fg4.tracing.re.pha],
               color.lineage),
    MyName2Col(src.fg4.tracing.re.pha$cluster.v06.26.re_mnn_umap_fta[cell_list], #[cellorder_cell_src.fg4.tracing.re.pha],
               cluster.endoderm.color.v5)
  ),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  #Rowv = "none",
  return.tree = "row",
  margins = c(10,10),
  graph = T)

src.fg4.tracing.re.pha.coltree  =  MyHeatmap(as.matrix(
  src.fg4.tracing.re.pha@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.fg4.tracing.re.pha$Time[cell_list], #[cellorder_cell_src.fg4.tracing.re.pha],
               colors.time.2),
    MyName2Col(src.fg4.tracing.re.pha$lineage[cell_list], #[cellorder_cell_src.fg4.tracing.re.pha],
               color.lineage),
    MyName2Col(src.fg4.tracing.re.pha$cluster.v06.26.re_mnn_umap_fta[cell_list], #[cellorder_cell_src.fg4.tracing.re.pha],
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
tree_src.fg4.tracing.re.pha.rowtree = as.dendrogram(src.fg4.tracing.re.pha.rowtree)
tree_src.fg4.tracing.re.pha.coltree = as.dendrogram(src.fg4.tracing.re.pha.coltree)

gene_order = function(gene, type){
  data_temp = src.fg4.tracing.re.pha@assays$RNA@data[
    gene,
    src.fg4.tracing.re.pha$cluster.v06.26.re_mnn_umap_fta%in%type]
  data_temp = rowSums(data_temp)
  gene = unlist(rev(names(data_temp)[order(data_temp)]))
  return(gene)
}

gene_src.fg4.tracing.re.pha.rowtree = c(
  rev(labels(tree_src.fg4.tracing.re.pha.rowtree[[1]])), 
  labels(tree_src.fg4.tracing.re.pha.rowtree[[2]][[1]])
)

names(gene_src.fg4.tracing.re.pha.rowtree) = c(
  rep(4, length(c(
    labels(tree_src.fg4.tracing.re.pha.rowtree[[1]])
  ))),
  rep(7, length(c(
    labels(tree_src.fg4.tracing.re.pha.rowtree[[2]][[1]])
  )))
)


pdf("figure.v08.07/organ_development_re_v240115/try.fg4.pha.pdf",10,10)
cell_list = cellorder_src.fg4.tracing.re.pha
gene_list = gene_src.fg4.tracing.re.pha.rowtree
#--------------------
src.fg4.tracing.re.pha.coltree1 = MyHeatmap(as.matrix(
  src.fg4.tracing.re.pha@assays$RNA@data[
    gene_list,
    cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.fg4.tracing.re.pha$Time[cell_list],
               colors.time.2),
    MyName2Col(src.fg4.tracing.re.pha$lineage[cell_list],
               color.lineage),
    MyName2Col(src.fg4.tracing.re.pha$cluster.v06.26.re_mnn_umap_fta[cell_list],
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

save(gene_src.fg4.tracing.re.pha.rowtree,
     src.fg4.tracing.re.pha,
     file = "figure.v08.07/organ_development_re_v240115/src.fg4.tracing.re.pha.parameter.Rdata")



# FG.4 liv
#-----------------

type_define = "cluster.v06.26.re_correct"
#---------------------
#-- mg3 --
cell_src.fg4.tracing.re.liv_mg3 = 
  rownames(src.fg4.tracing.re.liv@meta.data[unlist(src.fg4.tracing.re.liv[[type_define]]) %in%c("FG.4","FG.4-Liver"),])
coord_src.fg4.tracing.re.liv_mg3 = src.fg4.tracing.re.liv[["mnn_umap_fta"]]@cell.embeddings[cell_src.fg4.tracing.re.liv_mg3,c(2,1)]
pcurve_src.fg4.tracing.re.liv_mg3 = princurve::principal_curve(x = coord_src.fg4.tracing.re.liv_mg3, smoother = "smooth.spline")
src.fg4.tracing.re.liv$lambda_mg3 = pcurve_src.fg4.tracing.re.liv_mg3$lambda[cell_src.fg4.tracing.re.liv_mg3]
src.fg4.tracing.re.liv@meta.data[cell_src.fg4.tracing.re.liv_mg3,]$lambda_mg3 = 
  pcurve_src.fg4.tracing.re.liv_mg3$lambda[cell_src.fg4.tracing.re.liv_mg3]
src.fg4.tracing.re.liv@meta.data[!src.fg4.tracing.re.liv$lambda_mg3%in%NA,]$lambda_mg3 = 
  norm_range(src.fg4.tracing.re.liv@meta.data[!src.fg4.tracing.re.liv$lambda_mg3%in%NA,]$lambda_mg3)

#-- sto --
cell_src.fg4.tracing.re.liv_sto = 
  rownames(src.fg4.tracing.re.liv@meta.data[unlist(src.fg4.tracing.re.liv[[type_define]]) %in%c("AL.1/2-Liver","Liver"),])
coord_src.fg4.tracing.re.liv_sto = src.fg4.tracing.re.liv[["mnn_umap_fta"]]@cell.embeddings[cell_src.fg4.tracing.re.liv_sto, c(1,2)]
pcurve_src.fg4.tracing.re.liv_sto = princurve::principal_curve(x = coord_src.fg4.tracing.re.liv_sto, smoother = "smooth.spline")
src.fg4.tracing.re.liv$lambda_sto = pcurve_src.fg4.tracing.re.liv_sto$lambda[cell_src.fg4.tracing.re.liv_sto]
src.fg4.tracing.re.liv@meta.data[cell_src.fg4.tracing.re.liv_sto,]$lambda_sto = 
  pcurve_src.fg4.tracing.re.liv_sto$lambda[cell_src.fg4.tracing.re.liv_sto]
src.fg4.tracing.re.liv@meta.data[!src.fg4.tracing.re.liv$lambda_sto%in%NA,]$lambda_sto = 
  norm_range(src.fg4.tracing.re.liv@meta.data[!src.fg4.tracing.re.liv$lambda_sto%in%NA,]$lambda_sto)

cellorder_src.fg4.tracing.re.liv = c(
  rev(intersect(names(src.fg4.tracing.re.liv$lambda_mg3[order(src.fg4.tracing.re.liv$lambda_mg3)]), 
             colnames(src.fg4.tracing.re.liv[, unlist(src.fg4.tracing.re.liv[[type_define]]) %in%"FG.4"]))),
  (intersect(names(src.fg4.tracing.re.liv$lambda_mg3[order(src.fg4.tracing.re.liv$lambda_mg3)]), 
             colnames(src.fg4.tracing.re.liv[, unlist(src.fg4.tracing.re.liv[[type_define]]) %in%"FG.4-Liver"]))),
  (intersect(names(src.fg4.tracing.re.liv$lambda_sto[order(src.fg4.tracing.re.liv$lambda_sto)]), 
                colnames(src.fg4.tracing.re.liv[, unlist(src.fg4.tracing.re.liv[[type_define]]) %in%"AL.1/2-Liver"]))),
  (intersect(names(src.fg4.tracing.re.liv$lambda_sto[order(src.fg4.tracing.re.liv$lambda_sto)]), 
                colnames(src.fg4.tracing.re.liv[, unlist(src.fg4.tracing.re.liv[[type_define]]) %in%"Liver"])))
)
#---------------------


pdf("figure.v08.07/organ_development_re_v240115/try.fg4.liv.pdf",10,10)
gene_list = unique(c(markergene_src.fg4.tracing.re.liv
                     # src.fg4.tracing.re.liv.filtergene
                     ))
cell_list = colnames(src.fg4.tracing.re.liv)
#--------------------
src.fg4.tracing.re.liv.rowtree  =  MyHeatmap(as.matrix(
  src.fg4.tracing.re.liv@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.fg4.tracing.re.liv$Time[cell_list], #[cellorder_cell_src.fg4.tracing.re.liv],
               colors.time.2),
    MyName2Col(src.fg4.tracing.re.liv$lineage[cell_list], #[cellorder_cell_src.fg4.tracing.re.liv],
               color.lineage),
    MyName2Col(src.fg4.tracing.re.liv$cluster.v06.26.re_mnn_umap_fta[cell_list], #[cellorder_cell_src.fg4.tracing.re.liv],
               cluster.endoderm.color.v5)
  ),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  #Rowv = "none",
  return.tree = "row",
  margins = c(10,10),
  graph = T)

src.fg4.tracing.re.liv.coltree  =  MyHeatmap(as.matrix(
  src.fg4.tracing.re.liv@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.fg4.tracing.re.liv$Time[cell_list], #[cellorder_cell_src.fg4.tracing.re.liv],
               colors.time.2),
    MyName2Col(src.fg4.tracing.re.liv$lineage[cell_list], #[cellorder_cell_src.fg4.tracing.re.liv],
               color.lineage),
    MyName2Col(src.fg4.tracing.re.liv$cluster.v06.26.re_mnn_umap_fta[cell_list], #[cellorder_cell_src.fg4.tracing.re.liv],
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
tree_src.fg4.tracing.re.liv.rowtree = as.dendrogram(src.fg4.tracing.re.liv.rowtree)
tree_src.fg4.tracing.re.liv.coltree = as.dendrogram(src.fg4.tracing.re.liv.coltree)


gene_src.fg4.tracing.re.liv.rowtree = c(
  labels(tree_src.fg4.tracing.re.liv.rowtree[[2]][[2]][[2]][[1]]),
  
  labels(tree_src.fg4.tracing.re.liv.rowtree[[2]][[2]][[1]][[2]]),
  
  labels(tree_src.fg4.tracing.re.liv.rowtree[[2]][[1]][[2]][[2]][[1]]),
  
  labels(tree_src.fg4.tracing.re.liv.rowtree[[1]][[2]][[2]][[2]]),
  labels(tree_src.fg4.tracing.re.liv.rowtree[[1]][[1]])
)

names(gene_src.fg4.tracing.re.liv.rowtree) = c(
  rep(4, length(c(
    labels(tree_src.fg4.tracing.re.liv.rowtree[[2]][[2]][[2]][[1]])
  ))),
  rep(3, length(c(
    labels(tree_src.fg4.tracing.re.liv.rowtree[[2]][[2]][[1]][[2]])
  ))),
  rep(5, length(c(
    labels(tree_src.fg4.tracing.re.liv.rowtree[[2]][[1]][[2]][[2]][[1]])
  ))),
  rep(7, length(c(
    labels(tree_src.fg4.tracing.re.liv.rowtree[[1]][[2]][[2]][[2]]),
    labels(tree_src.fg4.tracing.re.liv.rowtree[[1]][[1]])
  )))
)

src.fg4.tracing.re.liv$tree = NA
for(i1 in c(1,2)){for(i2 in c(1,2)){for(i3 in c(1,2)){for(i4 in c(1,2)){
  src.fg4.tracing.re.liv@meta.data[
    labels(tree_src.fg4.tracing.re.liv.coltree[[i1]][[i2]][[i3]][[i4]]),]$tree = paste(i1,i2,i3,i4,sep = ".")
}}}}
DimPlot(src.fg4.tracing.re.liv, reduction = "mnn_umap_fta", group.by = "tree", label = T)


src.fg4.tracing.re.liv$cluster.v06.26.re_correct =
  src.fg4.tracing.re.liv$cluster.v06.26.re_mnn_umap_fta
src.fg4.tracing.re.liv@meta.data[
  src.fg4.tracing.re.liv$tree%in%"2.2.2.1",]$cluster.v06.26.re_correct = "FG.4-Liver"
src.fg4.tracing.re.liv@meta.data[
  src.fg4.tracing.re.liv$Time%in%"9ss",]$cluster.v06.26.re_correct = "FG.4"


gene_order = function(gene, type){
  data_temp = src.fg4.tracing.re.liv@assays$RNA@data[
    gene,
    src.fg4.tracing.re.liv$cluster.v06.26.re_correct%in%type]
  data_temp = rowSums(data_temp)
  gene = unlist(rev(names(data_temp)[order(data_temp)]))
  return(gene)
}

gene_src.fg4.tracing.re.liv.rowtree.fin = c(
  rev(gene_order(gene_src.fg4.tracing.re.liv.rowtree[names(gene_src.fg4.tracing.re.liv.rowtree)%in%4], type = "FG.4-Liver")),
  gene_src.fg4.tracing.re.liv.rowtree[names(gene_src.fg4.tracing.re.liv.rowtree)%in%3],
  gene_src.fg4.tracing.re.liv.rowtree[names(gene_src.fg4.tracing.re.liv.rowtree)%in%5],
  gene_order(gene_src.fg4.tracing.re.liv.rowtree[names(gene_src.fg4.tracing.re.liv.rowtree)%in%7], type = "Liver")
)
names(gene_src.fg4.tracing.re.liv.rowtree.fin) = c(
  rep(4, length(gene_src.fg4.tracing.re.liv.rowtree[names(gene_src.fg4.tracing.re.liv.rowtree)%in%4])),
  rep(3, length(gene_src.fg4.tracing.re.liv.rowtree[names(gene_src.fg4.tracing.re.liv.rowtree)%in%3])),
  rep(5, length(gene_src.fg4.tracing.re.liv.rowtree[names(gene_src.fg4.tracing.re.liv.rowtree)%in%5])),
  rep(7, length(gene_src.fg4.tracing.re.liv.rowtree[names(gene_src.fg4.tracing.re.liv.rowtree)%in%7]))
)

pdf("figure.v08.07/organ_development_re_v240115/try.fg4.liv.pdf",10,10)
cell_list = cellorder_src.fg4.tracing.re.liv
gene_list = gene_src.fg4.tracing.re.liv.rowtree.fin
#--------------------
src.fg4.tracing.re.liv.coltree1 = MyHeatmap(as.matrix(
  src.fg4.tracing.re.liv@assays$RNA@data[
    gene_list,
    cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.fg4.tracing.re.liv$Time[cell_list],
               colors.time.2),
    MyName2Col(src.fg4.tracing.re.liv$lineage[cell_list],
               color.lineage),
    MyName2Col(src.fg4.tracing.re.liv$cluster.v06.26.re_correct[cell_list],
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


save(gene_src.fg4.tracing.re.liv.rowtree,
     src.fg4.tracing.re.liv,
     file = "figure.v08.07/organ_development_re_v240115/src.fg4.tracing.re.liv.parameter.Rdata")



# FG.4 lun
#-----------------
DimPlot(src.fg4.tracing.re.lun, reduction = "mnn_umap_fta")

type_define = "cluster.v06.26.re_correct"
#---------------------
#-- mg3 --
cell_src.fg4.tracing.re.lun_mg3 = 
  rownames(src.fg4.tracing.re.lun@meta.data[unlist(src.fg4.tracing.re.lun[[type_define]]) %in%c("FG.4","FG.4-Lung/Stomach","Lung"),])
coord_src.fg4.tracing.re.lun_mg3 = src.fg4.tracing.re.lun[["mnn_umap_fta"]]@cell.embeddings[cell_src.fg4.tracing.re.lun_mg3,c(2,1)]
pcurve_src.fg4.tracing.re.lun_mg3 = princurve::principal_curve(x = coord_src.fg4.tracing.re.lun_mg3, smoother = "smooth.spline")
src.fg4.tracing.re.lun$lambda_mg3 = pcurve_src.fg4.tracing.re.lun_mg3$lambda[cell_src.fg4.tracing.re.lun_mg3]
src.fg4.tracing.re.lun@meta.data[cell_src.fg4.tracing.re.lun_mg3,]$lambda_mg3 = 
  pcurve_src.fg4.tracing.re.lun_mg3$lambda[cell_src.fg4.tracing.re.lun_mg3]
src.fg4.tracing.re.lun@meta.data[!src.fg4.tracing.re.lun$lambda_mg3%in%NA,]$lambda_mg3 = 
  norm_range(src.fg4.tracing.re.lun@meta.data[!src.fg4.tracing.re.lun$lambda_mg3%in%NA,]$lambda_mg3)

#-- sto --
cell_src.fg4.tracing.re.lun_sto = 
  rownames(src.fg4.tracing.re.lun@meta.data[unlist(src.fg4.tracing.re.lun[[type_define]]) %in%c("FG.4","FG.4-Lung/Stomach","Stomach"),])
coord_src.fg4.tracing.re.lun_sto = src.fg4.tracing.re.lun[["mnn_umap_fta"]]@cell.embeddings[cell_src.fg4.tracing.re.lun_sto, c(2,1)]
pcurve_src.fg4.tracing.re.lun_sto = princurve::principal_curve(x = coord_src.fg4.tracing.re.lun_sto, smoother = "smooth.spline")
src.fg4.tracing.re.lun$lambda_sto = pcurve_src.fg4.tracing.re.lun_sto$lambda[cell_src.fg4.tracing.re.lun_sto]
src.fg4.tracing.re.lun@meta.data[cell_src.fg4.tracing.re.lun_sto,]$lambda_sto = 
  pcurve_src.fg4.tracing.re.lun_sto$lambda[cell_src.fg4.tracing.re.lun_sto]
src.fg4.tracing.re.lun@meta.data[!src.fg4.tracing.re.lun$lambda_sto%in%NA,]$lambda_sto = 
  norm_range(src.fg4.tracing.re.lun@meta.data[!src.fg4.tracing.re.lun$lambda_sto%in%NA,]$lambda_sto)



cellorder_src.fg4.tracing.re.lun = c(
  (intersect(names(src.fg4.tracing.re.lun$lambda_mg3[order(src.fg4.tracing.re.lun$lambda_mg3)]), 
             colnames(src.fg4.tracing.re.lun[, unlist(src.fg4.tracing.re.lun[[type_define]]) %in%"FG.4"]))),
  (intersect(names(src.fg4.tracing.re.lun$lambda_mg3[order(src.fg4.tracing.re.lun$lambda_mg3)]), 
             colnames(src.fg4.tracing.re.lun[, unlist(src.fg4.tracing.re.lun[[type_define]]) %in%"FG.4-Lung/Stomach"]))),
  rev(intersect(rownames(
    src.fg4.tracing.re.lun@reductions$mnn_umap_fta@cell.embeddings[
      order(src.fg4.tracing.re.lun@reductions$mnn_umap_fta@cell.embeddings[,1]),]),
             colnames(src.fg4.tracing.re.lun[, unlist(src.fg4.tracing.re.lun[[type_define]]) %in%"Lung"]))),
  (intersect(names(src.fg4.tracing.re.lun$lambda_sto[order(src.fg4.tracing.re.lun$lambda_sto)]), 
                colnames(src.fg4.tracing.re.lun[, unlist(src.fg4.tracing.re.lun[[type_define]]) %in%"Stomach"])))
)
#---------------------

pdf("figure.v08.07/organ_development_re_v240115/try.fg4.lun.pdf",10,10)
gene_list = unique(c(markergene_src.fg4.tracing.re.lun
                     #src.fg4.tracing.re.lun.filtergene
                     ))
cell_list = colnames(src.fg4.tracing.re.lun)
#--------------------
src.fg4.tracing.re.lun.rowtree  =  MyHeatmap(as.matrix(
  src.fg4.tracing.re.lun@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.fg4.tracing.re.lun$Time[cell_list], #[cellorder_cell_src.fg4.tracing.re.lun],
               colors.time.2),
    MyName2Col(src.fg4.tracing.re.lun$lineage[cell_list], #[cellorder_cell_src.fg4.tracing.re.lun],
               color.lineage),
    MyName2Col(src.fg4.tracing.re.lun$cluster.v06.26.re_mnn_umap_fta[cell_list], #[cellorder_cell_src.fg4.tracing.re.lun],
               cluster.endoderm.color.v5)
  ),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  #Rowv = "none",
  return.tree = "row",
  margins = c(10,10),
  graph = T)

src.fg4.tracing.re.lun.coltree  =  MyHeatmap(as.matrix(
  src.fg4.tracing.re.lun@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.fg4.tracing.re.lun$Time[cell_list], #[cellorder_cell_src.fg4.tracing.re.lun],
               colors.time.2),
    MyName2Col(src.fg4.tracing.re.lun$lineage[cell_list], #[cellorder_cell_src.fg4.tracing.re.lun],
               color.lineage),
    MyName2Col(src.fg4.tracing.re.lun$cluster.v06.26.re_mnn_umap_fta[cell_list], #[cellorder_cell_src.fg4.tracing.re.lun],
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
tree_src.fg4.tracing.re.lun.rowtree = as.dendrogram(src.fg4.tracing.re.lun.rowtree)
tree_src.fg4.tracing.re.lun.coltree = as.dendrogram(src.fg4.tracing.re.lun.coltree)

gene_src.fg4.tracing.re.lun.rowtree = c(
  labels(tree_src.fg4.tracing.re.lun.rowtree[[1]][[2]][[2]][[2]]),
  labels(tree_src.fg4.tracing.re.lun.rowtree[[1]][[2]][[1]][[1]]),
  
  rev(labels(tree_src.fg4.tracing.re.lun.rowtree[[2]][[1]])),
  
  labels(tree_src.fg4.tracing.re.lun.rowtree[[2]][[2]][[2]][[2]][[1]][[2]])
)

names(gene_src.fg4.tracing.re.lun.rowtree) = c(
  rep(4, length(c(
    labels(tree_src.fg4.tracing.re.lun.rowtree[[1]][[2]][[2]][[2]]),
    labels(tree_src.fg4.tracing.re.lun.rowtree[[1]][[2]][[1]][[1]])
  ))),
  rep(3, length(c(
    labels(tree_src.fg4.tracing.re.lun.rowtree[[2]][[1]])
  ))),
  rep(7, length(c(
    labels(tree_src.fg4.tracing.re.lun.rowtree[[2]][[2]][[2]][[2]][[1]][[2]])
  )))
)



pdf("figure.v08.07/organ_development_re_v240115/try.fg4.lun.pdf",10,10)
gene_list = gene_src.fg4.tracing.re.lun.rowtree
cell_list = colnames(src.fg4.tracing.re.lun[,src.fg4.tracing.re.lun$cluster.v06.26.re_mnn_umap_fta%in%c("Lung","FG.4-Lung/Stomach")])
#--------------------
src.fg4.tracing.re.lun.rowtree1  =  MyHeatmap(as.matrix(
  src.fg4.tracing.re.lun@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.fg4.tracing.re.lun$Time[cell_list], #[cellorder_cell_src.fg4.tracing.re.lun],
               colors.time.2),
    MyName2Col(src.fg4.tracing.re.lun$lineage[cell_list], #[cellorder_cell_src.fg4.tracing.re.lun],
               color.lineage),
    MyName2Col(src.fg4.tracing.re.lun$tree[cell_list], #[cellorder_cell_src.fg4.tracing.re.lun],
               tree.hg1.sm2lar1),
    MyName2Col(src.fg4.tracing.re.lun$cluster.v06.26.re_mnn_umap_fta[cell_list], #[cellorder_cell_src.fg4.tracing.re.lun],
               cluster.endoderm.color.v5)
  ),
  RowSideColors = t(cbind(MyName2Col(
    names(gene_list), colors.geneset))),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  #Rowv = "none",
  return.tree = "row",
  margins = c(10,10),
  graph = T)

src.fg4.tracing.re.lun.coltree1  =  MyHeatmap(as.matrix(
  src.fg4.tracing.re.lun@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.fg4.tracing.re.lun$Time[cell_list], #[cellorder_cell_src.fg4.tracing.re.lun],
               colors.time.2),
    MyName2Col(src.fg4.tracing.re.lun$lineage[cell_list], #[cellorder_cell_src.fg4.tracing.re.lun],
               color.lineage),
    MyName2Col(src.fg4.tracing.re.lun$cluster.v06.26.re_mnn_umap_fta[cell_list], #[cellorder_cell_src.fg4.tracing.re.lun],
               cluster.endoderm.color.v5)
  ),
  RowSideColors = t(cbind(MyName2Col(
    names(gene_list), colors.geneset))),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  #Rowv = "none",
  return.tree = "col",
  margins = c(10,10),
  graph = T)
#--------------------
dev.off()
tree_src.fg4.tracing.re.lun.rowtree1 = as.dendrogram(src.fg4.tracing.re.lun.rowtree1)
tree_src.fg4.tracing.re.lun.coltree1 = as.dendrogram(src.fg4.tracing.re.lun.coltree1)

src.fg4.tracing.re.lun$tree = NA
for(i1 in c(1,2)){for(i2 in c(1,2)){for(i3 in c(1,2)){
  src.fg4.tracing.re.lun@meta.data[
    labels(tree_src.fg4.tracing.re.lun.coltree1[[i1]][[i2]][[i3]]), ]$tree = paste(i1,i2,i3,sep=".")
}}}
DimPlot(src.fg4.tracing.re.lun, reduction = "mnn_umap_fta", 
        group.by = "tree", label = T, cols = tree.hg1.sm2lar1)

src.fg4.tracing.re.lun$cluster.v06.26.re_correct = src.fg4.tracing.re.lun$cluster.v06.26.re_mnn_umap_fta
src.fg4.tracing.re.lun@meta.data[
  src.fg4.tracing.re.lun$tree%in%c("1.1.2","1.1.1","1.2.1","1.2.2","2.2.2","2.2.1"),]$cluster.v06.26.re_correct = "FG.4-Lung/Stomach"
src.fg4.tracing.re.lun@meta.data[
  src.fg4.tracing.re.lun$tree%in%c("2.1.1","2.1.2"),]$cluster.v06.26.re_correct = "Lung"
src.fg4.tracing.re.lun@meta.data[
  src.fg4.tracing.re.lun$tree%in%c("1.1.2","1.1.1","1.2.1","1.2.2","2.2.2","2.2.1") &
    src.fg4.tracing.re.lun$Time%in%c("24ss","27ss"),]$cluster.v06.26.re_correct = "Lung"

DimPlot(src.fg4.tracing.re.lun, reduction = "mnn_umap_fta", 
        group.by = "cluster.v06.26.re_correct", label = T)


pdf("figure.v08.07/organ_development_re_v240115/try.fg4.lun.pdf",10,10)
cell_list = cellorder_src.fg4.tracing.re.lun
gene_list = gene_src.fg4.tracing.re.lun.rowtree
#--------------------
src.fg4.tracing.re.lun.coltree1 = MyHeatmap(as.matrix(
  src.fg4.tracing.re.lun@assays$RNA@data[
    gene_list,
    cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.fg4.tracing.re.lun$Time[cell_list],
               colors.time.2),
    MyName2Col(src.fg4.tracing.re.lun$lineage[cell_list],
               color.lineage),
    MyName2Col(src.fg4.tracing.re.lun$cluster.v06.26.re_correct[cell_list],
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


save(gene_src.fg4.tracing.re.lun.rowtree,
     cellorder_src.fg4.tracing.re.lun,
     src.fg4.tracing.re.lun,
     file = "figure.v08.07/organ_development_re_v240115/src.fg4.tracing.re.lun.parameter.Rdata")



src.fg4.tracing.re$cluster.v06.26.re_correct = src.fg4.tracing.re$cluster.v06.26.re_mnn_umap_fta
src.fg4.tracing.re@meta.data[
  colnames(src.fg4.tracing.re.liv),]$cluster.v06.26.re_correct = src.fg4.tracing.re.liv$cluster.v06.26.re_correct
src.fg4.tracing.re@meta.data[
  colnames(src.fg4.tracing.re.lun),]$cluster.v06.26.re_correct = src.fg4.tracing.re.lun$cluster.v06.26.re_correct

save(src.fg4.tracing.re,
     file = "figure.v08.07/organ_development_re_v240115/src.fg4.tracing.re.Rdata")













