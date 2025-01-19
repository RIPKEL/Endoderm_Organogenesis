#-----------
# AL.3-Re
#-----------

src.al3.tracing_re@meta.data[
  intersect(rownames(src.9ss.integrated.merge@meta.data[
    src.9ss.integrated.merge$cluster.predict.umap_int.ext.v1.1%in%"AL.3",]),
    rownames(src.al3.tracing_re@meta.data)),]$cluster.v06.26.re_correct_mnn_umap_fta = "AL.3"

src.al3.tracing_re = FindVariableFeatures(src.al3.tracing_re, nfeatures = 2000)
src.al3.tracing_re.filtergene = 
  Myfilter(as.matrix(src.al3.tracing_re@assays$RNA@data),
           gene = src.al3.tracing_re@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.al3.tracing_re.filtergene.re = src.al3.tracing_re.filtergene; rm(src.al3.tracing_re.filtergene)
src.al3.tracing_re.filtergene = src.al3.tracing_re.filtergene.re; rm(src.al3.tracing_re.filtergene.re)

src.al3.tracing_re = SetIdent(src.al3.tracing_re, value = src.al3.tracing_re$cluster.v06.26.re_mnn_umap_fta)
marker_src.al3.tracing_re = FindAllMarkers(src.al3.tracing_re)
marker_src.al3.tracing_re$pct.ratio = marker_src.al3.tracing_re$pct.1 / marker_src.al3.tracing_re$pct.2
marker_src.al3.tracing_re$rank = marker_src.al3.tracing_re$pct.ratio * (-log(marker_src.al3.tracing_re$p_val_adj))
marker_src.al3.tracing_re = marker_src.al3.tracing_re[order(marker_src.al3.tracing_re$rank, decreasing = T),]
markergene_src.al3.tracing_re = unique(marker_src.al3.tracing_re$gene)

src.al3.tracing_re = RunPCA(src.al3.tracing_re, features = src.al3.tracing_re.filtergene)
src.al3.tracing_re = RunUMAP(src.al3.tracing_re, dims = 1:30, reduction = 'pca', 
                             n.neighbors = 100, n.components = 2)

src.al3.tracing_re[["mnn_umap_fta"]] = src.al3.tracing_re[["umap"]]
src.al3.tracing_re@reductions$mnn_umap_fta@key = "Coord_"
src.al3.tracing_re@reductions$mnn_umap_fta@cell.embeddings = 
  src.al3.integrated.merge.re@reductions$mnn_umap_fta@cell.embeddings[colnames(src.al3.tracing_re), c(1:2)]
colnames(src.al3.tracing_re@reductions$mnn_umap_fta@cell.embeddings) = c("Coord_1","Coord_2")
DimPlot(src.al3.tracing_re, reduction = "mnn_umap_fta", 
        group.by = 'cluster.v06.26.re_correct_mnn_umap_fta') +
  DimPlot(src.al3.tracing_re, reduction = "mnn_umap_fta", group.by = 'lineage')

src.al3.tracing_re.int = src.al3.tracing_re[,src.al3.tracing_re$cluster.v06.26.re_correct_mnn_umap_fta%in%c("AL.3","AL.3-Small.intestine.1",'Small.intestine.1')]
src.al3.tracing_re.pan = src.al3.tracing_re[,src.al3.tracing_re$cluster.v06.26.re_correct_mnn_umap_fta%in%c("AL.3","AL.3-EHBD/VP",'EHBD',"VP")]
src.al3.tracing_re.liv = src.al3.tracing_re[,src.al3.tracing_re$cluster.v06.26.re_correct_mnn_umap_fta%in%c("AL.3","AL.3-Liver",'Liver')]

for(i.seurat in c("src.al3.tracing_re.int", 
                  "src.al3.tracing_re.liv",
                  "src.al3.tracing_re.pan")){
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

# AL.3 pan
#-----------------

type_define = "cluster.v06.26.re_correct"
#---------------------
#-- mg3 --
cell_src.al3.tracing_re.pan_mg3 = 
  rownames(src.al3.tracing_re.pan@meta.data[unlist(src.al3.tracing_re.pan[[type_define]]) %in%c("AL.3"),])
coord_src.al3.tracing_re.pan_mg3 = src.al3.tracing_re.pan[["mnn_umap_fta"]]@cell.embeddings[cell_src.al3.tracing_re.pan_mg3,c(2,1)]
pcurve_src.al3.tracing_re.pan_mg3 = princurve::principal_curve(x = coord_src.al3.tracing_re.pan_mg3, smoother = "smooth.spline")
src.al3.tracing_re.pan$lambda_mg3 = pcurve_src.al3.tracing_re.pan_mg3$lambda[cell_src.al3.tracing_re.pan_mg3]
src.al3.tracing_re.pan@meta.data[cell_src.al3.tracing_re.pan_mg3,]$lambda_mg3 = 
  pcurve_src.al3.tracing_re.pan_mg3$lambda[cell_src.al3.tracing_re.pan_mg3]
src.al3.tracing_re.pan@meta.data[!src.al3.tracing_re.pan$lambda_mg3%in%NA,]$lambda_mg3 = 
  norm_range(src.al3.tracing_re.pan@meta.data[!src.al3.tracing_re.pan$lambda_mg3%in%NA,]$lambda_mg3)

#-- mg3m --
cell_src.al3.tracing_re.pan_mg3m = 
  rownames(src.al3.tracing_re.pan@meta.data[unlist(src.al3.tracing_re.pan[[type_define]]) %in%c("AL.3-EHBD/VP","AL.3"),])
coord_src.al3.tracing_re.pan_mg3m = src.al3.tracing_re.pan[["mnn_umap_fta"]]@cell.embeddings[cell_src.al3.tracing_re.pan_mg3m,c(2,1)]
pcurve_src.al3.tracing_re.pan_mg3m = princurve::principal_curve(x = coord_src.al3.tracing_re.pan_mg3m, smoother = "smooth.spline")
src.al3.tracing_re.pan$lambda_mg3m = pcurve_src.al3.tracing_re.pan_mg3m$lambda[cell_src.al3.tracing_re.pan_mg3m]
src.al3.tracing_re.pan@meta.data[cell_src.al3.tracing_re.pan_mg3m,]$lambda_mg3m = 
  pcurve_src.al3.tracing_re.pan_mg3m$lambda[cell_src.al3.tracing_re.pan_mg3m]
src.al3.tracing_re.pan@meta.data[!src.al3.tracing_re.pan$lambda_mg3m%in%NA,]$lambda_mg3m = 
  norm_range(src.al3.tracing_re.pan@meta.data[!src.al3.tracing_re.pan$lambda_mg3m%in%NA,]$lambda_mg3m)

#-- dp--
cell_src.al3.tracing_re.pan_dp = 
  rownames(src.al3.tracing_re.pan@meta.data[unlist(src.al3.tracing_re.pan[[type_define]]) %in%c("EHBD","AL.3-EHBD/VP"),])
coord_src.al3.tracing_re.pan_dp = src.al3.tracing_re.pan[["mnn_umap_fta"]]@cell.embeddings[cell_src.al3.tracing_re.pan_dp, c(2,1)]
pcurve_src.al3.tracing_re.pan_dp = princurve::principal_curve(x = coord_src.al3.tracing_re.pan_dp, smoother = "smooth.spline")
src.al3.tracing_re.pan$lambda_dp = pcurve_src.al3.tracing_re.pan_dp$lambda[cell_src.al3.tracing_re.pan_dp]
src.al3.tracing_re.pan@meta.data[cell_src.al3.tracing_re.pan_dp,]$lambda_dp = 
  pcurve_src.al3.tracing_re.pan_dp$lambda[cell_src.al3.tracing_re.pan_dp]
src.al3.tracing_re.pan@meta.data[!src.al3.tracing_re.pan$lambda_dp%in%NA,]$lambda_dp = 
  norm_range(src.al3.tracing_re.pan@meta.data[!src.al3.tracing_re.pan$lambda_dp%in%NA,]$lambda_dp)

#-- ep--
cell_src.al3.tracing_re.pan_ep = 
  rownames(src.al3.tracing_re.pan@meta.data[unlist(src.al3.tracing_re.pan[[type_define]]) %in%c("AL.3-EHBD/VP","VP"),])
coord_src.al3.tracing_re.pan_ep = src.al3.tracing_re.pan[["mnn_umap_fta"]]@cell.embeddings[cell_src.al3.tracing_re.pan_ep, c(2,1)]
pcurve_src.al3.tracing_re.pan_ep = princurve::principal_curve(x = coord_src.al3.tracing_re.pan_ep, smoother = "smooth.spline")
src.al3.tracing_re.pan$lambda_ep = pcurve_src.al3.tracing_re.pan_ep$lambda[cell_src.al3.tracing_re.pan_ep]
src.al3.tracing_re.pan@meta.data[cell_src.al3.tracing_re.pan_ep,]$lambda_ep = 
  pcurve_src.al3.tracing_re.pan_ep$lambda[cell_src.al3.tracing_re.pan_ep]
src.al3.tracing_re.pan@meta.data[!src.al3.tracing_re.pan$lambda_ep%in%NA,]$lambda_ep = 
  norm_range(src.al3.tracing_re.pan@meta.data[!src.al3.tracing_re.pan$lambda_ep%in%NA,]$lambda_ep)


cellorder_src.al3.tracing_re.pan = c(
  rev(intersect(names(src.al3.tracing_re.pan$lambda_mg3[order(src.al3.tracing_re.pan$lambda_mg3)]), 
                colnames(src.al3.tracing_re.pan[, unlist(src.al3.tracing_re.pan[[type_define]]) %in%"AL.3"]))),
  rev(intersect(names(src.al3.tracing_re.pan$lambda_dp[order(src.al3.tracing_re.pan$lambda_dp)]), 
                colnames(src.al3.tracing_re.pan[, unlist(src.al3.tracing_re.pan[[type_define]]) %in%"AL.3-EHBD/VP"]))),
  
  (intersect(names(src.al3.tracing_re.pan$lambda_dp[order(src.al3.tracing_re.pan$lambda_dp)]), 
                colnames(src.al3.tracing_re.pan[, unlist(src.al3.tracing_re.pan[[type_define]]) %in%"EHBD"]))),
  rev(intersect(names(src.al3.tracing_re.pan$lambda_ep[order(src.al3.tracing_re.pan$lambda_ep)]), 
             colnames(src.al3.tracing_re.pan[, unlist(src.al3.tracing_re.pan[[type_define]]) %in%"VP"])))
)
#---------------------


pdf("figure.v08.07/organ_development_re_v240115/try.al3.pan.pdf",10,10)
gene_list = unique(marker_src.al3.tracing_re.pan[
  marker_src.al3.tracing_re.pan$cluster%in%c("AL.3","AL.3-EHBD/VP"), "gene"])
cell_list = colnames(src.al3.tracing_re.pan)
#--------------------
src.al3.tracing_re.pan.rowtree  =  MyHeatmap(as.matrix(
  src.al3.tracing_re.pan@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.al3.tracing_re.pan$Time[cell_list], #[cellorder_cell_src.al3.tracing_re.pan],
               colors.time.2),
    MyName2Col(src.al3.tracing_re.pan$lineage[cell_list], #[cellorder_cell_src.al3.tracing_re.pan],
               color.lineage),
    MyName2Col(src.al3.tracing_re.pan$cluster.v06.26.re_correct_mnn_umap_fta[cell_list], #[cellorder_cell_src.al3.tracing_re.pan],
               cluster.endoderm.color.v5)
  ),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  #Rowv = "none",
  return.tree = "row",
  margins = c(10,10),
  graph = T)

src.al3.tracing_re.pan.coltree  =  MyHeatmap(as.matrix(
  src.al3.tracing_re.pan@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.al3.tracing_re.pan$Time[cell_list], #[cellorder_cell_src.al3.tracing_re.pan],
               colors.time.2),
    MyName2Col(src.al3.tracing_re.pan$lineage[cell_list], #[cellorder_cell_src.al3.tracing_re.pan],
               color.lineage),
    MyName2Col(src.al3.tracing_re.pan$cluster.v06.26.re_correct_mnn_umap_fta[cell_list], #[cellorder_cell_src.al3.tracing_re.pan],
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
tree_src.al3.tracing_re.pan.rowtree = as.dendrogram(src.al3.tracing_re.pan.rowtree)
tree_src.al3.tracing_re.pan.coltree = as.dendrogram(src.al3.tracing_re.pan.coltree)

gene_order = function(gene, type = "VP"){
  data_temp = src.al3.tracing_re.pan@assays$RNA@data[
    gene,
    src.al3.tracing_re.pan$cluster.v06.26.re_correct%in%type]
  data_temp = rowSums(data_temp)
  gene = unlist(rev(names(data_temp)[order(data_temp)]))
  return(gene)
}

gene_src.al3.tracing_re.pan.rowtree = c(
  gene_order(
    labels(tree_src.al3.tracing_re.pan.rowtree[[2]][[1]]), "AL.3"),
  # labels(tree_src.al3.tracing_re.pan.rowtree[[2]][[2]][[1]]),
  labels(tree_src.al3.tracing_re.pan.rowtree[[2]][[2]][[2]][[2]])
)
names(gene_src.al3.tracing_re.pan.rowtree) = c(
  rep(4, length(c(
    labels(tree_src.al3.tracing_re.pan.rowtree[[2]][[1]])
  ))),
  rep(3, length(c(
    # labels(tree_src.al3.tracing_re.pan.rowtree[[2]][[2]][[1]]),
    labels(tree_src.al3.tracing_re.pan.rowtree[[2]][[2]][[2]][[2]])
  )))
)

src.al3.tracing_re.pan$tree = NA
tree_src.al3.tracing_re.pan.coltree
for(i1 in c(1,2)){ for(i2 in c(1,2)){ for(i3 in c(1,2)){ for(i4 in c(1,2)){
        src.al3.tracing_re.pan@meta.data[
          labels(tree_src.al3.tracing_re.pan.coltree[[i1]][[i2]][[i3]][[i4]]), ]$tree = paste(i1,i2,i3,i4,sep=".")
  }}}
}
DimPlot(src.al3.tracing_re.pan, reduction = "mnn_umap_fta", group.by = "tree", label = T) +
  DimPlot(src.al3.tracing_re.pan, reduction = "mnn_umap_fta", group.by = "cluster.v06.26.re_correct_mnn_umap_fta")
src.al3.tracing_re.pan$cluster.v06.26.re_correct = 
  src.al3.tracing_re.pan$cluster.v06.26.re_correct_mnn_umap_fta
src.al3.tracing_re.pan@meta.data[
  src.al3.tracing_re.pan$tree%in%c("2.1.2.1","2.1.2.2","2.1.1.1","2.1.1.2","2.2.2.1") &
    src.al3.tracing_re.pan$Time%in%c("9ss","12ss"),]$cluster.v06.26.re_correct = "AL.3"


pdf("figure.v08.07/organ_development_re_v240115/try.al3.pan.pdf",10,10)
gene_list = setdiff(
  unique(marker_src.al3.tracing_re.pan[
    marker_src.al3.tracing_re.pan$cluster%in%c("EHBD0","VP"), "gene"]),
  labels(tree_src.al3.tracing_re.pan.rowtree))
cell_list = colnames(src.al3.tracing_re.pan)
#--------------------
src.al3.tracing_re.pan.rowtree1  =  MyHeatmap(as.matrix(
  src.al3.tracing_re.pan@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.al3.tracing_re.pan$Time[cell_list], #[cellorder_cell_src.al3.tracing_re.pan],
               colors.time.2),
    MyName2Col(src.al3.tracing_re.pan$lineage[cell_list], #[cellorder_cell_src.al3.tracing_re.pan],
               color.lineage),
    MyName2Col(src.al3.tracing_re.pan$cluster.v06.26.re_correct_mnn_umap_fta[cell_list], #[cellorder_cell_src.al3.tracing_re.pan],
               cluster.endoderm.color.v5),
    MyName2Col(src.al3.tracing_re.pan$cluster.v06.26.re_correct[cell_list], #[cellorder_cell_src.al3.tracing_re.pan],
               cluster.endoderm.color.v5)
  ),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  #Rowv = "none",
  return.tree = "row",
  margins = c(10,10),
  graph = T)

src.al3.tracing_re.pan.coltree1  =  MyHeatmap(as.matrix(
  src.al3.tracing_re.pan@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.al3.tracing_re.pan$Time[cell_list], #[cellorder_cell_src.al3.tracing_re.pan],
               colors.time.2),
    MyName2Col(src.al3.tracing_re.pan$lineage[cell_list], #[cellorder_cell_src.al3.tracing_re.pan],
               color.lineage),
    MyName2Col(src.al3.tracing_re.pan$cluster.v06.26.re_correct_mnn_umap_fta[cell_list], #[cellorder_cell_src.al3.tracing_re.pan],
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
tree_src.al3.tracing_re.pan.rowtree1 = as.dendrogram(src.al3.tracing_re.pan.rowtree1)
tree_src.al3.tracing_re.pan.coltree1 = as.dendrogram(src.al3.tracing_re.pan.coltree1)


gene_src.al3.tracing_re.pan.rowtree1 = c(
  # labels(tree_src.al3.tracing_re.pan.rowtree[[1]][[1]]),
  labels(tree_src.al3.tracing_re.pan.rowtree1[[1]]),
  gene_order(labels(tree_src.al3.tracing_re.pan.rowtree1[[2]]))
)
names(gene_src.al3.tracing_re.pan.rowtree1) = c(
  rep(5, length(c(
    labels(tree_src.al3.tracing_re.pan.rowtree1[[1]])
  ))),
  rep(7, length(c(
    labels(tree_src.al3.tracing_re.pan.rowtree1[[2]])
  )))
)

gene_src.al3.tracing_re.pan.rowtree2 = c(
  gene_src.al3.tracing_re.pan.rowtree,
  gene_src.al3.tracing_re.pan.rowtree1
)


#--- Based on new marker
src.al3.tracing_re.pan = SetIdent(src.al3.tracing_re.pan,
                                  value = src.al3.tracing_re.pan$cluster.v06.26.re_correct)
marker_src.al3.tracing_re.pan.re = FindAllMarkers(src.al3.tracing_re.pan)
pdf("figure.v08.07/organ_development_re_v240115/try.al3.pan.pdf",10,10)
gene_list = unique(marker_src.al3.tracing_re.pan.re[
  marker_src.al3.tracing_re.pan.re$cluster%in%c("AL.3","AL.3-EHBD/VP"), "gene"])
cell_list = colnames(src.al3.tracing_re.pan)
#--------------------
src.al3.tracing_re.pan.rowtree3  =  MyHeatmap(as.matrix(
  src.al3.tracing_re.pan@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.al3.tracing_re.pan$Time[cell_list], #[cellorder_cell_src.al3.tracing_re.pan],
               colors.time.2),
    MyName2Col(src.al3.tracing_re.pan$lineage[cell_list], #[cellorder_cell_src.al3.tracing_re.pan],
               color.lineage),
    MyName2Col(src.al3.tracing_re.pan$cluster.v06.26.re_correct_mnn_umap_fta[cell_list], #[cellorder_cell_src.al3.tracing_re.pan],
               cluster.endoderm.color.v5)
  ),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  #Rowv = "none",
  return.tree = "row",
  margins = c(10,10),
  graph = T)

src.al3.tracing_re.pan.coltree3  =  MyHeatmap(as.matrix(
  src.al3.tracing_re.pan@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.al3.tracing_re.pan$Time[cell_list], #[cellorder_cell_src.al3.tracing_re.pan],
               colors.time.2),
    MyName2Col(src.al3.tracing_re.pan$lineage[cell_list], #[cellorder_cell_src.al3.tracing_re.pan],
               color.lineage),
    MyName2Col(src.al3.tracing_re.pan$cluster.v06.26.re_correct_mnn_umap_fta[cell_list], #[cellorder_cell_src.al3.tracing_re.pan],
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
tree_src.al3.tracing_re.pan.rowtree3 = as.dendrogram(src.al3.tracing_re.pan.rowtree3)
tree_src.al3.tracing_re.pan.coltree3 = as.dendrogram(src.al3.tracing_re.pan.coltree3)

gene_src.al3.tracing_re.pan.rowtree3 = c(
  rev(labels(tree_src.al3.tracing_re.pan.rowtree3[[2]][[2]][[1]][[1]][[1]])),
  # labels(tree_src.al3.tracing_re.pan.rowtree[[2]][[2]][[1]]),
  rev(labels(tree_src.al3.tracing_re.pan.rowtree3[[2]][[1]]))
)
names(gene_src.al3.tracing_re.pan.rowtree3) = c(
  rep(4, length(c(
    labels(tree_src.al3.tracing_re.pan.rowtree3[[2]][[2]][[1]][[1]][[1]])
  ))),
  rep(3, length(c(
    # labels(tree_src.al3.tracing_re.pan.rowtree[[2]][[2]][[1]]),
    labels(tree_src.al3.tracing_re.pan.rowtree3[[2]][[1]])
  )))
)

gene_src.al3.tracing_re.pan.rowtree3.fin = c(
  gene_src.al3.tracing_re.pan.rowtree3,
  gene_src.al3.tracing_re.pan.rowtree1
)

pdf("figure.v08.07/organ_development_re_v240115/try.al3.pan.pdf",10,10)
gene_list = gene_src.al3.tracing_re.pan.rowtree3.fin
cell_list = cellorder_src.al3.tracing_re.pan
#--------------------
src.al3.tracing_re.pan.rowtree1  =  MyHeatmap(as.matrix(
  src.al3.tracing_re.pan@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.al3.tracing_re.pan$Time[cell_list], #[cellorder_cell_src.al3.tracing_re.pan],
               colors.time.2),
    MyName2Col(src.al3.tracing_re.pan$lineage[cell_list], #[cellorder_cell_src.al3.tracing_re.pan],
               color.lineage),
    MyName2Col(src.al3.tracing_re.pan$cluster.v06.26.re_correct[cell_list], #[cellorder_cell_src.al3.tracing_re.pan],
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

save(cellorder_src.al3.tracing_re.pan,
     src.al3.tracing_re.pan,
     markergene_src.al3.tracing_re.pan,
     marker_src.al3.tracing_re.pan,
     marker_src.al3.tracing_re.pan.re,
     gene_src.al3.tracing_re.pan.rowtree,
     gene_src.al3.tracing_re.pan.rowtree1,
     gene_src.al3.tracing_re.pan.rowtree2,
     gene_src.al3.tracing_re.pan.rowtree3,
     gene_src.al3.tracing_re.pan.rowtree3.fin,
     file = "figure.v08.07/organ_development_re_v240115/src.al3.tracing_re.pan.parameter.Rdata")


# for liv and sm1: after AL.3 correction
#========================================================
src.al3.tracing_re$cluster.v06.26.re_correct = 
  src.al3.tracing_re$cluster.v06.26.re_correct_mnn_umap_fta
src.al3.tracing_re@meta.data[colnames(src.al3.tracing_re.pan),]$cluster.v06.26.re_correct = 
  src.al3.tracing_re.pan$cluster.v06.26.re_correct

src.al3.tracing_re.int = src.al3.tracing_re[,src.al3.tracing_re$cluster.v06.26.re_correct%in%c("AL.3","AL.3-Small.intestine.1",'Small.intestine.1')]
src.al3.tracing_re.pan = src.al3.tracing_re[,src.al3.tracing_re$cluster.v06.26.re_correct_mnn_umap_fta%in%c("AL.3","AL.3-EHBD/VP",'EHBD',"VP")]
src.al3.tracing_re.liv = src.al3.tracing_re[,src.al3.tracing_re$cluster.v06.26.re_correct%in%c("AL.3","AL.3-Liver",'Liver')]

for(i.seurat in c(
  "src.al3.tracing_re.int",
  "src.al3.tracing_re.liv",
  "src.al3.tracing_re.pan"
  )){
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


type_define = "cluster.v06.26.re_correct.refine"
#---------------------
#-- mg3 --
cell_src.al3.tracing_re.int_mg3 = 
  rownames(src.al3.tracing_re.int@meta.data[unlist(src.al3.tracing_re.int[[type_define]]) %in%c("AL.3"),])
coord_src.al3.tracing_re.int_mg3 = src.al3.tracing_re.int[["mnn_umap_fta"]]@cell.embeddings[cell_src.al3.tracing_re.int_mg3,c(2,1)]
pcurve_src.al3.tracing_re.int_mg3 = princurve::principal_curve(x = coord_src.al3.tracing_re.int_mg3, smoother = "smooth.spline")
src.al3.tracing_re.int$lambda_mg3 = pcurve_src.al3.tracing_re.int_mg3$lambda[cell_src.al3.tracing_re.int_mg3]
src.al3.tracing_re.int@meta.data[cell_src.al3.tracing_re.int_mg3,]$lambda_mg3 = 
  pcurve_src.al3.tracing_re.int_mg3$lambda[cell_src.al3.tracing_re.int_mg3]
src.al3.tracing_re.int@meta.data[!src.al3.tracing_re.int$lambda_mg3%in%NA,]$lambda_mg3 = 
  norm_range(src.al3.tracing_re.int@meta.data[!src.al3.tracing_re.int$lambda_mg3%in%NA,]$lambda_mg3)

#-- mg3m --
cell_src.al3.tracing_re.int_mg3m = 
  rownames(src.al3.tracing_re.int@meta.data[unlist(src.al3.tracing_re.int[[type_define]]) %in%c("AL.3-Small.intestine.1","Small.intestine.1"),])
coord_src.al3.tracing_re.int_mg3m = src.al3.tracing_re.int[["mnn_umap_fta"]]@cell.embeddings[cell_src.al3.tracing_re.int_mg3m, c(2,1)]
pcurve_src.al3.tracing_re.int_mg3m = princurve::principal_curve(x = coord_src.al3.tracing_re.int_mg3m, smoother = "smooth.spline")
src.al3.tracing_re.int$lambda_mg3m = pcurve_src.al3.tracing_re.int_mg3m$lambda[cell_src.al3.tracing_re.int_mg3m]
src.al3.tracing_re.int@meta.data[cell_src.al3.tracing_re.int_mg3m,]$lambda_mg3m = 
  pcurve_src.al3.tracing_re.int_mg3m$lambda[cell_src.al3.tracing_re.int_mg3m]
src.al3.tracing_re.int@meta.data[!src.al3.tracing_re.int$lambda_mg3m%in%NA,]$lambda_mg3m = 
  norm_range(src.al3.tracing_re.int@meta.data[!src.al3.tracing_re.int$lambda_mg3m%in%NA,]$lambda_mg3m)

#-- dp--
cell_src.al3.tracing_re.int_dp = 
  rownames(src.al3.tracing_re.int@meta.data[unlist(src.al3.tracing_re.int[[type_define]]) %in%c("Small.intestine.1"),])
coord_src.al3.tracing_re.int_dp = src.al3.tracing_re.int[["mnn_umap_fta"]]@cell.embeddings[cell_src.al3.tracing_re.int_dp, c(2,1)]
pcurve_src.al3.tracing_re.int_dp = princurve::principal_curve(x = coord_src.al3.tracing_re.int_dp, smoother = "smooth.spline")
src.al3.tracing_re.int$lambda_dp = pcurve_src.al3.tracing_re.int_dp$lambda[cell_src.al3.tracing_re.int_dp]
src.al3.tracing_re.int@meta.data[cell_src.al3.tracing_re.int_dp,]$lambda_dp = 
  pcurve_src.al3.tracing_re.int_dp$lambda[cell_src.al3.tracing_re.int_dp]
src.al3.tracing_re.int@meta.data[!src.al3.tracing_re.int$lambda_dp%in%NA,]$lambda_dp = 
  norm_range(src.al3.tracing_re.int@meta.data[!src.al3.tracing_re.int$lambda_dp%in%NA,]$lambda_dp)



cellorder_src.al3.tracing_re.int = c(
  rev(intersect(names(src.al3.tracing_re.int$lambda_mg3[order(src.al3.tracing_re.int$lambda_mg3)]), 
                colnames(src.al3.tracing_re.int[, unlist(src.al3.tracing_re.int[[type_define]]) %in%"AL.3"]))),
  rev(intersect(names(src.al3.tracing_re.int$lambda_mg3m[order(src.al3.tracing_re.int$lambda_mg3m)]), 
             colnames(src.al3.tracing_re.int[, unlist(src.al3.tracing_re.int[[type_define]]) %in%"AL.3-Small.intestine.1"]))),
  (intersect(names(src.al3.tracing_re.int$lambda_dp[order(src.al3.tracing_re.int$lambda_dp)]), 
                colnames(src.al3.tracing_re.int[, unlist(src.al3.tracing_re.int[[type_define]]) %in%"Small.intestine.1"])))
)
#---------------------

pdf("figure.v08.07/organ_development_re_v240115/try.al3.int.pdf",10,10)
gene_list = unique(c(
  markergene_src.al3.tracing_re.int
  #src.al3.tracing_re.int.filtergene
))
cell_list = colnames(src.al3.tracing_re.int)
#--------------------
src.al3.tracing_re.int.rowtree  =  MyHeatmap(as.matrix(
  src.al3.tracing_re.int@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.al3.tracing_re.int$Time[cell_list], #[cellorder_cell_src.al3.tracing_re.int],
               colors.time.2),
    MyName2Col(src.al3.tracing_re.int$lineage[cell_list], #[cellorder_cell_src.al3.tracing_re.int],
               color.lineage),
    MyName2Col(src.al3.tracing_re.int$cluster.v06.26.re_correct_mnn_umap_fta[cell_list], #[cellorder_cell_src.al3.tracing_re.int],
               cluster.endoderm.color.v5)
  ),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  #Rowv = "none",
  return.tree = "row",
  margins = c(10,10),
  graph = T)

src.al3.tracing_re.int.coltree  =  MyHeatmap(as.matrix(
  src.al3.tracing_re.int@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.al3.tracing_re.int$Time[cell_list], #[cellorder_cell_src.al3.tracing_re.int],
               colors.time.2),
    MyName2Col(src.al3.tracing_re.int$lineage[cell_list], #[cellorder_cell_src.al3.tracing_re.int],
               color.lineage),
    MyName2Col(src.al3.tracing_re.int$cluster.v06.26.re_correct_mnn_umap_fta[cell_list], #[cellorder_cell_src.al3.tracing_re.int],
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
tree_src.al3.tracing_re.int.rowtree = as.dendrogram(src.al3.tracing_re.int.rowtree)
tree_src.al3.tracing_re.int.coltree = as.dendrogram(src.al3.tracing_re.int.coltree)

gene_src.al3.tracing_re.int.rowtree = c(
  labels(tree_src.al3.tracing_re.int.rowtree[[1]][[1]][[1]]),
  labels(tree_src.al3.tracing_re.int.rowtree[[1]][[1]][[2]][[1]]),
  labels(tree_src.al3.tracing_re.int.rowtree[[1]][[2]][[2]][[2]][[1]]),
  labels(tree_src.al3.tracing_re.int.rowtree[[2]][[2]][[2]][[2]][[2]][[2]]),
  
  labels(tree_src.al3.tracing_re.int.rowtree[[2]][[2]][[2]][[2]][[2]][[1]]),
  # labels(tree_src.al3.tracing_re.int.rowtree[[2]][[2]][[1]]),
  labels(tree_src.al3.tracing_re.int.rowtree[[2]][[2]][[1]]),
  labels(tree_src.al3.tracing_re.int.rowtree[[2]][[2]][[2]][[1]][[2]])
)
names(gene_src.al3.tracing_re.int.rowtree) = c(
  rep(4, length(c(
    labels(tree_src.al3.tracing_re.int.rowtree[[1]][[1]][[1]]),
    labels(tree_src.al3.tracing_re.int.rowtree[[1]][[1]][[2]][[1]])
  ))),
  rep(3, length(c(
    labels(tree_src.al3.tracing_re.int.rowtree[[1]][[2]][[2]][[2]][[1]])
  ))),
  rep(5, length(c(
    labels(tree_src.al3.tracing_re.int.rowtree[[2]][[2]][[2]][[2]][[2]][[2]])
  ))),
  rep(7, length(c(
    labels(tree_src.al3.tracing_re.int.rowtree[[2]][[2]][[2]][[2]][[2]][[1]]),
    labels(tree_src.al3.tracing_re.int.rowtree[[2]][[2]][[1]]),
    labels(tree_src.al3.tracing_re.int.rowtree[[2]][[2]][[2]][[1]][[2]])
  )))
)

src.al3.tracing_re.int$tree = NA
for(i1 in c(1,2)){ for(i2 in c(1,2)){ for(i3 in c(1,2)){ 
  src.al3.tracing_re.int@meta.data[
    labels(tree_src.al3.tracing_re.int.coltree[[i1]][[i2]][[i3]]), ]$tree = paste(i1,i2,i3,sep=".")
}}}

DimPlot(src.al3.tracing_re.int, reduction = "mnn_umap_fta", group.by = "tree", label = T) +
  DimPlot(src.al3.tracing_re.int, reduction = "mnn_umap_fta", group.by = "cluster.v06.26.re_correct.refine")


src.al3.tracing_re.int$cluster.v06.26.re_correct = src.al3.tracing_re$cluster.v06.26.re_correct[colnames(src.al3.tracing_re.int)]
src.al3.tracing_re.int$cluster.v06.26.re_correct.refine = src.al3.tracing_re.int$cluster.v06.26.re_correct
src.al3.tracing_re.int@meta.data[
  src.al3.tracing_re.int$tree%in%c("2.2.2"),]$cluster.v06.26.re_correct.refine = "AL.3-Small.intestine.1"

src.al3.tracing_re$cluster.v06.26.re_correct.refine = src.al3.tracing_re$cluster.v06.26.re_correct
src.al3.tracing_re$cluster.v06.26.re_correct.refine[colnames(src.al3.tracing_re.int)] = 
  src.al3.tracing_re.int$cluster.v06.26.re_correct.refine[colnames(src.al3.tracing_re.int)]


pdf("figure.v08.07/organ_development_re_v240115/try.al3.int.pdf",10,10)
gene_list = (gene_src.al3.tracing_re.int.rowtree)
cell_list = cellorder_src.al3.tracing_re.int
#--------------------
src.al3.tracing_re.int.rowtree1  =  MyHeatmap(as.matrix(
  src.al3.tracing_re.int@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.al3.tracing_re.int$Time[cell_list], #[cellorder_cell_src.al3.tracing_re.int],
               colors.time.2),
    MyName2Col(src.al3.tracing_re.int$lineage[cell_list], #[cellorder_cell_src.al3.tracing_re.int],
               color.lineage),
    MyName2Col(src.al3.tracing_re.int$cluster.v06.26.re_correct.refine[cell_list], #[cellorder_cell_src.al3.tracing_re.int],
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

save(cellorder_src.al3.tracing_re.int,
     src.al3.tracing_re.int,
     src.al3.tracing_re,
     markergene_src.al3.tracing_re.int,
     marker_src.al3.tracing_re.int,
     gene_src.al3.tracing_re.pan.rowtree,
     file = "figure.v08.07/organ_development_re_v240115/src.al3.tracing_re.int.parameter.Rdata")




# AL.3 int
#-----------------
# for liv and sm1: after AL.3 correction
#========================================================
src.al3.tracing_re.liv = src.al3.tracing_re[,src.al3.tracing_re$cluster.v06.26.re_correct%in%c("AL.3","AL.3-Liver","AL.1/2-Liver",'Liver')]

for(i.seurat in c(#"src.al3.tracing_re.int" #, 
                  "src.al3.tracing_re.liv"
                  #"src.al3.tracing_re.pan"
)){
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

type_define = "cluster.v06.26.re_correct.refine"
#---------------------
#-- mg3 --
cell_src.al3.tracing_re.liv_mg3 = 
  rownames(src.al3.tracing_re.liv@meta.data[unlist(src.al3.tracing_re.liv[[type_define]]) %in%c("AL.3"),])
coord_src.al3.tracing_re.liv_mg3 = src.al3.tracing_re.liv[["mnn_umap_fta"]]@cell.embeddings[cell_src.al3.tracing_re.liv_mg3,c(2,1)]
pcurve_src.al3.tracing_re.liv_mg3 = princurve::principal_curve(x = coord_src.al3.tracing_re.liv_mg3, smoother = "smooth.spline")
src.al3.tracing_re.liv$lambda_mg3 = pcurve_src.al3.tracing_re.liv_mg3$lambda[cell_src.al3.tracing_re.liv_mg3]
src.al3.tracing_re.liv@meta.data[cell_src.al3.tracing_re.liv_mg3,]$lambda_mg3 = 
  pcurve_src.al3.tracing_re.liv_mg3$lambda[cell_src.al3.tracing_re.liv_mg3]
src.al3.tracing_re.liv@meta.data[!src.al3.tracing_re.liv$lambda_mg3%in%NA,]$lambda_mg3 = 
  norm_range(src.al3.tracing_re.liv@meta.data[!src.al3.tracing_re.liv$lambda_mg3%in%NA,]$lambda_mg3)

#-- mg3m --
cell_src.al3.tracing_re.liv_mg3m = 
  rownames(src.al3.tracing_re.liv@meta.data[unlist(src.al3.tracing_re.liv[[type_define]]) %in%c("AL.3","AL.3-Liver","AL.1/2-Liver"),])
coord_src.al3.tracing_re.liv_mg3m = src.al3.tracing_re.liv[["mnn_umap_fta"]]@cell.embeddings[cell_src.al3.tracing_re.liv_mg3m, c(2,1)]
pcurve_src.al3.tracing_re.liv_mg3m = princurve::principal_curve(x = coord_src.al3.tracing_re.liv_mg3m, smoother = "smooth.spline")
src.al3.tracing_re.liv$lambda_mg3m = pcurve_src.al3.tracing_re.liv_mg3m$lambda[cell_src.al3.tracing_re.liv_mg3m]
src.al3.tracing_re.liv@meta.data[cell_src.al3.tracing_re.liv_mg3m,]$lambda_mg3m = 
  pcurve_src.al3.tracing_re.liv_mg3m$lambda[cell_src.al3.tracing_re.liv_mg3m]
src.al3.tracing_re.liv@meta.data[!src.al3.tracing_re.liv$lambda_mg3m%in%NA,]$lambda_mg3m = 
  norm_range(src.al3.tracing_re.liv@meta.data[!src.al3.tracing_re.liv$lambda_mg3m%in%NA,]$lambda_mg3m)

#-- dp--
cell_src.al3.tracing_re.liv_dp = 
  rownames(src.al3.tracing_re.liv@meta.data[unlist(src.al3.tracing_re.liv[[type_define]]) %in%c("AL.1/2-Liver","Liver"),])
coord_src.al3.tracing_re.liv_dp = src.al3.tracing_re.liv[["mnn_umap_fta"]]@cell.embeddings[cell_src.al3.tracing_re.liv_dp, c(2,1)]
pcurve_src.al3.tracing_re.liv_dp = princurve::principal_curve(x = coord_src.al3.tracing_re.liv_dp, smoother = "smooth.spline")
src.al3.tracing_re.liv$lambda_dp = pcurve_src.al3.tracing_re.liv_dp$lambda[cell_src.al3.tracing_re.liv_dp]
src.al3.tracing_re.liv@meta.data[cell_src.al3.tracing_re.liv_dp,]$lambda_dp = 
  pcurve_src.al3.tracing_re.liv_dp$lambda[cell_src.al3.tracing_re.liv_dp]
src.al3.tracing_re.liv@meta.data[!src.al3.tracing_re.liv$lambda_dp%in%NA,]$lambda_dp = 
  norm_range(src.al3.tracing_re.liv@meta.data[!src.al3.tracing_re.liv$lambda_dp%in%NA,]$lambda_dp)



cellorder_src.al3.tracing_re.liv = c(
  (intersect(names(src.al3.tracing_re.liv$lambda_mg3[order(src.al3.tracing_re.liv$lambda_mg3)]), 
                colnames(src.al3.tracing_re.liv[, unlist(src.al3.tracing_re.liv[[type_define]]) %in%"AL.3"]))),
  (intersect(names(src.al3.tracing_re.liv$lambda_mg3m[order(src.al3.tracing_re.liv$lambda_mg3m)]), 
                colnames(src.al3.tracing_re.liv[, unlist(src.al3.tracing_re.liv[[type_define]]) %in%"AL.3-Liver"]))),
  (intersect(names(src.al3.tracing_re.liv$lambda_dp[order(src.al3.tracing_re.liv$lambda_dp)]), 
                colnames(src.al3.tracing_re.liv[, unlist(src.al3.tracing_re.liv[[type_define]]) %in%"AL.1/2-Liver"]))),
  (intersect(names(src.al3.tracing_re.liv$lambda_dp[order(src.al3.tracing_re.liv$lambda_dp)]), 
             colnames(src.al3.tracing_re.liv[, unlist(src.al3.tracing_re.liv[[type_define]]) %in%"Liver"])))
)
#---------------------

pdf("figure.v08.07/organ_development_re_v240115/try.al3.liv.pdf",10,10)
gene_list = unique(
  marker_src.al3.tracing_re.liv[marker_src.al3.tracing_re.liv$cluster%in%c("AL.3","AL.3-Liver"),"gene"]
  #src.al3.tracing_re.liv.filtergene
)
cell_list = colnames(src.al3.tracing_re.liv) # [,src.al3.tracing_re.liv$cluster.v06.26.re_correct%in%c("AL.3","AL.3-Liver")])
#--------------------
src.al3.tracing_re.liv.rowtree  =  MyHeatmap(as.matrix(
  src.al3.tracing_re.liv@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.al3.tracing_re.liv$Time[cell_list], #[cellorder_cell_src.al3.tracing_re.liv],
               colors.time.2),
    MyName2Col(src.al3.tracing_re.liv$lineage[cell_list], #[cellorder_cell_src.al3.tracing_re.liv],
               color.lineage),
    MyName2Col(src.al3.tracing_re.liv$cluster.v06.26.re_correct[cell_list], #[cellorder_cell_src.al3.tracing_re.liv],
               cluster.endoderm.color.v5)
  ),
  ColSideColorsSize = 4,
  c.hc.method = "ward.D",
  r.hc.method = "ward.D2",
  #Rowv = "none",
  return.tree = "row",
  margins = c(10,10),
  graph = T)

src.al3.tracing_re.liv.coltree  =  MyHeatmap(as.matrix(
  src.al3.tracing_re.liv@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.al3.tracing_re.liv$Time[cell_list], #[cellorder_cell_src.al3.tracing_re.liv],
               colors.time.2),
    MyName2Col(src.al3.tracing_re.liv$lineage[cell_list], #[cellorder_cell_src.al3.tracing_re.liv],
               color.lineage),
    MyName2Col(src.al3.tracing_re.liv$cluster.v06.26.re_correct[cell_list], #[cellorder_cell_src.al3.tracing_re.liv],
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

tree_src.al3.tracing_re.liv.rowtree = as.dendrogram(src.al3.tracing_re.liv.rowtree)
tree_src.al3.tracing_re.liv.coltree = as.dendrogram(src.al3.tracing_re.liv.coltree)

gene_src.al3.tracing_re.liv.rowtree = c(
  labels(tree_src.al3.tracing_re.liv.rowtree[[2]][[2]][[1]][[2]]),
  
  rev(labels(tree_src.al3.tracing_re.liv.rowtree[[2]][[1]][[2]][[2]])),
  
  # labels(tree_src.al3.tracing_re.liv.rowtree[[1]][[1]]),
  labels(tree_src.al3.tracing_re.liv.rowtree[[1]][[2]][[1]][[1]]),
  labels(tree_src.al3.tracing_re.liv.rowtree[[1]][[2]][[2]][[2]][[1]][[2]]),
  labels(tree_src.al3.tracing_re.liv.rowtree[[1]][[2]][[2]][[2]][[2]])
)
names(gene_src.al3.tracing_re.liv.rowtree) = c(
  rep(4, length(c(
    labels(tree_src.al3.tracing_re.liv.rowtree[[2]][[2]][[1]][[2]])
  ))),
  rep(3, length(c(
    labels(tree_src.al3.tracing_re.liv.rowtree[[2]][[1]][[2]][[2]])
  ))),
  rep(7, length(c(
    # labels(tree_src.al3.tracing_re.liv.rowtree[[1]][[1]]),
    labels(tree_src.al3.tracing_re.liv.rowtree[[1]][[2]][[1]][[1]]),
    labels(tree_src.al3.tracing_re.liv.rowtree[[1]][[2]][[2]][[2]][[1]][[2]]),
    labels(tree_src.al3.tracing_re.liv.rowtree[[1]][[2]][[2]][[2]][[2]])
  )))
)

src.al3.tracing_re.liv$tree = NA
for(i1 in c(1,2)){ for(i2 in c(1,2)){ for(i3 in c(1,2)){ 
  src.al3.tracing_re.liv@meta.data[
    labels(tree_src.al3.tracing_re.liv.coltree[[i1]][[i2]][[i3]]), ]$tree = paste(i1,i2,i3,sep=".")
}}}

DimPlot(src.al3.tracing_re.liv, reduction = "mnn_umap_fta", group.by = "tree", label = T) +
  DimPlot(src.al3.tracing_re.liv, reduction = "mnn_umap_fta", group.by = "cluster.v06.26.re_correct.refine")

# src.al3.tracing_re.liv$cluster.v06.26.re_correct = src.al3.tracing_re$cluster.v06.26.re_correct[colnames(src.al3.tracing_re.liv)]
src.al3.tracing_re.liv$cluster.v06.26.re_correct.refine = src.al3.tracing_re.liv$cluster.v06.26.re_correct
src.al3.tracing_re.liv@meta.data[
  src.al3.tracing_re.liv$tree%in%c("2.1.1"),]$cluster.v06.26.re_correct.refine = "AL.3-Liver"
src.al3.tracing_re.liv@meta.data[
  src.al3.tracing_re.liv$Time%in%"9ss",]$cluster.v06.26.re_correct.refine = "AL.3"
src.al3.tracing_re.liv@meta.data[
  src.al3.tracing_re.liv$Time%in%c("24ss","27ss"),]$cluster.v06.26.re_correct.refine = "Liver"

# src.al3.tracing_re$cluster.v06.26.re_correct.refine = src.al3.tracing_re$cluster.v06.26.re_correct
src.al3.tracing_re$cluster.v06.26.re_correct.refine[colnames(src.al3.tracing_re.liv)] = 
  src.al3.tracing_re.liv$cluster.v06.26.re_correct.refine[colnames(src.al3.tracing_re.liv)]


pdf("figure.v08.07/organ_development_re_v240115/try.al3.liv.pdf",10,10)
gene_list = (gene_src.al3.tracing_re.liv.rowtree)
cell_list = cellorder_src.al3.tracing_re.liv
#--------------------
src.al3.tracing_re.liv.rowtree1  =  MyHeatmap(as.matrix(
  src.al3.tracing_re.liv@assays$RNA@data[
    gene_list, cell_list]),
  type = "row.relat",
  hc.c.data.type = "row.relat",
  hc.r.data.type = "row.relat",
  c.cov.method = "s",
  r.cov.method = "s",
  ColSideColors = cbind(
    MyName2Col(src.al3.tracing_re.liv$Time[cell_list], #[cellorder_cell_src.al3.tracing_re.liv],
               colors.time.2),
    MyName2Col(src.al3.tracing_re.liv$lineage[cell_list], #[cellorder_cell_src.al3.tracing_re.liv],
               color.lineage),
    MyName2Col(src.al3.tracing_re.liv$cluster.v06.26.re_correct.refine[cell_list], #[cellorder_cell_src.al3.tracing_re.liv],
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

save(cellorder_src.al3.tracing_re.liv,
     src.al3.tracing_re.liv,
     src.al3.tracing_re,
     markergene_src.al3.tracing_re.liv,
     marker_src.al3.tracing_re.liv,
     gene_src.al3.tracing_re.pan.rowtree,
     file = "figure.v08.07/organ_development_re_v240115/src.al3.tracing_re.liv.parameter.Rdata")












