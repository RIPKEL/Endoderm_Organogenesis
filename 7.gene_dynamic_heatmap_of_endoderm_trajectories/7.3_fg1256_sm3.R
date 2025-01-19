
#-----------
# FG.1
#-----------
# Correction for 9-SS
src.fg1.tracing@meta.data[
  src.fg1.tracing$cluster.extract.v1.1_define%in%"FG.1"&
    src.fg1.tracing$Time%in%"9ss",]$cluster.v06.26.re_mnn_umap_fta = "FG.1"
src.fg1.tracing.filtergene = 
  Myfilter(as.matrix(src.fg1.tracing@assays$RNA@data),
           gene = src.fg1.tracing@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)

pdf("figure.v08.07/organ_development_re_v240115/try.fg1.pdf",9,7)
src.fg1.tracing.rowtree =
  MyHeatmap(as.matrix(src.fg1.tracing@assays$RNA@data[
    src.fg1.tracing.filtergene,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg1.tracing$Time,
                 colors.time.2),
      MyName2Col(src.fg1.tracing$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg1.tracing$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg1.tracing$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()




#-----------
# FG.6
#-----------
# Correction for 9-SS
src.fg6.tracing@meta.data[
  src.fg6.tracing$cluster.extract.v1.1_define%in%"FG.6"&
    src.fg6.tracing$Time%in%"9ss",]$cluster.v06.26.re_mnn_umap_fta = "FG.6"
src.fg6.tracing.filtergene = 
  Myfilter(as.matrix(src.fg6.tracing@assays$RNA@data),
           gene = src.fg6.tracing@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)

pdf("figure.v08.07/organ_development_re_v240115/try.fg6.pdf",9,7)
src.fg6.tracing.rowtree =
  MyHeatmap(as.matrix(src.fg6.tracing@assays$RNA@data[
    src.fg6.tracing.filtergene,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg6.tracing$Time,
                 colors.time.2),
      MyName2Col(src.fg6.tracing$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg6.tracing$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg6.tracing$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()


#-----------
# FG.1-6
#-----------
src.fg1_6.tracing = merge(src.fg1.tracing,src.fg6.tracing)

src.fg1_6.tracing = FindVariableFeatures(src.fg1_6.tracing)
src.fg1_6.tracing = ScaleData(src.fg1_6.tracing,
                               rownames(src.fg1_6.tracing))
src.fg1_6.tracing.filtergene = 
  Myfilter(as.matrix(src.fg1_6.tracing@assays$RNA@data),
           gene = src.fg1_6.tracing@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)

pdf("figure.v08.07/organ_development_re_v240115/try.fg1_6.pdf",9,7)
src.fg1_6.tracing.colree =
  MyHeatmap(as.matrix(src.fg1_6.tracing@assays$RNA@data[
    src.fg1_6.tracing.filtergene,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg1_6.tracing$Time,
                 colors.time.2),
      MyName2Col(src.fg1_6.tracing$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg1_6.tracing$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg1_6.tracing$lineage,
                 color.lineage),
      MyName2Col(src.fg1_6.tracing$tree,
                 color.temp)
      
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "col",
    graph = T)
dev.off()

tree_src.fg1_6.tracing.colree =
  as.dendrogram(src.fg1_6.tracing.colree )

color.temp = c("#436B8D", "#807D80", "#C66A7A", "#C3F1E6", "#CFFE85",
               "#6AD4C2", "#D8B694", "#A8C2D5", "#82A843", "#807D08")
names(color.temp) = c(1:10)

src.fg1_6.tracing$tree = NA
src.fg1_6.tracing@meta.data[
  labels(tree_src.fg1_6.tracing.colree[[1]][[1]]),]$tree = 1
src.fg1_6.tracing@meta.data[
  labels(tree_src.fg1_6.tracing.colree[[1]][[2]]),]$tree = 2
src.fg1_6.tracing@meta.data[
  labels(tree_src.fg1_6.tracing.colree[[2]][[1]]),]$tree = 3
src.fg1_6.tracing@meta.data[
  labels(tree_src.fg1_6.tracing.colree[[2]][[2]]),]$tree = 4
src.fg1_6.tracing@meta.data[
  labels(tree_src.fg1_6.tracing.colree[[2]][[2]][[2]][[1]]),]$tree = 5
src.fg1_6.tracing@meta.data[
  labels(tree_src.fg1_6.tracing.colree[[2]][[2]][[2]][[2]][[2]]),]$tree = 6
src.fg1_6.tracing@meta.data[
  labels(tree_src.fg1_6.tracing.colree[[2]][[2]][[2]][[2]][[2]][[2]]),]$tree = 7


rownames(src.fg1_6.tracing@meta.data[src.fg1_6.tracing$tree%in%3,])


#-------------------------------------------------------------------------------
src.fg1.tracing.re = src.fg1_6.tracing[,src.fg1_6.tracing$tree%in%c(3,4,7)&
                                         src.fg1_6.tracing$lineage%in%t(
                                           values(endoderm_lineage_raw[["FG.1"]]["ss9"]))]
src.fg6.tracing.re = src.fg1_6.tracing[,src.fg1_6.tracing$tree%in%c(1,2,5,6)&
                                         src.fg1_6.tracing$lineage%in%t(
                                           values(endoderm_lineage_raw[["FG.6"]]["ss9"]))]

# Create :: src.fg1.integrated.merge.re 
#------------------------------------------
src.fg1.integrated.merge.re = 
  merge(src.fg1.integrated, src.fg1.tracing.re)
src.fg1.integrated.merge.re = 
  batch_process_integration_Re(
    cluster_names = "FG.1", set_src = T,
    seurat = src.fg1.integrated.merge.re,
    red_refer = values(endoderm_embedding, keys="FG.1"))
save(src.fg1.integrated.merge.re, 
     file = "figure.v08.07/organ_development_re_v240115/src.fg1.integrated.merge.re.Rdata")

src.fg6.integrated.merge.re = 
  merge(src.fg6.integrated, src.fg6.tracing.re)
src.fg6.integrated.merge.re = 
  batch_process_integration_Re(
    cluster_names = "FG.6",set_src = T,
    seurat = src.fg6.integrated.merge.re,
    red_refer = values(endoderm_embedding, keys="FG.6"))
save(src.fg6.integrated.merge.re, 
     file = "figure.v08.07/organ_development_re_v240115/src.fg6.integrated.merge.re.Rdata")
#------------------------------------------


colors.geneset = c(
  "#d76364", "#9dc3e7", "#f1d77e", "#b1ce46",
  "#63e398", "#9394e7", "#5f97d2", "#14517c",
  "#ef7a6d", "#f7e1ed", "#c497b2", "#f8f3f9")
names(colors.geneset) = c(1:12) # 4 3 9 7

#-----------
# FG.1-Re
#-----------
# Correction for 9-SS
src.fg1.tracing.re@meta.data[
  src.fg1.tracing.re$cluster.extract.v1.1_define%in%"FG.1"&
    src.fg1.tracing.re$Time%in%"9ss",]$cluster.v06.26.re_mnn_umap_fta = "FG.1"
src.fg1.tracing.re = FindVariableFeatures(src.fg1.tracing.re, nfeatures = 2000)
src.fg1.tracing.re.filtergene = 
  Myfilter(as.matrix(src.fg1.tracing.re@assays$RNA@data),
           gene = src.fg1.tracing.re@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.fg1.tracing.re = ScaleData(src.fg1.tracing.re, 
                               rownames(src.fg1.tracing.re))
src.fg1.tracing.re = RunPCA(src.fg1.tracing.re, 
                            features = src.fg1.tracing.re.filtergene)
src.fg1.tracing.re = RunUMAP(src.fg1.tracing.re,
                             dims = 1:30, reduction = "pca", assay = "RNA")
DimPlot(src.fg1.tracing.re, group.by = "Time")

src.fg1.tracing.re[["mnn_umap_fta"]] = src.fg1.tracing.re[["umap"]]
src.fg1.tracing.re@reductions$mnn_umap_fta@key = "Coord_"
src.fg1.tracing.re@reductions$mnn_umap_fta@cell.embeddings = 
  src.fg1.integrated.merge.re@reductions$mnn_umap_fta@cell.embeddings[colnames(src.fg1.tracing.re),]
colnames(src.fg1.tracing.re@reductions$mnn_umap_fta@cell.embeddings) = c("Coord_1","Coord_2")
DimPlot(src.fg1.tracing.re, reduction = "mnn_umap_fta",
        group.by = "Time", cols = colors.time.2)

coord_src.fg1.tracing.re = src.fg1.tracing.re[["mnn_umap_fta"]]@cell.embeddings
cell_src.fg1.tracing.re = rownames(src.fg1.tracing.re@meta.data)
pcurve_src.fg1.tracing.re = princurve::principal_curve(
  x = coord_src.fg1.tracing.re[cell_src.fg1.tracing.re,],  
  smoother = "smooth.spline")
src.fg1.tracing.re$lambda = pcurve_src.fg1.tracing.re$lambda[cell_src.fg1.tracing.re]
src.fg1.tracing.re$lambda = norm_range(src.fg1.tracing.re$lambda)
src.fg1.tracing.re$order = pcurve_src.fg1.tracing.re$ord[cell_src.fg1.tracing.re]
cellorder_src.fg1.tracing.re = c(
  intersect(names(src.fg1.tracing.re$lambda[order(src.fg1.tracing.re$lambda)]),
            rownames(src.fg1.tracing.re@meta.data[src.fg1.tracing.re$cluster.v06.26.re_correct%in%"FG.1",])),
  intersect(names(src.fg1.tracing.re$lambda[order(src.fg1.tracing.re$lambda)]),
            rownames(src.fg1.tracing.re@meta.data[src.fg1.tracing.re$cluster.v06.26.re_correct%in%"Pharynx.organ.2",]))
)
  


pdf("figure.v08.07/organ_development_re_v240115/try.fg1.pdf",9,7)
src.fg1.tracing.re.rowtree =
  MyHeatmap(as.matrix(src.fg1.tracing.re@assays$RNA@data[
    src.fg1.tracing.re.filtergene,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg1.tracing.re$Time,
                 colors.time.2),
      MyName2Col(src.fg1.tracing.re$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg1.tracing.re$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg1.tracing.re$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()
tree_src.fg1.tracing.re.rowtree = as.dendrogram(src.fg1.tracing.re.rowtree)
gene_src.fg1.tracing.re.rowtree = c(
  labels(tree_src.fg1.tracing.re.rowtree[[1]]),
  labels(tree_src.fg1.tracing.re.rowtree[[2]][[1]])) 

pdf("figure.v08.07/organ_development_re_v240115/try.fg1.pdf",9,7)
src.fg1.tracing.re.rowtree.1 =
  MyHeatmap(as.matrix(src.fg1.tracing.re@assays$RNA@data[
    gene_src.fg1.tracing.re.rowtree,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg1.tracing.re$Time,
                 colors.time.2),
      MyName2Col(src.fg1.tracing.re$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg1.tracing.re$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg1.tracing.re$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)

src.fg1.tracing.re.coltree.1 =
  MyHeatmap(as.matrix(src.fg1.tracing.re@assays$RNA@data[
    gene_src.fg1.tracing.re.rowtree,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg1.tracing.re$Time,
                 colors.time.2),
      MyName2Col(src.fg1.tracing.re$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg1.tracing.re$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg1.tracing.re$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "col",
    graph = T)
dev.off()

pdf("figure.v08.07/organ_development_re_v240115/try.fg1.tree.pdf",100,20)
plot(src.fg1.tracing.re.rowtree.1)
dev.off()

gene_order = function(gene, type){
  data_temp = src.fg1.tracing.re@assays$RNA@data[
    gene,
    src.fg1.tracing.re$cluster.v06.26.re_correct%in%type]
  data_temp = rowSums(data_temp)
  gene = unlist(rev(names(data_temp)[order(data_temp)]))
  return(gene)
}

tree_src.fg1.tracing.re.rowtree.1 = as.dendrogram(src.fg1.tracing.re.rowtree.1)
gene_src.fg1.tracing.re.rowtree.1 = c(
  gene_order(
    setdiff(labels(tree_src.fg1.tracing.re.rowtree.1[[1]]),
            labels(tree_src.fg1.tracing.re.rowtree.1[[1]][[2]][[1]]))
  ),
  rev(c(labels(tree_src.fg1.tracing.re.rowtree.1[[2]][[1]][[2]]),
        labels(tree_src.fg1.tracing.re.rowtree.1[[2]][[2]]))))
  
names(gene_src.fg1.tracing.re.rowtree.1) = c(
  rep(4, length(
    setdiff(labels(tree_src.fg1.tracing.re.rowtree.1[[1]]),
            labels(tree_src.fg1.tracing.re.rowtree.1[[1]][[2]][[1]])))),
  rep(7, length(c(
    labels(tree_src.fg1.tracing.re.rowtree.1[[2]][[2]]),
    labels(tree_src.fg1.tracing.re.rowtree.1[[2]][[1]][[2]])))))
save(gene_src.fg1.tracing.re.rowtree.1,
     file = "figure.v08.07/organ_development_re_v240115/gene_src.fg1.tracing.re.rowtree.1.Rdata")


gene_order = function(gene, type){
  data_temp = src.fg1.tracing.re@assays$RNA@data[
    gene,
    src.fg1.tracing.re$cluster.v06.26.re_correct%in%type]
  data_temp = rowSums(data_temp)
  gene = unlist(rev(names(data_temp)[order(data_temp)]))
  return(gene)
}
gene_src.fg1.tracing.re.rowtree.1.fin = c(
  rev(gene_order(gene_src.fg1.tracing.re.rowtree.1[names(gene_src.fg1.tracing.re.rowtree.1)%in%4], 'FG.1')),
  gene_order(gene_src.fg1.tracing.re.rowtree.1[names(gene_src.fg1.tracing.re.rowtree.1)%in%7], 'Pharynx.organ.2')
)
names(gene_src.fg1.tracing.re.rowtree.1.fin) = names(gene_src.fg1.tracing.re.rowtree.1)


tree_src.fg1.tracing.re.coltree.1 = as.dendrogram(src.fg1.tracing.re.coltree.1)
src.fg1.tracing.re$tree = NA
src.fg1.tracing.re@meta.data[labels(tree_src.fg1.tracing.re.coltree.1[[1]][[1]]),]$tree = 1
src.fg1.tracing.re@meta.data[labels(tree_src.fg1.tracing.re.coltree.1[[1]][[2]]),]$tree = 2
src.fg1.tracing.re@meta.data[labels(tree_src.fg1.tracing.re.coltree.1[[2]][[1]][[1]]),]$tree = 3
src.fg1.tracing.re@meta.data[labels(tree_src.fg1.tracing.re.coltree.1[[2]][[1]][[2]]),]$tree = 5
src.fg1.tracing.re@meta.data[labels(tree_src.fg1.tracing.re.coltree.1[[2]][[2]]),]$tree = 4
DimPlot(src.fg1.tracing.re, group.by = "tree", reduction = "umap")

src.fg1.tracing.re$cluster.v06.26.re_correct = src.fg1.tracing.re$cluster.v06.26.re_mnn_umap_fta
src.fg1.tracing.re@meta.data[src.fg1.tracing.re$tree%in%c(3,4),]$cluster.v06.26.re_correct = "FG.1"
src.fg1.tracing.re@meta.data[src.fg1.tracing.re$tree%in%c(1,2,5),]$cluster.v06.26.re_correct = "Pharynx.organ.2"


# pdf("figure.v08.07/organ_development_re_v240115/try.fg1.pdf",9,7)
pdf("figure.v08.07/organ_development_re_v240115//Heatmap_for_trajectory/try.fg1.pdf",10,10)
tree_src.fg1.tracing.re.rowtree.2 =
  MyHeatmap(
    as.matrix(src.fg1.tracing.re@assays$RNA@data[gene_src.fg1.tracing.re.rowtree.1.fin,
                                                 cellorder_src.fg1.tracing.re]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg1.tracing.re$Time[cellorder_src.fg1.tracing.re], 
                 colors.time.2),
      MyName2Col(src.fg1.tracing.re$lineage[cellorder_src.fg1.tracing.re], 
                 color.lineage),
      MyName2Col(src.fg1.tracing.re$cluster.v06.26.re_correct[cellorder_src.fg1.tracing.re], 
                 cluster.endoderm.color.v5)
      ),
    RowSideColors = t(cbind(
      MyName2Col(names(gene_src.fg1.tracing.re.rowtree.1),
                 colors.geneset))),
    ColSideColorsSize = 4.8,
    RowSideColorsSize = 1.2,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    Rowv = "none",
    Colv = "none",
    #return.tree = "row",
    margins = c(10,10), graph = T)
dev.off()


tf_gene_src.fg1.tracing.re.rowtree.1 = intersect(
  gene_src.fg1.tracing.re.rowtree.1, gi[gi$TF%in%T,]$SymbolDedu)

pdf("figure.v08.07/organ_development_re_v240115/try.fg1_TF.pdf",9,7)
tree_src.fg1.tracing.re.rowtree.3 =
  MyHeatmap(
    as.matrix(src.fg1.tracing.re@assays$RNA@data[tf_gene_src.fg1.tracing.re.rowtree.1,
                                                 cellorder_src.fg1.tracing.re]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg1.tracing.re$Time[cellorder_src.fg1.tracing.re], 
                 colors.time.2),
      MyName2Col(src.fg1.tracing.re$lineage[cellorder_src.fg1.tracing.re], 
                 color.lineage),
      MyName2Col(src.fg1.tracing.re$cluster.v06.26.re_correct[cellorder_src.fg1.tracing.re], 
                 cluster.endoderm.color.v5)
    ),
    ColSideColorsSize = 3,
    RowSideColorsSize = 1,
    labRow = tf_gene_src.fg1.tracing.re.rowtree.1,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    Rowv = "none",
    Colv = "none",
    #return.tree = "row",
    margins = c(10,10), graph = T)
dev.off()

save(src.fg1.tracing.re,
     file = "figure.v08.07/organ_development_re_v240115/src.fg1.tracing.re.Rdata")
save(coord_src.fg1.tracing.re,
     cellorder_src.fg1.tracing.re,
     pcurve_src.fg1.tracing.re,
     tree_src.fg1.tracing.re.rowtree,
     gene_src.fg1.tracing.re.rowtree,
     tree_src.fg1.tracing.re.rowtree.1,
     gene_src.fg1.tracing.re.rowtree.1,
     tree_src.fg1.tracing.re.rowtree.2,
     file = "figure.v08.07/organ_development_re_v240115/src.fg1.tracing.re.parameter.Rdata")


#-----------
# FG.6-Re
#-----------
# Correction for 9-SS
src.fg6.tracing.re@meta.data[
  src.fg6.tracing.re$cluster.extract.v1.1_define%in%"FG.6"&
    src.fg6.tracing.re$Time%in%"9ss",]$cluster.v06.26.re_mnn_umap_fta = "FG.6"
src.fg6.tracing.re = FindVariableFeatures(src.fg6.tracing.re, nfeatures = 2000)
src.fg6.tracing.re.filtergene = 
  Myfilter(as.matrix(src.fg6.tracing.re@assays$RNA@data),
           gene = src.fg6.tracing.re@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.fg6.tracing.re = ScaleData(src.fg6.tracing.re, 
                               rownames(src.fg6.tracing.re))
src.fg6.tracing.re = RunPCA(src.fg6.tracing.re, 
                            features = src.fg6.tracing.re.filtergene)
src.fg6.tracing.re = RunUMAP(src.fg6.tracing.re,
                             dims = 1:30, reduction = "pca", assay = "RNA")
DimPlot(src.fg6.tracing.re, group.by = "Time")

src.fg6.tracing.re[["mnn_umap_fta"]] = src.fg6.tracing.re[["umap"]]
src.fg6.tracing.re@reductions$mnn_umap_fta@key = "Coord_"
src.fg6.tracing.re@reductions$mnn_umap_fta@cell.embeddings = 
  src.fg6.integrated.merge.re@reductions$mnn_umap_fta@cell.embeddings[colnames(src.fg6.tracing.re),]
colnames(src.fg6.tracing.re@reductions$mnn_umap_fta@cell.embeddings) = c("Coord_1","Coord_2")
DimPlot(src.fg6.tracing.re, reduction = "mnn_umap_fta",
        group.by = "Time", cols = colors.time.2)

coord_src.fg6.tracing.re = src.fg6.tracing.re[["mnn_umap_fta"]]@cell.embeddings
cell_src.fg6.tracing.re = rownames(src.fg6.tracing.re@meta.data)
pcurve_src.fg6.tracing.re = princurve::principal_curve(
  x = coord_src.fg6.tracing.re[cell_src.fg6.tracing.re,],  
  smoother = "smooth.spline")


norm_range = function(x, start=0, end=1){
  x_min = min(x)
  x_max = max(x)
  
  x_re = (end - start) * (x - x_min) / (x_max - x_min) + start
  
  return(x_re)
}

src.fg6.tracing.re$lambda = pcurve_src.fg6.tracing.re$lambda[cell_src.fg6.tracing.re]
src.fg6.tracing.re$lambda = norm_range(src.fg6.tracing.re$lambda)
src.fg6.tracing.re$order = pcurve_src.fg6.tracing.re$ord[cell_src.fg6.tracing.re]

cellorder_src.fg6.tracing.re = 
  names(src.fg6.tracing.re$lambda[order(src.fg6.tracing.re$lambda)])


pdf("figure.v08.07/organ_development_re_v240115/try.fg6.pdf",9,7)
src.fg6.tracing.re.rowtree =
  MyHeatmap(as.matrix(src.fg6.tracing.re@assays$RNA@data[
    src.fg6.tracing.re.filtergene,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg6.tracing.re$Time,
                 colors.time.2),
      MyName2Col(src.fg6.tracing.re$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg6.tracing.re$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg6.tracing.re$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()

tree_src.fg6.tracing.re.rowtree = as.dendrogram(src.fg6.tracing.re.rowtree)
gene_src.fg6.tracing.re.rowtree = c(
  rev(labels(tree_src.fg6.tracing.re.rowtree[[1]])),
  labels(tree_src.fg6.tracing.re.rowtree[[2]][[1]])) 
names(gene_src.fg6.tracing.re.rowtree) = c(
  rep(4, length(labels(tree_src.fg6.tracing.re.rowtree[[1]])) ),
  rep(7, length(labels(tree_src.fg6.tracing.re.rowtree[[2]][[1]])) ))


gene_order = function(gene, type){
  data_temp = src.fg6.tracing.re@assays$RNA@data[
    gene,
    src.fg6.tracing.re$cluster.v06.26.re_correct%in%type]
  data_temp = rowSums(data_temp)
  gene = unlist(rev(names(data_temp)[order(data_temp)]))
  return(gene)
}
gene_src.fg6.tracing.re.rowtree.1.fin = c(
  rev(gene_order(gene_src.fg6.tracing.re.rowtree[names(gene_src.fg6.tracing.re.rowtree)%in%4], 'FG.6')),
  gene_order(gene_src.fg6.tracing.re.rowtree[names(gene_src.fg6.tracing.re.rowtree)%in%7], 'Esophagus')
)
names(gene_src.fg6.tracing.re.rowtree.1.fin) = names(gene_src.fg6.tracing.re.rowtree)




pdf("figure.v08.07/organ_development_re_v240115/try.fg6.pdf",9,7)
src.fg6.tracing.re.coltree =
  MyHeatmap(
    as.matrix(src.fg6.tracing.re@assays$RNA@data[gene_src.fg6.tracing.re.rowtree,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg6.tracing.re$Time,
                 colors.time.2),
      MyName2Col(src.fg6.tracing.re$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg6.tracing.re$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      # MyName2Col(src.fg6.tracing.re$tree,
      #            color.temp),
      MyName2Col(src.fg6.tracing.re$lineage,
                 color.lineage)),
    RowSideColors = t(cbind(
      MyName2Col(names(gene_src.fg6.tracing.re.rowtree),
                 colors.geneset))),
    ColSideColorsSize = 4,
    RowSideColorsSize = 1,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    # Rowv = "none",
    # Colv = "none",
    return.tree = "col",
    margins = c(10,10), graph = T)
dev.off()

tree_src.fg6.tracing.re.coltree = as.dendrogram(src.fg6.tracing.re.coltree)
src.fg6.tracing.re$tree = NA
src.fg6.tracing.re@meta.data[labels(tree_src.fg6.tracing.re.coltree[[1]]),]$tree = 1
src.fg6.tracing.re@meta.data[labels(tree_src.fg6.tracing.re.coltree[[2]][[1]]),]$tree = 2
src.fg6.tracing.re@meta.data[labels(tree_src.fg6.tracing.re.coltree[[2]][[2]][[1]]),]$tree = 3
src.fg6.tracing.re@meta.data[labels(tree_src.fg6.tracing.re.coltree[[2]][[2]][[2]]),]$tree = 4

src.fg6.tracing.re$cluster.v06.26.re_correct = 
  src.fg6.tracing.re$cluster.v06.26.re_mnn_umap_fta
src.fg6.tracing.re@meta.data[src.fg6.tracing.re$tree%in%c(1,4),]$cluster.v06.26.re_correct = "FG.6"
src.fg6.tracing.re@meta.data[src.fg6.tracing.re$tree%in%c(2,3),]$cluster.v06.26.re_correct = "Esophagus"

pdf("figure.v08.07/organ_development_re_v240115/Heatmap_for_trajectory/try.fg6.re.pdf", 10, 10)
src.fg6.tracing.re.tree =
  MyHeatmap(
    as.matrix(src.fg6.tracing.re@assays$RNA@data[gene_src.fg6.tracing.re.rowtree.1.fin,
                                                 cellorder_src.fg6.tracing.re]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg6.tracing.re$Time[cellorder_src.fg6.tracing.re],
                 colors.time.2),
      MyName2Col(src.fg6.tracing.re$lineage[cellorder_src.fg6.tracing.re],
                 color.lineage),
      MyName2Col(src.fg6.tracing.re$cluster.v06.26.re_correct[cellorder_src.fg6.tracing.re],
                 cluster.endoderm.color.v5)),
    RowSideColors = t(cbind(
      MyName2Col(names(gene_src.fg6.tracing.re.rowtree),
                 colors.geneset))),
    ColSideColorsSize = 4.8,
    RowSideColorsSize = 1.2,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    Rowv = "none",
    Colv = "none",
    #return.tree = "row",
    margins = c(10,10), graph = T)
dev.off()



pdf("figure.v08.07/organ_development_re_v240115/Heatmap_for_trajectory/try.fg6_TF.pdf",10,10)
src.fg6.tracing.re.tree =
  MyHeatmap(
    as.matrix(src.fg6.tracing.re@assays$RNA@data[intersect(gene_src.fg6.tracing.re.rowtree,
                                                           gi[gi$TF%in%T,]$SymbolDedu),
                                                 cellorder_src.fg6.tracing.re]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg6.tracing.re$Time[cellorder_src.fg6.tracing.re],
                 colors.time.2),
      MyName2Col(src.fg6.tracing.re$lineage[cellorder_src.fg6.tracing.re],
                 color.lineage),
      MyName2Col(src.fg6.tracing.re$cluster.v06.26.re_correct[cellorder_src.fg6.tracing.re],
                 cluster.endoderm.color.v5)),
    ColSideColorsSize = 4.8,
    RowSideColorsSize = 1.2,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    Rowv = "none",
    Colv = "none",
    labRow = intersect(gene_src.fg6.tracing.re.rowtree,
                       gi[gi$TF%in%T,]$SymbolDedu),
    #return.tree = "row",
    margins = c(10,10), graph = T)
dev.off()

save(src.fg6.tracing.re,
     file = "figure.v08.07/organ_development_re_v240115/src.fg6.tracing.re.Rdata")
save(coord_src.fg6.tracing.re,
     cellorder_src.fg6.tracing.re,
     pcurve_src.fg6.tracing.re,
     tree_src.fg6.tracing.re.rowtree,
     gene_src.fg6.tracing.re.rowtree,
     file = "figure.v08.07/organ_development_re_v240115/src.fg6.tracing.re.parameter.Rdata")




#-----------------------------------
# FG.2 - MNN re-Correction
#-----------------------------------

seurat = src.fg2.integrated.merge
seurat.selectgene = src.fg2.integrated.selectgene
cell_sort_raw = colnames(seurat)
cell_refer = colnames(seurat[,seurat$Source_tech%in%'refer'])
cell_query = colnames(seurat[,seurat$Source_tech%in%'query'])
assay = "mnnRNA"
src.10x.refer = src.fg2.integrated

MNN.res = 
  mnnCorrect(
    as.matrix(seurat@assays$RNA@data[seurat.selectgene, 
                                     seurat$Source_tech%in%"refer" & seurat$batch%in%1]),
    as.matrix(seurat@assays$RNA@data[seurat.selectgene, 
                                     seurat$Source_tech%in%"refer" & seurat$batch%in%2]),
    as.matrix(seurat@assays$RNA@data[seurat.selectgene,
                                     seurat$Source_tech%in%"query" & 
                                       (seurat$SeqDate%in%c("20230315","20230513","20230324")&
                                           seurat$Time%in%c("9ss","12ss"))]),
    as.matrix(seurat@assays$RNA@data[seurat.selectgene,
                                     seurat$Source_tech%in%"query" & 
                                       !(seurat$SeqDate%in%c("20230315","20230513","20230324")&
                                          seurat$Time%in%c("9ss","12ss"))]),
    k = 20, cos.norm.out=F)
seurat@assays$mnnRNA = CreateAssayObject(data = MNN.res@assays@data$corrected[,cell_sort_raw])
seurat@assays$mnnRNA@key = "mnn_"
seurat = ScaleData(seurat,  rownames(seurat@assays$mnnRNA), assay = "mnnRNA")

anchor.integrated = 
  FindTransferAnchors(reference = seurat[,cell_refer],
                      query = seurat[,cell_query],
                      reference.assay =  assay,
                      query.assay =  assay, scale = T,
                      features = seurat.selectgene)

red_refer = endoderm_embedding[["FG.2"]]
umap_embedding = src.10x.refer[[red_refer]]@cell.embeddings[cell_refer,c(1:2)]
umap.transfer = TransferData(anchor.integrated, t(umap_embedding))
assay_name = gsub('RNA',"", gsub("mnn","mnn_",assay))
seurat[[paste(assay_name, "umap_fta", sep="")]] =  seurat[["umap"]]
seurat[[paste(assay_name, "umap_fta", sep="")]]@cell.embeddings[cell_refer,] = umap_embedding
seurat[[paste(assay_name, "umap_fta", sep="")]]@cell.embeddings[cell_query,] = 
  as.matrix(t(umap.transfer@data))
src.fg2.integrated.merge = seurat

names = colnames(src.fg2.integrated.merge@reductions$mnn_umap_fta@cell.embeddings)
src.fg2.integrated.merge@reductions$mnn_umap_fta@cell.embeddings = cbind(
  -src.fg2.integrated.merge@reductions$mnn_umap_fta@cell.embeddings[,1],
  src.fg2.integrated.merge@reductions$mnn_umap_fta@cell.embeddings[,2])
colnames(src.fg2.integrated.merge@reductions$mnn_umap_fta@cell.embeddings) = names


#------------------
#  FG.2 Re
#------------------
src.fg2.tracing@meta.data[
  src.fg2.tracing$cluster.extract.v1.1_define%in%"FG.2"&
    src.fg2.tracing$Time%in%"9ss",]$cluster.v06.26.re_mnn_umap_fta = "FG.2"
src.fg2.tracing.filtergene = 
  Myfilter(as.matrix(src.fg2.tracing@assays$RNA@data),
           gene = src.fg2.tracing@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)

src.fg2.tracing[["mnn_umap_fta"]] = src.fg2.tracing[["umap"]]
src.fg2.tracing@reductions$mnn_umap_fta@key = "Coord_"
src.fg2.tracing@reductions$mnn_umap_fta@cell.embeddings = 
  src.fg2.integrated.merge@reductions$mnn_umap_fta@cell.embeddings[colnames(src.fg2.tracing),]
colnames(src.fg2.tracing@reductions$mnn_umap_fta@cell.embeddings) = c("Coord_1","Coord_2")
DimPlot(src.fg2.tracing, reduction = "mnn_umap_fta",
        group.by = "Time", cols = colors.time.2)

coord_src.fg2.tracing = src.fg2.tracing[["mnn_umap_fta"]]@cell.embeddings
cell_src.fg2.tracing = rownames(src.fg2.tracing@meta.data)
pcurve_src.fg2.tracing = princurve::principal_curve(
  x = coord_src.fg2.tracing[cell_src.fg2.tracing,],  
  smoother = "smooth.spline")

src.fg2.tracing$lambda = pcurve_src.fg2.tracing$lambda[cell_src.fg2.tracing]
src.fg2.tracing$lambda = norm_range(src.fg2.tracing$lambda)
src.fg2.tracing$order = pcurve_src.fg2.tracing$ord[cell_src.fg2.tracing]
DimPlot(src.fg2.tracing, reduction = "umap", group.by = "Time", label = T)

mat_1 = src.fg2.tracing@meta.data[src.fg2.tracing$cluster.v06.26.re_mnn_umap_fta%in%"FG.2",]
mat_2 = src.fg2.tracing@meta.data[src.fg2.tracing$cluster.v06.26.re_mnn_umap_fta%in%"Pharynx.organ.1",]
cellorder_src.fg2.tracing = c(
  rev(rownames(mat_1[order(mat_1$lambda),])),
  rev(rownames(mat_2[order(mat_2$lambda),])))


src.fg2.tracing = SetIdent(src.fg2.tracing, value = src.fg2.tracing$cluster.v06.26.re_mnn_umap_fta)
marker_src.fg2.tracing = FindAllMarkers(src.fg2.tracing)
marker_src.fg2.tracing$pct.ratio = marker_src.fg2.tracing$pct.1 / marker_src.fg2.tracing$pct.2
marker_src.fg2.tracing$rank = marker_src.fg2.tracing$pct.ratio * (-log(marker_src.fg2.tracing$p_val_adj))
marker_src.fg2.tracing = marker_src.fg2.tracing[order(marker_src.fg2.tracing$rank, decreasing = T),]
markergene_src.fg2.tracing = unique(marker_src.fg2.tracing$gene)

pdf("figure.v08.07/organ_development_re_v240115/try.fg2.pdf",9,7)
src.fg2.tracing.rowtree =
  MyHeatmap(as.matrix(src.fg2.tracing@assays$RNA@data[
    #unique(c(markergene_src.fg2.tracing, src.fg2.tracing.filtergene)),
    markergene_src.fg2.tracing,
    cellorder_src.fg2.tracing]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg2.tracing$Time[cellorder_src.fg2.tracing],
                 colors.time.2),
      MyName2Col(src.fg2.tracing$cluster.extract.v1.1_define[cellorder_src.fg2.tracing],
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg2.tracing$cluster.v06.26.re_mnn_umap_fta[cellorder_src.fg2.tracing],
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg2.tracing$lineage[cellorder_src.fg2.tracing],
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    # Colv = "none",
    # Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()

tree_src.fg2.tracing.rowtree = as.dendrogram(src.fg2.tracing.rowtree)
gene_src.fg2.tracing.rowtree = c(
  setdiff(labels(tree_src.fg2.tracing.rowtree[[1]]),
          labels(tree_src.fg2.tracing.rowtree[[1]][[2]][[1]])),
  labels(tree_src.fg2.tracing.rowtree[[2]][[2]][[1]]),
  labels(tree_src.fg2.tracing.rowtree[[2]][[1]]))
names(gene_src.fg2.tracing.rowtree) = c(
  rep(4, length(setdiff(labels(tree_src.fg2.tracing.rowtree[[1]]),
                        labels(tree_src.fg2.tracing.rowtree[[1]][[2]][[1]])))),
  rep(7, length(setdiff(labels(tree_src.fg2.tracing.rowtree[[2]]),
                        labels(tree_src.fg2.tracing.rowtree[[2]][[2]][[2]])))))


pdf("figure.v08.07/organ_development_re_v240115/try.fg2.pdf",9,7)
src.fg2.tracing.coltree =
  MyHeatmap(as.matrix(src.fg2.tracing@assays$RNA@data[
    gene_src.fg2.tracing.rowtree,
    cellorder_src.fg2.tracing]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg2.tracing$Time[cellorder_src.fg2.tracing],
                 colors.time.2),
      MyName2Col(src.fg2.tracing$cluster.extract.v1.1_define[cellorder_src.fg2.tracing],
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg2.tracing$cluster.v06.26.re_mnn_umap_fta[cellorder_src.fg2.tracing],
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg2.tracing$lineage[cellorder_src.fg2.tracing],
                 color.lineage)
    ),
    RowSideColors = t(cbind(
      MyName2Col(names(gene_src.fg2.tracing.rowtree),
                 colors.geneset))),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    # Colv = "none",
    # Rowv = "none",
    return.tree = "col",
    graph = T)
dev.off()

tree_src.fg2.tracing.coltree = as.dendrogram(src.fg2.tracing.coltree)
cell_src.fg2.tracing.coltree = c(
  labels(tree_src.fg2.tracing.coltree[[1]]),
  labels(tree_src.fg2.tracing.coltree[[2]]))
names(cell_src.fg2.tracing.coltree) = c(
  rep(1, length(labels(tree_src.fg2.tracing.coltree[[1]]))),
  rep(2, length(labels(tree_src.fg2.tracing.coltree[[1]]))))


src.fg2.tracing$tree = NA
src.fg2.tracing$tree[labels(tree_src.fg2.tracing.coltree[[1]][[1]][[1]])] = 1
src.fg2.tracing$tree[labels(tree_src.fg2.tracing.coltree[[1]][[1]][[2]])] = 2
src.fg2.tracing$tree[labels(tree_src.fg2.tracing.coltree[[1]][[2]][[1]])] = 3
src.fg2.tracing$tree[labels(tree_src.fg2.tracing.coltree[[1]][[2]][[2]])] = 4
src.fg2.tracing$tree[labels(tree_src.fg2.tracing.coltree[[2]][[1]][[1]])] = 5
src.fg2.tracing$tree[labels(tree_src.fg2.tracing.coltree[[2]][[1]][[2]])] = 6
src.fg2.tracing$tree[labels(tree_src.fg2.tracing.coltree[[2]][[2]][[1]])] = 7
src.fg2.tracing$tree[labels(tree_src.fg2.tracing.coltree[[2]][[2]][[2]])] = 8
DimPlot(src.fg2.tracing, reduction = "umap_fta", group.by = "tree")
DimPlot(src.fg2.tracing, reduction = "umap_fta", group.by = "Time", cols = colors.time.2)


src.fg2.tracing$cluster.v06.26.re_correct = src.fg2.tracing$cluster.v06.26.re_mnn_umap_fta
src.fg2.tracing@meta.data[labels(tree_src.fg2.tracing.coltree[[1]]),]$cluster.v06.26.re_correct  = "FG.2"
src.fg2.tracing@meta.data[labels(tree_src.fg2.tracing.coltree[[2]]),]$cluster.v06.26.re_correct  = "Pharynx.organ.1"
src.fg2.tracing@meta.data[
  src.fg2.tracing$Time%in%c("9ss","12ss","15ss") &
    src.fg2.tracing$tree%in%c(7),]$cluster.v06.26.re_correct  = "FG.2"
src.fg2.tracing@meta.data[
  src.fg2.tracing$Time%in%c("9ss") &
    src.fg2.tracing$tree%in%c(5,6,7,8),]$cluster.v06.26.re_correct  = "FG.2"
save(src.fg2.tracing, file ="figure.v08.07/organ_development_re_v240115/src.fg2.tracing.Rdata")

mat_1 = src.fg2.tracing@meta.data[src.fg2.tracing$cluster.v06.26.re_correct%in%"FG.2",]
mat_2 = src.fg2.tracing@meta.data[src.fg2.tracing$cluster.v06.26.re_correct%in%"Pharynx.organ.1",]
cellorder_src.fg2.tracing = c(
  rev(rownames(mat_1[order(mat_1$lambda),])),
  rev(rownames(mat_2[order(mat_2$lambda),])))


pdf("figure.v08.07/organ_development_re_v240115/Heatmap_for_trajectory/try.fg2.pdf",10,10)
src.fg2.tracing.rowtree.1 =
  MyHeatmap(
    as.matrix(src.fg2.tracing@assays$RNA@data[
      gene_src.fg2.tracing.rowtree,
      cellorder_src.fg2.tracing]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg2.tracing$Time[cellorder_src.fg2.tracing],
                 colors.time.2),
      MyName2Col(src.fg2.tracing$lineage[cellorder_src.fg2.tracing],
                 color.lineage),
      MyName2Col(src.fg2.tracing$cluster.v06.26.re_correct[cellorder_src.fg2.tracing],
                 cluster.endoderm.color.v5)),
    RowSideColors = t(cbind(
      MyName2Col(names(gene_src.fg2.tracing.rowtree),
                 colors.geneset))),
    ColSideColorsSize = 4.8,
    RowSideColorsSize = 1.2,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    Colv = "none",
    Rowv = "none",
    # return.tree = "row",
    margins = c(10,10), graph = T)
dev.off()


save(src.fg2.integrated.merge,
     file = "figure.v08.07/organ_development_re_v240115/src.fg2.integrated.merge.Rdata")
save(src.fg2.tracing,
     file = "figure.v08.07/organ_development_re_v240115/src.fg2.tracing.Rdata")
save(coord_src.fg2.tracing,
     cellorder_src.fg2.tracing,
     pcurve_src.fg2.tracing,
     markergene_src.fg2.tracing,
     
     tree_src.fg2.tracing.rowtree,
     gene_src.fg2.tracing.rowtree,
     tree_src.fg2.tracing.coltree,
     cell_src.fg2.tracing.coltree,
     
     file = "figure.v08.07/organ_development_re_v240115/src.fg2.tracing.parameter.Rdata")



#-----------
# FG.5-Re
#-----------
# Correction for 9-SS
src.fg5.tracing@meta.data[
  src.fg5.tracing$cluster.extract.v1.1_define%in%"FG.5"&
    src.fg5.tracing$Time%in%"9ss",]$cluster.v06.26.re_mnn_umap_fta = "FG.5"
src.fg5.tracing = FindVariableFeatures(src.fg5.tracing, nfeatures = 2000)
src.fg5.tracing.filtergene = 
  Myfilter(as.matrix(src.fg5.tracing@assays$RNA@data),
           gene = src.fg5.tracing@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.fg5.tracing = ScaleData(src.fg5.tracing, 
                            rownames(src.fg5.tracing))
src.fg5.tracing = RunPCA(src.fg5.tracing, 
                         features = src.fg5.tracing.filtergene)
src.fg5.tracing = RunUMAP(src.fg5.tracing,
                          dims = 1:30, reduction = "pca", assay = "RNA")
DimPlot(src.fg5.tracing, group.by = "Time")

src.fg5.tracing[["mnn_umap_fta"]] = src.fg5.tracing[["umap"]]
src.fg5.tracing@reductions$mnn_umap_fta@key = "Coord_"
src.fg5.tracing@reductions$mnn_umap_fta@cell.embeddings = 
  src.fg5.integrated.merge@reductions$mnn_umap_fta@cell.embeddings[colnames(src.fg5.tracing),]
colnames(src.fg5.tracing@reductions$mnn_umap_fta@cell.embeddings) = c("Coord_1","Coord_2")
DimPlot(src.fg5.tracing, reduction = "mnn_umap_fta",
        group.by = "Time", cols = colors.time.2)

# Pseudo-time :: cluster.v06.26.re_mnn_umap_fta
#---------------------
#-- pha3 --
cell_src.fg5.tracing_pha3 = 
  rownames(src.fg5.tracing@meta.data[src.fg5.tracing$cluster.v06.26.re_mnn_umap_fta%in%c("FG.5",'Pharynx.organ.3'),])
coord_src.fg5.tracing_pha3 = src.fg5.tracing[["mnn_umap_fta"]]@cell.embeddings[cell_src.fg5.tracing_pha3,]
pcurve_src.fg5.tracing_pha3 = princurve::principal_curve(x = coord_src.fg5.tracing_pha3, smoother = "smooth.spline")
src.fg5.tracing$lambda_pha3 = pcurve_src.fg5.tracing_pha3$lambda[cell_src.fg5.tracing_pha3]
src.fg5.tracing@meta.data[cell_src.fg5.tracing_pha3,]$lambda_pha3 = 
  pcurve_src.fg5.tracing_pha3$lambda[cell_src.fg5.tracing_pha3]
src.fg5.tracing@meta.data[!src.fg5.tracing$lambda_pha3%in%NA,]$lambda_pha3 = 
  norm_range(src.fg5.tracing@meta.data[!src.fg5.tracing$lambda_pha3%in%NA,]$lambda_pha3)

#-- th --
cell_src.fg5.tracing_th = 
  rownames(src.fg5.tracing@meta.data[src.fg5.tracing$cluster.v06.26.re_mnn_umap_fta%in%c("FG.5",'Thyroid'),])
coord_src.fg5.tracing_th = src.fg5.tracing[["mnn_umap_fta"]]@cell.embeddings[cell_src.fg5.tracing_th,]
pcurve_src.fg5.tracing_th = princurve::principal_curve(x = coord_src.fg5.tracing_th, smoother = "smooth.spline")
src.fg5.tracing$lambda_th = pcurve_src.fg5.tracing_th$lambda[cell_src.fg5.tracing_th]
src.fg5.tracing@meta.data[cell_src.fg5.tracing_th,]$lambda_th = 
  pcurve_src.fg5.tracing_th$lambda[cell_src.fg5.tracing_th]
src.fg5.tracing@meta.data[!src.fg5.tracing$lambda_th%in%NA,]$lambda_th = 
  norm_range(src.fg5.tracing@meta.data[!src.fg5.tracing$lambda_th%in%NA,]$lambda_th)

cell_src.fg5.tracing_fg5 = 
  rownames(src.fg5.tracing@meta.data[src.fg5.tracing$cluster.v06.26.re_mnn_umap_fta%in%c("FG.5"),])
src.fg5.tracing$lambda_fg5 = NA
w1 = length(cell_src.fg5.tracing_pha3) / (length(cell_src.fg5.tracing_th) + length(cell_src.fg5.tracing_pha3))
w2 = length(cell_src.fg5.tracing_th) / (length(cell_src.fg5.tracing_th) + length(cell_src.fg5.tracing_pha3))
src.fg5.tracing$lambda_fg5 = 
  w1 * src.fg5.tracing$lambda_pha3[cell_src.fg5.tracing_fg5] + 
  w2 * src.fg5.tracing$lambda_th[cell_src.fg5.tracing_fg5]


cellorder_src.fg5.tracing = c(
  rev(intersect(names(src.fg5.tracing$lambda_fg5[order(src.fg5.tracing$lambda_fg5)]), 
                colnames(src.fg5.tracing[,src.fg5.tracing$cluster.v06.26.re_mnn_umap_fta%in%"FG.5"]))),
  rev(intersect(names(src.fg5.tracing$lambda_pha3[order(src.fg5.tracing$lambda_pha3)]), 
                colnames(src.fg5.tracing[,src.fg5.tracing$cluster.v06.26.re_mnn_umap_fta%in%"Pharynx.organ.3"]))),
  intersect(names(src.fg5.tracing$lambda_th[order(src.fg5.tracing$lambda_th)]), 
            colnames(src.fg5.tracing[,src.fg5.tracing$cluster.v06.26.re_mnn_umap_fta%in%"Thyroid"])))
#---------------------

src.fg5.tracing = SetIdent(src.fg5.tracing, value = src.fg5.tracing$cluster.v06.26.re_mnn_umap_fta)
marker_src.fg5.tracing = FindAllMarkers(src.fg5.tracing)
marker_src.fg5.tracing$pct.ratio = marker_src.fg5.tracing$pct.1 / marker_src.fg5.tracing$pct.2
marker_src.fg5.tracing$rank = marker_src.fg5.tracing$pct.ratio * (-log(marker_src.fg5.tracing$p_val_adj))
marker_src.fg5.tracing = marker_src.fg5.tracing[order(marker_src.fg5.tracing$rank, decreasing = T),]
markergene_src.fg5.tracing = unique(marker_src.fg5.tracing$gene)
markergene_src.fg5.tracing.raw = markergene_src.fg5.tracing


pdf("figure.v08.07/organ_development_re_v240115/try.fg5.pdf",9,7)
src.fg5.tracing.rowtree =
  MyHeatmap(as.matrix(src.fg5.tracing@assays$RNA@data[
    unique(c(markergene_src.fg5.tracing,
             src.fg5.tracing.filtergene,
             c()
             )),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg5.tracing$Time,
                 colors.time.2),
      MyName2Col(src.fg5.tracing$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg5.tracing$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg5.tracing$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()

pdf("figure.v08.07/organ_development_re_v240115/try.fg5.tree.pdf",150,20)
plot(src.fg5.tracing.rowtree)
dev.off()

tree_src.fg5.tracing.rowtree = as.dendrogram(src.fg5.tracing.rowtree)
gene_src.fg5.tracing.rowtree = setdiff(
  unique(c(markergene_src.fg5.tracing, src.fg5.tracing.filtergene)),
  labels(tree_src.fg5.tracing.rowtree[[2]][[2]][[1]]))
  

pdf("figure.v08.07/organ_development_re_v240115/try.fg5.pdf",9,7)
src.fg5.tracing.rowtree.1 =
  MyHeatmap(as.matrix(src.fg5.tracing@assays$RNA@data[
    gene_src.fg5.tracing.rowtree,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg5.tracing$Time,
                 colors.time.2),
      MyName2Col(src.fg5.tracing$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg5.tracing$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg5.tracing$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()
tree_src.fg5.tracing.rowtree.1 = as.dendrogram(src.fg5.tracing.rowtree.1)
gene_src.fg5.tracing.rowtree.1 = setdiff(
  gene_src.fg5.tracing.rowtree,
  c(labels(tree_src.fg5.tracing.rowtree.1[[2]][[2]][[1]][[2]]),
    labels(tree_src.fg5.tracing.rowtree.1[[2]][[2]][[2]][[1]])))


pdf("figure.v08.07/organ_development_re_v240115/try.fg5.pdf",9,7)
src.fg5.tracing.rowtree.2 =
  MyHeatmap(as.matrix(src.fg5.tracing@assays$RNA@data[
    gene_src.fg5.tracing.rowtree.1,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg5.tracing$Time,
                 colors.time.2),
      MyName2Col(src.fg5.tracing$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg5.tracing$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      #MyName2Col(src.fg5.tracing$cluster.v06.26.re_correct,
      #           cluster.endoderm.color.v5),
      MyName2Col(src.fg5.tracing$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()
pdf("figure.v08.07/organ_development_re_v240115/try.fg5.tree.pdf",150,20)
plot(src.fg5.tracing.rowtree.2)
dev.off()
tree_src.fg5.tracing.rowtree.2 = as.dendrogram(src.fg5.tracing.rowtree.2)
gene_src.fg5.tracing.rowtree.2 = setdiff(
  gene_src.fg5.tracing.rowtree.1,
  c(labels(tree_src.fg5.tracing.rowtree.2[[2]][[1]][[1]][[1]]),
    labels(tree_src.fg5.tracing.rowtree.2[[2]][[2]][[1]][[1]])))


pdf("figure.v08.07/organ_development_re_v240115/try.fg5.pdf",9,7)
src.fg5.tracing.rowtree.3 =
  MyHeatmap(as.matrix(src.fg5.tracing@assays$RNA@data[
    gene_src.fg5.tracing.rowtree.2,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg5.tracing$Time,
                 colors.time.2),
      MyName2Col(src.fg5.tracing$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg5.tracing$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg5.tracing$lineage,
                 color.lineage)
      # MyName2Col(src.fg5.tracing$tree,
      #            color.temp)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()


tree_src.fg5.tracing.rowtree.3 = as.dendrogram(src.fg5.tracing.rowtree.3)
gene_src.fg5.tracing.rowtree.3 = c(
  labels(tree_src.fg5.tracing.rowtree.3[[2]][[2]][[2]]),

  labels(tree_src.fg5.tracing.rowtree.3[[1]][[1]]),
  labels(tree_src.fg5.tracing.rowtree.3[[1]][[2]][[1]]),
  rev(labels(tree_src.fg5.tracing.rowtree.3[[1]][[2]][[2]])),
  
  labels(tree_src.fg5.tracing.rowtree.3[[2]][[2]][[1]]),
  labels(tree_src.fg5.tracing.rowtree.3[[2]][[1]][[2]]),
  labels(tree_src.fg5.tracing.rowtree.3[[2]][[1]][[1]]))

names(gene_src.fg5.tracing.rowtree.3) = c(
  rep(4, length(labels(tree_src.fg5.tracing.rowtree.3[[2]][[2]][[2]]))),

  rep(3, length(c(labels(tree_src.fg5.tracing.rowtree.3[[1]][[1]]),
                  labels(tree_src.fg5.tracing.rowtree.3[[1]][[2]][[1]]),
                  labels(tree_src.fg5.tracing.rowtree.3[[1]][[2]][[2]])))),
  
  rep(7, length(labels(tree_src.fg5.tracing.rowtree.3[[2]][[2]][[1]]))),
  rep(7, length(c(labels(tree_src.fg5.tracing.rowtree.3[[2]][[1]][[2]]),
                  labels(tree_src.fg5.tracing.rowtree.3[[2]][[1]][[1]])))))



tree_src.fg5.tracing.rowtree.4 = as.dendrogram(src.fg5.tracing.rowtree.3)
gene_src.fg5.tracing.rowtree.4 = c(
  labels(tree_src.fg5.tracing.rowtree.4[[2]][[2]][[2]]),
  
  rev(labels(tree_src.fg5.tracing.rowtree.4[[1]][[1]])),
  labels(tree_src.fg5.tracing.rowtree.4[[1]][[2]][[2]][[2]][[1]][[1]]),
  labels(tree_src.fg5.tracing.rowtree.4[[1]][[2]][[2]][[2]][[1]][[2]][[1]]),
  labels(tree_src.fg5.tracing.rowtree.4[[1]][[2]][[2]][[2]][[1]][[2]][[2]][[1]]),
  labels(tree_src.fg5.tracing.rowtree.4[[1]][[2]][[2]][[1]][[2]]),
  setdiff(labels(tree_src.fg5.tracing.rowtree.4[[1]][[2]][[2]][[2]]),
    c(labels(tree_src.fg5.tracing.rowtree.4[[1]][[2]][[2]][[2]][[1]][[1]]),
      labels(tree_src.fg5.tracing.rowtree.4[[1]][[2]][[2]][[2]][[1]][[2]][[1]]),
      labels(tree_src.fg5.tracing.rowtree.4[[1]][[2]][[2]][[2]][[1]][[2]][[2]][[1]]))),
  labels(tree_src.fg5.tracing.rowtree.4[[1]][[2]][[2]][[1]][[1]]),
  
  # labels(tree_src.fg5.tracing.rowtree.4[[2]][[2]][[1]]),
  rev(
    setdiff(labels(tree_src.fg5.tracing.rowtree.4[[2]][[1]][[2]]),
          labels(tree_src.fg5.tracing.rowtree.4[[2]][[1]][[2]][[1]]))),
  setdiff(labels(tree_src.fg5.tracing.rowtree.4[[2]][[1]][[1]]),
          labels(tree_src.fg5.tracing.rowtree.4[[2]][[1]][[1]][[2]][[2]])))

names(gene_src.fg5.tracing.rowtree.4) = c(
  rep(4, length(labels(tree_src.fg5.tracing.rowtree.4[[2]][[2]][[2]]))),
  
  rep(3, length(c(labels(tree_src.fg5.tracing.rowtree.4[[1]][[1]]),
                  # labels(tree_src.fg5.tracing.rowtree.4[[1]][[2]][[1]][[1]]),
                  labels(tree_src.fg5.tracing.rowtree.4[[1]][[2]][[2]])))),
  
  # rep(7, length(labels(tree_src.fg5.tracing.rowtree.4[[2]][[2]][[1]]))),
  rep(7, length(c(setdiff(labels(tree_src.fg5.tracing.rowtree.4[[2]][[1]][[2]]),
                          labels(tree_src.fg5.tracing.rowtree.4[[2]][[1]][[2]][[1]])),
                  setdiff(labels(tree_src.fg5.tracing.rowtree.4[[2]][[1]][[1]]),
                          labels(tree_src.fg5.tracing.rowtree.4[[2]][[1]][[1]][[2]][[2]]))))))



pdf("figure.v08.07/organ_development_re_v240115/try.fg5.pdf",10,10)
src.fg5.tracing.coltree =
  MyHeatmap(as.matrix(src.fg5.tracing@assays$RNA@data[
    gene_src.fg5.tracing.rowtree.3,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg5.tracing$Time,
                 colors.time.2),
      MyName2Col(src.fg5.tracing$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg5.tracing$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      # MyName2Col(src.fg5.tracing$cluster.v06.26.re_correct,
      #            cluster.endoderm.color.v5),
      MyName2Col(src.fg5.tracing$lineage,
                 color.lineage)
      # MyName2Col(src.fg5.tracing$tree,
      #            color.temp)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "col",
    graph = T)
dev.off()

tree_src.fg5.tracing.coltree = as.dendrogram(src.fg5.tracing.coltree)
src.fg5.tracing$tree = NA
src.fg5.tracing@meta.data[labels(tree_src.fg5.tracing.coltree[[1]][[1]]),]$tree = 1
src.fg5.tracing@meta.data[labels(tree_src.fg5.tracing.coltree[[1]][[2]]),]$tree = 2
src.fg5.tracing@meta.data[labels(tree_src.fg5.tracing.coltree[[2]][[1]]),]$tree = 3
src.fg5.tracing@meta.data[labels(tree_src.fg5.tracing.coltree[[2]][[2]][[1]]),]$tree = 4
src.fg5.tracing@meta.data[labels(tree_src.fg5.tracing.coltree[[2]][[2]][[2]][[1]]),]$tree = 5
src.fg5.tracing@meta.data[labels(tree_src.fg5.tracing.coltree[[2]][[2]][[2]][[2]][[1]]),]$tree = 6
src.fg5.tracing@meta.data[labels(tree_src.fg5.tracing.coltree[[2]][[2]][[2]][[2]][[2]]),]$tree = 7
DimPlot(src.fg5.tracing, group.by = "tree", reduction = 'mnn_umap_fta')


pdf("figure.v08.07/organ_development_re_v240115/try.fg5.pdf",10,10)
src.fg5.tracing.coltree.2 =
  MyHeatmap(as.matrix(src.fg5.tracing@assays$RNA@data[
    gene_src.fg5.tracing.rowtree.4,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg5.tracing$Time,
                 colors.time.2),
      MyName2Col(src.fg5.tracing$cluster.extract.v1.1_define, 
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg5.tracing$cluster.v06.26.re_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg5.tracing$cluster.v06.26.re_correct,
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg5.tracing$lineage,
                 color.lineage)
      # MyName2Col(src.fg5.tracing$tree,
      #            color.temp),
      # MyName2Col(src.fg5.tracing$tree.4,
      #            color.temp)
      ),
    # RowSideColors = t(cbind(
    #   MyName2Col(names(gene_src.fg5.tracing.rowtree.4),
    #              colors.geneset))),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "col",
    graph = T)
dev.off()

tree_src.fg5.tracing.coltree.2 = as.dendrogram(src.fg5.tracing.coltree.2)
src.fg5.tracing$tree.2 = NA
src.fg5.tracing@meta.data[labels(tree_src.fg5.tracing.coltree.2[[1]][[1]]),]$tree.2 = 1
src.fg5.tracing@meta.data[labels(tree_src.fg5.tracing.coltree.2[[1]][[2]]),]$tree.2 = 2
src.fg5.tracing@meta.data[labels(tree_src.fg5.tracing.coltree.2[[2]][[1]]),]$tree.2 = 3
src.fg5.tracing@meta.data[labels(tree_src.fg5.tracing.coltree.2[[2]][[2]][[1]]),]$tree.2 = 4
src.fg5.tracing@meta.data[labels(tree_src.fg5.tracing.coltree.2[[2]][[2]][[2]][[1]]),]$tree.2 = 5
src.fg5.tracing@meta.data[labels(tree_src.fg5.tracing.coltree.2[[2]][[2]][[2]][[2]][[1]]),]$tree.2 = 6
src.fg5.tracing@meta.data[labels(tree_src.fg5.tracing.coltree.2[[2]][[2]][[2]][[2]][[2]]),]$tree.2 = 7

src.fg5.tracing = FindNeighbors(src.fg5.tracing, dims = 1:30)
src.fg5.tracing = FindClusters(src.fg5.tracing, resolution = 1)
DimPlot(src.fg5.tracing, group.by = "tree.2", reduction = 'mnn_umap_fta') +
  DimPlot(src.fg5.tracing, group.by = "tree", reduction = 'mnn_umap_fta') +
  DimPlot(src.fg5.tracing, reduction = "mnn_umap_fta", group.by = "RNA_snn_res.1")

src.fg5.tracing$cluster.v06.26.re_correct = src.fg5.tracing$cluster.v06.26.re_mnn_umap_fta
src.fg5.tracing@meta.data[src.fg5.tracing$tree%in%c(3),]$cluster.v06.26.re_correct = "FG.5"
src.fg5.tracing@meta.data[src.fg5.tracing$tree%in%c(4,5,6,7),]$cluster.v06.26.re_correct = "Thyroid"
src.fg5.tracing@meta.data[src.fg5.tracing$tree%in%c(1,2),]$cluster.v06.26.re_correct = "Pharynx.organ.3"
src.fg5.tracing@meta.data[src.fg5.tracing$cluster.v06.26.re_mnn_umap_fta%in%c("FG.5")&
                            src.fg5.tracing$cluster.v06.26.re_correct%in%c("Thyroid","Pharynx.organ.3"),]$cluster.v06.26.re_correct = "FG.5"
src.fg5.tracing@meta.data[src.fg5.tracing$RNA_snn_res.1%in%c(2,4),]$cluster.v06.26.re_correct = "Pharynx.organ.3"

DimPlot(src.fg5.tracing, group.by = "cluster.v06.26.re_correct", reduction = 'mnn_umap_fta') +
  DimPlot(src.fg5.tracing, group.by = "cluster.v06.26.re_mnn_umap_fta", reduction = 'mnn_umap_fta') +
  DimPlot(src.fg5.tracing, group.by = "tree", reduction = 'mnn_umap_fta')

# Pseudo-time :: cluster.v06.26.re_correct
#---------------------
#-- pha3 --
cell_src.fg5.tracing_pha3 = 
  rownames(src.fg5.tracing@meta.data[src.fg5.tracing$cluster.v06.26.re_correct%in%c("FG.5",'Pharynx.organ.3'),])
coord_src.fg5.tracing_pha3 = src.fg5.tracing[["mnn_umap_fta"]]@cell.embeddings[cell_src.fg5.tracing_pha3,]
pcurve_src.fg5.tracing_pha3 = princurve::principal_curve(x = coord_src.fg5.tracing_pha3, smoother = "smooth.spline")
src.fg5.tracing$lambda_pha3 = pcurve_src.fg5.tracing_pha3$lambda[cell_src.fg5.tracing_pha3]
src.fg5.tracing@meta.data[cell_src.fg5.tracing_pha3,]$lambda_pha3 = 
  pcurve_src.fg5.tracing_pha3$lambda[cell_src.fg5.tracing_pha3]
src.fg5.tracing@meta.data[!src.fg5.tracing$lambda_pha3%in%NA,]$lambda_pha3 = 
  norm_range(src.fg5.tracing@meta.data[!src.fg5.tracing$lambda_pha3%in%NA,]$lambda_pha3)

#-- th --
cell_src.fg5.tracing_th = 
  rownames(src.fg5.tracing@meta.data[src.fg5.tracing$cluster.v06.26.re_correct%in%c("FG.5",'Thyroid'),])
coord_src.fg5.tracing_th = src.fg5.tracing[["mnn_umap_fta"]]@cell.embeddings[cell_src.fg5.tracing_th,]
pcurve_src.fg5.tracing_th = princurve::principal_curve(x = coord_src.fg5.tracing_th, smoother = "smooth.spline")
src.fg5.tracing$lambda_th = pcurve_src.fg5.tracing_th$lambda[cell_src.fg5.tracing_th]
src.fg5.tracing@meta.data[cell_src.fg5.tracing_th,]$lambda_th = 
  pcurve_src.fg5.tracing_th$lambda[cell_src.fg5.tracing_th]
src.fg5.tracing@meta.data[!src.fg5.tracing$lambda_th%in%NA,]$lambda_th = 
  norm_range(src.fg5.tracing@meta.data[!src.fg5.tracing$lambda_th%in%NA,]$lambda_th)

cell_src.fg5.tracing_fg5 = rownames(src.fg5.tracing@meta.data[src.fg5.tracing$cluster.v06.26.re_correct%in%c("FG.5"),])
src.fg5.tracing$lambda_fg5 = NA
w1 = length(cell_src.fg5.tracing_pha3) / (length(cell_src.fg5.tracing_th) + length(cell_src.fg5.tracing_pha3))
w2 = length(cell_src.fg5.tracing_th) / (length(cell_src.fg5.tracing_th) + length(cell_src.fg5.tracing_pha3))
src.fg5.tracing$lambda_fg5 = 
  w1 * src.fg5.tracing$lambda_pha3[cell_src.fg5.tracing_fg5] + 
  w2 * src.fg5.tracing$lambda_th[cell_src.fg5.tracing_fg5]
# src.fg5.tracing$lambda_fg5 = src.fg5.tracing$lambda_pha3[cell_src.fg5.tracing_fg5] 

cellorder_src.fg5.tracing = c(
  intersect(names(src.fg5.tracing$lambda_th[order(src.fg5.tracing$lambda_th)]), 
            colnames(src.fg5.tracing[,src.fg5.tracing$cluster.v06.26.re_correct%in%"FG.5"])),
  intersect(names(src.fg5.tracing$lambda_th[order(src.fg5.tracing$lambda_th)]), 
            colnames(src.fg5.tracing[,src.fg5.tracing$cluster.v06.26.re_correct%in%"Thyroid"])),
  rev(intersect(names(src.fg5.tracing$lambda_pha3[order(src.fg5.tracing$lambda_pha3)]), 
                colnames(src.fg5.tracing[,src.fg5.tracing$cluster.v06.26.re_correct%in%"Pharynx.organ.3"]))))


#---------------------


pdf("figure.v08.07/organ_development_re_v240115/try.fg5.pdf",10,10)
src.fg5.tracing.coltree.5 =
  MyHeatmap(
    as.matrix(src.fg5.tracing@assays$RNA@data[
      unique(c(gene_src.fg5.tracing.rowtree.3,
               gene_src.fg5.tracing.rowtree.4)),
      cellorder_src.fg5.tracing]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg5.tracing$Time[cellorder_src.fg5.tracing],
                 colors.time.2),
      MyName2Col(src.fg5.tracing$lineage[cellorder_src.fg5.tracing],
                 color.lineage),
      MyName2Col(src.fg5.tracing$cluster.v06.26.re_correct[cellorder_src.fg5.tracing],
                 cluster.endoderm.color.v5)
      # MyName2Col(src.fg5.tracing$tree,
      #            color.temp)
    ),
    # RowSideColors = t(cbind(
    #   MyName2Col(names(gene_src.fg5.tracing.rowtree.4),
    #              colors.geneset))),
    ColSideColorsSize = 4.8,
    RowSideColorsSize = 1.2,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    # Colv = "none",
    # Rowv = "none",
    return.tree = "col",
    margins = c(10,10), graph = T)
dev.off()
tree_src.fg5.tracing.rowtree.5 = as.dendrogram(src.fg5.tracing.rowtree.5)
tree_src.fg5.tracing.coltree.5 = as.dendrogram(src.fg5.tracing.coltree.5)

gene_src.fg5.tracing.rowtree.5 = c(
  labels(tree_src.fg5.tracing.rowtree.5[[2]][[1]]),
  
  labels(tree_src.fg5.tracing.rowtree.5[[2]][[2]][[1]][[1]]),
  labels(tree_src.fg5.tracing.rowtree.5[[2]][[2]][[1]][[2]]),
  labels(tree_src.fg5.tracing.rowtree.5[[2]][[2]][[2]][[2]]),
  
  labels(tree_src.fg5.tracing.rowtree.5[[1]][[1]]),
  labels(tree_src.fg5.tracing.rowtree.5[[1]][[2]][[1]]),
  labels(tree_src.fg5.tracing.rowtree.5[[1]][[2]][[2]][[2]]))

names(gene_src.fg5.tracing.rowtree.5) = c(
  rep(4, length(
    labels(tree_src.fg5.tracing.rowtree.5[[2]][[1]])
  )),
  rep(3, length(c(
    labels(tree_src.fg5.tracing.rowtree.5[[2]][[2]][[1]][[1]]),
    labels(tree_src.fg5.tracing.rowtree.5[[2]][[2]][[1]][[2]]),
    labels(tree_src.fg5.tracing.rowtree.5[[2]][[2]][[2]][[2]])
  ))),
  
  rep(7, length(c(
    labels(tree_src.fg5.tracing.rowtree.5[[1]][[1]]),
    labels(tree_src.fg5.tracing.rowtree.5[[1]][[2]][[1]]),
    labels(tree_src.fg5.tracing.rowtree.5[[1]][[2]][[2]][[2]])
  ))))



gene_order = function(gene, type){
  data_temp = src.fg5.tracing@assays$RNA@data[
    gene,
    src.fg5.tracing$cluster.v06.26.re_correct%in%type]
  data_temp = rowSums(data_temp)
  gene = unlist(rev(names(data_temp)[order(data_temp)]))
  return(gene)
}

gene_src.fg5.tracing.rowtree.5.correct = c(
  gene_order(gene_src.fg5.tracing.rowtree.5[
    names(gene_src.fg5.tracing.rowtree.5)%in%4], type = "Pharynx.organ.3"),
  gene_order(gene_src.fg5.tracing.rowtree.5[
    names(gene_src.fg5.tracing.rowtree.5)%in%3], type = "Pharynx.organ.3"),
  gene_order(gene_src.fg5.tracing.rowtree.5[
    names(gene_src.fg5.tracing.rowtree.5)%in%7], type = "FG.5"))

names(gene_src.fg5.tracing.rowtree.5.correct) = c(
  rep(4, length(gene_src.fg5.tracing.rowtree.5[
    names(gene_src.fg5.tracing.rowtree.5)%in%4])),
  rep(3, length(gene_src.fg5.tracing.rowtree.5[
    names(gene_src.fg5.tracing.rowtree.5)%in%3])),
  rep(7, length(gene_src.fg5.tracing.rowtree.5[
    names(gene_src.fg5.tracing.rowtree.5)%in%7]))
)


pdf("figure.v08.07/organ_development_re_v240115/Heatmap_for_trajectory/try.fg5.re.pdf",10,10)
#-----------
src.fg5.tracing.rowtree.5 =
  MyHeatmap(
    as.matrix(src.fg5.tracing@assays$RNA@data[
      # gene_src.fg5.tracing.rowtree.5,
      gene_src.fg5.tracing.rowtree.5.correct,
      cellorder_src.fg5.tracing]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg5.tracing$Time[cellorder_src.fg5.tracing],
                 colors.time.2),
      MyName2Col(src.fg5.tracing$lineage[cellorder_src.fg5.tracing],
                 color.lineage),
      MyName2Col(src.fg5.tracing$cluster.v06.26.re_correct[cellorder_src.fg5.tracing],
                 cluster.endoderm.color.v5)
      # MyName2Col(src.fg5.tracing$tree,
      #            color.temp)
    ),
    RowSideColors = t(cbind(
      MyName2Col(names(gene_src.fg5.tracing.rowtree.5.correct),
                 colors.geneset))),
    ColSideColorsSize = 4.8,
    RowSideColorsSize = 1.2,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    Colv = "none",
    Rowv = "none",
    # return.tree = "row",
    margins = c(10,10), graph = T)
#-----------
dev.off()

save(src.fg5.tracing,
     cell_src.fg5.tracing_fg5,
     gene_src.fg5.tracing.rowtree.5,
     file = "figure.v08.07/organ_development_re_v240115/src.fg5.tracing.parameter.Rdata")















