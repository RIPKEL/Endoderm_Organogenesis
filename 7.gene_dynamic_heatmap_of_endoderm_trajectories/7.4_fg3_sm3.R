#-----------
# FG.3-Re
#-----------
src.fg3.tracing.re = src.fg3.integrated.merge[,!src.fg3.integrated.merge$lineage%in%NA]
src.fg3.tracing.re@meta.data[
  src.fg3.tracing.re$cluster.extract.v1.1_define%in%"FG.3"&
    src.fg3.tracing.re$Time%in%"9ss",]$cluster.v06.26.re_correct_refine_mnn_umap_fta = "FG.3"

src.fg3.tracing.re = FindVariableFeatures(src.fg3.tracing.re, nfeatures = 2000)
src.fg3.tracing.re.filtergene = 
  Myfilter(as.matrix(src.fg3.tracing.re@assays$RNA@data),
           gene = src.fg3.tracing.re@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)

src.fg3.tracing.re = ScaleData(src.fg3.tracing.re, rownames(src.fg3.tracing.re))
src.fg3.tracing.re = RunPCA(src.fg3.tracing.re, features = src.fg3.tracing.re.filtergene)
src.fg3.tracing.re = RunUMAP(src.fg3.tracing.re, dims = 1:30, reduction = "pca", assay = "RNA", n.components = 2)

src.fg3.tracing.re[["mnn_umap_fta"]] = src.fg3.tracing.re[["umap"]]
src.fg3.tracing.re@reductions$mnn_umap_fta@key = "Coord_"
src.fg3.tracing.re@reductions$mnn_umap_fta@cell.embeddings = 
  src.fg3.integrated.merge@reductions$mnn_umap_fta@cell.embeddings[colnames(src.fg3.tracing.re),]
colnames(src.fg3.tracing.re@reductions$mnn_umap_fta@cell.embeddings) = c("Coord_1","Coord_2")

DimPlot(src.fg3.tracing.re, reduction = "mnn_umap_fta",
        group.by = "cluster.v06.26.re_correct_refine_mnn_umap_fta", 
        cols = cluster.endoderm.color.v5) 
DimPlot(src.fg3.tracing.re, reduction = "mnn_umap_fta",
        group.by = "Time", cols = colors.time.2)

# == FG.3 Basic process 
#===============================================================================
src.fg3.tracing.re = SetIdent(src.fg3.tracing.re, value = src.fg3.tracing.re$cluster.v06.26.re_correct_refine_mnn_umap_fta)
marker_src.fg3.tracing.re = FindAllMarkers(src.fg3.tracing.re)
marker_src.fg3.tracing.re$pct.ratio = marker_src.fg3.tracing.re$pct.1 / marker_src.fg3.tracing.re$pct.2
marker_src.fg3.tracing.re$rank = marker_src.fg3.tracing.re$pct.ratio * (-log(marker_src.fg3.tracing.re$p_val_adj))
marker_src.fg3.tracing.re = marker_src.fg3.tracing.re[order(marker_src.fg3.tracing.re$rank, decreasing = T),]
markergene_src.fg3.tracing.re = unique(marker_src.fg3.tracing.re$gene)
# markergene_src.fg3.tracing.re.raw = markergene_src.fg3.tracing.re

pdf("figure.v08.07/organ_development_re_v240115/try.fg3.pdf",9,7)
src.fg3.tracing.re.rowtree =
  MyHeatmap(as.matrix(src.fg3.tracing.re@assays$RNA@data[
    unique(c(markergene_src.fg3.tracing.re,
             # src.fg3.tracing.re.filtergene,
             c()
    )),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg3.tracing.re$Time,
                 colors.time.2),
      MyName2Col(src.fg3.tracing.re$cluster.v06.26.re_correct_refine_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg3.tracing.re$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()

pdf("figure.v08.07/organ_development_re_v240115/try.fg3.tree.pdf",150,20)
plot(src.fg3.tracing.re.rowtree)
dev.off()

tree_src.fg3.tracing.re.rowtree = as.dendrogram(src.fg3.tracing.re.rowtree)
gene_src.fg3.tracing.re.rowtree = c(
  labels(tree_src.fg3.tracing.re.rowtree[[2]][[1]]),
  
  # labels(tree_src.fg3.tracing.re.rowtree[[2]][[2]][[2]][[2]][[2]][[2]][[1]]),
  labels(tree_src.fg3.tracing.re.rowtree[[2]][[2]][[2]][[2]][[1]]),
  labels(tree_src.fg3.tracing.re.rowtree[[2]][[2]][[2]][[2]][[2]][[2]][[2]]),
  labels(tree_src.fg3.tracing.re.rowtree[[2]][[2]][[2]][[2]][[2]][[1]]))
names(gene_src.fg3.tracing.re.rowtree) = c(
  rep(5, length(labels(tree_src.fg3.tracing.re.rowtree[[2]][[1]]))),
  rep(7, length(
    setdiff(labels(tree_src.fg3.tracing.re.rowtree[[2]][[2]][[2]][[2]]),
            c(
              labels(tree_src.fg3.tracing.re.rowtree[[2]][[2]][[2]][[2]][[2]][[2]][[1]])
            )))))

gene_src.fg3.tracing.re.rowtree.con = setdiff(
  unique(marker_src.fg3.tracing.re[marker_src.fg3.tracing.re$cluster%in%c("FG.3", "Pharynx.organ.4"), "gene"]),
  gene_src.fg3.tracing.re.rowtree)

pdf("figure.v08.07/organ_development_re_v240115/try.fg3.pdf",9,7)
src.fg3.tracing.re.rowtree.1 =
  MyHeatmap(as.matrix(src.fg3.tracing.re@assays$RNA@data[
    gene_src.fg3.tracing.re.rowtree.con,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg3.tracing.re$Time,
                 colors.time.2),
      MyName2Col(src.fg3.tracing.re$cluster.v06.26.re_correct_refine_mnn_umap_fta,
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg3.tracing.re$lineage,
                 color.lineage)
    ),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()
tree_src.fg3.tracing.re.rowtree.1 = as.dendrogram(src.fg3.tracing.re.rowtree.1)

gene_src.fg3.tracing.re.rowtree.1 = c(
  labels(tree_src.fg3.tracing.re.rowtree.1[[2]][[2]][[1]][[2]]),
  
  labels(tree_src.fg3.tracing.re.rowtree.1[[2]][[2]][[2]][[2]][[2]]),
  rev(c(
    labels(tree_src.fg3.tracing.re.rowtree.1[[2]][[2]][[2]][[1]]),
    labels(tree_src.fg3.tracing.re.rowtree.1[[2]][[2]][[2]][[2]][[1]]),
    labels(tree_src.fg3.tracing.re.rowtree.1[[2]][[1]])
  )))
names(gene_src.fg3.tracing.re.rowtree.1) = c(
  rep(4, length(labels(tree_src.fg3.tracing.re.rowtree.1[[2]][[2]][[1]][[2]]))),
  rep(3, length(c(
    labels(tree_src.fg3.tracing.re.rowtree.1[[2]][[2]][[2]]),
    labels(tree_src.fg3.tracing.re.rowtree.1[[2]][[1]])
  ))))

gene_src.fg3.tracing.re.rowtree.2 = c(
  gene_src.fg3.tracing.re.rowtree.1,
  gene_src.fg3.tracing.re.rowtree)


# Pseudo-time :: cluster.v06.26.re_correct_refine_mnn_umap_fta
#---------------------
#-- pha4 --
cell_src.fg3.tracing.re_pha4 = 
  rownames(src.fg3.tracing.re@meta.data[src.fg3.tracing.re$cluster.v06.26.re_correct_refine_mnn_umap_fta%in%c("FG.3",'Pharynx.organ.4'),])
coord_src.fg3.tracing.re_pha4 = src.fg3.tracing.re[["mnn_umap_fta"]]@cell.embeddings[cell_src.fg3.tracing.re_pha4,]
pcurve_src.fg3.tracing.re_pha4 = princurve::principal_curve(x = coord_src.fg3.tracing.re_pha4, smoother = "smooth.spline")
src.fg3.tracing.re$lambda_pha4 = pcurve_src.fg3.tracing.re_pha4$lambda[cell_src.fg3.tracing.re_pha4]
src.fg3.tracing.re@meta.data[cell_src.fg3.tracing.re_pha4,]$lambda_pha4 = 
  pcurve_src.fg3.tracing.re_pha4$lambda[cell_src.fg3.tracing.re_pha4]
src.fg3.tracing.re@meta.data[!src.fg3.tracing.re$lambda_pha4%in%NA,]$lambda_pha4 = 
  norm_range(src.fg3.tracing.re@meta.data[!src.fg3.tracing.re$lambda_pha4%in%NA,]$lambda_pha4)

#-- pha5 --
cell_src.fg3.tracing.re_pha5 = 
  rownames(src.fg3.tracing.re@meta.data[src.fg3.tracing.re$cluster.v06.26.re_correct_refine_mnn_umap_fta%in%c('Pharynx.organ.5'),])
coord_src.fg3.tracing.re_pha5 = src.fg3.tracing.re[["mnn_umap_fta"]]@cell.embeddings[cell_src.fg3.tracing.re_pha5, c(2,1)]
pcurve_src.fg3.tracing.re_pha5 = princurve::principal_curve(x = coord_src.fg3.tracing.re_pha5, smoother = "smooth.spline")
src.fg3.tracing.re$lambda_pha5 = pcurve_src.fg3.tracing.re_pha5$lambda[cell_src.fg3.tracing.re_pha5]
src.fg3.tracing.re@meta.data[cell_src.fg3.tracing.re_pha5,]$lambda_pha5 = 
  pcurve_src.fg3.tracing.re_pha5$lambda[cell_src.fg3.tracing.re_pha5]
src.fg3.tracing.re@meta.data[!src.fg3.tracing.re$lambda_pha5%in%NA,]$lambda_pha5 = 
  norm_range(src.fg3.tracing.re@meta.data[!src.fg3.tracing.re$lambda_pha5%in%NA,]$lambda_pha5)

#-- lung --
cell_src.fg3.tracing.re_lung = 
  rownames(src.fg3.tracing.re@meta.data[src.fg3.tracing.re$cluster.v06.26.re_correct_refine_mnn_umap_fta%in%c('Lung'),])
coord_src.fg3.tracing.re_lung = src.fg3.tracing.re[["mnn_umap_fta"]]@cell.embeddings[cell_src.fg3.tracing.re_lung,]
pcurve_src.fg3.tracing.re_lung = princurve::principal_curve(x = coord_src.fg3.tracing.re_lung, smoother = "smooth.spline")
src.fg3.tracing.re$lambda_lung = pcurve_src.fg3.tracing.re_lung$lambda[cell_src.fg3.tracing.re_lung]
src.fg3.tracing.re@meta.data[cell_src.fg3.tracing.re_lung,]$lambda_lung = 
  pcurve_src.fg3.tracing.re_lung$lambda[cell_src.fg3.tracing.re_lung]
src.fg3.tracing.re@meta.data[!src.fg3.tracing.re$lambda_lung%in%NA,]$lambda_lung = 
  norm_range(src.fg3.tracing.re@meta.data[!src.fg3.tracing.re$lambda_lung%in%NA,]$lambda_lung)


cellorder_src.fg3.tracing.re = c(
  intersect(names(src.fg3.tracing.re$lambda_pha4[order(src.fg3.tracing.re$lambda_pha4)]), 
            colnames(src.fg3.tracing.re[,src.fg3.tracing.re$cluster.v06.26.re_correct_refine_mnn_umap_fta%in%"FG.3"])),
  rev(intersect(names(src.fg3.tracing.re$lambda_pha4[order(src.fg3.tracing.re$lambda_pha4)]), 
            colnames(src.fg3.tracing.re[,src.fg3.tracing.re$cluster.v06.26.re_correct_refine_mnn_umap_fta%in%"Pharynx.organ.4"]))),
  intersect(names(src.fg3.tracing.re$lambda_pha5[order(src.fg3.tracing.re$lambda_pha5)]), 
            colnames(src.fg3.tracing.re[,src.fg3.tracing.re$cluster.v06.26.re_correct_refine_mnn_umap_fta%in%"Pharynx.organ.5"])),
  intersect(names(src.fg3.tracing.re$lambda_lung[order(src.fg3.tracing.re$lambda_lung)]), 
            colnames(src.fg3.tracing.re[,src.fg3.tracing.re$cluster.v06.26.re_correct_refine_mnn_umap_fta%in%"Lung"])))
#---------------------

save(src.fg3.tracing.re,
     gene_src.fg3.tracing.re.rowtree,
     gene_src.fg3.tracing.re.rowtree.1,
     gene_src.fg3.tracing.re.rowtree.2,
     gene_src.fg3.tracing.re.rowtree.con,
     file = "figure.v08.07/organ_development_re_v240115/src.fg3.tracing.re.parameter.Rdata")

pdf("figure.v08.07/organ_development_re_v240115//Heatmap_for_trajectory/try.fg3.pdf",10,10)
src.fg3.tracing.re.rowtree.2 =
  MyHeatmap(as.matrix(src.fg3.tracing.re@assays$RNA@data[
    gene_src.fg3.tracing.re.rowtree.2,
    cellorder_src.fg3.tracing.re]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg3.tracing.re$Time[cellorder_src.fg3.tracing.re],
                 colors.time.2),
      MyName2Col(src.fg3.tracing.re$lineage[cellorder_src.fg3.tracing.re],
                 color.lineage),
      MyName2Col(src.fg3.tracing.re$cluster.v06.26.re_correct_refine_mnn_umap_fta[cellorder_src.fg3.tracing.re],
                 cluster.endoderm.color.v5)
    ),
    RowSideColors = t(cbind(
      MyName2Col(names(gene_src.fg3.tracing.re.rowtree.2),
                 colors.geneset))),
    ColSideColorsSize = 4.8,
    RowSideColorsSize = 1.2,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    Rowv = "none",
    Colv = "none",
    # return.tree = "row",
    margins = c(10,10), 
    graph = T)
dev.off()




