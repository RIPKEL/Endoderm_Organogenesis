#---------------------------------
# -- AL.3 Intermeditae Stage
#=================================================================================
src.al3.integrated_re_al3 = 
  src.al3.integrated[, src.al3.integrated$cluster.v06.26.re_correct%in%c("AL.3",
                                                                         "AL.3-EHBD/VP",
                                                                         "AL.3-Liver",
                                                                         "AL.3-Small.intestine.1")]

src.al3.integrated_re_al3 = NormalizeData(src.al3.integrated_re_al3, scale.factor = 10^5)
src.al3.integrated_re_al3 = ScaleData(src.al3.integrated_re_al3, 
                                      features = rownames(src.al3.integrated_re_al3),
                                      split.by = "batch_phase")
src.al3.integrated_re_al3 = FindVariableFeatures(src.al3.integrated_re_al3, assay = "RNA", nfeatures = 2000)
src.al3.integrated_re_al3_filtergene = 
  Myfilter(as.matrix(src.al3.integrated_re_al3@assays$RNA@data),
           gene = src.al3.integrated_re_al3@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.al3.integrated_re_al3_filtergene_re = src.al3.integrated_re_al3_filtergene
rm(src.al3.integrated_re_al3_filtergene)
src.al3.integrated_re_al3_filtergene = src.al3.integrated_re_al3_filtergene_re

# -- reduction
src.al3.integrated_re_al3 = RunPCA(src.al3.integrated_re_al3, 
                                   features = src.al3.integrated_re_al3_filtergene)
src.al3.integrated_re_al3 = RunUMAP(src.al3.integrated_re_al3, reduction = "pca",
                                    dims = 1:25, n.components = 3)

pdf("figure.v08.07/organ_development_re_v240115/try.al3_al3_pca.pdf",10,10)
DimPlot(src.al3.integrated_re_al3, group.by = "Time", 
        cols = colors.time, pt.size = 1.25,
        reduction = "pca", dims = c(1,2)) +
  theme_bw() + p_add
DimPlot(src.al3.integrated_re_al3, group.by = "cluster.v06.26.re_correct",
        cols = cluster.endoderm.color.v5, pt.size = 1.25,
        reduction = "pca", dims = c(1,3)) +
  theme_bw() + p_add
DimPlot(src.al3.integrated_re_al3, group.by = "cluster.v06.26.re_correct",
        cols = cluster.endoderm.color.v5, pt.size = 1.25,
        reduction = "umap", dims = c(1,2)) +
  theme_bw() + p_add
dev.off()


load_src.al3.integrated_re_al3_pca1 = abs(src.al3.integrated_re_al3@reductions$pca@feature.loadings[,1])
gene_src.al3.integrated_re_al3_pca1 = 
  load_src.al3.integrated_re_al3_pca1[order(load_src.al3.integrated_re_al3_pca1,decreasing=T)][1:150]

load_src.al3.integrated_re_al3_pca2 = abs(src.al3.integrated_re_al3@reductions$pca@feature.loadings[,3])
gene_src.al3.integrated_re_al3_pca2 = 
  load_src.al3.integrated_re_al3_pca2[order(load_src.al3.integrated_re_al3_pca2,decreasing=T)][1:150]

load_src.al3.integrated_re_al3_pca3 = abs(src.al3.integrated_re_al3@reductions$pca@feature.loadings[,3])
gene_src.al3.integrated_re_al3_pca3 = 
  load_src.al3.integrated_re_al3_pca3[order(load_src.al3.integrated_re_al3_pca3,decreasing=T)][1:150]


# -- marker
src.al3.integrated_re_al3 = SetIdent(src.al3.integrated_re_al3, 
                                     value = src.al3.integrated_re_al3$cluster.v06.26.re_correct)
src.al3.integrated_re_al3_marker = FindAllMarkers(src.al3.integrated_re_al3)
src.al3.integrated_re_al3_marker$pct.ratio =
  src.al3.integrated_re_al3_marker$pct.1 / src.al3.integrated_re_al3_marker$pct.2
src.al3.integrated_re_al3_markergene = unique(
  src.al3.integrated_re_al3_marker[
    src.al3.integrated_re_al3_marker$pct.ratio>2,]$gene)


#-- Basic processing :: Filter & Marker
#-------------------------------------------------------------------------------
pdf("figure.v08.07/organ_development_re_v240115/try.al3_al3.pdf",9,7)
src.al3.integrated_re_al3.rowtree =
  MyHeatmap(as.matrix(src.al3.integrated_re_al3@assays$RNA@data[
    unique(c(
      src.al3.integrated_re_al3_markergene,
      src.al3.integrated_re_al3_filtergene)),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.al3.integrated_re_al3$Time, colors.time),
      MyName2Col(src.al3.integrated_re_al3$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src.al3.integrated_re_al3$cluster.v06.26.re, cluster.endoderm.color.v5)),
    ColSideColorsSize = 5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()
pdf("figure.v08.07/organ_development_re_v240115/try.al3_al3_tree.pdf",100,20)
plot(src.al3.integrated_re_al3.rowtree)
dev.off()

tree_src.al3.integrated_re_al3.rowtree = as.dendrogram(src.al3.integrated_re_al3.rowtree)
gene_src.al3.integrated_re_al3.rowtree = setdiff(
  labels(tree_src.al3.integrated_re_al3.rowtree),
  c(labels(tree_src.al3.integrated_re_al3.rowtree[[2]][[1]]),
    labels(tree_src.al3.integrated_re_al3.rowtree[[2]][[2]][[2]][[2]])))


pdf("figure.v08.07/organ_development_re_v240115/try.al3_al3.pdf",9,7)
src.al3.integrated_re_al3.rowtree.1 =
  MyHeatmap(as.matrix(src.al3.integrated_re_al3@assays$RNA@data[
    gene_src.al3.integrated_re_al3.rowtree, ]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.al3.integrated_re_al3$Time, colors.time),
      MyName2Col(src.al3.integrated_re_al3$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src.al3.integrated_re_al3$cluster.v06.26.re, cluster.endoderm.color.v5)),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    return.tree = "row",
    graph = T)
src.al3.integrated_re_al3.coltree.1 =
  MyHeatmap(as.matrix(src.al3.integrated_re_al3@assays$RNA@data[
    gene_src.al3.integrated_re_al3.rowtree, ]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    return.tree = "col",
    graph = F)
dev.off()

tree_src.al3.integrated_re_al3.rowtree.1 = as.dendrogram(src.al3.integrated_re_al3.rowtree.1)
tree_src.al3.integrated_re_al3.coltree.1 = as.dendrogram(src.al3.integrated_re_al3.coltree.1)
gene_src.al3.integrated_re_al3.rowtree.1 = setdiff(
  labels(tree_src.al3.integrated_re_al3.rowtree.1),
  c(labels(tree_src.al3.integrated_re_al3.rowtree.1[[2]][[2]][[1]][[2]][[2]][[2]][[2]]),
    labels(tree_src.al3.integrated_re_al3.rowtree.1[[2]][[1]][[1]][[2]][[2]][[2]][[2]]),
    labels(tree_src.al3.integrated_re_al3.rowtree.1[[2]][[1]][[1]][[1]][[1]][[1]]),
    labels(tree_src.al3.integrated_re_al3.rowtree.1[[1]][[2]][[2]][[2]][[1]][[1]]),
    labels(tree_src.al3.integrated_re_al3.rowtree.1[[1]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]])))



pdf("figure.v08.07/organ_development_re_v240115/try.al3_al3.pdf",9,7)
src.al3.integrated_re_al3.rowtree.2 =
  MyHeatmap(as.matrix(src.al3.integrated_re_al3@assays$RNA@data[
    gene_src.al3.integrated_re_al3.rowtree.1, ]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.al3.integrated_re_al3$Time, colors.time),
      # MyName2Col(src.al3.integrated_re_al3$tree.2, color.temp),
      MyName2Col(src.al3.integrated_re_al3$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src.al3.integrated_re_al3$cluster.v06.26.re, cluster.endoderm.color.v5),
      MyName2Col(src.al3.integrated_re_al3$cluster.v06.26.re_correct, cluster.endoderm.color.v5)),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    return.tree = "row",
    graph = T)
src.al3.integrated_re_al3.coltree.2 =
  MyHeatmap(as.matrix(src.al3.integrated_re_al3@assays$RNA@data[
    gene_src.al3.integrated_re_al3.rowtree.1, ]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    return.tree = "col",
    graph = F)
dev.off()

# -- Row 2
tree_src.al3.integrated_re_al3.rowtree.2 = as.dendrogram(src.al3.integrated_re_al3.rowtree.2)
gene_src.al3.integrated_re_al3.rowtree.2 = c(
  # AL.3
  rev(labels(tree_src.al3.integrated_re_al3.rowtree.2[[2]][[2]][[1]])),
  # labels(tree_src.al3.integrated_re_al3.rowtree.2[[2]][[2]][[2]][[2]][[2]]),
  
  # AL.3-Small.intestine.1
  labels(tree_src.al3.integrated_re_al3.rowtree.2[[2]][[1]][[1]]),
  labels(tree_src.al3.integrated_re_al3.rowtree.2[[2]][[1]][[2]][[2]][[2]]),
  # labels(tree_src.al3.integrated_re_al3.rowtree.2[[2]][[1]][[2]][[2]][[1]]),
  
  # AL.3-EHBD/VP
  labels(tree_src.al3.integrated_re_al3.rowtree.2[[2]][[2]][[2]][[1]]),
  
  # AL.3-Liver
  rev(c(labels(tree_src.al3.integrated_re_al3.rowtree.2[[1]][[2]][[1]]),
        labels(tree_src.al3.integrated_re_al3.rowtree.2[[1]][[2]][[2]][[2]][[2]][[2]][[2]])
        # labels(tree_src.al3.integrated_re_al3.rowtree.2[[1]][[1]])
        )))

names(gene_src.al3.integrated_re_al3.rowtree.2) = c(
  rep(4, length(c(labels(tree_src.al3.integrated_re_al3.rowtree.2[[2]][[2]][[1]])))),
  
  rep(3, length(c(labels(tree_src.al3.integrated_re_al3.rowtree.2[[2]][[1]][[1]]),
                  labels(tree_src.al3.integrated_re_al3.rowtree.2[[2]][[1]][[2]][[2]][[2]])))),
  
  rep(9, length(c(labels(tree_src.al3.integrated_re_al3.rowtree.2[[2]][[2]][[2]][[1]])))),
  
  rep(7, length(c(labels(tree_src.al3.integrated_re_al3.rowtree.2[[1]][[2]][[1]]),
                  labels(tree_src.al3.integrated_re_al3.rowtree.2[[1]][[2]][[2]][[2]][[2]][[2]][[2]])
                  #labels(tree_src.al3.integrated_re_al3.rowtree.2[[1]][[1]])
                  ))))

# -- Col 2
tree_src.al3.integrated_re_al3.coltree.2 = as.dendrogram(src.al3.integrated_re_al3.coltree.2)
src.al3.integrated_re_al3$tree.2 =NA
src.al3.integrated_re_al3@meta.data[labels(tree_src.al3.integrated_re_al3.coltree.2[[1]][[1]][[1]]), ]$tree.2 = 1
src.al3.integrated_re_al3@meta.data[labels(tree_src.al3.integrated_re_al3.coltree.2[[1]][[1]][[2]]), ]$tree.2 = 2
src.al3.integrated_re_al3@meta.data[labels(tree_src.al3.integrated_re_al3.coltree.2[[1]][[2]][[1]]), ]$tree.2 = 3
src.al3.integrated_re_al3@meta.data[labels(tree_src.al3.integrated_re_al3.coltree.2[[1]][[2]][[2]]), ]$tree.2 = 4
src.al3.integrated_re_al3@meta.data[labels(tree_src.al3.integrated_re_al3.coltree.2[[2]][[1]][[1]]), ]$tree.2 = 5
src.al3.integrated_re_al3@meta.data[labels(tree_src.al3.integrated_re_al3.coltree.2[[2]][[1]][[2]]), ]$tree.2 = 6
src.al3.integrated_re_al3@meta.data[labels(tree_src.al3.integrated_re_al3.coltree.2[[2]][[2]][[1]]), ]$tree.2 = 7
src.al3.integrated_re_al3@meta.data[labels(tree_src.al3.integrated_re_al3.coltree.2[[2]][[2]][[2]]), ]$tree.2 = 8
DimPlot(src.al3.integrated_re_al3, reduction = 'umap_integrated', group.by = "tree.2", label = T, label.size = 5)

src.al3.integrated_re_al3$cluster.v06.26.re_correct_re = src.al3.integrated_re_al3$cluster.v06.26.re_correct
src.al3.integrated_re_al3@meta.data[src.al3.integrated_re_al3$tree.2%in%c(7), ]$cluster.v06.26.re_correct_re = "AL.3"
src.al3.integrated_re_al3@meta.data[src.al3.integrated_re_al3$tree.2%in%c(1,2,3,4), ]$cluster.v06.26.re_correct_re = "AL.3-Liver"
src.al3.integrated_re_al3@meta.data[src.al3.integrated_re_al3$tree.2%in%c(5,6), ]$cluster.v06.26.re_correct_re = "AL.3-Small.intestine.1"
src.al3.integrated_re_al3@meta.data[src.al3.integrated_re_al3$tree.2%in%c(8), ]$cluster.v06.26.re_correct_re = "AL.3-EHBD/VP"
src.al3.integrated_re_al3@meta.data[src.al3.integrated_re_al3$Time%in%"ss9",]$cluster.v06.26.re_correct_re = "AL.3"

src.al3.integrated_re_al3 = RunPCA(src.al3.integrated_re_al3, 
                                   features = gene_src.al3.integrated_re_al3.rowtree.2)
src.al3.integrated_re_al3 = RunUMAP(src.al3.integrated_re_al3, reduction = "pca",
                                    dims = 1:20, n.components = 2)

DimPlot(src.al3.integrated_re_al3, reduction = "umap", dims = c(1,2),
        group.by = "Time", cols = colors.time)
DimPlot(src.al3.integrated_re_al3, reduction = "umap", dims = c(1,2),
        group.by = "cluster.v06.26.re_correct_re", cols = cluster.endoderm.color.v5)

DimPlot(src.al3.integrated_re_al3, reduction = "umap_integrated", dims = c(1,2),
        group.by = "Time", cols = colors.time)
DimPlot(src.al3.integrated_re_al3, reduction = "umap_integrated", dims = c(1,2),
        group.by = "Time", cols = colors.time)
DimPlot(src.al3.integrated_re_al3, reduction = "umap_mnn", dims = c(1,2),
        group.by = "cluster.v06.26.re_correct_re", cols = cluster.endoderm.color.v5)
DimPlot(src.al3.integrated_re_al3, reduction = "umap_integrated", dims = c(1,2),
        group.by = "cluster.v06.26.re_correct_re", cols = cluster.endoderm.color.v5)

#-- Pesudo-time :: cluster.v06.26.re_correct_re
#------------------------------------------------
#-- EHBD/VP --
cell_src.al3.integrated_re_al3_ev = 
  rownames(src.al3.integrated_re_al3@meta.data[src.al3.integrated_re_al3$cluster.v06.26.re_correct_re%in%c("AL.3-EHBD/VP"),])
coord_src.al3.integrated_re_al3_ev = cbind(
  src.al3.integrated_re_al3[["umap_integrated"]]@cell.embeddings[cell_src.al3.integrated_re_al3_ev, 2],
  src.al3.integrated_re_al3[["umap_integrated"]]@cell.embeddings[cell_src.al3.integrated_re_al3_ev, 2])
pcurve_src.al3.integrated_re_al3_ev = princurve::principal_curve(x = coord_src.al3.integrated_re_al3_ev, smoother = "smooth.spline")
src.al3.integrated_re_al3$lambda_ev = NA
src.al3.integrated_re_al3@meta.data[cell_src.al3.integrated_re_al3_ev,]$lambda_ev = 
  pcurve_src.al3.integrated_re_al3_ev$lambda[cell_src.al3.integrated_re_al3_ev]
src.al3.integrated_re_al3@meta.data[!src.al3.integrated_re_al3$lambda_ev%in%NA,]$lambda_ev = 
  norm_range(src.al3.integrated_re_al3@meta.data[!src.al3.integrated_re_al3$lambda_ev%in%NA,]$lambda_ev)

#-- Liver --
cell_src.al3.integrated_re_al3_lv = 
  rownames(src.al3.integrated_re_al3@meta.data[src.al3.integrated_re_al3$cluster.v06.26.re_correct_re%in%c("AL.3", "AL.3-Liver"),])
coord_src.al3.integrated_re_al3_lv = cbind(
  src.al3.integrated_re_al3[["umap_integrated"]]@cell.embeddings[cell_src.al3.integrated_re_al3_lv, 1],
  src.al3.integrated_re_al3[["umap_integrated"]]@cell.embeddings[cell_src.al3.integrated_re_al3_lv, 2])
pcurve_src.al3.integrated_re_al3_lv = princurve::principal_curve(x = coord_src.al3.integrated_re_al3_lv, smoother = "smooth.spline")
src.al3.integrated_re_al3$lambda_lv = NA
src.al3.integrated_re_al3@meta.data[cell_src.al3.integrated_re_al3_lv,]$lambda_lv = 
  pcurve_src.al3.integrated_re_al3_lv$lambda[cell_src.al3.integrated_re_al3_lv]
src.al3.integrated_re_al3@meta.data[!src.al3.integrated_re_al3$lambda_lv%in%NA,]$lambda_lv = 
  norm_range(src.al3.integrated_re_al3@meta.data[!src.al3.integrated_re_al3$lambda_lv%in%NA,]$lambda_lv)

#-- Small intestine 1 --
cell_src.al3.integrated_re_al3_sm = 
  rownames(src.al3.integrated_re_al3@meta.data[src.al3.integrated_re_al3$cluster.v06.26.re_correct_re%in%c("AL.3", "AL.3-Small.intestine.1"),])
coord_src.al3.integrated_re_al3_sm = cbind(
  -src.al3.integrated_re_al3[["umap_integrated"]]@cell.embeddings[cell_src.al3.integrated_re_al3_sm, 1],
  src.al3.integrated_re_al3[["umap_integrated"]]@cell.embeddings[cell_src.al3.integrated_re_al3_sm, 2])
pcurve_src.al3.integrated_re_al3_sm = princurve::principal_curve(x = coord_src.al3.integrated_re_al3_sm, smoother = "smooth.spline")
src.al3.integrated_re_al3$lambda_sm = NA
src.al3.integrated_re_al3@meta.data[cell_src.al3.integrated_re_al3_sm,]$lambda_sm = 
  pcurve_src.al3.integrated_re_al3_sm$lambda[cell_src.al3.integrated_re_al3_sm]
src.al3.integrated_re_al3@meta.data[!src.al3.integrated_re_al3$lambda_sm%in%NA,]$lambda_sm = 
  norm_range(src.al3.integrated_re_al3@meta.data[!src.al3.integrated_re_al3$lambda_sm%in%NA,]$lambda_sm)

#-- AL.3 --
cell_src.al3.integrated_re_al3_al3 = 
  rownames(src.al3.integrated_re_al3@meta.data[src.al3.integrated_re_al3$cluster.v06.26.re_correct_re%in%c("AL.3"),])
w1 = length(cell_src.al3.integrated_re_al3_lv)
w2 = length(cell_src.al3.integrated_re_al3_sm)
src.al3.integrated_re_al3$lambda_al3 = NA
src.al3.integrated_re_al3$lambda_al3[cell_src.al3.integrated_re_al3_al3] = 
  (w1 * src.al3.integrated_re_al3$lambda_lv[cell_src.al3.integrated_re_al3_al3] +
   w2 * src.al3.integrated_re_al3$lambda_sm[cell_src.al3.integrated_re_al3_al3]) / (w1 + w2)
src.al3.integrated_re_al3@meta.data[!src.al3.integrated_re_al3$lambda_al3%in%NA,]$lambda_al3 = 
  norm_range(src.al3.integrated_re_al3@meta.data[!src.al3.integrated_re_al3$lambda_al3%in%NA,]$lambda_al3)


cellorder_src.al3.integrated_re_al3 = c(
  intersect(names(src.al3.integrated_re_al3$lambda_al3[order(src.al3.integrated_re_al3$lambda_al3)]), 
            colnames(src.al3.integrated_re_al3[,src.al3.integrated_re_al3$cluster.v06.26.re_correct_re%in%c("AL.3")])),
  intersect(names(src.al3.integrated_re_al3$lambda_sm[order(src.al3.integrated_re_al3$lambda_sm)]), 
            colnames(src.al3.integrated_re_al3[,src.al3.integrated_re_al3$cluster.v06.26.re_correct_re%in%"AL.3-Small.intestine.1"])),
  intersect(names(src.al3.integrated_re_al3$lambda_ev[order(src.al3.integrated_re_al3$lambda_ev)]), 
            colnames(src.al3.integrated_re_al3[,src.al3.integrated_re_al3$cluster.v06.26.re_correct_re%in%c("AL.3-EHBD/VP"),])),
  intersect(names(src.al3.integrated_re_al3$lambda_lv[order(src.al3.integrated_re_al3$lambda_lv)]), 
            colnames(src.al3.integrated_re_al3[,src.al3.integrated_re_al3$cluster.v06.26.re_correct_re%in%c('AL.3-Liver'),])))
#------------------------------------------------

pdf("figure.v08.07/organ_development_re_v240115/try.al3_al3.pdf",9,7)
src.al3.integrated_re_al3.rowtree.3 =
  MyHeatmap(as.matrix(src.al3.integrated_re_al3@assays$RNA@data[
    gene_src.al3.integrated_re_al3.rowtree.2, 
    cellorder_src.al3.integrated_re_al3]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.al3.integrated_re_al3$Time[cellorder_src.al3.integrated_re_al3], 
                 colors.time),
      MyName2Col(src.al3.integrated_re_al3$cluster.extract.v1.1[cellorder_src.al3.integrated_re_al3], 
                 cluster.endoderm.color.v5),
      MyName2Col(src.al3.integrated_re_al3$cluster.v06.26.re_correct_re[cellorder_src.al3.integrated_re_al3], 
                 cluster.endoderm.color.v5)),
    RowSideColors = t(cbind(
      MyName2Col(names(gene_src.al3.integrated_re_al3.rowtree.2), colors.geneset))),
    ColSideColorsSize = 3,
    RowSideColorsSize = 1.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    # return.tree = "row",
    Colv="none", Rowv = "none",
    graph = T)
dev.off()


# -- GCN for cluster
#----------------------
pearson.gene = WGCNA::cor(
  t(src.al3.integrated_re_al3@assays$RNA@data[
    gene_src.al3.integrated_re_al3.rowtree.2,]), method = "p")
pearson.gene = logistic(pearson.gene, threshold = 0.25)
pearson.gene[is.na(pearson.gene)]=0
pearson.gene[pearson.gene<0.53]=0

graph.gene = graph.adjacency(pearson.gene,mode = "undirected",weighted = T)
graph.gene = igraph::simplify(graph.gene,remove.multiple = T,remove.loops = T)
graph.gene = delete.vertices(
  graph.gene,
  v = V(graph.gene)[clusters(graph.gene)$membership%in%which(clusters(graph.gene)$csize<=10)])
degree.Foregut.gene=igraph::degree(graph.gene)
layout.gene = layout_nicely(graph.gene)

graph.cluster <- cluster_walktrap(graph.gene,steps = 5)
gene.cluster <- graph.cluster$membership
names(gene.cluster) <- graph.cluster$names
set.seed(3)

rev_gene_src.al3.integrated_re_al3.rowtree.2 = names(gene_src.al3.integrated_re_al3.rowtree.2)
names(rev_gene_src.al3.integrated_re_al3.rowtree.2) = gene_src.al3.integrated_re_al3.rowtree.2

pdf("figure.v08.07/organ_development_re_v240115/try.src.al3.integrated_re_al3_gcn_tf.pdf",5,5)
plot.igraph(graph.gene,
            layout = layout.gene,
            vertex.size=7,
            label.cex =8,
            vertex.label = ifelse(names(gene.cluster)%in%gi[gi$TF%in%T,]$SymbolDedu,
                                  names(gene.cluster), NA),
            vertex.label.font = 3,
            vertex.label.cex = 1.2,
            vertex.label.color = "black",
            vertex.color = colors.geneset[
              rev_gene_src.al3.integrated_re_al3.rowtree.2[
                names(gene.cluster)]])
dev.off()
#------------------------


save(src.al3.integrated_re_al3, file = 'figure.v08.07/organ_development_re_v240115/src.al3.integrated_re_al3.Rdata')
save(src.al3.integrated_re_al3_filtergene,
     load_src.al3.integrated_re_al3_pca1,
     load_src.al3.integrated_re_al3_pca2,
     load_src.al3.integrated_re_al3_pca3,
     src.al3.integrated_re_al3_marker,
     src.al3.integrated_re_al3_markergene,
     tree_src.al3.integrated_re_al3.rowtree,
     tree_src.al3.integrated_re_al3.rowtree.1,
     tree_src.al3.integrated_re_al3.rowtree.2,
     tree_src.al3.integrated_re_al3.coltree.1,
     tree_src.al3.integrated_re_al3.coltree.2,
     
     gene_src.al3.integrated_re_al3.rowtree,
     gene_src.al3.integrated_re_al3.rowtree.1,
     gene_src.al3.integrated_re_al3.rowtree.2,
     
     cellorder_src.al3.integrated_re_al3,
     file = 'figure.v08.07/organ_development_re_v240115/src.al3.integrated_re_al3.parameter.Rdata')









