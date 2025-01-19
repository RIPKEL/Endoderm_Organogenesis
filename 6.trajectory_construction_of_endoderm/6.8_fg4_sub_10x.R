#---------------------------------------------------------------------------
# --  FG.4_sub_intermediate 
#---------------------------------------------------------------------------
src.fg4.integrated_re = src.endoderm.fg4.ext_v1.1
DimPlot(src.endoderm.fg4.ext_v1.1, reduction = "umap_mnn")

src.fg4.integrated_re$cluster.v06.26.re_correct = src.fg4.integrated_re$cluster.v06.26.re

pdf("figure.v08.07/organ_development_re_v240115/try.src.fg4.integrated_re_umap.pdf",9,9)
DimPlot(src.fg4.integrated_re, group.by = "cluster.v06.26.re_correct", 
        reduction = "umap_mnn", pt.size = 1.2,
        cols = cluster.endoderm.color.v5) + 
  theme_bw() + p_add
DimPlot(src.fg4.integrated_re, group.by = "cluster.extract.v1.1", 
        reduction = "umap_mnn", pt.size = 1.2,
        cols = cluster.endoderm.color.v5)+
  theme_bw() + p_add
DimPlot(src.fg4.integrated_re, group.by = "Time", 
        reduction = "umap_mnn", pt.size = 1.2,
        cols = c(cluster.endoderm.color.v5, colors.time))+
  theme_bw() + p_add
dev.off()

p.3d=plot_ly(x=src.fg4.integrated_re[['umap_mnn']]@cell.embeddings[,1], 
             y=src.fg4.integrated_re[['umap_mnn']]@cell.embeddings[,2],
             z=src.fg4.integrated_re[['umap_mnn']]@cell.embeddings[,3],
             type = "scatter3d", mode="markers",
             color = src.fg4.integrated_re$cluster.v06.26.re_correct,
             colors= cluster.endoderm.color.v5[unique(src.fg4.integrated_re$cluster.v06.26.re_correct)],
             #color = src.fg4.integrated_re$Time,
             #colors= colors.time[unique(src.fg4.integrated_re$Time)],
             xlab = "", ylab = "", zlab = "",
             box=T, axes = T, size=5, alpha = 1.2)%>%
  layout(scene = list(xaxis = list(title = 'Coord_1', autorange = T, showgrid = T, zeroline = F, showline = F, autotick = F, ticks = '', showticklabels = F),
                      yaxis = list(title = 'Coord_2', autorange = T, showgrid = T, zeroline = F, showline = F, autotick = F, ticks = '', showticklabels = F),
                      zaxis = list(title = 'Coord_3', autorange = T, showgrid = T, zeroline = F, showline = F, autotick = F, ticks = '', showticklabels = F),
                      aspectmode = "manual", 
                      aspectratio = list(x=1, y=1, z=1),
                      camera = list(eye = list(x = 1.4, y = -1, z = 0.3))))
p.3d
htmlwidgets::saveWidget(as_widget(p.3d), "figure.v08.07/organ_development_re_v240115/html/celltype_src.fg4.integrated_re.html")
htmlwidgets::saveWidget(as_widget(p.3d), "figure.v08.07/organ_development_re_v240115/html/time_src.fg4.integrated_re.html")




src.fg4.integrated_re = NormalizeData(src.fg4.integrated_re, scale.factor = 10^5)
src.fg4.integrated_re = ScaleData(src.fg4.integrated_re, 
                                  features = rownames(src.fg4.integrated_re),
                                  split.by = "batch_phase")
src.fg4.integrated_re = FindVariableFeatures(src.fg4.integrated_re, assay = "RNA", nfeatures = 2000)
src.fg4.integrated_re_filtergene = 
  Myfilter(as.matrix(src.fg4.integrated_re@assays$RNA@data),
           gene = src.fg4.integrated_re@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.fg4.integrated_re_filtergene_re = src.fg4.integrated_re_filtergene
rm(src.fg4.integrated_re_filtergene)
src.fg4.integrated_re_filtergene = src.fg4.integrated_re_filtergene_re

src.fg4.integrated_re = RunPCA(src.fg4.integrated_re, 
                               features = src.endoderm.fg4.ext_v1.1.gene.fin)
src.fg4.integrated_re = RunUMAP(src.fg4.integrated_re, reduction = "pca",
                                n.neighbors = 100, dims = c(1:30), n.components = 3)

#---------------------------------
# -- FG.4 Intermeditae Stage
#=================================================================================
src.fg4.integrated_re_fg4 = 
  src.fg4.integrated_re[, src.fg4.integrated_re$cluster.v06.26.re_correct%in%c("FG.4","FG.4-Liver","FG.4-Lung/Stomach")]
DimPlot(src.fg4.integrated_re_fg4, group.by = "Time")

src.fg4.integrated_re_fg4 = NormalizeData(src.fg4.integrated_re_fg4, scale.factor = 10^5)
src.fg4.integrated_re_fg4 = ScaleData(src.fg4.integrated_re_fg4, 
                                      features = rownames(src.fg4.integrated_re_fg4),
                                      split.by = "batch_phase")
src.fg4.integrated_re_fg4 = FindVariableFeatures(src.fg4.integrated_re_fg4, assay = "RNA", nfeatures = 2000)
src.fg4.integrated_re_fg4_filtergene = 
  Myfilter(as.matrix(src.fg4.integrated_re_fg4@assays$RNA@data),
           gene = src.fg4.integrated_re_fg4@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.fg4.integrated_re_fg4_filtergene_re = src.fg4.integrated_re_fg4_filtergene
rm(src.fg4.integrated_re_fg4_filtergene)
src.fg4.integrated_re_fg4_filtergene = src.fg4.integrated_re_fg4_filtergene_re

# -- reduction
src.fg4.integrated_re_fg4 = RunPCA(src.fg4.integrated_re_fg4, 
                                   features = src.fg4.integrated_re_fg4_filtergene)
src.fg4.integrated_re_fg4 = RunUMAP(src.fg4.integrated_re_fg4, reduction = "pca",
                                    dims = 1:25, n.components = 2)

pdf("figure.v08.07/organ_development_re_v240115/try.fg4_fg4_pca.pdf",10,10)
src.fg4.integrated_re_fg4 = RunPCA(src.fg4.integrated_re_fg4, 
                                   features = src.fg4.integrated_re_fg4_filtergene)
DimPlot(src.fg4.integrated_re_fg4, group.by = "Time", 
        cols = colors.time, pt.size = 1.25,
        reduction = "pca", dims = c(1,3)) +
  theme_bw() + p_add
DimPlot(src.fg4.integrated_re_fg4, group.by = "cluster.v06.26.re_correct_re", 
        cols = cluster.endoderm.color.v5, pt.size = 1.25,
        reduction = "pca", dims = c(1,3)) +
  theme_bw() + p_add


# src.fg4.integrated_re_fg4 = RunPCA(src.fg4.integrated_re_fg4, 
#                                    features = gene_src.fg4.integrated_re_fg4.rowtree.3)
DimPlot(src.fg4.integrated_re_fg4, group.by = "Time", 
        cols = colors.time, pt.size = 1.25,
        reduction = "pca", dims = c(1,2)) +
  theme_bw() + p_add
DimPlot(src.fg4.integrated_re_fg4, group.by = "cluster.v06.26.re_correct_re", 
        cols = cluster.endoderm.color.v5, pt.size = 1.25,
        reduction = "pca", dims = c(1,2)) +
  theme_bw() + p_add
DimPlot(src.fg4.integrated_re_fg4, group.by = "cluster.v06.26.re_correct_re",
        cols = cluster.endoderm.color.v5, pt.size = 1.25,
        reduction = "umap", dims = c(1,2)) +
  theme_bw() + p_add
dev.off()


load_src.fg4.integrated_re_fg4_pca3 = abs(src.fg4.integrated_re_fg4@reductions$pca@feature.loadings[,3])
gene_src.fg4.integrated_re_fg4_pca3 = 
  load_src.fg4.integrated_re_fg4_pca3[order(load_src.fg4.integrated_re_fg4_pca3,decreasing=T)][1:50]

load_src.fg4.integrated_re_fg4_pca1 = abs(src.fg4.integrated_re_fg4@reductions$pca@feature.loadings[,1])
gene_src.fg4.integrated_re_fg4_pca1 = 
  load_src.fg4.integrated_re_fg4_pca1[order(load_src.fg4.integrated_re_fg4_pca1,decreasing=T)][1:150]

# -- marker
src.fg4.integrated_re_fg4 = SetIdent(src.fg4.integrated_re_fg4, 
                                     value = src.fg4.integrated_re_fg4$cluster.v06.26.re_correct)
src.fg4.integrated_re_fg4_marker = FindAllMarkers(src.fg4.integrated_re_fg4)
src.fg4.integrated_re_fg4_marker$pct.ratio =
  src.fg4.integrated_re_fg4_marker$pct.1 / src.fg4.integrated_re_fg4_marker$pct.2
src.fg4.integrated_re_fg4_markergene = unique(
  src.fg4.integrated_re_fg4_marker[
    src.fg4.integrated_re_fg4_marker$pct.ratio>2,]$gene)



pdf("figure.v08.07/organ_development_re_v240115/try.fg4_fg4.pdf",9,7)
src.fg4.integrated_re_fg4.rowtree =
  MyHeatmap(as.matrix(src.fg4.integrated_re_fg4@assays$RNA@data[
    unique(c(
      src.fg4.integrated_re_fg4_markergene,
      src.fg4.integrated_re_fg4_filtergene)),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg4.integrated_re_fg4$Time, colors.time),
      MyName2Col(src.fg4.integrated_re_fg4$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src.fg4.integrated_re_fg4$cluster.v06.26.re, cluster.endoderm.color.v5),
      MyName2Col(src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree1, cluster.endoderm.color.v5)),
    ColSideColorsSize = 5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()
pdf("figure.v08.07/organ_development_re_v240115/try.fg4_fg4_tree.pdf",100,20)
plot(src.fg4.integrated_re_fg4.rowtree)
dev.off()

tree_src.fg4.integrated_re_fg4.rowtree = as.dendrogram(src.fg4.integrated_re_fg4.rowtree)
gene_src.fg4.integrated_re_fg4.rowtree = setdiff(
  labels(tree_src.fg4.integrated_re_fg4.rowtree),
  c(labels(tree_src.fg4.integrated_re_fg4.rowtree[[1]][[1]]),
    labels(tree_src.fg4.integrated_re_fg4.rowtree[[2]][[1]]),
    labels(tree_src.fg4.integrated_re_fg4.rowtree[[2]][[2]][[1]][[1]]),
    labels(tree_src.fg4.integrated_re_fg4.rowtree[[2]][[2]][[1]][[2]][[1]])))



pdf("figure.v08.07/organ_development_re_v240115/try.fg4_fg4.pdf",9,7)
# -- Cell type
src.fg4.integrated_re_fg4.coltree.1 =
  MyHeatmap(as.matrix(src.fg4.integrated_re_fg4@assays$RNA@data[
    names(gene_src.fg4.integrated_re_fg4_pca3),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg4.integrated_re_fg4$Time, colors.time),
      MyName2Col(src.fg4.integrated_re_fg4$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src.fg4.integrated_re_fg4$cluster.v06.26.re, cluster.endoderm.color.v5),
      MyName2Col(src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree3, cluster.endoderm.color.v5)),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "col",
    graph = T)
src.fg4.integrated_re_fg4.rowtree.1 =
  MyHeatmap(as.matrix(src.fg4.integrated_re_fg4@assays$RNA@data[
    names(gene_src.fg4.integrated_re_fg4_pca3),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg4.integrated_re_fg4$Time, colors.time),
      MyName2Col(src.fg4.integrated_re_fg4$tree_1, color.temp),
      MyName2Col(src.fg4.integrated_re_fg4$tree_2, color.temp),
      MyName2Col(src.fg4.integrated_re_fg4$tree_3, color.temp),
      MyName2Col(src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree3, cluster.endoderm.color.v5)),
    ColSideColorsSize = 5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)

# -- Time order
src.fg4.integrated_re_fg4.coltree.2 =
  MyHeatmap(as.matrix(src.fg4.integrated_re_fg4@assays$RNA@data[
    names(gene_src.fg4.integrated_re_fg4_pca1),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg4.integrated_re_fg4$Time, colors.time),
      MyName2Col(src.fg4.integrated_re_fg4$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src.fg4.integrated_re_fg4$cluster.v06.26.re, cluster.endoderm.color.v5),
      MyName2Col(src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree3, cluster.endoderm.color.v5)),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "col",
    graph = T)
src.fg4.integrated_re_fg4.rowtree.2 =
  MyHeatmap(as.matrix(src.fg4.integrated_re_fg4@assays$RNA@data[
    names(gene_src.fg4.integrated_re_fg4_pca1),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg4.integrated_re_fg4$Time, colors.time),
      MyName2Col(src.fg4.integrated_re_fg4$tree_2, color.temp),
      MyName2Col(src.fg4.integrated_re_fg4$tree_3, color.temp),
      MyName2Col(src.fg4.integrated_re_fg4$cluster.extract.v1.1, cluster.endoderm.color.v5),
      #MyName2Col(src.fg4.integrated_re_fg4$cluster.v06.26.re, cluster.endoderm.color.v5),
      #MyName2Col(src.fg4.integrated_re_fg4$cluster.v06.26.re_correct, cluster.endoderm.color.v5),
      MyName2Col(src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree3, cluster.endoderm.color.v5)),
    ColSideColorsSize = 1.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)

dev.off()

#=== R2
tree_src.fg4.integrated_re_fg4.rowtree.2 = as.dendrogram(src.fg4.integrated_re_fg4.rowtree.2)
gene_src.fg4.integrated_re_fg4.rowtree.2 = c(labels(tree_src.fg4.integrated_re_fg4.rowtree.2[[1]]))

#=== C2
tree_src.fg4.integrated_re_fg4.coltree.2 = as.dendrogram(src.fg4.integrated_re_fg4.coltree.2)
src.fg4.integrated_re_fg4$tree_2 = NA
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.2[[2]][[1]]),]$tree_2 = 3
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.2[[2]][[2]][[1]]),]$tree_2 = 1
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.2[[2]][[2]][[2]]),]$tree_2 = 2

src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree2 = src.fg4.integrated_re_fg4$cluster.v06.26.re_correct
src.fg4.integrated_re_fg4@meta.data[src.fg4.integrated_re_fg4$tree_2%in%c(1,2,3),]$cluster.v06.26.re_correct_re.tree2 = "FG.4"

#=== R1
tree_src.fg4.integrated_re_fg4.rowtree.1 = as.dendrogram(src.fg4.integrated_re_fg4.rowtree.1)
gene_src.fg4.integrated_re_fg4.rowtree.1 = c(
  labels(tree_src.fg4.integrated_re_fg4.rowtree.1[[2]][[2]]),
  labels(tree_src.fg4.integrated_re_fg4.rowtree.1[[1]]))

#=== C1
tree_src.fg4.integrated_re_fg4.coltree.1 = as.dendrogram(src.fg4.integrated_re_fg4.coltree.1)
src.fg4.integrated_re_fg4$tree_1 = NA
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.1[[1]][[1]]),]$tree_1 = 1
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.1[[1]][[2]]),]$tree_1 = 2
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.1[[2]][[1]]),]$tree_1 = 3
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.1[[2]][[2]][[1]]),]$tree_1 = 4
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.1[[2]][[2]][[2]]),]$tree_1 = 5

src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree1 = src.fg4.integrated_re_fg4$cluster.v06.26.re_correct
src.fg4.integrated_re_fg4@meta.data[src.fg4.integrated_re_fg4$tree_1%in%c(4),]$cluster.v06.26.re_correct_re.tree1 = "FG.4"
src.fg4.integrated_re_fg4@meta.data[src.fg4.integrated_re_fg4$tree_1%in%c(1,2),]$cluster.v06.26.re_correct_re.tree1 = "FG.4-Liver"
src.fg4.integrated_re_fg4@meta.data[src.fg4.integrated_re_fg4$tree_1%in%c(3,5),]$cluster.v06.26.re_correct_re.tree1 = "FG.4-Lung/Stomach"



# Pseudo-time :: cluster.v06.26.re_correct_re
#---------------------------------------------------------
cell_src.fg4.integrated_re_fg4_lu = 
  rownames(src.fg4.integrated_re_fg4@meta.data[src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re%in%c("FG.4",'FG.4-Liver'),])
coord_src.fg4.integrated_re_fg4_lu = src.fg4.integrated_re_fg4[["pca"]]@cell.embeddings[cell_src.fg4.integrated_re_fg4_lu, c(1,2)]
pcurve_src.fg4.integrated_re_fg4_lu = princurve::principal_curve(x = coord_src.fg4.integrated_re_fg4_lu, smoother = "smooth.spline")
src.fg4.integrated_re_fg4$lambda_lu = pcurve_src.fg4.integrated_re_fg4_lu$lambda[cell_src.fg4.integrated_re_fg4_lu]
src.fg4.integrated_re_fg4@meta.data[cell_src.fg4.integrated_re_fg4_lu,]$lambda_lu = 
  pcurve_src.fg4.integrated_re_fg4_lu$lambda[cell_src.fg4.integrated_re_fg4_lu]
src.fg4.integrated_re_fg4@meta.data[!src.fg4.integrated_re_fg4$lambda_lu%in%NA,]$lambda_lu = 
  norm_range(src.fg4.integrated_re_fg4@meta.data[!src.fg4.integrated_re_fg4$lambda_lu%in%NA,]$lambda_lu)

cell_src.fg4.integrated_re_fg4_lv = 
  rownames(src.fg4.integrated_re_fg4@meta.data[src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re%in%c("FG.4",'FG.4-Lung/Stomach'),])
coord_src.fg4.integrated_re_fg4_lv = src.fg4.integrated_re_fg4[["pca"]]@cell.embeddings[cell_src.fg4.integrated_re_fg4_lv, c(1,2)]
pcurve_src.fg4.integrated_re_fg4_lv = princurve::principal_curve(x = coord_src.fg4.integrated_re_fg4_lv, smoother = "smooth.spline")
src.fg4.integrated_re_fg4$lambda_lv = pcurve_src.fg4.integrated_re_fg4_lv$lambda[cell_src.fg4.integrated_re_fg4_lv]
src.fg4.integrated_re_fg4@meta.data[cell_src.fg4.integrated_re_fg4_lv,]$lambda_lv = 
  pcurve_src.fg4.integrated_re_fg4_lv$lambda[cell_src.fg4.integrated_re_fg4_lv]
src.fg4.integrated_re_fg4@meta.data[!src.fg4.integrated_re_fg4$lambda_lv%in%NA,]$lambda_lv = 
  norm_range(src.fg4.integrated_re_fg4@meta.data[!src.fg4.integrated_re_fg4$lambda_lv%in%NA,]$lambda_lv)


cellorder_src.fg4.integrated_re_fg4 = c(
  (intersect(names(src.fg4.integrated_re_fg4$lambda_lu[order(src.fg4.integrated_re_fg4$lambda_lu)]), 
             colnames(src.fg4.integrated_re_fg4[,src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re%in%c("FG.4"),]))),
  
  (intersect(names(src.fg4.integrated_re_fg4$lambda_lv[order(src.fg4.integrated_re_fg4$lambda_lv)]), 
             colnames(src.fg4.integrated_re_fg4[,src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re%in%c('FG.4-Lung/Stomach'),]))),
  
  (intersect(names(src.fg4.integrated_re_fg4$lambda_lu[order(src.fg4.integrated_re_fg4$lambda_lu)]), 
             colnames(src.fg4.integrated_re_fg4[,src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re%in%"FG.4-Liver"]))))
#------------------------------------------------


pdf("figure.v08.07/organ_development_re_v240115/try.fg4_fg4.pdf",9,7)
src.fg4.integrated_re_fg4.rowtree.3 =
  MyHeatmap(as.matrix(src.fg4.integrated_re_fg4@assays$RNA@data[
    unique(c(
      gene_src.fg4.integrated_re_fg4.rowtree.1,
      gene_src.fg4.integrated_re_fg4.rowtree.2)),
    cellorder_src.fg4.integrated_re_fg4]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg4.integrated_re_fg4$Time[cellorder_src.fg4.integrated_re_fg4], 
                 colors.time),
      MyName2Col(src.fg4.integrated_re_fg4$cluster.v06.26.re[cellorder_src.fg4.integrated_re_fg4], 
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg4.integrated_re_fg4$tree_3[cellorder_src.fg4.integrated_re_fg4], 
                 color.temp),
      MyName2Col(src.fg4.integrated_re_fg4$cluster.extract.v1.1[cellorder_src.fg4.integrated_re_fg4], 
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree3[cellorder_src.fg4.integrated_re_fg4], 
                 cluster.endoderm.color.v5)),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none", Colv = "none",
    return.tree = "row",
    graph = T)

src.fg4.integrated_re_fg4.coltree.3 =
  MyHeatmap(as.matrix(src.fg4.integrated_re_fg4@assays$RNA@data[
    unique(c(
      gene_src.fg4.integrated_re_fg4.rowtree.1,
      gene_src.fg4.integrated_re_fg4.rowtree.2)),
    cellorder_src.fg4.integrated_re_fg4]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg4.integrated_re_fg4$Time[cellorder_src.fg4.integrated_re_fg4], 
                 colors.time),
      MyName2Col(src.fg4.integrated_re_fg4$cluster.extract.v1.1[cellorder_src.fg4.integrated_re_fg4], 
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg4.integrated_re_fg4$tree_3[cellorder_src.fg4.integrated_re_fg4],
                 color.temp),
      MyName2Col(src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re[cellorder_src.fg4.integrated_re_fg4],  
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree3[cellorder_src.fg4.integrated_re_fg4], 
                 cluster.endoderm.color.v5)),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none", Colv = "none",
    return.tree = "col",
    graph = T)
dev.off()

#== R3
tree_src.fg4.integrated_re_fg4.rowtree.3 = as.dendrogram(src.fg4.integrated_re_fg4.rowtree.3)
gene_src.fg4.integrated_re_fg4.rowtree.3 = c(
  labels(tree_src.fg4.integrated_re_fg4.rowtree.3[[1]][[1]]),
  rev(labels(tree_src.fg4.integrated_re_fg4.rowtree.3[[1]][[2]])),
  
  labels(tree_src.fg4.integrated_re_fg4.rowtree.3[[2]][[1]]),
  rev(labels(tree_src.fg4.integrated_re_fg4.rowtree.3[[2]][[2]])))
names(gene_src.fg4.integrated_re_fg4.rowtree.3) = c(
  rep(4, length(labels(tree_src.fg4.integrated_re_fg4.rowtree.3[[1]]))),
  rep(3, length(labels(tree_src.fg4.integrated_re_fg4.rowtree.3[[2]][[1]]))),
  rep(7, length(labels(tree_src.fg4.integrated_re_fg4.rowtree.3[[2]][[2]]))))


#=== C3
tree_src.fg4.integrated_re_fg4.coltree.3 = as.dendrogram(src.fg4.integrated_re_fg4.coltree.3)
src.fg4.integrated_re_fg4$tree_3 = NA
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.3[[1]][[1]]),]$tree_3 = 1
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.3[[1]][[2]]),]$tree_3 = 2
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.3[[2]][[1]]),]$tree_3 = 3
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.3[[2]][[2]][[1]][[1]]),]$tree_3 = 4
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.3[[2]][[2]][[1]][[2]]),]$tree_3 = 5
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.3[[2]][[2]][[2]][[1]]),]$tree_3 = 6
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.3[[2]][[2]][[2]][[2]][[1]]),]$tree_3 = 7
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.3[[2]][[2]][[2]][[2]][[2]][[1]]),]$tree_3 = 8
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.3[[2]][[2]][[2]][[2]][[2]][[2]]),]$tree_3 = 9


src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree3 = 
  src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re
src.fg4.integrated_re_fg4@meta.data[
  src.fg4.integrated_re_fg4$tree_3%in%c(1,2,4,5,8),]$cluster.v06.26.re_correct_re.tree3 = "FG.4"
src.fg4.integrated_re_fg4@meta.data[
  src.fg4.integrated_re_fg4$tree_3%in%c(3,7),]$cluster.v06.26.re_correct_re.tree3 = "FG.4-Lung/Stomach"
src.fg4.integrated_re_fg4@meta.data[
  src.fg4.integrated_re_fg4$tree_3%in%c(6,9),]$cluster.v06.26.re_correct_re.tree3 = "FG.4-Liver"

DimPlot(src.fg4.integrated_re_fg4, group.by = "cluster.v06.26.re_correct_re.tree3")


# -- Refined
# src.fg4.integrated_re_fg4@meta.data[src.fg4.integrated_re_fg4$tree_2%in%c(1,2,3),]$cluster.v06.26.re_correct_re.tree3 = "FG.4"

# Pseudo-time :: cluster.v06.26.re_correct_re.tree3
#---------------------------------------------------------
cell_src.fg4.integrated_re_fg4_lu = 
  rownames(src.fg4.integrated_re_fg4@meta.data[src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree3%in%c("FG.4",'FG.4-Liver'),])
coord_src.fg4.integrated_re_fg4_lu = src.fg4.integrated_re_fg4[["pca"]]@cell.embeddings[cell_src.fg4.integrated_re_fg4_lu, c(1,2)]
pcurve_src.fg4.integrated_re_fg4_lu = princurve::principal_curve(x = coord_src.fg4.integrated_re_fg4_lu, smoother = "smooth.spline")
src.fg4.integrated_re_fg4$lambda_lu = pcurve_src.fg4.integrated_re_fg4_lu$lambda[cell_src.fg4.integrated_re_fg4_lu]
src.fg4.integrated_re_fg4@meta.data[cell_src.fg4.integrated_re_fg4_lu,]$lambda_lu = 
  pcurve_src.fg4.integrated_re_fg4_lu$lambda[cell_src.fg4.integrated_re_fg4_lu]
src.fg4.integrated_re_fg4@meta.data[!src.fg4.integrated_re_fg4$lambda_lu%in%NA,]$lambda_lu = 
  norm_range(src.fg4.integrated_re_fg4@meta.data[!src.fg4.integrated_re_fg4$lambda_lu%in%NA,]$lambda_lu)

cell_src.fg4.integrated_re_fg4_lv = 
  rownames(src.fg4.integrated_re_fg4@meta.data[src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree3%in%c("FG.4",'FG.4-Lung/Stomach'),])
coord_src.fg4.integrated_re_fg4_lv = src.fg4.integrated_re_fg4[["pca"]]@cell.embeddings[cell_src.fg4.integrated_re_fg4_lv, c(1,2)]
pcurve_src.fg4.integrated_re_fg4_lv = princurve::principal_curve(x = coord_src.fg4.integrated_re_fg4_lv, smoother = "smooth.spline")
src.fg4.integrated_re_fg4$lambda_lv = pcurve_src.fg4.integrated_re_fg4_lv$lambda[cell_src.fg4.integrated_re_fg4_lv]
src.fg4.integrated_re_fg4@meta.data[cell_src.fg4.integrated_re_fg4_lv,]$lambda_lv = 
  pcurve_src.fg4.integrated_re_fg4_lv$lambda[cell_src.fg4.integrated_re_fg4_lv]
src.fg4.integrated_re_fg4@meta.data[!src.fg4.integrated_re_fg4$lambda_lv%in%NA,]$lambda_lv = 
  norm_range(src.fg4.integrated_re_fg4@meta.data[!src.fg4.integrated_re_fg4$lambda_lv%in%NA,]$lambda_lv)


cellorder_src.fg4.integrated_re_fg4 = c(
  (intersect(names(src.fg4.integrated_re_fg4$lambda_lv[order(src.fg4.integrated_re_fg4$lambda_lv)]), 
             colnames(src.fg4.integrated_re_fg4[,src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree3%in%c("FG.4"),]))),
  
  (intersect(names(src.fg4.integrated_re_fg4$lambda_lv[order(src.fg4.integrated_re_fg4$lambda_lv)]), 
             colnames(src.fg4.integrated_re_fg4[,src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree3%in%c('FG.4-Lung/Stomach'),]))),
  
  (intersect(names(src.fg4.integrated_re_fg4$lambda_lu[order(src.fg4.integrated_re_fg4$lambda_lu)]), 
             colnames(src.fg4.integrated_re_fg4[,src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree3%in%"FG.4-Liver"]))))
#------------------------------------------------

pdf("figure.v08.07/organ_development_re_v240115/try.fg4_fg4_final.pdf",9,7)
src.fg4.integrated_re_fg4.coltree.4 =
  MyHeatmap(as.matrix(src.fg4.integrated_re_fg4@assays$RNA@data[
    gene_src.fg4.integrated_re_fg4.rowtree.3,
    cellorder_src.fg4.integrated_re_fg4]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg4.integrated_re_fg4$Time[cellorder_src.fg4.integrated_re_fg4], 
                 colors.time),
      MyName2Col(src.fg4.integrated_re_fg4$tree_3[cellorder_src.fg4.integrated_re_fg4], 
                 color.temp),
      MyName2Col(src.fg4.integrated_re_fg4$tree_4[cellorder_src.fg4.integrated_re_fg4], 
                 color.temp),
      MyName2Col(src.fg4.integrated_re_fg4$cluster.extract.v1.1[cellorder_src.fg4.integrated_re_fg4], 
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree3[cellorder_src.fg4.integrated_re_fg4], 
                 cluster.endoderm.color.v5)),
    RowSideColors = t(cbind(
      MyName2Col(names(gene_src.fg4.integrated_re_fg4.rowtree.3), colors.geneset))),
    ColSideColorsSize = 3,
    RowSideColorsSize = 1.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none", Colv = "none",
    return.tree = "col",
    graph = T)
dev.off()

#=== C4
tree_src.fg4.integrated_re_fg4.coltree.4 = as.dendrogram(src.fg4.integrated_re_fg4.coltree.4)
src.fg4.integrated_re_fg4$tree_4 = NA
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.4[[1]][[1]]),]$tree_4 = 1
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.4[[1]][[2]]),]$tree_4 = 2
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.4[[2]][[1]]),]$tree_4 = 3
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.4[[2]][[2]][[1]][[1]]),]$tree_4 = 4
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.4[[2]][[2]][[1]][[2]][[1]]),]$tree_4 = 5
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.4[[2]][[2]][[1]][[2]][[2]]),]$tree_4 = 10
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.4[[2]][[2]][[2]][[1]][[1]]),]$tree_4 = 6
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.4[[2]][[2]][[2]][[1]][[2]]),]$tree_4 = 11
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.4[[2]][[2]][[2]][[2]][[1]]),]$tree_4 = 7
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.4[[2]][[2]][[2]][[2]][[2]][[1]]),]$tree_4 = 8
src.fg4.integrated_re_fg4@meta.data[labels(tree_src.fg4.integrated_re_fg4.coltree.4[[2]][[2]][[2]][[2]][[2]][[2]]),]$tree_4 = 9



src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree4 = src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re
src.fg4.integrated_re_fg4@meta.data[src.fg4.integrated_re_fg4$tree_4%in%c(1,2,7,9,11,6,10),]$cluster.v06.26.re_correct_re.tree4 = "FG.4"
src.fg4.integrated_re_fg4@meta.data[src.fg4.integrated_re_fg4$tree_4%in%c(3,8),]$cluster.v06.26.re_correct_re.tree4 = "FG.4-Lung/Stomach"
src.fg4.integrated_re_fg4@meta.data[src.fg4.integrated_re_fg4$tree_4%in%c(4,5),]$cluster.v06.26.re_correct_re.tree4 = "FG.4-Liver"


# Pseudo-time :: cluster.v06.26.re_correct_re.tree4
#---------------------------------------------------------
cell_src.fg4.integrated_re_fg4_lu = 
  rownames(src.fg4.integrated_re_fg4@meta.data[src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree4%in%c("FG.4",'FG.4-Liver'),])
coord_src.fg4.integrated_re_fg4_lu = src.fg4.integrated_re_fg4[["pca"]]@cell.embeddings[cell_src.fg4.integrated_re_fg4_lu, c(1,2)]
pcurve_src.fg4.integrated_re_fg4_lu = princurve::principal_curve(x = coord_src.fg4.integrated_re_fg4_lu, smoother = "smooth.spline")
src.fg4.integrated_re_fg4$lambda_lu = pcurve_src.fg4.integrated_re_fg4_lu$lambda[cell_src.fg4.integrated_re_fg4_lu]
src.fg4.integrated_re_fg4@meta.data[cell_src.fg4.integrated_re_fg4_lu,]$lambda_lu = 
  pcurve_src.fg4.integrated_re_fg4_lu$lambda[cell_src.fg4.integrated_re_fg4_lu]
src.fg4.integrated_re_fg4@meta.data[!src.fg4.integrated_re_fg4$lambda_lu%in%NA,]$lambda_lu = 
  norm_range(src.fg4.integrated_re_fg4@meta.data[!src.fg4.integrated_re_fg4$lambda_lu%in%NA,]$lambda_lu)

cell_src.fg4.integrated_re_fg4_lv = 
  rownames(src.fg4.integrated_re_fg4@meta.data[src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree4%in%c("FG.4",'FG.4-Lung/Stomach'),])
coord_src.fg4.integrated_re_fg4_lv = src.fg4.integrated_re_fg4[["pca"]]@cell.embeddings[cell_src.fg4.integrated_re_fg4_lv, c(1,2)]
pcurve_src.fg4.integrated_re_fg4_lv = princurve::principal_curve(x = coord_src.fg4.integrated_re_fg4_lv, smoother = "smooth.spline")
src.fg4.integrated_re_fg4$lambda_lv = pcurve_src.fg4.integrated_re_fg4_lv$lambda[cell_src.fg4.integrated_re_fg4_lv]
src.fg4.integrated_re_fg4@meta.data[cell_src.fg4.integrated_re_fg4_lv,]$lambda_lv = 
  pcurve_src.fg4.integrated_re_fg4_lv$lambda[cell_src.fg4.integrated_re_fg4_lv]
src.fg4.integrated_re_fg4@meta.data[!src.fg4.integrated_re_fg4$lambda_lv%in%NA,]$lambda_lv = 
  norm_range(src.fg4.integrated_re_fg4@meta.data[!src.fg4.integrated_re_fg4$lambda_lv%in%NA,]$lambda_lv)


cellorder_src.fg4.integrated_re_fg4 = c(
  (intersect(names(src.fg4.integrated_re_fg4$lambda_lv[order(src.fg4.integrated_re_fg4$lambda_lv)]), 
             colnames(src.fg4.integrated_re_fg4[,src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree4%in%c("FG.4"),]))),
  
  (intersect(names(src.fg4.integrated_re_fg4$lambda_lv[order(src.fg4.integrated_re_fg4$lambda_lv)]), 
             colnames(src.fg4.integrated_re_fg4[,src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree4%in%c('FG.4-Lung/Stomach'),]))),
  
  (intersect(names(src.fg4.integrated_re_fg4$lambda_lu[order(src.fg4.integrated_re_fg4$lambda_lu)]), 
             colnames(src.fg4.integrated_re_fg4[,src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree4%in%"FG.4-Liver"]))))
#------------------------------------------------

pdf("figure.v08.07/organ_development_re_v240115/try.fg4_fg4_final.pdf",9,7)
src.fg4.integrated_re_fg4.coltree.5 =
  MyHeatmap(as.matrix(src.fg4.integrated_re_fg4@assays$RNA@data[
    gene_src.fg4.integrated_re_fg4.rowtree.3,
    cellorder_src.fg4.integrated_re_fg4]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg4.integrated_re_fg4$Time[cellorder_src.fg4.integrated_re_fg4], 
                 colors.time),
      MyName2Col(src.fg4.integrated_re_fg4$tree_3[cellorder_src.fg4.integrated_re_fg4], 
                 color.temp),
      MyName2Col(src.fg4.integrated_re_fg4$tree_4[cellorder_src.fg4.integrated_re_fg4], 
                 color.temp),
      MyName2Col(src.fg4.integrated_re_fg4$cluster.extract.v1.1[cellorder_src.fg4.integrated_re_fg4], 
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg4.integrated_re_fg4$cluster.v06.26.re_correct_re.tree4[cellorder_src.fg4.integrated_re_fg4], 
                 cluster.endoderm.color.v5)),
    RowSideColors = t(cbind(
      MyName2Col(names(gene_src.fg4.integrated_re_fg4.rowtree.3), colors.geneset))),
    ColSideColorsSize = 3,
    RowSideColorsSize = 1.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    Rowv = "none", Colv = "none",
    #return.tree = "col",
    graph = T)
dev.off()


# -- GCN for cluster
#----------------------
pearson.gene = WGCNA::cor(
  t(src.fg4.integrated_re_fg4@assays$RNA@data[
    gene_src.fg4.integrated_re_fg4.rowtree.3,
  ]), method = "p")
pearson.gene = logistic(pearson.gene, threshold = 0.25)
pearson.gene[is.na(pearson.gene)]=0
pearson.gene[pearson.gene<0.55]=0

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

rev_gene_src.fg4.integrated_re_fg4.rowtree.3 = names(gene_src.fg4.integrated_re_fg4.rowtree.3)
names(rev_gene_src.fg4.integrated_re_fg4.rowtree.3) = gene_src.fg4.integrated_re_fg4.rowtree.3

pdf("figure.v08.07/organ_development_re_v240115/try.src.fg4.integrated_re_fg4_gcn_tf.pdf",5,5)
plot.igraph(graph.gene,
            layout = layout.gene,
            vertex.size=7,
            label.cex =8,
            vertex.label = NA,
            vertex.label.font = 4,
            vertex.label.cex = 1.2,
            vertex.label.color = "black",
            vertex.color = colors.geneset[
              rev_gene_src.fg4.integrated_re_fg4.rowtree.3[
                names(gene.cluster)]])
dev.off()
#------------------------

# -- Save
#------------------------
save(src.fg4.integrated_re_fg4, 
     file = "figure.v08.07/organ_development_re_v240115/src.fg4.integrated_re_fg4.Rdata")
save(src.fg4.integrated_re_fg4, 
     src.fg4.integrated_re_fg4_filtergene,
     src.fg4.integrated_re_fg4_markergene,
     
     load_src.fg4.integrated_re_fg4_pca3,
     load_src.fg4.integrated_re_fg4_pca1,
     gene_src.fg4.integrated_re_fg4_pca3,
     gene_src.fg4.integrated_re_fg4_pca1,
     
     tree_src.fg4.integrated_re_fg4.rowtree,
     gene_src.fg4.integrated_re_fg4.rowtree,
     
     tree_src.fg4.integrated_re_fg4.rowtree.2,
     gene_src.fg4.integrated_re_fg4.rowtree.2,
     tree_src.fg4.integrated_re_fg4.coltree.2,
     
     tree_src.fg4.integrated_re_fg4.rowtree.1,
     gene_src.fg4.integrated_re_fg4.rowtree.1,
     tree_src.fg4.integrated_re_fg4.coltree.1,
     
     tree_src.fg4.integrated_re_fg4.rowtree.3,
     gene_src.fg4.integrated_re_fg4.rowtree.3,
     tree_src.fg4.integrated_re_fg4.coltree.3,
     cell_src.fg4.integrated_re_fg4_lu, #cluster.v06.26.re_correct_re.tree3
     
     file = "figure.v08.07/organ_development_re_v240115/src.fg4.integrated_re_fg4_parameter.Rdata")

save(cell_src.fg4.integrated_re_fg4_lv, #cluster.v06.26.re_correct_re.tree3
     cell_src.fg4.integrated_re_fg4_lu, #cluster.v06.26.re_correct_re.tree3
     cellorder_src.fg4.integrated_re_fg4,
     file = "~/Bioinformatic/project_20221224_endoderm.refine/figure_v6.05/figure.v08.07/organ_development_re_v240115/src.fg4.integrated_re_fg4_parameter.cell.Rdata")
