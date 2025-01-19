#---------------------------------
# -- MG.3 Intermeditae Stage
#=================================================================================
src.mg3.integrated_re_mg3 = 
  src.mg3.integrated[, src.mg3.integrated$cluster.v06.26.re%in%c("MG.3",'MG.3.A/M',"MG.3.P")]

src.mg3.integrated_re_mg3 = NormalizeData(src.mg3.integrated_re_mg3, scale.factor = 10^5)
src.mg3.integrated_re_mg3 = ScaleData(src.mg3.integrated_re_mg3, 
                                      features = rownames(src.mg3.integrated_re_mg3),
                                      split.by = "batch_phase")
src.mg3.integrated_re_mg3 = FindVariableFeatures(src.mg3.integrated_re_mg3, assay = "RNA", nfeatures = 2000)
src.mg3.integrated_re_mg3_filtergene = 
  Myfilter(as.matrix(src.mg3.integrated_re_mg3@assays$RNA@data),
           gene = src.mg3.integrated_re_mg3@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.mg3.integrated_re_mg3_filtergene_re = src.mg3.integrated_re_mg3_filtergene
rm(src.mg3.integrated_re_mg3_filtergene)
src.mg3.integrated_re_mg3_filtergene = src.mg3.integrated_re_mg3_filtergene_re

# -- reduction
src.mg3.integrated_re_mg3 = RunPCA(src.mg3.integrated_re_mg3, 
                                   features = src.mg3.integrated_re_mg3_filtergene)
src.mg3.integrated_re_mg3 = RunUMAP(src.mg3.integrated_re_mg3, reduction = "pca",
                                    dims = 1:25, n.components = 3)

pdf("figure.v08.07/organ_development_re_v240115/try.mg3_mg3_pca.pdf",10,10)
DimPlot(src.mg3.integrated_re_mg3, group.by = "Time", 
        cols = colors.time, pt.size = 1.25,
        reduction = "pca", dims = c(1,2)) +
  theme_bw() + p_add
DimPlot(src.mg3.integrated_re_mg3, group.by = "cluster.v06.26.re",
        cols = cluster.endoderm.color.v5, pt.size = 1.25,
        reduction = "pca", dims = c(1,2)) +
  theme_bw() + p_add
dev.off()


# -- marker
src.mg3.integrated_re_mg3 = SetIdent(src.mg3.integrated_re_mg3, 
                                     value = src.mg3.integrated_re_mg3$cluster.v06.26.re)
src.mg3.integrated_re_mg3_marker = FindAllMarkers(src.mg3.integrated_re_mg3)
src.mg3.integrated_re_mg3_marker$pct.ratio =
  src.mg3.integrated_re_mg3_marker$pct.1 / src.mg3.integrated_re_mg3_marker$pct.2
src.mg3.integrated_re_mg3_markergene = unique(
  src.mg3.integrated_re_mg3_marker[
    src.mg3.integrated_re_mg3_marker$pct.ratio>2,]$gene)


#-- Basic processing :: Filter & Marker
#-------------------------------------------------------------------------------
pdf("figure.v08.07/organ_development_re_v240115/try.mg3_mg3.pdf",9,7)
src.mg3.integrated_re_mg3.rowtree =
  MyHeatmap(as.matrix(src.mg3.integrated_re_mg3@assays$RNA@data[
    unique(c(
      src.mg3.integrated_re_mg3_markergene,
      src.mg3.integrated_re_mg3_filtergene)),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.mg3.integrated_re_mg3$Time, colors.time),
      MyName2Col(src.mg3.integrated_re_mg3$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src.mg3.integrated_re_mg3$cluster.v06.26.re, cluster.endoderm.color.v5)),
    ColSideColorsSize = 5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()
pdf("figure.v08.07/organ_development_re_v240115/try.mg3_mg3_tree.pdf",100,20)
plot(src.mg3.integrated_re_mg3.rowtree)
dev.off()

tree_src.mg3.integrated_re_mg3.rowtree = as.dendrogram(src.mg3.integrated_re_mg3.rowtree)
gene_src.mg3.integrated_re_mg3.rowtree = setdiff(
  labels(tree_src.mg3.integrated_re_mg3.rowtree),
  c(labels(tree_src.mg3.integrated_re_mg3.rowtree[[2]][[2]][[1]]),
    labels(tree_src.mg3.integrated_re_mg3.rowtree[[2]][[2]][[2]][[2]])))


pdf("figure.v08.07/organ_development_re_v240115/try.mg3_mg3.pdf",9,7)
src.mg3.integrated_re_mg3.rowtree.1 =
  MyHeatmap(as.matrix(src.mg3.integrated_re_mg3@assays$RNA@data[
    gene_src.mg3.integrated_re_mg3.rowtree, ]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.mg3.integrated_re_mg3$Time, colors.time),
      MyName2Col(src.mg3.integrated_re_mg3$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src.mg3.integrated_re_mg3$cluster.v06.26.re, cluster.endoderm.color.v5)),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    return.tree = "row",
    graph = T)
src.mg3.integrated_re_mg3.coltree.1 =
  MyHeatmap(as.matrix(src.mg3.integrated_re_mg3@assays$RNA@data[
    gene_src.mg3.integrated_re_mg3.rowtree, ]),
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
tree_src.mg3.integrated_re_mg3.rowtree.1 = as.dendrogram(src.mg3.integrated_re_mg3.rowtree.1)
tree_src.mg3.integrated_re_mg3.coltree.1 = as.dendrogram(src.mg3.integrated_re_mg3.coltree.1)
gene_src.mg3.integrated_re_mg3.rowtree.1 = c(
  labels(tree_src.mg3.integrated_re_mg3.rowtree.1[[2]][[1]][[2]]),
  labels(tree_src.mg3.integrated_re_mg3.rowtree.1[[2]][[2]][[2]][[2]]),
  labels(tree_src.mg3.integrated_re_mg3.rowtree.1[[2]][[2]][[1]]),
  labels(tree_src.mg3.integrated_re_mg3.rowtree.1[[1]][[1]]),
  labels(tree_src.mg3.integrated_re_mg3.rowtree.1[[1]][[2]]))


pdf("figure.v08.07/organ_development_re_v240115/try.mg3_mg3.pdf",9,7)
src.mg3.integrated_re_mg3.rowtree.2 =
  MyHeatmap(as.matrix(src.mg3.integrated_re_mg3@assays$RNA@data[
    gene_src.mg3.integrated_re_mg3.rowtree.1, ]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.mg3.integrated_re_mg3$Time, colors.time),
      #MyName2Col(src.mg3.integrated_re_mg3$tree.2, color.temp),
      MyName2Col(src.mg3.integrated_re_mg3$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src.mg3.integrated_re_mg3$cluster.v06.26.re, cluster.endoderm.color.v5)),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    return.tree = "row",
    graph = T)
src.mg3.integrated_re_mg3.coltree.2 =
  MyHeatmap(as.matrix(src.mg3.integrated_re_mg3@assays$RNA@data[
    gene_src.mg3.integrated_re_mg3.rowtree.1, ]),
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

tree_src.mg3.integrated_re_mg3.rowtree.2 = as.dendrogram(src.mg3.integrated_re_mg3.rowtree.2)
gene_src.mg3.integrated_re_mg3.rowtree.2 = c(
  labels(tree_src.mg3.integrated_re_mg3.rowtree.2[[2]][[1]][[1]]),
  # labels(tree_src.mg3.integrated_re_mg3.rowtree.2[[2]][[1]][[2]][[2]]),
  
  labels(tree_src.mg3.integrated_re_mg3.rowtree.2[[2]][[2]][[1]][[2]]),
  labels(tree_src.mg3.integrated_re_mg3.rowtree.2[[2]][[2]][[1]][[1]]),
  labels(tree_src.mg3.integrated_re_mg3.rowtree.2[[2]][[2]][[2]][[1]]),
  
  labels(tree_src.mg3.integrated_re_mg3.rowtree.2[[1]][[2]][[1]]),
  labels(tree_src.mg3.integrated_re_mg3.rowtree.2[[1]][[2]][[2]][[1]]),
  
  labels(tree_src.mg3.integrated_re_mg3.rowtree.2[[1]][[1]][[1]][[1]]),
  labels(tree_src.mg3.integrated_re_mg3.rowtree.2[[1]][[1]][[1]][[2]][[1]]),
  labels(tree_src.mg3.integrated_re_mg3.rowtree.2[[1]][[1]][[1]][[2]][[2]][[1]]),
  labels(tree_src.mg3.integrated_re_mg3.rowtree.2[[1]][[1]][[2]]))
names(gene_src.mg3.integrated_re_mg3.rowtree.2) = c(
  rep(4, length(c(labels(tree_src.mg3.integrated_re_mg3.rowtree.2[[2]][[1]][[1]])
                  # labels(tree_src.mg3.integrated_re_mg3.rowtree.2[[2]][[1]][[2]][[2]])
  ))),
  rep(3, length(c(labels(tree_src.mg3.integrated_re_mg3.rowtree.2[[2]][[2]][[1]][[2]])))),
  rep(3, length(c(labels(tree_src.mg3.integrated_re_mg3.rowtree.2[[2]][[2]][[1]][[1]])))),
  rep(3, length(c(labels(tree_src.mg3.integrated_re_mg3.rowtree.2[[2]][[2]][[2]][[1]])))),
  
  rep(7, length(c(labels(tree_src.mg3.integrated_re_mg3.rowtree.2[[1]][[2]][[1]]),
                  labels(tree_src.mg3.integrated_re_mg3.rowtree.2[[1]][[2]][[2]][[1]]),
                  labels(tree_src.mg3.integrated_re_mg3.rowtree.2[[1]][[1]][[1]][[1]]),
                  labels(tree_src.mg3.integrated_re_mg3.rowtree.2[[1]][[1]][[1]][[2]][[1]]),
                  labels(tree_src.mg3.integrated_re_mg3.rowtree.2[[1]][[1]][[1]][[2]][[2]][[1]]),
                  labels(tree_src.mg3.integrated_re_mg3.rowtree.2[[1]][[1]][[2]])))))


tree_src.mg3.integrated_re_mg3.coltree.2 = as.dendrogram(src.mg3.integrated_re_mg3.coltree.2)
src.mg3.integrated_re_mg3$tree.2 =NA
src.mg3.integrated_re_mg3@meta.data[labels(tree_src.mg3.integrated_re_mg3.coltree.2[[1]][[1]][[1]]), ]$tree.2 = 1
src.mg3.integrated_re_mg3@meta.data[labels(tree_src.mg3.integrated_re_mg3.coltree.2[[1]][[1]][[2]]), ]$tree.2 = 2
src.mg3.integrated_re_mg3@meta.data[labels(tree_src.mg3.integrated_re_mg3.coltree.2[[1]][[2]][[1]]), ]$tree.2 = 3
src.mg3.integrated_re_mg3@meta.data[labels(tree_src.mg3.integrated_re_mg3.coltree.2[[1]][[2]][[2]]), ]$tree.2 = 4
src.mg3.integrated_re_mg3@meta.data[labels(tree_src.mg3.integrated_re_mg3.coltree.2[[2]][[1]][[1]]), ]$tree.2 = 5
src.mg3.integrated_re_mg3@meta.data[labels(tree_src.mg3.integrated_re_mg3.coltree.2[[2]][[1]][[2]]), ]$tree.2 = 6
src.mg3.integrated_re_mg3@meta.data[labels(tree_src.mg3.integrated_re_mg3.coltree.2[[2]][[2]][[1]]), ]$tree.2 = 7
src.mg3.integrated_re_mg3@meta.data[labels(tree_src.mg3.integrated_re_mg3.coltree.2[[2]][[2]][[2]][[1]]), ]$tree.2 = 8
src.mg3.integrated_re_mg3@meta.data[labels(tree_src.mg3.integrated_re_mg3.coltree.2[[2]][[2]][[2]][[2]]), ]$tree.2 = 9

src.mg3.integrated_re_mg3$cluster.v06.26.re_correct = src.mg3.integrated_re_mg3$cluster.v06.26.re
src.mg3.integrated_re_mg3@meta.data[src.mg3.integrated_re_mg3$tree.2%in%c(7,8), ]$cluster.v06.26.re_correct = "MG.3"
src.mg3.integrated_re_mg3@meta.data[src.mg3.integrated_re_mg3$tree.2%in%c(5,6,9), ]$cluster.v06.26.re_correct = "MG.3.P"
src.mg3.integrated_re_mg3@meta.data[src.mg3.integrated_re_mg3$tree.2%in%c(1,2), ]$cluster.v06.26.re_correct = "MG.3.A/M"
src.mg3.integrated_re_mg3@meta.data[src.mg3.integrated_re_mg3$tree.2%in%c(3,4), ]$cluster.v06.26.re_correct = "MG.3.A/M"
src.mg3.integrated_re_mg3@meta.data[src.mg3.integrated_re_mg3$Time%in%c("ss9"), ]$cluster.v06.26.re_correct = "MG.3"

DimPlot(src.mg3.integrated_re_mg3, reduction = "umap_mnn", group.by = "cluster.v06.26.re_correct")
#------------------------------------------------
#-- MG.3.A/M --
cell_src.mg3.integrated_re_mg3_m3a = 
  rownames(src.mg3.integrated_re_mg3@meta.data[src.mg3.integrated_re_mg3$cluster.v06.26.re_correct%in%c("MG.3","MG.3.A/M"),])
coord_src.mg3.integrated_re_mg3_m3a = cbind(
  src.mg3.integrated_re_mg3[["umap_mnn"]]@cell.embeddings[cell_src.mg3.integrated_re_mg3_m3a, 2],
  src.mg3.integrated_re_mg3[["umap_mnn"]]@cell.embeddings[cell_src.mg3.integrated_re_mg3_m3a, 1])
pcurve_src.mg3.integrated_re_mg3_m3a = princurve::principal_curve(x = coord_src.mg3.integrated_re_mg3_m3a, smoother = "smooth.spline")
src.mg3.integrated_re_mg3$lambda_m3a = NA
src.mg3.integrated_re_mg3@meta.data[cell_src.mg3.integrated_re_mg3_m3a,]$lambda_m3a = 
  pcurve_src.mg3.integrated_re_mg3_m3a$lambda[cell_src.mg3.integrated_re_mg3_m3a]
src.mg3.integrated_re_mg3@meta.data[!src.mg3.integrated_re_mg3$lambda_m3a%in%NA,]$lambda_m3a = 
  norm_range(src.mg3.integrated_re_mg3@meta.data[!src.mg3.integrated_re_mg3$lambda_m3a%in%NA,]$lambda_m3a)

#-- MG.3.P --
cell_src.mg3.integrated_re_mg3_m3p = 
  rownames(src.mg3.integrated_re_mg3@meta.data[src.mg3.integrated_re_mg3$cluster.v06.26.re_correct%in%c("MG.3", "MG.3.P"),])
coord_src.mg3.integrated_re_mg3_m3p = cbind(
  src.mg3.integrated_re_mg3[["umap_mnn"]]@cell.embeddings[cell_src.mg3.integrated_re_mg3_m3p, 2],
  src.mg3.integrated_re_mg3[["umap_mnn"]]@cell.embeddings[cell_src.mg3.integrated_re_mg3_m3p, 1])
pcurve_src.mg3.integrated_re_mg3_m3p = princurve::principal_curve(x = coord_src.mg3.integrated_re_mg3_m3p, smoother = "smooth.spline")
src.mg3.integrated_re_mg3$lambda_m3p = NA
src.mg3.integrated_re_mg3@meta.data[cell_src.mg3.integrated_re_mg3_m3p,]$lambda_m3p = 
  pcurve_src.mg3.integrated_re_mg3_m3p$lambda[cell_src.mg3.integrated_re_mg3_m3p]
src.mg3.integrated_re_mg3@meta.data[!src.mg3.integrated_re_mg3$lambda_m3p%in%NA,]$lambda_m3p = 
  norm_range(src.mg3.integrated_re_mg3@meta.data[!src.mg3.integrated_re_mg3$lambda_m3p%in%NA,]$lambda_m3p)

#-- MG.3 --
cell_src.mg3.integrated_re_mg3_mg3 = 
  rownames(src.mg3.integrated_re_mg3@meta.data[src.mg3.integrated_re_mg3$cluster.v06.26.re_correct%in%c("MG.3"),])

w1 = length(cell_src.mg3.integrated_re_mg3_m3a)
# w2 = length(cell_src.mg3.integrated_re_mg3_m3m)
w3 = length(cell_src.mg3.integrated_re_mg3_m3p)

src.mg3.integrated_re_mg3$lambda_mg3 = NA
src.mg3.integrated_re_mg3$lambda_mg3[cell_src.mg3.integrated_re_mg3_mg3] = 
  (w1 * src.mg3.integrated_re_mg3$lambda_m3a[cell_src.mg3.integrated_re_mg3_mg3] +
     # w2 * src.mg3.integrated_re_mg3$lambda_m3m[cell_src.mg3.integrated_re_mg3_mg3] + 
     w3 * src.mg3.integrated_re_mg3$lambda_m3p[cell_src.mg3.integrated_re_mg3_mg3]) / (w1 +# w2+
                                                                                         w3)
src.mg3.integrated_re_mg3@meta.data[!src.mg3.integrated_re_mg3$lambda_mg3%in%NA,]$lambda_mg3 = 
  norm_range(src.mg3.integrated_re_mg3@meta.data[!src.mg3.integrated_re_mg3$lambda_mg3%in%NA,]$lambda_mg3)


cellorder_src.mg3.integrated_re_mg3 = c(
  intersect(names(src.mg3.integrated_re_mg3$lambda_mg3[order(src.mg3.integrated_re_mg3$lambda_mg3)]), 
            colnames(src.mg3.integrated_re_mg3[,src.mg3.integrated_re_mg3$cluster.v06.26.re_correct%in%c("MG.3")])),
  intersect(names(src.mg3.integrated_re_mg3$lambda_m3a[order(src.mg3.integrated_re_mg3$lambda_m3a)]), 
            colnames(src.mg3.integrated_re_mg3[,src.mg3.integrated_re_mg3$cluster.v06.26.re_correct%in%c("MG.3.A/M"),])),
  # intersect(names(src.mg3.integrated_re_mg3$lambda_m3m[order(src.mg3.integrated_re_mg3$lambda_m3m)]), 
  #           colnames(src.mg3.integrated_re_mg3[,src.mg3.integrated_re_mg3$cluster.v06.26.re_correct%in%c("MG.3.A/M"),])),
  intersect(names(src.mg3.integrated_re_mg3$lambda_m3p[order(src.mg3.integrated_re_mg3$lambda_m3p)]), 
            colnames(src.mg3.integrated_re_mg3[,src.mg3.integrated_re_mg3$cluster.v06.26.re_correct%in%c('MG.3.P'),])))
#------------------------------------------------


pdf("figure.v08.07/organ_development_re_v240115/try.mg3_mg3.pdf",9,7)
src.mg3.integrated_re_mg3.rowtree.3 =
  MyHeatmap(as.matrix(src.mg3.integrated_re_mg3@assays$RNA@data[
    gene_src.mg3.integrated_re_mg3.rowtree.2,
    cellorder_src.mg3.integrated_re_mg3]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.mg3.integrated_re_mg3$Time[cellorder_src.mg3.integrated_re_mg3], 
                 colors.time),
      MyName2Col(src.mg3.integrated_re_mg3$cluster.extract.v1.1[cellorder_src.mg3.integrated_re_mg3], 
                 cluster.endoderm.color.v5),
      MyName2Col(src.mg3.integrated_re_mg3$cluster.v06.26.re_correct[cellorder_src.mg3.integrated_re_mg3], 
                 cluster.endoderm.color.v5)),
    RowSideColors = t(cbind(
      MyName2Col(names(gene_src.mg3.integrated_re_mg3.rowtree.2), colors.geneset))),
    ColSideColorsSize = 3,
    RowSideColorsSize = 1.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    Colv = "none", Rowv = "none",
    # return.tree = "row",
    graph = T)
dev.off()

save(tree_src.mg3.integrated_re_mg3.rowtree.1,
     tree_src.mg3.integrated_re_mg3.rowtree.2,
     gene_src.mg3.integrated_re_mg3.rowtree.1,
     gene_src.mg3.integrated_re_mg3.rowtree.2,
     src.mg3.integrated_re_mg3,
     file = "figure.v08.07/organ_development_re_v240115/src.mg3.integrated_re_mg3.parameter.Rdata")

# -- GCN for cluster
#----------------------
pearson.gene = WGCNA::cor(
  t(src.mg3.integrated_re_mg3@assays$RNA@data[
    gene_src.mg3.integrated_re_mg3.rowtree.2,]), method = "p")
pearson.gene = logistic(pearson.gene, threshold = 0)
pearson.gene[is.na(pearson.gene)]=0
pearson.gene[pearson.gene<0.15]=0

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
set.seed(6)

rev_gene_src.mg3.integrated_re_mg3.rowtree.2 = names(gene_src.mg3.integrated_re_mg3.rowtree.2)
names(rev_gene_src.mg3.integrated_re_mg3.rowtree.2) = gene_src.mg3.integrated_re_mg3.rowtree.2

pdf("figure.v08.07/organ_development_re_v240115/try.src.mg3.integrated_re_mg3_gcn_tf.pdf",5,5)
plot.igraph(graph.gene,
            layout = layout.gene,
            vertex.size=7,
            label.cex =8,
            vertex.label = NA,
            vertex.label.font = 4,
            vertex.label.cex = 1.2,
            vertex.label.color = "black",
            vertex.color = colors.geneset[
              rev_gene_src.mg3.integrated_re_mg3.rowtree.2[
                names(gene.cluster)]])
dev.off()
#------------------------