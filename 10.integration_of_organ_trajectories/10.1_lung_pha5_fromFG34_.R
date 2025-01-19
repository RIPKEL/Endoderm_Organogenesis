
#-- Lung
#-----------------------
cell_name_lung_fg3 = rownames(src.fg3.integrated_re@meta.data[
  src.fg3.integrated_re$cluster.v06.26.re_correct_refine%in%c("FG.3","Lung"),])
cell_name_lung_fg4 = rownames(src.fg4.integrated_re@meta.data[
  src.fg4.integrated_re$cluster.v06.26.re_correct%in%c("FG.4","FG.4-Lung/Stomach","Lung"),])

src.lung.integrated = merge(
  src.fg4.integrated_re[,cell_name_lung_fg4],
  src.fg3.integrated_re[,setdiff(cell_name_lung_fg3, cell_name_lung_fg4)])


src.lung.integrated = NormalizeData(src.lung.integrated, assay = "RNA", scale.factor = 10^5)
src.lung.integrated$batch_phase
src.lung.integrated = ScaleData(src.lung.integrated, split.by = "batch_phase",
                                features = rownames(src.lung.integrated))
src.lung.integrated = FindVariableFeatures(src.lung.integrated, nfeatures = 2000, assay = "RNA")

src.lung.integrated_filtergene = 
  Myfilter(as.matrix(src.lung.integrated@assays$RNA@data),
           gene = src.lung.integrated@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.lung.integrated_filtergene_re = src.lung.integrated_filtergene
rm(src.lung.integrated_filtergene)
src.lung.integrated_filtergene = src.lung.integrated_filtergene_re

src.lung.integrated = RunPCA(src.lung.integrated, features = src.lung.integrated_filtergene)
src.lung.integrated = RunUMAP(src.lung.integrated, reduction = "pca",
                              dims = 1:25, n.components = 3)

#-- update
src.lung.integrated$cluster.v06.26.re_correct = NA
src.lung.integrated$cluster.v06.26.re_correct[cell_name_lung_fg4] = 
  src.fg4.integrated_re$cluster.v06.26.re_correct[cell_name_lung_fg4]
src.lung.integrated$cluster.v06.26.re_correct[cell_name_lung_fg3] = 
  src.fg3.integrated_re$cluster.v06.26.re_correct_refine[cell_name_lung_fg3]

src.lung.integrated$cluster.v06.26.re_correct.fg4 = NA
src.lung.integrated$cluster.v06.26.re_correct.fg4[cell_name_lung_fg4] = 
  src.fg4.integrated_re$cluster.v06.26.re_correct[cell_name_lung_fg4]

src.lung.integrated$cluster.v06.26.re_correct.fg3 = NA
src.lung.integrated$cluster.v06.26.re_correct.fg3[cell_name_lung_fg3] = 
  src.fg3.integrated_re$cluster.v06.26.re_correct_refine[cell_name_lung_fg3]


DimPlot(src.lung.integrated, group.by = "Time")
DimPlot(src.lung.integrated, group.by = "cluster.v06.26.re_correct",
        cols = cluster.endoderm.color.v5)


pdf("organ_development_re_v240115/try.lung.pdf",9,7)
src.lung.integrated.rowtree =
  MyHeatmap(as.matrix(src.lung.integrated@assays$RNA@data[
    unique(src.lung.integrated_filtergene),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.lung.integrated$Time,
                 colors.time.2),
      MyName2Col(src.lung.integrated$cluster.extract.v1.1.re,
                 cluster.endoderm.color.v5),
      MyName2Col(src.lung.integrated$cluster.v06.26.re_correct,
                 cluster.endoderm.color.v5)),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()
pdf("organ_development_re_v240115/try.lung.tree.pdf",150,20)
plot(src.lung.integrated.rowtree)
dev.off()
tree_src.lung.integrated.rowtree = as.dendrogram(src.lung.integrated.rowtree)
gene_src.lung.integrated.rowtree = 
  setdiff(src.lung.integrated_filtergene, 
          c(labels(tree_src.lung.integrated.rowtree[[1]][[2]][[1]][[2]]),
            labels(tree_src.lung.integrated.rowtree[[2]][[1]]),
            # labels(tree_src.lung.integrated.rowtree[[2]][[2]][[1]][[1]]),
            labels(tree_src.lung.integrated.rowtree[[2]][[2]][[1]][[2]][[1]]),
            labels(tree_src.lung.integrated.rowtree[[2]][[2]][[2]][[1]])))

src.lung.integrated = RunPCA(src.lung.integrated, features = gene_src.lung.integrated.rowtree)
src.lung.integrated = RunUMAP(src.lung.integrated, reduction = "pca",
                              dims = 1:30, n.components = 3)

DimPlot(src.lung.integrated, group.by = "Time", dims = c(1,3))
DimPlot(src.lung.integrated, group.by = "cluster.v06.26.re_correct", dims = c(1,3))

src.lung.integrated = FindNeighbors(src.lung.integrated, reduction = "pca")
src.lung.integrated = FindClusters(src.lung.integrated, resolution = 2)
DimPlot(src.lung.integrated, group.by = "RNA_snn_res.2", dims = c(1,2),label = T)

data_src.lung.integrated = cbind(
  src.lung.integrated@meta.data,
  src.lung.integrated@reductions$umap@cell.embeddings)
colnames(data_src.lung.integrated) = gsub("UMAP","Coord",colnames(data_src.lung.integrated)) 
save(data_src.lung.integrated, file = 'organ_development_re_v240115/data_src.lung.integrated.Rdata')
#-----------------------

#-- Pha.5
#-----------------------
cell_name_pha5_fg3 = rownames(src.fg3.integrated_re@meta.data[
  src.fg3.integrated_re$cluster.v06.26.re_correct_refine%in%c("FG.3","Pharynx.organ.5"),])
cell_name_pha5_fg4 = rownames(src.fg4.integrated_re@meta.data[
  src.fg4.integrated_re$cluster.v06.26.re_correct%in%c("FG.4","Pharynx.organ.5"),])
src.pha5.integrated = merge(
  src.fg4.integrated_re[,cell_name_pha5_fg4],
  src.fg3.integrated_re[,setdiff(cell_name_pha5_fg3, cell_name_pha5_fg4)])


src.pha5.integrated = NormalizeData(src.pha5.integrated, assay = "RNA", scale.factor = 10^5)
src.pha5.integrated = ScaleData(src.pha5.integrated, features = rownames(src.pha5.integrated))
src.pha5.integrated = FindVariableFeatures(src.pha5.integrated, nfeatures = 2000, assay = "RNA")

#-- update
src.pha5.integrated$cluster.v06.26.re_correct = NA
src.pha5.integrated$cluster.v06.26.re_correct[cell_name_pha5_fg4] = 
  src.fg4.integrated_re$cluster.v06.26.re_correct[cell_name_pha5_fg4]
src.pha5.integrated$cluster.v06.26.re_correct[cell_name_pha5_fg3] = 
  src.fg3.integrated_re$cluster.v06.26.re_correct_refine[cell_name_pha5_fg3]

src.pha5.integrated_filtergene = 
  Myfilter(as.matrix(src.pha5.integrated@assays$RNA@data),
           gene = src.pha5.integrated@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.pha5.integrated_filtergene_re = src.pha5.integrated_filtergene
rm(src.pha5.integrated_filtergene)
src.pha5.integrated_filtergene = src.pha5.integrated_filtergene_re

pdf("figure.v08.07/organ_development_re_v240115/try.pha5.pdf",9,7)
src.pha5.integrated.rowtree =
  MyHeatmap(as.matrix(src.pha5.integrated@assays$RNA@data[
    unique(src.pha5.integrated_filtergene),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.pha5.integrated$Time,
                 colors.time.2),
      MyName2Col(src.pha5.integrated$cluster.extract.v1.1.re,
                 cluster.endoderm.color.v5),
      MyName2Col(src.pha5.integrated$cluster.v06.26.re_correct,
                 cluster.endoderm.color.v5)),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()
pdf("figure.v08.07/organ_development_re_v240115/try.pha5.tree.pdf",150,20)
plot(src.pha5.integrated.rowtree)
dev.off()

tree_src.pha5.integrated.rowtree = as.dendrogram(src.pha5.integrated.rowtree)
gene_src.pha5.integrated.rowtree = setdiff(
  src.pha5.integrated_filtergene,
  c(labels(tree_src.pha5.integrated.rowtree[[1]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]]),
    labels(tree_src.pha5.integrated.rowtree[[2]][[2]][[1]]),
    labels(tree_src.pha5.integrated.rowtree[[2]][[1]])))

src.pha5.integrated = RunPCA(src.pha5.integrated, features = gene_src.pha5.integrated.rowtree)
src.pha5.integrated = RunUMAP(src.pha5.integrated, reduction = "pca",
                              dims = 1:30, n.components = 3, seed.use = 3)

src.pha5.integrated = FindNeighbors(src.pha5.integrated, reduction = "pca", dims = 1:30)
src.pha5.integrated = FindClusters(src.pha5.integrated, resolution = 2)

DimPlot(src.pha5.integrated, group.by = "Time")
DimPlot(src.pha5.integrated, group.by = "RNA_snn_res.2", label = T, label.size = 6)
DimPlot(src.pha5.integrated, group.by = "cluster.v06.26.re_correct")
DimPlot(src.pha5.integrated, group.by = "cluster.extract.v1.1")

src.pha5.integrated$cluster.extract.v1.1.re = src.pha5.integrated$cluster.extract.v1.1
src.pha5.integrated@meta.data[src.pha5.integrated$RNA_snn_res.2%in%c(2,3,10,12),]$cluster.extract.v1.1.re = "FG.3/4"
src.pha5.integrated@meta.data[src.pha5.integrated$RNA_snn_res.2%in%c(0,4,14,6,7),]$cluster.extract.v1.1.re = "FG.4"
src.pha5.integrated@meta.data[src.pha5.integrated$RNA_snn_res.2%in%c(8,9,15,11,5,13,1),]$cluster.extract.v1.1.re = "FG.3"
DimPlot(src.pha5.integrated, group.by = "cluster.extract.v1.1.re")

src.pha5.integrated$check = NA
src.pha5.integrated@meta.data[intersect(colnames(src.pha5.integrated.re_v0), colnames(src.pha5.integrated)),]$check = "check"
DimPlot(src.pha5.integrated, group.by = "check")

src.pha5.integrated.re$check = NA
src.pha5.integrated.re@meta.data[intersect(colnames(src.pha5.integrated.re_v0), colnames(src.pha5.integrated)),]$check = "check"
DimPlot(src.pha5.integrated.re, group.by = "check", reduction = "umap.rot")

save(src.pha5.integrated, file = "figure.v08.07/organ_development_re_v240115/src.pha5.integrated.Rdata")
save(tree_src.pha5.integrated.rowtree, 
     gene_src.pha5.integrated.rowtree,
     src.pha5.integrated_filtergene,
     file = "figure.v08.07/organ_development_re_v240115/src.pha5.integrated.parameter.Rdata")

# data_src.pha5.integrated.re = cbind(
#   src.pha5.integrated.re@meta.data,
#   src.pha5.integrated.re@reductions$umap@cell.embeddings)
# colnames(data_src.pha5.integrated.re) = gsub("UMAP","Coord",colnames(data_src.pha5.integrated.re)) 
# save(data_src.pha5.integrated.re_v0, file = 'figure.v08.07/organ_development_re_v240115/data_src.pha5.integrated.re_v0.Rdata')
#
# src.pha5.integrated.re_v0 = src.pha5.integrated.re
# save(src.pha5.integrated.re_v0, file = "figure.v08.07/organ_development_re_v240115/src.pha5.integrated.re_v0.Rdata")

#-- Pha.5-Re
src.pha5.integrated.re = src.pha5.integrated[, !src.pha5.integrated$RNA_snn_res.2%in%c(5,14,15)]
src.pha5.integrated.re = ScaleData(src.pha5.integrated.re, features = rownames(src.pha5.integrated.re), split.by = "batch_phase")
src.pha5.integrated.re = RunPCA(src.pha5.integrated.re, features = gene_src.pha5.integrated.rowtree)
src.pha5.integrated.re = RunUMAP(src.pha5.integrated.re, reduction = "pca",
                                 dims = 1:30, n.components = 3, seed.use = 3)

DimPlot(src.pha5.integrated.re, group.by = "Time")
DimPlot(src.pha5.integrated.re, group.by = "RNA_snn_res.2", label = T, label.size = 6)
DimPlot(src.pha5.integrated.re, group.by = "cluster.v06.26.re_correct")
DimPlot(src.pha5.integrated.re, group.by = "cluster.extract.v1.1.re")
save(src.pha5.integrated.re, file = "figure.v08.07/organ_development_re_v240115/src.pha5.integrated.re.Rdata")
#-----------------------

#-- pha5 * lung --
#-----------------------
src.lung.integrated.pha5 = 
  merge(src.lung.integrated,
        src.pha5.integrated.re[, setdiff(colnames(src.pha5.integrated.re), colnames(src.lung.integrated))])

src.lung.integrated.pha5 = ScaleData(src.lung.integrated.pha5, features = rownames(src.lung.integrated.pha5))
src.lung.integrated.pha5 = FindVariableFeatures(src.lung.integrated.pha5, nfeatures = 2000, assay = "RNA")

src.lung.integrated.pha5_filtergene = 
  Myfilter(as.matrix(src.lung.integrated.pha5@assays$RNA@data),
           gene = src.lung.integrated.pha5@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.lung.integrated.pha5_filtergene_re = src.lung.integrated.pha5_filtergene
rm(src.lung.integrated.pha5_filtergene)
src.lung.integrated.pha5_filtergene = src.lung.integrated.pha5_filtergene_re


#-- update
src.lung.integrated.pha5$cluster.v06.26.re_correct = NA
src.lung.integrated.pha5$cluster.v06.26.re_correct[union(cell_name_pha5_fg4, cell_name_lung_fg4)] = 
  src.fg4.integrated_re$cluster.v06.26.re_correct[union(cell_name_pha5_fg4, cell_name_lung_fg4)]
src.lung.integrated.pha5$cluster.v06.26.re_correct[union(cell_name_pha5_fg3, cell_name_lung_fg3)] = 
  src.fg3.integrated_re$cluster.v06.26.re_correct_refine[union(cell_name_pha5_fg3, cell_name_lung_fg3)]


src.lung.integrated.pha5 = RunPCA(src.lung.integrated.pha5, features = src.lung.integrated.pha5_filtergene)
src.lung.integrated.pha5 = RunUMAP(src.lung.integrated.pha5, reduction = "pca",
                                   dims = 1:30, n.components = 3)


pdf("figure.v08.07/organ_development_re_v240115/try.lung_pha5.pdf",9,7)
src.lung.integrated.pha5.rowtree =
  MyHeatmap(as.matrix(src.lung.integrated.pha5@assays$RNA@data[
    unique(src.lung.integrated.pha5_filtergene),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.lung.integrated.pha5$Time,
                 colors.time.2),
      MyName2Col(src.lung.integrated.pha5$cluster.extract.v1.1.re,
                 cluster.endoderm.color.v5),
      MyName2Col(src.lung.integrated.pha5$cluster.v06.26.re_correct,
                 cluster.endoderm.color.v5)),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()

pdf("figure.v08.07/organ_development_re_v240115/try.lung_pha5.tree.pdf",150,20)
plot(src.lung.integrated.pha5.rowtree)
dev.off()

tree_src.lung.integrated.pha5.rowtree = as.dendrogram(src.lung.integrated.pha5.rowtree)
gene_src.lung.integrated.pha5.rowtree = 
  setdiff(src.lung.integrated.pha5_filtergene, 
          c(labels(tree_src.lung.integrated.pha5.rowtree[[2]][[2]][[2]][[1]]),
            labels(tree_src.lung.integrated.pha5.rowtree[[2]][[2]][[2]][[2]][[1]]),
            labels(tree_src.lung.integrated.pha5.rowtree[[2]][[2]][[2]][[2]][[2]][[1]]),
            labels(tree_src.lung.integrated.pha5.rowtree[[2]][[1]])))

src.lung.integrated.pha5 = RunPCA(src.lung.integrated.pha5, features = gene_src.lung.integrated.pha5.rowtree)
src.lung.integrated.pha5 = RunUMAP(src.lung.integrated.pha5, reduction = "pca",
                                   dims = 1:30, n.components = 3)


src.lung.integrated.pha5@reductions$umap.rot = src.lung.integrated.pha5@reductions$umap
src.lung.integrated.pha5@reductions$umap.rot@cell.embeddings = 
  t(view_src.lung.integrated.pha5[1:3,1:3] %*% 
      t(as.matrix(src.lung.integrated.pha5@meta.data[,c("Coord_1","Coord_2",'Coord_3')])))
colnames(src.lung.integrated.pha5@reductions$umap.rot@cell.embeddings) =  paste("Coord_",c(1:3),sep="")
src.lung.integrated.pha5@reductions$umap.rot@key = "Coord_"

DimPlot(src.lung.integrated.pha5, group.by = "Time", reduction = 'umap.rot')
DimPlot(src.lung.integrated.pha5, group.by = "cluster.v06.26.re_correct")
save(src.lung.integrated.pha5, file = "figure.v08.07/organ_development_re_v240115/src.lung.integrated.pha5.Rdata")

# data_src.lung.integrated.pha5 = cbind(
#   src.lung.integrated.pha5@meta.data,
#   src.lung.integrated.pha5@reductions$umap@cell.embeddings)
# colnames(data_src.lung.integrated.pha5) = gsub("UMAP","Coord",colnames(data_src.lung.integrated.pha5)) 
# data_src.lung.integrated.pha5_v0 = data_src.lung.integrated.pha5
# save(data_src.lung.integrated.pha5_v0, file = 'figure.v08.07/organ_development_re_v240115/data_src.lung.integrated.pha5_v0.Rdata')
#-----------------------


#-- lung * pha5 Add FG.5
#-----------------------
src.lung.integrated.pha5_fg5 = merge(src.lung.integrated.pha5, src.fg5.integrated)
src.lung.integrated.pha5_fg5 = ScaleData(src.lung.integrated.pha5_fg5, features = rownames(src.lung.integrated.pha5_fg5))
src.lung.integrated.pha5_fg5 = FindVariableFeatures(src.lung.integrated.pha5_fg5, nfeatures = 2000, assay = "RNA")

src.lung.integrated.pha5_fg5_filtergene = 
  Myfilter(as.matrix(src.lung.integrated.pha5_fg5@assays$RNA@data),
           gene = src.lung.integrated.pha5_fg5@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.lung.integrated.pha5_fg5_filtergene_re = src.lung.integrated.pha5_fg5_filtergene
rm(src.lung.integrated.pha5_fg5_filtergene)
src.lung.integrated.pha5_fg5_filtergene = src.lung.integrated.pha5_fg5_filtergene_re

src.lung.integrated.pha5_fg5 = RunPCA(src.lung.integrated.pha5_fg5, features = src.lung.integrated.pha5_fg5_filtergene)
src.lung.integrated.pha5_fg5 = RunUMAP(src.lung.integrated.pha5_fg5, reduction = "pca",
                                       dims = 1:30, n.components = 3)


pdf("organ_development_re_v240115/try.lung_pha5_fg5.pdf",9,7)
src.lung.integrated.pha5_fg5.rowtree =
  MyHeatmap(as.matrix(src.lung.integrated.pha5_fg5@assays$RNA@data[
    unique(src.lung.integrated.pha5_fg5_filtergene),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.lung.integrated.pha5_fg5$Time,
                 colors.time.2),
      MyName2Col(src.lung.integrated.pha5_fg5$cluster.extract.v1.1.re,
                 cluster.endoderm.color.v5),
      MyName2Col(src.lung.integrated.pha5_fg5$cluster.v06.26.re_correct,
                 cluster.endoderm.color.v5)),
    ColSideColorsSize = 4,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()
pdf("organ_development_re_v240115/try.lung_pha5_fg5.tree.pdf",150,20)
plot(src.lung.integrated.pha5_fg5.rowtree)
dev.off()
tree_src.lung.integrated.pha5_fg5.rowtree = as.dendrogram(src.lung.integrated.pha5_fg5.rowtree)
gene_src.lung.integrated.pha5_fg5.rowtree = 
  setdiff(src.lung.integrated.pha5_fg5_filtergene, 
          c(labels(tree_src.lung.integrated.pha5_fg5.rowtree[[1]][[1]][[2]][[2]][[2]][[2]]),
            labels(tree_src.lung.integrated.pha5_fg5.rowtree[[2]][[2]])))

src.lung.integrated.pha5_fg5 = RunPCA(src.lung.integrated.pha5_fg5, 
                                      features = gene_src.lung.integrated.pha5_fg5.rowtree)
src.lung.integrated.pha5_fg5 = RunUMAP(src.lung.integrated.pha5_fg5, 
                                       reduction = "pca",
                                       dims = 1:30, n.components = 3)
DimPlot(src.lung.integrated.pha5_fg5, group.by = "Time",
        cols = colors.time)
DimPlot(src.lung.integrated.pha5_fg5, group.by = "cluster.v06.26.re_correct",
        cols = cluster.endoderm.color.v5)

data_src.lung.integrated.pha5_fg5 = cbind(
  src.lung.integrated.pha5_fg5@meta.data,
  src.lung.integrated.pha5_fg5@reductions$umap@cell.embeddings)
colnames(data_src.lung.integrated.pha5_fg5) = gsub("UMAP","Coord",colnames(data_src.lung.integrated.pha5_fg5)) 
save(data_src.lung.integrated.pha5_fg5, file = 'organ_development_re_v240115/data_src.lung.integrated.pha5_fg5.Rdata')
# rm(src.lung.integrated.pha5_fg5)
#-----------------------





