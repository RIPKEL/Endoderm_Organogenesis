#-------------------------------
#  Part1: FG.3
#-------------------------------
src.endoderm.fg3.ext_v1.1 = FindVariableFeatures(src.endoderm.fg3.ext_v1.1, nfeatures = 2000)
src.endoderm.fg3.ext_v1.1 = ScaleData(src.endoderm.fg3.ext_v1.1, split.by = "batch_phase",
                                      features = rownames(src.endoderm.fg3.ext_v1.1))
src.endoderm.fg3.ext_v1.1.filtergene = 
  Myfilter(as.matrix(src.endoderm.fg3.ext_v1.1@assays$RNA@data),
           gene = src.endoderm.fg3.ext_v1.1@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
#rm(src.endoderm.fg3.ext_v1.1.filtergene)
src.endoderm.fg3.ext.v1.1.selectgene = unique(c(
  select_gene_fg3,src.endoderm.fg3.ext_v1.1.filtergene))

#=============================
# Row
pdf("figure.v08.07/try.pdf")
#-----------------
fg3.ext.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.fg3.ext_v1.1@assays$RNA@data[
    src.endoderm.fg3.ext_v1.1.gene.fin,]),
            type = "row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.fg3.ext_v1.1$Time,colors.time),
              MyName2Col(src.endoderm.fg3.ext_v1.1$cluster.extract.v1.1,cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.fg3.ext_v1.1$cluster.v06.26.re,cluster.endoderm.color.v5)
            ),
            ColSideColorsSize = 1.5,
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            #Rowv = "none",
            return.tree = "row",
            graph = T)
#-------------------
dev.off()

#fg3.ext.re.rowtree_1 = as.dendrogram(fg3.ext.re.rowtree)
#fg3.ext.re.rowtree_2 = as.dendrogram(fg3.ext.re.rowtree)
#fg3.ext.re.rowtree_3 = as.dendrogram(fg3.ext.re.rowtree)
fg3.ext.re.rowtree_4 = as.dendrogram(fg3.ext.re.rowtree)

src.endoderm.fg3.ext_v1.1.gene.fin = setdiff(
  unique(src.endoderm.fg3.ext.re_v1.1.gene.fin,
         src.endoderm.fg3.ext_v1.1.gene.fin),
  c(labels(fg3.ext.re.rowtree_1[[2]][[1]]),
    labels(fg3.ext.re.rowtree_1[[2]][[2]][[1]][[2]][[2]]),
    labels(fg3.ext.re.rowtree_2[[2]][[1]]),
    labels(fg3.ext.re.rowtree_3[[2]][[2]][[2]][[1]][[1]]),
    labels(fg3.ext.re.rowtree_4[[2]][[2]])))

#=============================
# Col
pdf("figure.v08.07/try.pdf")
#-----------------
fg3.ext.re.coltree =
  MyHeatmap(as.matrix(src.endoderm.fg3.ext_v1.1@assays$RNA@data[
    src.endoderm.fg3.ext_v1.1.gene.fin,]),
            type = "row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.fg3.ext_v1.1$Time,colors.time),
              MyName2Col(src.endoderm.fg3.ext_v1.1$cluster.extract.v1.1,cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.fg3.ext_v1.1$cluster.v06.26.re,cluster.endoderm.color.v5)
            ),
            ColSideColorsSize = 1.5,
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            #Rowv = "none",
            return.tree = "col",
            graph = T)
#-------------------
dev.off()

fg3.ext.re.coltree_1 = as.dendrogram(fg3.ext.re.coltree)

src.endoderm.fg3.ext_v1.1$cluster.v06.26.re_hc = NA
src.endoderm.fg3.ext_v1.1@meta.data[c(
  labels(fg3.ext.re.coltree_1[[1]][[2]][[2]][[1]])), 
  "cluster.v06.26.re_hc"] = "FG.3"

src.endoderm.fg3.ext_v1.1@meta.data[c(
  labels(fg3.ext.re.coltree_1[[1]][[2]][[2]][[2]][[2]])), 
  "cluster.v06.26.re_hc"] = "FG.3"

src.endoderm.fg3.ext_v1.1@meta.data[setdiff(
  labels(fg3.ext.re.coltree_1[[1]]),
  c(  
    labels(fg3.ext.re.coltree_1[[1]][[2]][[2]][[1]]),
    labels(fg3.ext.re.coltree_1[[1]][[2]][[2]][[2]][[2]])
  )), "cluster.v06.26.re_hc"] = "Pharynx.organ.4"

src.endoderm.fg3.ext_v1.1@meta.data[c(
  labels(fg3.ext.re.coltree_1[[2]][[2]][[2]][[1]])), 
  "cluster.v06.26.re_hc"] = "FG.3-Lung"

src.endoderm.fg3.ext_v1.1@meta.data[setdiff(
  labels(fg3.ext.re.coltree_1[[2]]),
  c(  
    labels(fg3.ext.re.coltree_1[[2]][[2]][[2]][[1]])
  )), "cluster.v06.26.re_hc"] = "Lung"


#--------------------------------------------------------------------------------
src.endoderm.fg3.ext.re_v1.1.gene.fin =
  rownames(src.endoderm.fg3.ext.re@reductions$pca@feature.loadings)

src.endoderm.fg3.ext_v1.1 = 
  seurat_mnn(seurat_name = "src.endoderm.fg3.ext_v1.1", dims = 1:30, 
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.fg3.ext_v1.1.gene.fin)
src.endoderm.fg3.ext_v1.1 = 
  seurat_int(seurat_name = "src.endoderm.fg3.ext_v1.1", dims = 1:30,
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.fg3.ext_v1.1.gene.fin)


# UMAP Integrated Rotation
#-------------------------------------------------
src.endoderm.fg3.ext_v1.1@reductions$umap_int_rotated = 
  src.endoderm.fg3.ext_v1.1@reductions$umap
src.endoderm.fg3.ext_v1.1@reductions$umap_int_rotated@cell.embeddings = 
  fg3_umap_integrated_rotatedEmbedding[colnames(src.endoderm.fg3.ext_v1.1),1:3]
src.endoderm.fg3.ext_v1.1@reductions$umap_int_rotated@feature.loadings.projected = 
  fg3_umap_integrated_rotation
colnames(src.endoderm.fg3.ext_v1.1@reductions$umap_int_rotated@cell.embeddings) = 
  c("UMAP_1","UMAP_2","UMAP_3")
#-------------------------------------------------

src.endoderm.fg3.ext_v1.1 = 
  FindNeighbors(src.endoderm.fg3.ext_v1.1, dims=1:30, 
                reduction = "pca_integrated", assay = "integrated")
src.endoderm.fg3.ext_v1.1 = 
  FindNeighbors(src.endoderm.fg3.ext_v1.1, dims=1:30, 
                reduction = "pca", assay = "RNA")
src.endoderm.fg3.ext_v1.1 = 
  FindClusters(src.endoderm.fg3.ext_v1.1, 
               resolution = 2, graph.name = "RNA_snn",)


DimPlot(src.endoderm.fg3.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.2', reduction = "fdl_pca_integrated")
DimPlot(src.endoderm.fg3.ext_v1.1, group.by = 'Time', 
        reduction = "umap_mnn", cols = colors.time)
DimPlot(src.endoderm.fg3.ext_v1.1, group.by = 'Time', 
        reduction = "umap_integrated", cols = colors.time)
DimPlot(src.endoderm.fg3.ext_v1.1, group.by = 'batch', 
        reduction = "umap_integrated") #, cols = colors.time)

DimPlot(src.endoderm.fg3.ext_v1.1, group.by = 'cluster.extract.v1.1', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.fg3.ext_v1.1, group.by = 'cluster.v06.26.re', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.fg3.ext_v1.1, group.by = 'cluster.v06.26.re_hc',
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.fg3.ext_v1.1, group.by = 'cluster.v06.26.re_hc',
        reduction = "umap_int_rotated", cols = cluster.endoderm.color.v5)


src.endoderm.fg3.ext_v1.1$cluster.v06.26.re = NA
src.endoderm.fg3.ext_v1.1@meta.data[
  intersect(colnames(src.endoderm.fg3.ext_v1.1), colnames(src.endoderm.fg3.ext.re)),]$cluster.v06.26.re =
  src.endoderm.fg3.ext.re@meta.data[
    intersect(colnames(src.endoderm.fg3.ext_v1.1), colnames(src.endoderm.fg3.ext.re)),]$cluster.v06.26
a = FNN::knn(
  src.endoderm.fg3.ext_v1.1@reductions$umap_integrated@cell.embeddings[
    rownames(src.endoderm.fg3.ext_v1.1@meta.data[!src.endoderm.fg3.ext_v1.1$cluster.v06.26.re%in%NA,]),],
  src.endoderm.fg3.ext_v1.1@reductions$umap_integrated@cell.embeddings[
    rownames(src.endoderm.fg3.ext_v1.1@meta.data[src.endoderm.fg3.ext_v1.1$cluster.v06.26.re%in%NA,]),],
  src.endoderm.fg3.ext_v1.1@meta.data[!src.endoderm.fg3.ext_v1.1$cluster.v06.26.re%in%NA,]$cluster.v06.26.re, k = 10)
src.endoderm.fg3.ext_v1.1@meta.data[src.endoderm.fg3.ext_v1.1$cluster.v06.26.re%in%NA,]$cluster.v06.26.re = as.character(a)



# FDL
src.endoderm.fg3.ext_v1.1 = 
  seurat_fdl("src.endoderm.fg3.ext_v1.1", "pca","RNA")
src.endoderm.fg3.ext_v1.1 = 
  seurat_fdl("src.endoderm.fg3.ext_v1.1", "mnn","RNA")
src.endoderm.fg3.ext_v1.1 = 
  seurat_fdl("src.endoderm.fg3.ext_v1.1", "pca_integrated","integrated")

# Graph
graph.fg3.ext_v1.1 = 
  seurat_GCN(seurat_name = "src.endoderm.fg3.ext_v1.1",
             src.endoderm.fg3.ext_v1.1.gene.fin,
             reduction = "pca_integrated", assay = "integrated")
graph.fg3_cluster = cluster_walktrap(graph.fg3.ext_v1.1, steps = 1)
#----------------------------


pdf("figure.v08.07/organ_development_re/fg3_summary_graph.re_snn.pdf",12,7)
# SNN
for(color in c("cluster.extract.v1.1","cluster.extract.v1.1.re",
               "cluster.v06.26.re","cluster.v06.26.re_hc","Time")){
  for(edge.color in c(NA,"lightgray")){
    for(edge.color in c(NA,"lightgray")){
      set.seed(1)
      plot.igraph(graph.fg3.ext_v1.1,
                  layout = layout_with_fr(graph.fg3.ext_v1.1),
                  edge.color = edge.color,
                  vertex.size=2.5, vertex.label=NA,
                  vertex.label.color = "black", 
                  vertex.frame.color = NA, vertex.frame.width = 0.5,
                  #vertex.color =  colors.time[src.endoderm.fg3.ext_v1.1[["Time"]][graph.fg3_cluster$names,]],
                  #vertex.color =  cluster.endoderm.color.v5[src.endoderm.fg3.ext_v1.1[["cluster.v06.26.re"]][graph.fg3_cluster$names,]],
                  vertex.color =  colors_bg[src.endoderm.fg3.ext_v1.1[[color]][graph.fg3_cluster$names,]])
    }}}
dev.off()


src.endoderm.fg3.ext_v1.1$cluster.v06.26.re = src.endoderm.fg3.ext_v1.1$cluster.v06.26.re_hc
src.endoderm.fg3.ext_v1.1@meta.data[
  src.endoderm.fg3.ext_v1.1$cluster.v06.26.re_hc%in%"FG.3-Lung",]$cluster.v06.26.re_hc = "Lung"

pdf("figure.v08.07/organ_development_re/fg3_summary_graph.re_other.pdf",12,7)
# Reduction
for(red in c("umap_integrated","umap_int_rotated","umap_mnn",
             "fdl_pca",'fdl_mnn',"fdl_pca_integrated")){
  p1 = DimPlot(src.endoderm.fg3.ext_v1.1,group.by = 'cluster.v06.26.re', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p2 = DimPlot(src.endoderm.fg3.ext_v1.1,group.by = 'cluster.extract.v1.1', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p21 = DimPlot(src.endoderm.fg3.ext_v1.1,group.by = 'cluster.extract.v1.1.re', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p3 = DimPlot(src.endoderm.fg3.ext_v1.1,group.by = 'Time', 
               reduction = red , cols = colors.time,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p4 = DimPlot(src.endoderm.fg3.ext_v1.1,group.by = 'cluster.v06.26.re_hc', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  print(p1); print(p2); print(p21);print(p3);print(p4)
}
dev.off()


pdf("figure.v08.07/organ_development_re/fg3_summary_proportion.re.pdf",12,7)
#-----------------------------
cell_type = c("FG.3","FG.3-Lung","Lung","Pharynx.organ.4")
meta_fg3 = src.endoderm.fg3.ext_v1.1@meta.data
meta_fg3$Time = factor(meta_fg3$Time, levels = names(colors.time))
meta_fg3$cluster.extract.v1.1 = factor(meta_fg3$cluster.extract.v1.1, levels = c("FG.3","FG.3/4"))
meta_fg3$cluster.extract.v1.1 = factor(meta_fg3$cluster.extract.v1.1.re, levels = c("FG.3","FG.3/4"))
meta_fg3$cluster.v06.26.re = factor(meta_fg3$cluster.v06.26.re, levels = rev(cell_type))
meta_fg3$cluster.v06.26.re_hc = factor(meta_fg3$cluster.v06.26.re_hc, levels = rev(cell_type))

meta_fg3_B0 = meta_fg3
meta_fg3_B1 = meta_fg3[meta_fg3$batch%in%1&!meta_fg3$Time%in%"ss9",]
meta_fg3_B2 = meta_fg3[meta_fg3$batch%in%2,]

for(data in c("meta_fg3_B0","meta_fg3_B1","meta_fg3_B2")){
  p = ggplot() +
    scale_fill_manual(values = cluster.endoderm.color.v5)+
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle=90, size = 13),
          axis.text.y = element_text(size = 13),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          #legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          plot.title = element_blank(),
          aspect.ratio=1/0.5) +
    guides(color = guide_legend(override.aes = list(size = 5))) 
  
  print(p + 
          geom_bar(data = get(data), 
                   mapping = aes(x = Time,
                                 group = cluster.extract.v1.1, fill = cluster.extract.v1.1),
                   stat = "count", position = 'fill') + 
          xlab("Time") + ylab("Cluster proportion of
 tracing code")+ ggtitle(data))
  
  print(p + 
          geom_bar(data = get(data), 
                   mapping = aes(x = Time,
                                 group = cluster.extract.v1.1, fill = cluster.extract.v1.1.re),
                   stat = "count", position = 'fill') + 
          xlab("Time") + ylab("Cluster proportion of
 tracing code")+ ggtitle(data))
  
  print(p + 
          geom_bar(data = get(data), 
                   mapping = aes(x = Time,
                                 group = cluster.v06.26.re, fill = cluster.v06.26.re),
                   stat = "count", position = 'fill') + 
          xlab("Time") + ylab("Cluster proportion of
 cell type")+ ggtitle(data))
  
  print(p + 
          geom_bar(data = get(data), 
                   mapping = aes(x = Time,
                                 group = cluster.v06.26.re_hc, fill = cluster.v06.26.re_hc),
                   stat = "count", position = 'fill') + 
          xlab("Time") + ylab("Cluster proportion of
 cell type")+ ggtitle(data))
  
}
#-----------------------------
dev.off()

save(src.endoderm.fg3.ext_v1.1,
     file = "figure.v08.07/organ_development_re/src.endoderm.fg3.ext_v1.1.re.Rdata")


# Marker
#------------------------
src.endoderm.fg3.ext_v1.1 = 
  SetIdent(src.endoderm.fg3.ext_v1.1,
           value = src.endoderm.fg3.ext_v1.1$cluster.v06.26.re_hc)
markergene.fg3.ext_v1.1 = 
  FindAllMarkers(src.endoderm.fg3.ext_v1.1, assay = "RNA")
markergene.fg3.ext_v1.1 = markergene.fg3.ext_v1.1[
  markergene.fg3.ext_v1.1$avg_log2FC>0.2&
    markergene.fg3.ext_v1.1$p_val_adj<0.1, "gene"]


fg3.umap.embedding = src.endoderm.fg3.ext_v1.1[["fdl_pca_integrated"]]@cell.embeddings[,1:3]
fg3.umap.embedding = src.endoderm.fg3.ext_v1.1[["umap_int_rotated"]]@cell.embeddings[,1:3]
fg3.umap.embedding = src.endoderm.fg3.ext_v1.1[["fdl_mnn"]]@cell.embeddings[,1:3]
# Cell order :: Princurve
#---------------------------
#  FG.3
cellorder.fg3_fg3 = 
  rownames(src.endoderm.fg3.ext_v1.1@meta.data[
    src.endoderm.fg3.ext_v1.1$cluster.v06.26.re_hc%in%"FG.3",])
princurve.fg3_fg3 = princurve::principal_curve(
  x = fg3.umap.embedding[cellorder.fg3_fg3,],  smoother = "smooth.spline")
cellorder.fg3_fg3 =  names(
  princurve.fg3_fg3$lambda[cellorder.fg3_fg3][order(princurve.fg3_fg3$lambda[cellorder.fg3_fg3])])

# FG.3 - Pha.2
cellorder.fg3_pha4 = 
  rownames(src.endoderm.fg3.ext_v1.1@meta.data[
    src.endoderm.fg3.ext_v1.1$cluster.v06.26.re_hc%in%"Pharynx.organ.4",])
princurve.fg3_pha4 = princurve::principal_curve(
  x = fg3.umap.embedding[cellorder.fg3_pha4,],  smoother = "smooth.spline")
cellorder.fg3_pha4 =  names(
  princurve.fg3_pha4$lambda[cellorder.fg3_pha4][order(princurve.fg3_pha4$lambda[cellorder.fg3_pha4])])

# FG.3-Lung
cellorder.fg3_fg3lu = 
  rownames(src.endoderm.fg3.ext_v1.1@meta.data[
    src.endoderm.fg3.ext_v1.1$cluster.v06.26.re_hc%in%"FG.3-Lung",])
princurve.fg3_fg3lu = princurve::principal_curve(
  x = fg3.umap.embedding[cellorder.fg3_fg3lu,],  smoother = "smooth.spline")
cellorder.fg3_fg3lu =  names(
  princurve.fg3_fg3lu$lambda[cellorder.fg3_fg3lu][order(princurve.fg3_fg3lu$lambda[cellorder.fg3_fg3lu])])

# FG.3-Lung - Eso
cellorder.fg3_lung = 
  rownames(src.endoderm.fg3.ext_v1.1@meta.data[
    src.endoderm.fg3.ext_v1.1$cluster.v06.26.re_hc%in%"Lung",])
princurve.fg3_lung = princurve::principal_curve(
  x = fg3.umap.embedding[cellorder.fg3_lung,],  smoother = "smooth.spline")
cellorder.fg3_lung =  names(
  princurve.fg3_lung$lambda[cellorder.fg3_lung][order(princurve.fg3_lung$lambda[cellorder.fg3_lung])])
#---------------------------

#----------------------------------------------
plot(c(1:length(cellorder.fg3_fg3)),
     c(1:length(cellorder.fg3_fg3)),
     col = colors.time[src.endoderm.fg3.ext_v1.1[,cellorder.fg3_fg3]$Time])
plot(c(1:length(cellorder.fg3_fg3lu)),
     c(1:length(cellorder.fg3_fg3lu)),
     col = colors.time[src.endoderm.fg3.ext_v1.1[,cellorder.fg3_fg3lu]$Time])
plot(c(1:length(cellorder.fg3_lung)),
     c(1:length(cellorder.fg3_lung)),
     col = colors.time[src.endoderm.fg3.ext_v1.1[,cellorder.fg3_lung]$Time])
plot(c(1:length(cellorder.fg3_pha4)),
     c(1:length(cellorder.fg3_pha4)),
     col = colors.time[src.endoderm.fg3.ext_v1.1[,cellorder.fg3_pha4]$Time])
#----------------------------------------------


cellorder.fg3_fg3lu_FPI = cellorder.fg3_fg3lu # 1:3
cellorder.fg3_lung_FPI = cellorder.fg3_lung # 1:3
cellorder.fg3_fg3_FPI = cellorder.fg3_fg3
cellorder.fg3_pha4_FPI = cellorder.fg3_pha4 # 1:3

pdf("figure.v08.07/organ_development_re/fg3_heatmap_marker.re.pdf",6,9)
cellorder.fg3 = c(cellorder.fg3_fg3_PI, cellorder.fg3_fg3lu_FPI,
                  cellorder.fg3_lung_FPI, cellorder.fg3_pha4_FPI)
# cellorder.fg3 = colnames(src.endoderm.fg3.ext_v1.1)
# selectgene = markergene.fg3.ext_v1.1
selectgene = src.endoderm.fg3.ext_v1.1.filtergene_fin

#-------------------
fg3.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.fg3.ext_v1.1@assays$RNA@data[selectgene, cellorder.fg3]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.fg3.ext_v1.1$Time[cellorder.fg3], colors.time),
              MyName2Col(src.endoderm.fg3.ext_v1.1$cluster.v06.26.re_hc[cellorder.fg3], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.fg3.ext_v1.1$cluster.extract.v1.1.re[cellorder.fg3], cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            #return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.fg3.ext_v1.1@meta.data[cell_fg3_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()


fg3.ext_v1.1.re.rowtree.1 = as.dendrogram(fg3.ext_v1.1.re.rowtree)
src.endoderm.fg3.ext_v1.1.filtergene_fin = c(
  labels(fg3.ext_v1.1.re.rowtree.1[[2]][[2]]),
  labels(fg3.ext_v1.1.re.rowtree.1[[2]][[1]]),
  labels(fg3.ext_v1.1.re.rowtree.1[[1]][[1]]),
  labels(fg3.ext_v1.1.re.rowtree.1[[1]][[2]][[1]][[1]]),
  labels(fg3.ext_v1.1.re.rowtree.1[[1]][[2]][[1]][[2]][[2]]),
  labels(fg3.ext_v1.1.re.rowtree.1[[1]][[2]][[2]][[2]]),
  labels(fg3.ext_v1.1.re.rowtree.1[[1]][[2]][[1]][[2]][[1]]),
  labels(fg3.ext_v1.1.re.rowtree.1[[1]][[2]][[2]][[1]])
)

names(src.endoderm.fg3.ext_v1.1.filtergene_fin) = c(
  rep(4, length(labels(fg3.ext_v1.1.re.rowtree.1[[2]][[2]]))),
  rep(3, length(labels(fg3.ext_v1.1.re.rowtree.1[[2]][[1]]))),

  rep(9, length(c(labels(fg3.ext_v1.1.re.rowtree.1[[1]][[1]]),
                  labels(fg3.ext_v1.1.re.rowtree.1[[1]][[2]][[1]][[1]]),
                  labels(fg3.ext_v1.1.re.rowtree.1[[1]][[2]][[1]][[2]][[2]]),
                  labels(fg3.ext_v1.1.re.rowtree.1[[1]][[2]][[2]][[2]])))),
  rep(7, length(c(labels(fg3.ext_v1.1.re.rowtree.1[[1]][[2]][[1]][[2]][[1]]),
                  labels(fg3.ext_v1.1.re.rowtree.1[[1]][[2]][[2]][[1]]))))
)

save(src.endoderm.fg3.ext_v1.1.filtergene_fin,
     #src.endoderm.fg3.ext_v1.1.filtergene,
     src.endoderm.fg3.ext_v1.1.selectgene,
     src.endoderm.fg3.ext.re_v1.1.gene.fin,
     cellorder.fg3_fg3lu_FPI,cellorder.fg3_lung_FPI,
     cellorder.fg3_fg3_FPI,cellorder.fg3_pha4_FPI,
     file= "figure.v08.07/organ_development_re/fg3_heatmap_parameter.Rdata")



#-------------------------------
#  Part2: FG.3
#-------------------------------
src.fg3.integrated_re = src.endoderm.fg3.ext_v1.1
src.fg3.integrated_re$tree = 
  src.fg34.integrated_phlu$tree[colnames(src.fg3.integrated_re)]

src.fg3.integrated_re$cluster.v06.26.re_correct =
  src.fg3.integrated_re$cluster.v06.26.re_hc
src.fg3.integrated_re@meta.data[colnames(src.fg34.integrated_phlu),]$cluster.v06.26.re_correct =
  src.fg34.integrated_phlu$cluster.v06.26.re_correct[colnames(src.fg34.integrated_phlu)]


src.fg3.integrated_re = NormalizeData(src.fg3.integrated_re, scale.factor = 10^5)
src.fg3.integrated_re = ScaleData(src.fg3.integrated_re, rownames(src.fg3.integrated_re))
src.fg3.integrated_re = FindVariableFeatures(src.fg3.integrated_re, assay = "RNA", nfeatures = 2000)
src.fg3.integrated_re_filtergene = 
  Myfilter(as.matrix(src.fg3.integrated_re@assays$RNA@data),
           gene = src.fg3.integrated_re@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)

src.fg3.integrated_re_filtergene_re = src.fg3.integrated_re_filtergene
rm(src.fg3.integrated_re_filtergene)
src.fg3.integrated_re_filtergene = src.fg3.integrated_re_filtergene_re


pdf("figure.v08.07/organ_development_re_v240115/try.fg3_phlu_re.pdf")
src.fg3.integrated_re.rowtree =
  MyHeatmap(as.matrix(src.fg3.integrated_re@assays$RNA@data[
    src.fg3.integrated_re_filtergene,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg3.integrated_re$Time, colors.time),
      MyName2Col(src.fg3.integrated_re$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src.fg3.integrated_re$cluster.v06.26.re_hc, cluster.endoderm.color.v5),
      MyName2Col(src.fg3.integrated_re$tree, c(cluster.endoderm.color.v5, color.temp))),
    ColSideColorsSize = 1.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()
pdf("figure.v08.07/organ_development_re_v240115/try.fg3_phlu_re.tree.pdf",100,20)
plot(src.fg3.integrated_re.rowtree)
dev.off()

tree_src.fg3.integrated_re.rowtree = as.dendrogram(src.fg3.integrated_re.rowtree)
gene_src.fg3.integrated_re.rowtree = setdiff(
  src.fg3.integrated_re_filtergene,
  c(labels(tree_src.fg3.integrated_re.rowtree[[2]][[1]]),
    labels(tree_src.fg3.integrated_re.rowtree[[2]][[2]][[2]][[2]][[1]]),
    labels(tree_src.fg3.integrated_re.rowtree[[2]][[2]][[2]][[2]][[2]][[1]])))


pdf("figure.v08.07/organ_development_re_v240115/try.fg3_phlu_re.pdf")
src.fg3.integrated_re.rowtree.1 =
  MyHeatmap(as.matrix(src.fg3.integrated_re@assays$RNA@data[
    gene_src.fg3.integrated_re.rowtree,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg3.integrated_re$Time, colors.time),
      MyName2Col(src.fg3.integrated_re$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src.fg3.integrated_re$cluster.v06.26.re_hc, cluster.endoderm.color.v5),
      MyName2Col(src.fg3.integrated_re$tree, c(cluster.endoderm.color.v5, color.temp))),
    ColSideColorsSize = 1.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()

pdf("figure.v08.07/organ_development_re_v240115/try.fg3_phlu_re.tree.pdf",100,20)
plot(src.fg3.integrated_re.rowtree.1)
dev.off()

tree_src.fg3.integrated_re.rowtree.1 = as.dendrogram(src.fg3.integrated_re.rowtree.1)
gene_src.fg3.integrated_re.rowtree.1 = setdiff(
  gene_src.fg3.integrated_re.rowtree,
  c(labels(tree_src.fg3.integrated_re.rowtree.1[[1]][[1]]),
    labels(tree_src.fg3.integrated_re.rowtree.1[[2]][[1]])))


pdf("figure.v08.07/organ_development_re_v240115/try.fg3_phlu_re.pdf",9,7)
src.fg3.integrated_re.rowtree.2 =
  MyHeatmap(as.matrix(src.fg3.integrated_re@assays$RNA@data[
    gene_src.fg3.integrated_re.rowtree.1,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg3.integrated_re$Time, colors.time),
      MyName2Col(src.fg3.integrated_re$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src.fg3.integrated_re$cluster.v06.26.re_hc, cluster.endoderm.color.v5),
      MyName2Col(src.fg3.integrated_re$tree, c(cluster.endoderm.color.v5, color.temp))),
    ColSideColorsSize = 1.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()
tree_src.fg3.integrated_re.rowtree.2 = as.dendrogram(src.fg3.integrated_re.rowtree.2)
gene_src.fg3.integrated_re.rowtree.2 = c(
  labels(tree_src.fg3.integrated_re.rowtree.2[[2]][[2]]),
  labels(tree_src.fg3.integrated_re.rowtree.2[[2]][[1]]),
  rev(labels(tree_src.fg3.integrated_re.rowtree.2[[1]])))
names(gene_src.fg3.integrated_re.rowtree.2) = c(
  rep(4, length(labels(tree_src.fg3.integrated_re.rowtree.2[[2]][[2]]))),
  rep(3, length(labels(tree_src.fg3.integrated_re.rowtree.2[[2]][[1]]))),
  rep(7, length(labels(tree_src.fg3.integrated_re.rowtree.2[[1]]))))


pdf("figure.v08.07/organ_development_re_v240115/try.fg3_phlu_re.pdf",9,7)
src.fg3.integrated_re@meta.data = 
  src.fg3.integrated_re@meta.data[!src.fg3.integrated_re$Time%in%NA,]
src.fg3.integrated_re.coltree.2 =
  MyHeatmap(as.matrix(src.fg3.integrated_re@assays$RNA@data[
    gene_src.fg3.integrated_re.rowtree.1,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg3.integrated_re$Time, colors.time),
      MyName2Col(src.fg3.integrated_re$Phase, cluster.endoderm.color.v5),
      MyName2Col(src.fg3.integrated_re$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src.fg3.integrated_re$cluster.v06.26.re_hc, cluster.endoderm.color.v5),
      MyName2Col(src.fg3.integrated_re$tree, c(cluster.endoderm.color.v5, color.temp)),
      MyName2Col(src.fg3.integrated_re$tree_re, c(cluster.endoderm.color.v5, color.temp))
    ),
    ColSideColorsSize = 1.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "col",
    graph = T)
dev.off()

tree_src.fg3.integrated_re.coltree.2 = as.dendrogram(src.fg3.integrated_re.coltree.2)
src.fg3.integrated_re$tree_re = NA
src.fg3.integrated_re@meta.data[labels(tree_src.fg3.integrated_re.coltree.2[[1]][[1]][[1]]),]$tree_re = 1
src.fg3.integrated_re@meta.data[labels(tree_src.fg3.integrated_re.coltree.2[[1]][[1]][[2]]),]$tree_re = 2
src.fg3.integrated_re@meta.data[labels(tree_src.fg3.integrated_re.coltree.2[[1]][[2]][[1]]),]$tree_re = 3

src.fg3.integrated_re@meta.data[labels(tree_src.fg3.integrated_re.coltree.2[[1]][[2]][[2]][[1]]),]$tree_re = 4
src.fg3.integrated_re@meta.data[labels(tree_src.fg3.integrated_re.coltree.2[[1]][[2]][[2]][[2]][[1]]),]$tree_re = 9
src.fg3.integrated_re@meta.data[labels(tree_src.fg3.integrated_re.coltree.2[[1]][[2]][[2]][[2]][[2]]),]$tree_re = 10

src.fg3.integrated_re@meta.data[labels(tree_src.fg3.integrated_re.coltree.2[[2]][[1]][[1]]),]$tree_re = 5
src.fg3.integrated_re@meta.data[labels(tree_src.fg3.integrated_re.coltree.2[[2]][[1]][[2]]),]$tree_re = 6
src.fg3.integrated_re@meta.data[labels(tree_src.fg3.integrated_re.coltree.2[[2]][[2]][[1]]),]$tree_re = 7
src.fg3.integrated_re@meta.data[labels(tree_src.fg3.integrated_re.coltree.2[[2]][[2]][[2]]),]$tree_re = 8

src.fg3.integrated_re$cluster.v06.26.re_correct_re = src.fg3.integrated_re$tree_re
src.fg3.integrated_re@meta.data[src.fg3.integrated_re$tree_re%in%c(4,10),]$cluster.v06.26.re_correct_re = "FG.3"
src.fg3.integrated_re@meta.data[src.fg3.integrated_re$tree_re%in%c(3,9),]$cluster.v06.26.re_correct_re = "Pharynx.organ.4"
src.fg3.integrated_re@meta.data[src.fg3.integrated_re$tree_re%in%c(1,2),]$cluster.v06.26.re_correct_re = "Pharynx.organ.5"
src.fg3.integrated_re@meta.data[src.fg3.integrated_re$tree_re%in%c(5,6,7,8),]$cluster.v06.26.re_correct_re = "Lung"


src.fg3.integrated_re = RunPCA(src.fg3.integrated_re, 
                               features = unique(c(# gene_src.fg34.integrated_phlu.rowtree.2,
                                 gene_src.fg3.integrated_re.rowtree.2)))
src.fg3.integrated_re = RunUMAP(src.fg3.integrated_re, reduction = "pca", assay = "RNA",
                                n.neighbors = 150,
                                dims = c(1:25), n.components = 3)


src.fg3.integrated_re@reductions$umap_monocle3 = 
  src.fg3.integrated_re@reductions$umap
src.fg3.integrated_re@reductions$umap_monocle3@cell.embeddings = 
  umap_embedding_fg3_mon3[colnames(src.fg3.integrated_re),]
colnames(src.fg3.integrated_re@reductions$umap_monocle3@cell.embeddings) = c("UMAP_1", "UMAP_2") #, "UMAP_3")

pdf("figure.v08.07/organ_development_re_v240115/try.fg3.integrated_re_umap.pdf",9,9)
#-----------------------
DimPlot(src.fg3.integrated_re, group.by = "cluster.v06.26.re_correct_re", pt.size = 1.5,
        dims = c(1,2), reduction = 'umap_monocle3', cols = cluster.endoderm.color.v5) +
  theme_bw() + p_add + ggtitle('')+ theme(legend.position = "none")
DimPlot(src.fg3.integrated_re, group.by = "Time",pt.size = 1.5,
        dims = c(1,2), reduction = 'umap_monocle3', cols = colors.time) +
  theme_bw() + p_add + ggtitle('')+ theme(legend.position = "none")
DimPlot(src.fg3.integrated_re, group.by = "Phase",pt.size = 1.5,
        dims = c(1,2), reduction = 'umap_monocle3') +
  theme_bw() + p_add + ggtitle('')+ theme(legend.position = "none")
DimPlot(src.fg3.integrated_re, group.by = "cluster.extract.v1.1",pt.size = 1.5,
        dims = c(1,2), reduction = 'umap_monocle3', cols = cluster.endoderm.color.v5) +
  theme_bw() + p_add + ggtitle('')+ theme(legend.position = "none")
#-----------------------
dev.off()


pdf("figure.v08.07/organ_development_re_v240115/try.fg3.integrated_re_pca_umap.pdf",9,9)
for(red in c("pca","umap")){
  if(red == "pca"){
    i_list = c(1:5)
    j_list = c(1:5)
  }else{
    i_list = c(1:2)
    j_list = c(1:2)
  }
  
  
  for(i in i_list){
    for(j in j_list){
      p1 = DimPlot(src.fg3.integrated_re, group.by = "Time", dims = c(i,j),
                   cols = colors.time, reduction = red) + theme_bw() + p_add_leg 
      p2 = DimPlot(src.fg3.integrated_re, group.by = "tree", dims = c(i,j),
                   cols = color.temp, reduction = red)+ theme_bw() + p_add_leg 
      p3 = DimPlot(src.fg3.integrated_re, group.by = "cluster.v06.26.re_hc", dims = c(i,j),
                   cols = cluster.endoderm.color.v5, reduction = red)+ theme_bw() + p_add_leg 
      print((p1/p3)|p2)
    }}
}
dev.off()


src.fg3.integrated_re@reductions$umap_monocle3_3d = 
  src.fg3.integrated_re@reductions$umap
src.fg3.integrated_re@reductions$umap_monocle3_3d@cell.embeddings = 
  umap_embedding_fg3_mon3_3d[colnames(src.fg3.integrated_re),]
colnames(src.fg3.integrated_re@reductions$umap_monocle3_3d@cell.embeddings) = c("UMAP_1", "UMAP_2", "UMAP_3")

data_src.fg3.integrated_re = cbind(
  src.fg3.integrated_re@meta.data[colnames(src.fg3.integrated_re),],
  src.fg3.integrated_re@reductions$umap_monocle3_3d@cell.embeddings[
    colnames(src.fg3.integrated_re), c(1:3)])
colnames(data_src.fg3.integrated_re) = c(
  colnames(data_src.fg3.integrated_re)[1:(length(colnames(data_src.fg3.integrated_re))-3)],
  'Coord_1','Coord_2','Coord_3')
p.3d=plot_ly(x=data_src.fg3.integrated_re[,"Coord_1"], 
             y=-data_src.fg3.integrated_re[,"Coord_3"],
             z=data_src.fg3.integrated_re[,"Coord_2"],
             type = "scatter3d", mode="markers",
             color = data_src.fg3.integrated_re$cluster.v06.26.re_correct_re,
             colors= cluster.endoderm.color.v5[unique(data_src.fg3.integrated_re$cluster.v06.26.re_correct_re)],
             # color = data_src.fg3.integrated_re$Time,
             # colors= colors.time[unique(data_src.fg3.integrated_re$Time)],
             xlab = "", ylab = "", zlab = "",
             box=T, axes = T, size=8, alpha = 1)
p.3d=p.3d%>%
  layout(scene = list(xaxis = list(title = 'Coord_1', autorange = T, showgrid = F, zeroline = T, showline = T, autotick = F, ticks = '', showticklabels = F),
                      yaxis = list(title = 'Coord_2', autorange = T, showgrid = F, zeroline = T, showline = T, autotick = F, ticks = '', showticklabels = F),
                      zaxis = list(title = 'Coord_3', autorange = T, showgrid = F, zeroline = T, showline = T, autotick = F, ticks = '', showticklabels = F),
                      aspectmode = "manual", aspectratio = list(x=1, y=1, z=1),
                      camera = list(eye = list(x = -1, y = 1, z = 0.5))))
p.3d

# htmlwidgets::saveWidget(as_widget(p.3d), "figure.v08.07/organ_development_re_v240115/html/time_src.fg3.integrated_re.html")
htmlwidgets::saveWidget(as_widget(p.3d), "figure.v08.07/organ_development_re_v240115/html/typee_src.fg3.integrated_re.html")


src.fg3.integrated_re = SetIdent(src.fg3.integrated_re,
                                 value = src.fg3.integrated_re$cluster.v06.26.re_correct_re)
src.fg3.integrated_re_marker = FindAllMarkers(src.fg3.integrated_re) 
src.fg3.integrated_re_marker$pct.ratio = 
  src.fg3.integrated_re_marker$pct.1 / src.fg3.integrated_re_marker$pct.2
src.fg3.integrated_re_markergene = unique(
  src.fg3.integrated_re_marker[src.fg3.integrated_re_marker$pct.ratio>2,]$gene)



cell_src.fg3.integrated_re = rownames(src.fg3.integrated_re@meta.data[src.fg3.integrated_re$cluster.v06.26.re_correct_re%in%c("FG.3","Lung"),])
pcurve_src.fg3.integrated_re = 
  princurve::principal_curve(
    x = src.fg3.integrated_re@reductions$umap@cell.embeddings[cell_src.fg3.integrated_re, c(1,2,3)], smoother = "smooth.spline")
src.fg3.integrated_re$lambda = NA
src.fg3.integrated_re@meta.data[cell_src.fg3.integrated_re,]$lambda = norm_range(pcurve_src.fg3.integrated_re$lambda[cell_src.fg3.integrated_re])

cell_src.fg3.integrated_re.rowtree.2 = c(
  intersect(cell_src.fg3.integrated_re,
            rownames(src.fg3.integrated_re@meta.data[src.fg3.integrated_re$cluster.v06.26.re_correct_re%in%c("FG.3"),])),
  cell_src.fg34.integrated_phlu_re.rowtree)


pdf("figure.v08.07/organ_development_re_v240115/try.fg3_phlu_re.pdf",9,7)
src.fg3.integrated_re.rowtree.3 =
  MyHeatmap(as.matrix(src.fg3.integrated_re@assays$RNA@data[
    unique(c(gene_src.fg3.integrated_re.rowtree.2,
             src.fg3.integrated_re_marker[
               src.fg3.integrated_re_marker$cluster%in%"FG.3"&
                 src.fg3.integrated_re_marker$pct.ratio>1.5,]$gene)),
    cell_src.fg3.integrated_re.rowtree.2]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg3.integrated_re$Time[cell_src.fg3.integrated_re.rowtree.2], 
                 colors.time),
      MyName2Col(src.fg3.integrated_re$cluster.extract.v1.1[cell_src.fg3.integrated_re.rowtree.2],
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg3.integrated_re$cluster.v06.26.re_correct_re[cell_src.fg3.integrated_re.rowtree.2], 
                 cluster.endoderm.color.v5)),
    #RowSideColors = t(cbind(
    #  MyName2Col(names(gene_src.fg3.integrated_re.rowtree.2),
    #             colors.geneset))),
    ColSideColorsSize = 3,
    RowSideColorsSize = 1.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    # Colv = "none",
    return.tree = "row",
    graph = T)
dev.off()

tree_src.fg3.integrated_re.rowtree.3 = as.dendrogram(src.fg3.integrated_re.rowtree.3)
gene_src.fg3.integrated_re.rowtree.3 = c(
  labels(tree_src.fg3.integrated_re.rowtree.3[[2]][[1]][[1]]),
  labels(tree_src.fg3.integrated_re.rowtree.3[[2]][[1]][[2]]),
  
  labels(tree_src.fg3.integrated_re.rowtree.3[[2]][[2]][[2]][[2]][[2]]),
  
  labels(tree_src.fg3.integrated_re.rowtree.3[[2]][[2]][[1]]),
  
  labels(tree_src.fg3.integrated_re.rowtree.3[[1]][[1]]),
  labels(tree_src.fg3.integrated_re.rowtree.3[[1]][[2]][[2]][[1]]),
  labels(tree_src.fg3.integrated_re.rowtree.3[[1]][[2]][[2]][[2]]),
  labels(tree_src.fg3.integrated_re.rowtree.3[[1]][[2]][[1]])
)

pdf("figure.v08.07/organ_development_re_v240115/try.fg3_phlu_re.pdf",9,7)
src.fg3.integrated_re.rowtree.4 =
  MyHeatmap(as.matrix(src.fg3.integrated_re@assays$RNA@data[
    gene_src.fg3.integrated_re.rowtree.3,
    cell_src.fg3.integrated_re.rowtree.2]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg3.integrated_re$Time[cell_src.fg3.integrated_re.rowtree.2], 
                 colors.time),
      MyName2Col(src.fg3.integrated_re$cluster.extract.v1.1[cell_src.fg3.integrated_re.rowtree.2],
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg3.integrated_re$cluster.v06.26.re_correct_re[cell_src.fg3.integrated_re.rowtree.2], 
                 cluster.endoderm.color.v5)),
    #RowSideColors = t(cbind(
    #  MyName2Col(names(gene_src.fg3.integrated_re.rowtree.2),
    #             colors.geneset))),
    ColSideColorsSize = 3,
    RowSideColorsSize = 1.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    Colv = "none",Rowv = "none",
    # return.tree = "row",
    graph = T)
dev.off()



save(src.fg3.integrated_re, file = "figure.v08.07/organ_development_re_v240115/src.fg3.integrated_re.Rdata")
save(gene_src.fg3.integrated_re.rowtree.2, 
     gene_src.fg34.integrated_phlu.rowtree.2,
     file = "figure.v08.07/organ_development_re_v240115/src.fg3.integrated_re_parameter.Rdata")


src.fg3.integrated@reductions$umap_monocle3 = 
  src.fg3.integrated_re@reductions$umap_monocle3
DimPlot(src.fg3.integrated_re, reduction = "umap_monocle3")


