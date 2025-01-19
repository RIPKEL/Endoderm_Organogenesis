#----------------
#>>  Part1: FG.5
#----------------
src.endoderm.fg5.ext_v1.1 = FindVariableFeatures(src.endoderm.fg5.ext_v1.1, nfeatures = 2000)
src.endoderm.fg5.ext_v1.1 = ScaleData(src.endoderm.fg5.ext_v1.1, split.by = "batch_phase",
                                      features = rownames(src.endoderm.fg5.ext_v1.1))
src.endoderm.fg5.ext_v1.1.filtergene = 
  Myfilter(as.matrix(src.endoderm.fg5.ext_v1.1@assays$RNA@data),
           gene = src.endoderm.fg5.ext_v1.1@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
#rm(src.endoderm.fg5.ext_v1.1.filtergene)
src.endoderm.fg5.ext.v1.1.selectgene = unique(c(
  select_gene_fg5,src.endoderm.fg5.ext_v1.1.filtergene))

#=============================
# Row
pdf("figure.v08.07/try.pdf")
#-----------------
fg5.ext.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.fg5.ext_v1.1@assays$RNA@data[
    src.endoderm.fg5.ext.v1.1.selectgene,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.endoderm.fg5.ext_v1.1$Time,colors.time),
      MyName2Col(src.endoderm.fg5.ext_v1.1$cluster.extract.v1.1,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.fg5.ext_v1.1$cluster.v06.26.re,cluster.endoderm.color.v5)
    ),
    ColSideColorsSize = 1.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
#-------------------
dev.off()

pdf("figure.v08.07/try.tree.pdf",100,20)
plot(fg5.ext.re.rowtree_1)
dev.off()

fg5.ext.re.rowtree_1 = as.dendrogram(fg5.ext.re.rowtree)
src.endoderm.fg5.ext_v1.1.gene.fin = setdiff(
  src.endoderm.fg5.ext.v1.1.selectgene,
  c(labels(fg5.ext.re.rowtree_1[[2]][[1]]),
    labels(fg5.ext.re.rowtree_1[[2]][[2]][[2]][[2]][[1]][[2]])))

#=============================
# Col
pdf("figure.v08.07/try.pdf")
#-----------------
fg5.ext.re.coltree =
  MyHeatmap(as.matrix(src.endoderm.fg5.ext_v1.1@assays$RNA@data[
    src.endoderm.fg5.ext_v1.1.gene.fin,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.endoderm.fg5.ext_v1.1$Time,colors.time),
      MyName2Col(src.endoderm.fg5.ext_v1.1$cluster.extract.v1.1,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.fg5.ext_v1.1$cluster.v06.26.re,cluster.endoderm.color.v5)
    ),
    ColSideColorsSize = 1.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "col",
    graph = T)
#-------------------
dev.off()

fg5.ext.re.coltree_1 = as.dendrogram(fg5.ext.re.coltree)

src.endoderm.fg5.ext_v1.1$cluster.v06.26.re_hc = NA
src.endoderm.fg5.ext_v1.1@meta.data[c(
  labels(fg5.ext.re.coltree_1[[1]][[2]][[2]][[1]])), 
  "cluster.v06.26.re_hc"] = "FG.3"

src.endoderm.fg5.ext_v1.1@meta.data[c(
  labels(fg5.ext.re.coltree_1[[1]][[2]][[2]][[2]][[2]])), 
  "cluster.v06.26.re_hc"] = "FG.3"

src.endoderm.fg5.ext_v1.1@meta.data[setdiff(
  labels(fg5.ext.re.coltree_1[[1]]),
  c(  
    labels(fg5.ext.re.coltree_1[[1]][[2]][[2]][[1]]),
    labels(fg5.ext.re.coltree_1[[1]][[2]][[2]][[2]][[2]])
  )), "cluster.v06.26.re_hc"] = "Pharynx.organ.4"

src.endoderm.fg5.ext_v1.1@meta.data[c(
  labels(fg5.ext.re.coltree_1[[2]][[2]][[2]][[1]])), 
  "cluster.v06.26.re_hc"] = "FG.3-Lung"

src.endoderm.fg5.ext_v1.1@meta.data[setdiff(
  labels(fg5.ext.re.coltree_1[[2]]),
  c(  
    labels(fg5.ext.re.coltree_1[[2]][[2]][[2]][[1]])
  )), "cluster.v06.26.re_hc"] = "Lung"


#--------------------------------------------------------------------------------
# src.endoderm.fg5.ext.re_v1.1.gene.fin =
#   rownames(src.endoderm.fg5.ext.re@reductions$pca@feature.loadings)

src.endoderm.fg5.ext_v1.1 = 
  seurat_mnn(seurat_name = "src.endoderm.fg5.ext_v1.1", dims = 1:30, 
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.fg5.ext_v1.1.gene.fin)
src.endoderm.fg5.ext_v1.1 = 
  seurat_int(seurat_name = "src.endoderm.fg5.ext_v1.1", dims = 1:30,
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.fg5.ext_v1.1.gene.fin)


#---------------
src.endoderm.fg5.ext_v1.1 = 
  FindNeighbors(src.endoderm.fg5.ext_v1.1, dims=1:30, 
                reduction = "pca_integrated", assay = "integrated")
# src.endoderm.fg5.ext_v1.1 = 
#   FindNeighbors(src.endoderm.fg5.ext_v1.1, dims=1:30, 
#                 reduction = "pca", assay = "RNA")
src.endoderm.fg5.ext_v1.1 = 
  FindClusters(src.endoderm.fg5.ext_v1.1, 
               resolution = 2, graph.name = "RNA_snn",)
#---------------

DimPlot(src.endoderm.fg5.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.2', reduction = "umap_integrated")

DimPlot(src.endoderm.fg5.ext_v1.1, group.by = 'Time',
        reduction = "umap_mnn", cols = colors.time)
DimPlot(src.endoderm.fg5.ext_v1.1, group.by = 'Time', 
        reduction = "umap_integrated", cols = colors.time)
DimPlot(src.endoderm.fg5.ext_v1.1, group.by = 'batch', 
        reduction = "umap_integrated") #, cols = colors.time)

DimPlot(src.endoderm.fg5.ext_v1.1, group.by = 'cluster.extract.v1.1', 
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.fg5.ext_v1.1, group.by = 'cluster.extract.v1.1', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.fg5.ext_v1.1, group.by = 'cluster.v06.26.re', 
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.fg5.ext_v1.1, group.by = 'cluster.v06.26.re', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)

src.endoderm.fg5.ext_v1.1@reductions$pca_integrated@cell.embeddings = 
  src.endoderm.fg5.ext_v1.1@reductions$pca_integrated@cell.embeddings[colnames(src.endoderm.fg5.ext_v1.1),]

src.endoderm.fg5.ext_v1.1$cluster.v06.26.re = NA
src.endoderm.fg5.ext_v1.1@meta.data[
  intersect(colnames(src.endoderm.fg5.ext_v1.1),colnames(src.endoderm.fg5.ext.re)),]$cluster.v06.26.re = 
  src.endoderm.fg5.ext.re@meta.data[
    intersect(colnames(src.endoderm.fg5.ext_v1.1),colnames(src.endoderm.fg5.ext.re)),]$cluster.v06.26
a = FNN::knn(
  src.endoderm.fg5.ext_v1.1@reductions$pca_integrated@cell.embeddings[
    rownames(src.endoderm.fg5.ext_v1.1@meta.data[!src.endoderm.fg5.ext_v1.1$cluster.v06.26.re%in%NA,]),],
  src.endoderm.fg5.ext_v1.1@reductions$pca_integrated@cell.embeddings[
    rownames(src.endoderm.fg5.ext_v1.1@meta.data[src.endoderm.fg5.ext_v1.1$cluster.v06.26.re%in%NA,]),],
  src.endoderm.fg5.ext_v1.1@meta.data[!src.endoderm.fg5.ext_v1.1$cluster.v06.26.re%in%NA,]$cluster.v06.26.re, k = 10)
src.endoderm.fg5.ext_v1.1@meta.data[
  src.endoderm.fg5.ext_v1.1$cluster.v06.26.re%in%NA,]$cluster.v06.26.re = as.character(a)

src.endoderm.fg5.ext_v1.1@meta.data[
  src.endoderm.fg5.ext_v1.1$RNA_snn_res.2%in%c(0,11,14),]$cluster.v06.26.re = "Pharynx.organ.3"

src.endoderm.fg5.ext_v1.1@meta.data[
  src.endoderm.fg5.ext_v1.1$Time%in%c("ss24","ss27")&
    src.endoderm.fg5.ext_v1.1$cluster.revise.re.v1.30.re.fin%in%c("Pharynx.organ.3"),]$cluster.v06.26.re = "Pharynx.organ.3"
src.endoderm.fg5.ext_v1.1@meta.data[
  src.endoderm.fg5.ext_v1.1$Time%in%c("ss24","ss27")&
    src.endoderm.fg5.ext_v1.1$cluster.revise.re.v1.30.re.fin%in%c("Thyroid"),]$cluster.v06.26.re = "Thyroid"
  
DimPlot(src.endoderm.fg5.ext_v1.1, group.by = 'cluster.v06.26.re', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)


# FDL
src.endoderm.fg5.ext_v1.1 = 
  seurat_fdl("src.endoderm.fg5.ext_v1.1", "pca","RNA")
src.endoderm.fg5.ext_v1.1 = 
  seurat_fdl("src.endoderm.fg5.ext_v1.1", "mnn","RNA")
src.endoderm.fg5.ext_v1.1 = 
  seurat_fdl("src.endoderm.fg5.ext_v1.1", "pca_integrated","integrated")


# Graph
pdf("figure.v08.07/organ_development_re/fg5_summary_graph.re_snn_PI.pdf",12,7)
# SNN
for(red in c('pca',"mnn","pca_integrated")){
  
  assay = ifelse(red=="mnn","RNA","integrated")
  
  graph.fg5.ext_v1.1 = 
    seurat_GCN(seurat_name = "src.endoderm.fg5.ext_v1.1",
               src.endoderm.fg5.ext_v1.1.gene.fin,
               reduction = red, assay =  assay)
  graph.fg5_cluster = cluster_walktrap(graph.fg5.ext_v1.1, steps = 1)
  
  for(color in c("cluster.extract.v1.1","cluster.v06.26.re",
                 "Time")){
    for(edge.color in c(NA,"lightgray")){
      for(edge.color in c(NA,"lightgray")){
        set.seed(1)
        plot.igraph(graph.fg5.ext_v1.1,
                    layout = layout_with_fr(graph.fg5.ext_v1.1),
                    edge.color = edge.color,
                    vertex.size=2.5, vertex.label=NA,
                    vertex.label.color = "black", 
                    vertex.frame.color = NA, vertex.frame.width = 0.5,
                    #vertex.color =  colors.time[src.endoderm.fg5.ext_v1.1[["Time"]][graph.fg5_cluster$names,]],
                    #vertex.color =  cluster.endoderm.color.v5[src.endoderm.fg5.ext_v1.1[["cluster.v06.26.re"]][graph.fg5_cluster$names,]],
                    vertex.color =  colors_bg[src.endoderm.fg5.ext_v1.1[[color]][graph.fg5_cluster$names,]])
      }}}
  
}
dev.off()

pdf("figure.v08.07/organ_development_re/fg5_summary_graph.re_other.pdf",12,7)
# Reduction
for(red in c("umap_integrated","umap_mnn",
             "fdl_pca",'fdl_mnn',"fdl_pca_integrated")){
  p1 = DimPlot(src.endoderm.fg5.ext_v1.1,group.by = 'cluster.v06.26.re', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p2 = DimPlot(src.endoderm.fg5.ext_v1.1,group.by = 'cluster.extract.v1.1', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p3 = DimPlot(src.endoderm.fg5.ext_v1.1,group.by = 'Time', 
               reduction = red , cols = colors.time,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  # p4 = DimPlot(src.endoderm.fg5.ext_v1.1,group.by = 'cluster.v06.26.re_hc', 
  #              reduction = red, cols = cluster.endoderm.color.v5,
  #              pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  print(p1); print(p2); print(p3);# print(p4)
}
dev.off()


pdf("figure.v08.07/organ_development_re/fg5_summary_proportion.re.pdf",12,7)
#-----------------------------
cell_type = c("FG.5","Pharynx.organ.3","Thyroid")
meta_fg5 = src.endoderm.fg5.ext_v1.1@meta.data
meta_fg5$Time = factor(meta_fg5$Time, levels = names(colors.time))
meta_fg5$cluster.extract.v1.1 = factor(meta_fg5$cluster.extract.v1.1, levels = c("FG.5"))
meta_fg5$cluster.v06.26.re = factor(meta_fg5$cluster.v06.26.re, levels = rev(cell_type))

meta_fg5_B0 = meta_fg5
meta_fg5_B1 = meta_fg5[meta_fg5$batch%in%1&!meta_fg5$Time%in%"ss9",]
meta_fg5_B2 = meta_fg5[meta_fg5$batch%in%2,]

for(data in c("meta_fg5_B0","meta_fg5_B1","meta_fg5_B2")){
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
                                 group = cluster.v06.26.re, fill = cluster.v06.26.re),
                   stat = "count", position = 'fill') + 
          xlab("Time") + ylab("Cluster proportion of
 cell type")+ ggtitle(data))
  
  #  print(p + 
  #          geom_bar(data = get(data), 
  #                   mapping = aes(x = Time,
  #                                 group = cluster.v06.26.re_hc, fill = cluster.v06.26.re_hc),
  #                   stat = "count", position = 'fill') + 
  #          xlab("Time") + ylab("Cluster proportion of
  # cell type")+ ggtitle(data))
  
}
#-----------------------------
dev.off()

save(src.endoderm.fg5.ext_v1.1,
     file = "figure.v08.07/organ_development_re/src.endoderm.fg5.ext_v1.1.re.Rdata")

src.endoderm.fg5.ext_v1.1@reductions$umap_integrated@cell.embeddings = cbind(
  src.endoderm.fg5.ext_v1.1@reductions$umap_integrated@cell.embeddings[,2],
  -src.endoderm.fg5.ext_v1.1@reductions$umap_integrated@cell.embeddings[,1],
  src.endoderm.fg5.ext_v1.1@reductions$umap_integrated@cell.embeddings[,3])
colnames(src.endoderm.fg5.ext_v1.1@reductions$umap_integrated@cell.embeddings) = 
  c("UMAP_1","UMAP_2","UMAP_3")


# Marker
#------------------------
src.endoderm.fg5.ext_v1.1 = 
  SetIdent(src.endoderm.fg5.ext_v1.1,
           value = src.endoderm.fg5.ext_v1.1$cluster.v06.26.re)
markergene.fg5.ext_v1.1 = 
  FindAllMarkers(src.endoderm.fg5.ext_v1.1, assay = "RNA")
markergene.fg5.ext_v1.1 = markergene.fg5.ext_v1.1[
  markergene.fg5.ext_v1.1$avg_log2FC>0.2&
    markergene.fg5.ext_v1.1$p_val_adj<0.1, "gene"]

#------------------------
fg5.umap.embedding = src.endoderm.fg5.ext_v1.1[["fdl_pca_integrated"]]@cell.embeddings[,1:2]
fg5.umap.embedding = src.endoderm.fg5.ext_v1.1[["fdl_mnn"]]@cell.embeddings[,1:2]
fg5.umap.embedding = src.endoderm.fg5.ext_v1.1[["umap_mnn"]]@cell.embeddings[,2:3]
fg5.umap.embedding = src.endoderm.fg5.ext_v1.1[["umap_integrated"]]@cell.embeddings[,2:3]

# Cell order :: Princurve
#---------------------------
#  FG.5
cellorder.fg5_fg5 = 
  rownames(src.endoderm.fg5.ext_v1.1@meta.data[
    src.endoderm.fg5.ext_v1.1$cluster.v06.26.re%in%"FG.5",])
princurve.fg5_fg5 = princurve::principal_curve(
  x = fg5.umap.embedding[cellorder.fg5_fg5,],  smoother = "smooth.spline")
cellorder.fg5_fg5 =  names(
  princurve.fg5_fg5$lambda[cellorder.fg5_fg5][order(princurve.fg5_fg5$lambda[cellorder.fg5_fg5])])

# Pha
cellorder.fg5_pha3 = 
  rownames(src.endoderm.fg5.ext_v1.1@meta.data[
    src.endoderm.fg5.ext_v1.1$cluster.v06.26.re%in%"Pharynx.organ.3",])
princurve.fg5_pha3 = princurve::principal_curve(
  x = fg5.umap.embedding[cellorder.fg5_pha3,],  smoother = "smooth.spline")
cellorder.fg5_pha3 =  names(
  princurve.fg5_pha3$lambda[cellorder.fg5_pha3][order(princurve.fg5_pha3$lambda[cellorder.fg5_pha3])])

# Th
cellorder.fg5_th = 
  rownames(src.endoderm.fg5.ext_v1.1@meta.data[
    src.endoderm.fg5.ext_v1.1$cluster.v06.26.re%in%"Thyroid",])
princurve.fg5_th = princurve::principal_curve(
  x = fg5.umap.embedding[cellorder.fg5_th,],  smoother = "smooth.spline")
cellorder.fg5_th =  names(
  princurve.fg5_th$lambda[cellorder.fg5_th][order(princurve.fg5_th$lambda[cellorder.fg5_th])])
#---------------------------

pdf("figure.v08.07/try.pdf")
plot(c(1:length(cellorder.fg5_fg5)),
     c(1:length(cellorder.fg5_fg5)),
     col = colors.time[src.endoderm.fg5.ext_v1.1[,cellorder.fg5_fg5]$Time])
plot(c(1:length(cellorder.fg5_pha3)),
     c(1:length(cellorder.fg5_pha3)),
     col = colors.time[src.endoderm.fg5.ext_v1.1[,cellorder.fg5_pha3]$Time])
plot(c(1:length(cellorder.fg5_th)),
     c(1:length(cellorder.fg5_th)),
     col = colors.time[src.endoderm.fg5.ext_v1.1[,cellorder.fg5_th]$Time])
dev.off()


cellorder.fg5_fg5_FPI = cellorder.fg5_fg5
cellorder.fg5_pha3_FPI = cellorder.fg5_pha3
cellorder.fg5_th_FPI = cellorder.fg5_th


pdf("figure.v08.07/organ_development_re/fg5_heatmap_marker.re.pdf",6,9)
cellorder.fg5 = c(cellorder.fg5_fg5_FPI,
                  rev(cellorder.fg5_pha3_FPI), 
                  cellorder.fg5_th_FPI)
selectgene = src.endoderm.fg5.ext_v1.1.filtergene_fin
#-------------------
fg5.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.fg5.ext_v1.1@assays$RNA@data[selectgene, cellorder.fg5]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.fg5.ext_v1.1$Time[cellorder.fg5], colors.time),
              MyName2Col(src.endoderm.fg5.ext_v1.1$cluster.v06.26.re[cellorder.fg5], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.fg5.ext_v1.1$cluster.extract.v1.1[cellorder.fg5], cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            #return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.fg5.ext_v1.1@meta.data[cell_fg5_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

fg5.ext_v1.1.re.rowtree.1 = as.dendrogram(fg5.ext_v1.1.re.rowtree)

src.endoderm.fg5.ext_v1.1.filtergene_fin = c(
  labels(fg5.ext_v1.1.re.rowtree.1[[1]][[2]]),
  rev(labels(fg5.ext_v1.1.re.rowtree.1[[2]][[2]][[2]][[1]])),
  labels(fg5.ext_v1.1.re.rowtree.1[[2]][[2]][[2]][[2]]),
  labels(fg5.ext_v1.1.re.rowtree.1[[2]][[2]][[1]]),
  labels(fg5.ext_v1.1.re.rowtree.1[[2]][[1]]))
names(src.endoderm.fg5.ext_v1.1.filtergene_fin) = c(
  rep(4,length(c(
    labels(fg5.ext_v1.1.re.rowtree.1[[1]][[2]]),
    labels(fg5.ext_v1.1.re.rowtree.1[[2]][[2]][[2]][[1]][[2]][[2]])
  ))),
  rep(3,length(setdiff(c(
    rev(labels(fg5.ext_v1.1.re.rowtree.1[[2]][[2]][[2]][[1]])),
    labels(fg5.ext_v1.1.re.rowtree.1[[2]][[2]][[2]][[2]]),
    labels(fg5.ext_v1.1.re.rowtree.1[[2]][[2]][[1]])
  ),c(
    labels(fg5.ext_v1.1.re.rowtree.1[[2]][[2]][[2]][[1]][[2]][[2]])
  )))),
  rep(7,length(c(
    labels(fg5.ext_v1.1.re.rowtree.1[[2]][[1]])
  )))
)



save(
  # Gene 
  src.endoderm.fg5.ext.v1.1.selectgene,
  src.endoderm.fg5.ext_v1.1.gene.fin,
  src.endoderm.fg5.ext_v1.1.filtergene_fin,
  # Cell type
  cellorder.fg5_fg5_FPI, cellorder.fg5_pha3_FPI, cellorder.fg5_th_FPI,  
  file= "figure.v08.07/organ_development_re/fg5_heatmap_parameter.Rdata")




#----------------
#>>  Part2: FG5, remake
#----------------
# src.fg5.integrated = src.endoderm.fg5.ext_v1.1

src = src.fg5.integrated
src = NormalizeData(src, assay = "RNA", scale.factor = 10^5)
src = FindVariableFeatures(src, nfeatures = 2000)
src = ScaleData(src, features = rownames(src), assay = "RNA")

src.filtergene = 
  Myfilter(as.matrix(src@assays$RNA@data),
           gene = src@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)

src = SetIdent(src, value = src$cluster.v06.26.re)
marker_src = FindAllMarkers(src)
marker_src$pct.ratio = marker_src$pct.1 / marker_src$pct.2
markergene_src = unique(marker_src[(marker_src$pct.ratio)>2 &
                                     (marker_src$p_val_adj)<0.05,]$gene)


coverage_rate(unique(c(src.filtergene, markergene_src)), 
              gene_src.fg5.tracing.rowtree.3)
coverage_rate(unique(c(src.filtergene, markergene_src)), 
              src.endoderm.fg5.ext_v1.1.filtergene_fin)
coverage_rate(src.endoderm.fg5.ext_v1.1.filtergene_fin, 
              gene_src.fg5.tracing.rowtree.3)

# for(names in names(src.fg5.integrated@reductions)){
#   print(DimPlot(src.fg5.integrated, reduction = names,
#                 group.by = "Time", cols = colors.time)+
#           ggtitle(names))
# }

# coord_src = src[["fdl_mnn"]]@cell.embeddings
# cell_src = rownames(src@meta.data)
# pcurve_src = princurve::principal_curve(
#   x = coord_src[cell_src,], smoother = "smooth.spline")
# src$lambda = pcurve_src$lambda[cell_src]
# src$lambda = norm_range(src$lambda)
# src$order = pcurve_src$ord[cell_src]
# cellorder_src.fg5 = names(src$lambda[order(src$lambda)])
# #cellorder_src.fg5 = names(src$order[order(src$order)])

cellorder_src.fg5 = c(
  cellorder.fg5_fg5_FPI, 
  cellorder.fg5_pha3_FPI, cellorder.fg5_th_FPI)

# Basic::Marker+Filter
#----------------------------
# --- Marker
pdf("figure.v08.07/organ_development_re_v240115/try.src.pdf",9,7)
src.rowtree =
  MyHeatmap(as.matrix(src@assays$RNA@data[markergene_src,]),
            type = "row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            ColSideColors = cbind(
              MyName2Col(src$Time, colors.time),
              MyName2Col(src$cluster.extract.v1.1, cluster.endoderm.color.v5),
              MyName2Col(src$cluster.v06.26.re, cluster.endoderm.color.v5),
              MyName2Col(src$batch, color.temp)),
            ColSideColorsSize = 3,
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            #Rowv = "none",
            return.tree = "row",
            graph = T)
dev.off()

tree_src.rowtree = as.dendrogram(src.rowtree)
gene_src.rowtree = c(
  labels(tree_src.rowtree[[2]][[1]]),
  labels(tree_src.rowtree[[1]][[2]][[2]]),
  labels(tree_src.rowtree[[2]][[2]][[1]]),
  labels(tree_src.rowtree[[2]][[2]][[2]][[2]]))
coverage_rate(gene_src.rowtree, gene_src.fg5.tracing.rowtree.3)


# --- Filter
pdf("figure.v08.07/organ_development_re_v240115/try.src.pdf",9,7)
src.rowtree.1 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    src.filtergene,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time, colors.time),
      MyName2Col(src$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re, cluster.endoderm.color.v5)),
    ColSideColorsSize = 3,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()

tree_src.rowtree.1 = as.dendrogram(src.rowtree.1)
gene_src.rowtree.1 = c(
  labels(tree_src.rowtree.1[[1]]))
coverage_rate(gene_src.rowtree.1, gene_src.fg5.tracing.rowtree.3)
#----------------------------

# Intersect
#----------------------------
pdf("figure.v08.07/organ_development_re_v240115/try.src.pdf",9,7)
src.coltree.2 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    intersect(unique(c(gene_src.rowtree.1, gene_src.rowtree)), 
              gene_src.fg5.tracing.rowtree.3),
    cellorder_src.fg5]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time[cellorder_src.fg5], 
                 colors.time),
      MyName2Col(src$cluster.extract.v1.1[cellorder_src.fg5],  
                 cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re[cellorder_src.fg5],
                 cluster.endoderm.color.v5)),
    ColSideColorsSize = 3,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    # Rowv = "none", Colv = "none",
    return.tree = "col",
    graph = T)

src.rowtree.2 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    intersect(unique(c(gene_src.rowtree.1, gene_src.rowtree)), 
              gene_src.fg5.tracing.rowtree.3),
    cellorder_src.fg5]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time[cellorder_src.fg5], 
                 colors.time),
      MyName2Col(src$cluster.extract.v1.1[cellorder_src.fg5],  
                 cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re[cellorder_src.fg5],
                 cluster.endoderm.color.v5)),
    ColSideColorsSize = 3,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    # Rowv = "none", Colv = "none",
    return.tree = "row",
    graph = T)
dev.off()

tree_src.coltree.2 = as.dendrogram(src.coltree.2)
tree_src.rowtree.2 = as.dendrogram(src.rowtree.2)
gene_src.rowtree.2 = c(
  labels(tree_src.rowtree.2[[2]][[1]]),
  
  labels(tree_src.rowtree.2[[2]][[2]][[1]][[2]]),
  labels(tree_src.rowtree.2[[2]][[2]][[2]][[2]][[2]]),
  
  setdiff(labels(tree_src.rowtree.2[[1]][[2]]),
          labels(tree_src.rowtree.2[[1]][[2]][[2]][[2]][[2]])),
  # labels(tree_src.rowtree.2[[1]][[2]][[2]][[2]][[2]]),
  labels(tree_src.rowtree.2[[1]][[1]][[2]]),
  labels(tree_src.rowtree.2[[1]][[1]][[1]]))

coverage_rate(gene_src.rowtree.2,
              gene_src.fg5.tracing.rowtree.3)
#----------------------------

# Setdiff
#----------------------------
pdf("figure.v08.07/organ_development_re_v240115/try.src.pdf",9,7)
src.rowtree.3 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    setdiff(unique(c(gene_src.rowtree.1, gene_src.rowtree)),
            gene_src.fg5.tracing.rowtree.3),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time, colors.time),
      MyName2Col(src$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re, cluster.endoderm.color.v5)),
    ColSideColorsSize = 3,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()

tree_src.rowtree.3 = as.dendrogram(src.rowtree.3)
gene_src.rowtree.3 = c(
  labels(tree_src.rowtree.3[[1]]))
#-----------------------------

# Finally! 
#-----------------------------
pdf("figure.v08.07/organ_development_re_v240115/try.src.pdf",9,7)
src.rowtree.5 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    c(gene_src.rowtree.2,
      gene_src.rowtree.3),]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time, colors.time),
      MyName2Col(src$cluster.extract.v1.1, cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re, cluster.endoderm.color.v5)),
    ColSideColorsSize = 3,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
dev.off()

tree_src.rowtree.5 = as.dendrogram(src.rowtree.5)
gene_src.rowtree.5 = c(
  labels(tree_src.rowtree.5[[2]][[1]]),
  
  rev(labels(tree_src.rowtree.5[[1]])),
  
  rev(c(labels(tree_src.rowtree.5[[2]][[2]][[1]]),
        labels(tree_src.rowtree.5[[2]][[2]][[2]][[2]][[2]][[1]]))))

names(gene_src.rowtree.5) = c(
  rep(4, length(labels(tree_src.rowtree.5[[2]][[1]]))),
  rep(3, length(rev(labels(tree_src.rowtree.5[[1]])))),
  rep(7, length(c(labels(tree_src.rowtree.5[[2]][[2]][[1]]),
                  labels(tree_src.rowtree.5[[2]][[2]][[2]][[2]][[2]][[1]])))))
table(names(gene_src.rowtree.5))


pdf("figure.v08.07/organ_development_re_v240115/try.src.pdf",9,7)
# pdf("figure.v08.07/organ_development_re_v240115/try.src.fg5.integrated.pdf",9,7)
src.coltree.5 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    gene_src.rowtree.5,
    cellorder_src.fg5]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time[cellorder_src.fg5], 
                 colors.time),
      MyName2Col(src$tree[cellorder_src.fg5], 
                 color.temp),
      MyName2Col(src$cluster.extract.v1.1[cellorder_src.fg5],  
                 cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re[cellorder_src.fg5],
                 cluster.endoderm.color.v5)),
    ColSideColorsSize = 3,
    RowSideColorsSize = 1,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    # Rowv = "none",
    # Colv = "none",
    return.tree = "col",
    margins = c(10,10), graph = T)
dev.off()

# tree_src.coltree.5 = as.dendrogram(src.coltree.5)
src$tree = NA
src@meta.data[labels(tree_src.coltree.5[[2]][[1]]),]$tree = 1
src@meta.data[labels(tree_src.coltree.5[[2]][[2]]),]$tree = 2
src@meta.data[labels(tree_src.coltree.5[[1]][[1]]),]$tree = 3
src@meta.data[labels(tree_src.coltree.5[[1]][[2]][[1]]),]$tree = 4
src@meta.data[labels(tree_src.coltree.5[[1]][[2]][[2]][[1]]),]$tree = 5
src@meta.data[labels(tree_src.coltree.5[[1]][[2]][[2]][[2]]),]$tree = 6

src$cluster.v06.26.re_correct = src$cluster.v06.26.re
src@meta.data[src$tree%in%c(1,6),]$cluster.v06.26.re_correct = "FG.5"
src@meta.data[src$tree%in%c(2),]$cluster.v06.26.re_correct = "Pharynx.organ.3"
src@meta.data[src$tree%in%c(3,4,5),]$cluster.v06.26.re_correct = "Thyroid"


for(names in names(src.fg5.integrated@reductions)){
  print(DimPlot(src.fg5.integrated, reduction = names,
                group.by = "Time", cols = colors.time)+
          ggtitle(names))
}

coord_src = src[["umap_integrated"]]@cell.embeddings
cell_src = rownames(src@meta.data)

# Pseudo-time :: cluster.v06.26.re_correct
#---------------------
#-- pha3 --
cell_src_pha3 = 
  rownames(src@meta.data[src$cluster.v06.26.re_correct%in%c('Pharynx.organ.3'),])
coord_src_pha3 = src[["umap_integrated"]]@cell.embeddings[cell_src_pha3,]
pcurve_src_pha3 = princurve::principal_curve(x = coord_src_pha3, smoother = "smooth.spline")
src$lambda_pha3 = pcurve_src_pha3$lambda[cell_src_pha3]
src@meta.data[cell_src_pha3,]$lambda_pha3 = 
  pcurve_src_pha3$lambda[cell_src_pha3]
src@meta.data[!src$lambda_pha3%in%NA,]$lambda_pha3 = 
  norm_range(src@meta.data[!src$lambda_pha3%in%NA,]$lambda_pha3)

#-- th --
cell_src_th = 
  rownames(src@meta.data[src$cluster.v06.26.re_correct%in%c("FG.5",'Thyroid'),])
coord_src_th = src[["umap_integrated"]]@cell.embeddings[cell_src_th,]
pcurve_src_th = princurve::principal_curve(x = coord_src_th, smoother = "smooth.spline")
src$lambda_th = pcurve_src_th$lambda[cell_src_th]
src@meta.data[cell_src_th,]$lambda_th = 
  pcurve_src_th$lambda[cell_src_th]
src@meta.data[!src$lambda_th%in%NA,]$lambda_th = 
  norm_range(src@meta.data[!src$lambda_th%in%NA,]$lambda_th)

#-- fg5 --
cell_src_fg5 = 
  rownames(src@meta.data[src$cluster.v06.26.re_correct%in%c("FG.5"),])
coord_src_fg5 = src[["umap_integrated"]]@cell.embeddings[cell_src_fg5,]
pcurve_src_fg5 = princurve::principal_curve(x = coord_src_fg5, smoother = "smooth.spline")
src$lambda_fg5 = pcurve_src_fg5$lambda[cell_src_fg5]
src@meta.data[cell_src_fg5,]$lambda_fg5 = 
  pcurve_src_fg5$lambda[cell_src_fg5]
src@meta.data[!src$lambda_fg5%in%NA,]$lambda_fg5 = 
  norm_range(src@meta.data[!src$lambda_fg5%in%NA,]$lambda_fg5)


cellorder_src.fg5 = c(
  (intersect(names(src$lambda_fg5[order(src$lambda_fg5)]), 
             colnames(src[,src$cluster.v06.26.re_correct%in%"FG.5"]))),
  intersect(names(src$lambda_th[order(src$lambda_th)]), 
            colnames(src[,src$cluster.v06.26.re_correct%in%"Thyroid"])),
  (intersect(names(src$lambda_pha3[order(src$lambda_pha3)]), 
             colnames(src[,src$cluster.v06.26.re_correct%in%"Pharynx.organ.3"]))))
#---------------------


# pdf("figure.v08.07/organ_development_re_v240115/try.src.pdf",9,7)
pdf("figure.v08.07/organ_development_re_v240115/try.src.fg5.integrated.pdf",9,7)
src.rowtree.7 =
  MyHeatmap(as.matrix(src@assays$RNA@data[
    gene_src.rowtree.5,
    cellorder_src.fg5]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src$Time[cellorder_src.fg5], 
                 colors.time),
      MyName2Col(src$cluster.extract.v1.1[cellorder_src.fg5],  
                 cluster.endoderm.color.v5),
      MyName2Col(src$cluster.v06.26.re_correct[cellorder_src.fg5],
                 cluster.endoderm.color.v5)),
    RowSideColors = t(cbind(
      MyName2Col(names(gene_src.rowtree.5),
                 colors.geneset))),
    ColSideColorsSize = 3,
    RowSideColorsSize = 1,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    Rowv = "none",
    Colv = "none",
    # return.tree = "row",
    margins = c(10,10), graph = T)
dev.off()

tf_gene_src.fg5.integrated.rowtree.5 = intersect(
  gene_src.fg5.integrated.rowtree.5, gi[gi$TF%in%T,]$SymbolDedu)
pdf("figure.v08.07/organ_development_re_v240115/try.src.fg5.integrated.fg5.integrated_tf.pdf",9,7)
src.fg5.integrated.rowtree.7 =
  MyHeatmap(as.matrix(src.fg5.integrated@assays$RNA@data[
    tf_gene_src.fg5.integrated.rowtree.5,
    cellorder_src.fg5.integrated.fg5]),
    type = "log.row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.fg5.integrated$Time[cellorder_src.fg5.integrated.fg5], 
                 colors.time),
      MyName2Col(src.fg5.integrated$cluster.extract.v1.1[cellorder_src.fg5.integrated.fg5],  
                 cluster.endoderm.color.v5),
      MyName2Col(src.fg5.integrated$cluster.v06.26.re_correct[cellorder_src.fg5.integrated.fg5],
                 cluster.endoderm.color.v5)),
    # RowSideColors = t(cbind(
    #   MyName2Col(names(gene_src.fg5.integrated.rowtree.5),
    #              colors.geneset))),
    ColSideColorsSize = 3,
    RowSideColorsSize = 1,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    labRow = tf_gene_src.fg5.integrated.rowtree.5,
    Rowv = "none",
    Colv = "none",
    #return.tree = "row",
    margins = c(10,10), graph = T)
dev.off()

coverage_rate(tf_gene_src.rowtree.5, 
              intersect(gene_src.fg5.tracing.rowtree.3, 
                        gi[gi$TF%in%T,]$SymbolDedu))
#-----------------------------

# Save and Rename
#-----------------------------
src.fg5.integrated$tree = src$tree
src.fg5.integrated$cluster.v06.26.re_correct = src$cluster.v06.26.re_correct

for(i in c("markergene_src", "tree_src.rowtree", "gene_src.rowtree",
           "src.filtergene", "tree_src.rowtree.1", "gene_src.rowtree.1",
           "tree_src.coltree.2", "tree_src.rowtree.2", "gene_src.rowtree.2",
           "tree_src.rowtree.3", "gene_src.rowtree.3",
           "tree_src.rowtree.4", "gene_src.rowtree.4",
           "tree_src.rowtree.5", "gene_src.rowtree.5", 
           "tree_src.coltree.5",
           #"tree_src.rowtree.6", "gene_src.rowtree.6",
           "tf_gene_src.rowtree.5", 
           "cellorder_src.fg5")){
  assign(gsub("src","src.fg5.integrated",i), get(i))
}

save(src.fg5.integrated,
     file = "figure.v08.07/organ_development_re_v240115/src.fg5.integrated.Rdata")
save(markergene_src.fg5.integrated, 
     tree_src.fg5.integrated.rowtree, gene_src.fg5.integrated.rowtree,
     src.fg5.integrated.filtergene, 
     tree_src.fg5.integrated.rowtree.1, gene_src.fg5.integrated.rowtree.1,
     tree_src.fg5.integrated.coltree.2, tree_src.fg5.integrated.rowtree.2, gene_src.fg5.integrated.rowtree.2,
     tree_src.fg5.integrated.rowtree.3, gene_src.fg5.integrated.rowtree.3,
     tree_src.fg5.integrated.rowtree.4, gene_src.fg5.integrated.rowtree.4,
     tree_src.fg5.integrated.rowtree.5, gene_src.fg5.integrated.rowtree.5,
     tree_src.fg5.integrated.coltree.5,
     # tree_src.fg5.integrated.rowtree.6, gene_src.fg5.integrated.rowtree.6,
     tf_gene_src.fg5.integrated.rowtree.5,
     cellorder_src.fg5.integrated.fg5,
     file = "figure.v08.07/organ_development_re_v240115/fg5_10x_parmeter.Rdata")

rm(markergene_src, tree_src.rowtree, gene_src.rowtree,
   src.filtergene, tree_src.rowtree.1, gene_src.rowtree.1,
   tree_src.coltree.2, tree_src.rowtree.2, gene_src.rowtree.2,
   tree_src.rowtree.3, gene_src.rowtree.3,
   tree_src.rowtree.4, gene_src.rowtree.4,
   tree_src.rowtree.5, gene_src.rowtree.5,
   tree_src.coltree.5,
   #tree_src.rowtree.6, gene_src.rowtree.6,
   tf_gene_src.rowtree.5,
   cellorder_src.fg5)
#================================================================================




