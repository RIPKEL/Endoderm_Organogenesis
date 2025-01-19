#----------------
#   HG.2
#----------------
src.endoderm.hg2.ext_v1.1 = FindVariableFeatures(src.endoderm.hg2.ext_v1.1, nfeatures = 2000)
src.endoderm.hg2.ext_v1.1 = ScaleData(src.endoderm.hg2.ext_v1.1, split.by = "batch_phase",
                                      features = rownames(src.endoderm.hg2.ext_v1.1))
src.endoderm.hg2.ext_v1.1.filtergene = 
  Myfilter(as.matrix(src.endoderm.hg2.ext_v1.1@assays$RNA@data),
           gene = src.endoderm.hg2.ext_v1.1@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
#rm(src.endoderm.hg2.ext_v1.1.filtergene)
src.endoderm.hg2.ext.v1.1.selectgene = unique(c(
  select_gene_hg2, src.endoderm.hg2.ext_v1.1.filtergene))

#=============================
# Row
pdf("figure.v08.07/try.pdf")
#-----------------
hg2.ext.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.hg2.ext_v1.1@assays$RNA@data[
    src.endoderm.hg2.ext.v1.1.selectgene,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.endoderm.hg2.ext_v1.1$Time,colors.time),
      MyName2Col(src.endoderm.hg2.ext_v1.1$cluster.extract.v1.1,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.hg2.ext_v1.1$cluster.v06.26.re,cluster.endoderm.color.v5)
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
plot(hg2.ext.re.rowtree)
dev.off()

hg2.ext.re.rowtree_1 = as.dendrogram(hg2.ext.re.rowtree)
src.endoderm.hg2.ext_v1.1.gene.fin = setdiff(
  src.endoderm.hg2.ext.v1.1.selectgene,
  c(labels(hg2.ext.re.rowtree_1[[2]][[1]]),
    labels(hg2.ext.re.rowtree_1[[1]][[2]][[2]])))

#--------------------------------------------------------------------------------
src.endoderm.hg2.ext_v1.1 = RunPCA(src.endoderm.hg2.ext_v1.1,
                                   features = src.endoderm.hg2.ext_v1.1.gene.fin)
src.endoderm.hg2.ext_v1.1 = 
  seurat_mnn(seurat_name = "src.endoderm.hg2.ext_v1.1", dims = 1:30, 
             n.neighbors = 100,n.components = 2,
             seurat.selectgene = src.endoderm.hg2.ext_v1.1.gene.fin)
src.endoderm.hg2.ext_v1.1 = 
  seurat_int(seurat_name = "src.endoderm.hg2.ext_v1.1", dims = 1:30,
             n.neighbors = 100,n.components = 2,
             seurat.selectgene = src.endoderm.hg2.ext_v1.1.gene.fin)

#---------------
src.endoderm.hg2.ext_v1.1@reductions$mnn@cell.embeddings = 
  src.endoderm.hg2.ext_v1.1@reductions$mnn@cell.embeddings[colnames(src.endoderm.hg2.ext_v1.1),]
src.endoderm.hg2.ext_v1.1 = 
  FindNeighbors(src.endoderm.hg2.ext_v1.1, dims=1:30, 
                reduction = "mnn", assay = "RNA")
src.endoderm.hg2.ext_v1.1 = 
  FindClusters(src.endoderm.hg2.ext_v1.1, resolution = 2)
src.endoderm.hg2.ext_v1.1 = 
  FindClusters(src.endoderm.hg2.ext_v1.1, resolution = 4)
#---------------

pdf("figure.v08.07/try.pdf")
DimPlot(src.endoderm.hg2.ext_v1.1, label = T, label.size = 4, cols = colors.num, pt.size = 1.2,
        group.by = 'RNA_snn_res.2', reduction = "umap_integrated")
DimPlot(src.endoderm.hg2.ext_v1.1, label = T, label.size = 4, cols = colors.num, pt.size = 1.2,
        group.by = 'RNA_snn_res.4', reduction = "umap_integrated")

DimPlot(src.endoderm.hg2.ext_v1.1, group.by = 'Time',pt.size = 1.2,
        reduction = "umap_mnn", cols = colors.time)
DimPlot(src.endoderm.hg2.ext_v1.1, group.by = 'Time', pt.size = 1.2,
        reduction = "umap_integrated", cols = colors.time)
DimPlot(src.endoderm.hg2.ext_v1.1, group.by = 'batch', 
        reduction = "umap_integrated") #, cols = colors.time)

DimPlot(src.endoderm.hg2.ext_v1.1, group.by = 'cluster.extract.v1.1', 
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.hg2.ext_v1.1, group.by = 'cluster.extract.v1.1', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.hg2.ext_v1.1, group.by = 'cluster.v06.26.re', pt.size = 1.2,
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.hg2.ext_v1.1, group.by = 'cluster.v06.26.re', pt.size = 1.2,
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
dev.off()

src.endoderm.hg2.ext_v1.1@reductions$pca_integrated@cell.embeddings = 
  src.endoderm.hg2.ext_v1.1@reductions$pca_integrated@cell.embeddings[colnames(src.endoderm.hg2.ext_v1.1),]
#--------------------
src.endoderm.hg2.ext_v1.1$cluster.v06.26.re_raw = NA
src.endoderm.hg2.ext_v1.1@meta.data[
  intersect(colnames(src.endoderm.hg2.ext_v1.1),colnames(src.endoderm.hg2.ext.re)),]$cluster.v06.26.re_raw = 
  src.endoderm.hg2.ext.re@meta.data[
    intersect(colnames(src.endoderm.hg2.ext_v1.1),colnames(src.endoderm.hg2.ext.re)),]$cluster.v06.26
# src.endoderm.hg2.ext.re$cluster.v06.26.re = NA
# src.endoderm.hg2.ext.re@meta.data[colnames(src.endoderm.hg2.ext.re.sma.1_lar.1),]$cluster.v06.26.re =
#   src.endoderm.hg2.ext.re.sma.1_lar.1$cluster.v06.26.re
# src.endoderm.hg2.ext.re@meta.data[colnames(src.endoderm.hg2.ext.re.lar.2),]$cluster.v06.26.re =
#   src.endoderm.hg2.ext.re.lar.2$cluster.v06.26
#---------------------
src.endoderm.hg2.ext_v1.1$cluster.v06.26.re = NA

src.endoderm.hg2.ext_v1.1@meta.data[
  src.endoderm.hg2.ext_v1.1$RNA_snn_res.2%in%c(9,11,10,6),]$cluster.v06.26.re = "Large.intestine.1"
src.endoderm.hg2.ext_v1.1@meta.data[
  src.endoderm.hg2.ext_v1.1$RNA_snn_res.2%in%c(12,16,14,5),]$cluster.v06.26.re = "Large.intestine.3"
src.endoderm.hg2.ext_v1.1@meta.data[
  src.endoderm.hg2.ext_v1.1$RNA_snn_res.2%in%c(13,2,8,3,15,0),]$cluster.v06.26.re = "HG.2" 

src.endoderm.hg2.ext_v1.1@meta.data[
  src.endoderm.hg2.ext_v1.1$Time%in%c("ss24","ss27")&
    src.endoderm.hg2.ext_v1.1$cluster.revise.re.v1.30.re.fin%in%c("Large.intestine.3"),]$cluster.v06.26.re = "Large.intestine.3"
src.endoderm.hg2.ext_v1.1@meta.data[
  src.endoderm.hg2.ext_v1.1$Time%in%c("ss24","ss27")&
    src.endoderm.hg2.ext_v1.1$cluster.revise.re.v1.30.re.fin%in%c("Large.intestine.1"),]$cluster.v06.26.re = "Large.intestine.1"
# src.endoderm.hg2.ext_v1.1@meta.data[
#   src.endoderm.hg2.ext_v1.1$Time%in%c("ss24","ss27")&
#     src.endoderm.hg2.ext_v1.1$cluster.revise.re.v1.30.re.fin%in%c("Large.intestine.2"),]$cluster.v06.26.re = "Large.intestine.2"

# src.endoderm.hg2.ext_v1.1.re = src.endoderm.hg2.ext_v1.1 #Save Lar.2
# src.endoderm.hg2.ext_v1.1 = src.endoderm.hg2.ext_v1.1.re[,!src.endoderm.hg2.ext_v1.1$cluster.v06.26.re%in%"Large.intestine.2"]


a = FNN::knn(
  src.endoderm.hg2.ext_v1.1@reductions$pca_integrated@cell.embeddings[
    rownames(src.endoderm.hg2.ext_v1.1@meta.data[!src.endoderm.hg2.ext_v1.1$cluster.v06.26.re%in%NA,]),],
  src.endoderm.hg2.ext_v1.1@reductions$pca_integrated@cell.embeddings[
    rownames(src.endoderm.hg2.ext_v1.1@meta.data[src.endoderm.hg2.ext_v1.1$cluster.v06.26.re%in%NA,]),],
  src.endoderm.hg2.ext_v1.1@meta.data[!src.endoderm.hg2.ext_v1.1$cluster.v06.26.re%in%NA,]$cluster.v06.26.re, k = 10)
src.endoderm.hg2.ext_v1.1@meta.data[
  src.endoderm.hg2.ext_v1.1$cluster.v06.26.re%in%NA,]$cluster.v06.26.re = as.character(a)

DimPlot(src.endoderm.hg2.ext_v1.1, group.by = 'cluster.v06.26.re', pt.size = 1.2,
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)

# FDL
src.endoderm.hg2.ext_v1.1 = 
  seurat_fdl("src.endoderm.hg2.ext_v1.1", "pca","RNA")
src.endoderm.hg2.ext_v1.1 = 
  seurat_fdl("src.endoderm.hg2.ext_v1.1", "mnn","RNA")
src.endoderm.hg2.ext_v1.1 = 
  seurat_fdl("src.endoderm.hg2.ext_v1.1", "pca_integrated","integrated")

# Graph
pdf("figure.v08.07/organ_development_re/hg2_summary_graph.re_snn_PI.pdf",12,7)
# SNN
for(red in c('pca',"mnn","pca_integrated")){
  
  assay = ifelse(red=="mnn","RNA","integrated")
  
  graph.hg2.ext_v1.1 = 
    seurat_GCN(seurat_name = "src.endoderm.hg2.ext_v1.1",
               src.endoderm.hg2.ext_v1.1.gene.fin,
               reduction = red, assay =  assay)
  graph.hg2_cluster = cluster_walktrap(graph.hg2.ext_v1.1, steps = 1)
  
  for(color in c("cluster.extract.v1.1","cluster.v06.26.re",#"cluster.v06.26.re_hc",
                 "Time")){
    for(edge.color in c(NA,"lightgray")){
      for(edge.color in c(NA,"lightgray")){
        set.seed(1)
        plot.igraph(graph.hg2.ext_v1.1,
                    layout = layout_with_fr(graph.hg2.ext_v1.1),
                    edge.color = edge.color,
                    vertex.size=2.5, vertex.label=NA,
                    vertex.label.color = "black", 
                    vertex.frame.color = NA, vertex.frame.width = 0.5,
                    #vertex.color =  colors.time[src.endoderm.hg2.ext_v1.1[["Time"]][graph.hg2_cluster$names,]],
                    #vertex.color =  cluster.endoderm.color.v5[src.endoderm.hg2.ext_v1.1[["cluster.v06.26.re"]][graph.hg2_cluster$names,]],
                    vertex.color =  colors_bg[src.endoderm.hg2.ext_v1.1[[color]][graph.hg2_cluster$names,]])
      }}}
  
}
dev.off()

pdf("figure.v08.07/organ_development_re/hg2_summary_graph.re_other.pdf",12,7)
# Reduction
for(red in c("umap_integrated","umap_mnn",
             "fdl_pca",'fdl_mnn',"fdl_pca_integrated")){
  p1 = DimPlot(src.endoderm.hg2.ext_v1.1,group.by = 'cluster.v06.26.re', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p2 = DimPlot(src.endoderm.hg2.ext_v1.1,group.by = 'cluster.extract.v1.1', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p3 = DimPlot(src.endoderm.hg2.ext_v1.1,group.by = 'Time', 
               reduction = red , cols = colors.time,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  # p4 = DimPlot(src.endoderm.hg2.ext_v1.1,group.by = 'cluster.v06.26.re_hc', 
  #              reduction = red, cols = cluster.endoderm.color.v5,
  #              pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  print(p1); print(p2); print(p3);#  print(p4)
}
dev.off()


pdf("figure.v08.07/organ_development_re/hg2_summary_proportion.re.pdf",12,7)
#-----------------------------
cell_type = c("HG.2","Large.intestine.1","Large.intestine.3")
meta_hg2 = src.endoderm.hg2.ext_v1.1@meta.data
meta_hg2$Time = factor(meta_hg2$Time, levels = names(colors.time))
meta_hg2$cluster.extract.v1.1 = factor(meta_hg2$cluster.extract.v1.1, 
                                       levels = c("HG.2"))
meta_hg2$cluster.v06.26.re = factor(meta_hg2$cluster.v06.26.re, levels = rev(cell_type))

meta_hg2_B0 = meta_hg2
meta_hg2_B1 = meta_hg2[meta_hg2$batch%in%1&!meta_hg2$Time%in%"ss9",]
meta_hg2_B2 = meta_hg2[meta_hg2$batch%in%2,]

for(data in c("meta_hg2_B0","meta_hg2_B1","meta_hg2_B2")){
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
  
  # print(p + 
  #         geom_bar(data = get(data), 
  #                  mapping = aes(x = Time,
  #                                group = cluster.v06.26.re_hc, fill = cluster.v06.26.re_hc),
  #                  stat = "count", position = 'fill') + 
  #         xlab("Time") + ylab("Cluster proportion of
  # cell type")+ ggtitle(data))
  
}
#-----------------------------
dev.off()

save(src.endoderm.hg2.ext_v1.1,
     src.endoderm.hg2.ext_v1.1.re,
     file = "figure.v08.07/organ_development_re/src.endoderm.hg2.ext_v1.1.re.Rdata")


# Marker
#------------------------
src.endoderm.hg2.ext_v1.1 = 
  SetIdent(src.endoderm.hg2.ext_v1.1,
           value = src.endoderm.hg2.ext_v1.1$cluster.v06.26.re)
markergene.hg2.ext_v1.1 = 
  FindAllMarkers(src.endoderm.hg2.ext_v1.1, assay = "RNA")
markergene.hg2.ext_v1.1 = markergene.hg2.ext_v1.1[
  markergene.hg2.ext_v1.1$avg_log2FC>0.2&
    markergene.hg2.ext_v1.1$p_val_adj<0.1, "gene"]

#------------------------
hg2.umap.embedding = src.endoderm.hg2.ext_v1.1[["fdl_pca_integrated"]]@cell.embeddings[,1:3]
hg2.umap.embedding = src.endoderm.hg2.ext_v1.1[["fdl_mnn"]]@cell.embeddings[,1:3]
hg2.umap.embedding = src.endoderm.hg2.ext_v1.1[["umap_mnn"]]@cell.embeddings[,1:3]
hg2.umap.embedding = src.endoderm.hg2.ext_v1.1[["umap_integrated"]]@cell.embeddings[,1:3]

# Cell order :: Princurve
#---------------------------
# HG.2
cellorder.hg2_hg2 = 
  rownames(src.endoderm.hg2.ext_v1.1@meta.data[
    src.endoderm.hg2.ext_v1.1$cluster.v06.26.re%in%"HG.2",])
princurve.hg2_hg2 = princurve::principal_curve(
  x = hg2.umap.embedding[cellorder.hg2_hg2,],  smoother = "smooth.spline")
cellorder.hg2_hg2 =  names(
  princurve.hg2_hg2$lambda[cellorder.hg2_hg2][order(princurve.hg2_hg2$lambda[cellorder.hg2_hg2])])
# Large.intestine.1
cellorder.hg2_l1 = 
  rownames(src.endoderm.hg2.ext_v1.1@meta.data[
    src.endoderm.hg2.ext_v1.1$cluster.v06.26.re%in%"Large.intestine.1",])
princurve.hg2_l1 = princurve::principal_curve(
  x = hg2.umap.embedding[cellorder.hg2_l1,],  smoother = "smooth.spline")
cellorder.hg2_l1 =  names(
  princurve.hg2_l1$lambda[cellorder.hg2_l1][order(princurve.hg2_l1$lambda[cellorder.hg2_l1])])
# Large.intestine.3
cellorder.hg2_l3 = 
  rownames(src.endoderm.hg2.ext_v1.1@meta.data[
    src.endoderm.hg2.ext_v1.1$cluster.v06.26.re%in%"Large.intestine.3",])
princurve.hg2_l3 = princurve::principal_curve(
  x = hg2.umap.embedding[cellorder.hg2_l3,],  smoother = "smooth.spline")
cellorder.hg2_l3 =  names(
  princurve.hg2_l3$lambda[cellorder.hg2_l3][order(princurve.hg2_l3$lambda[cellorder.hg2_l3])])
#---------------------------

pdf("figure.v08.07/try.pdf")
#---------------------------
plot(c(1:length(cellorder.hg2_hg2)),
     c(1:length(cellorder.hg2_hg2)),
     col = colors.time[src.endoderm.hg2.ext_v1.1[,cellorder.hg2_hg2]$Time])
plot(c(1:length(cellorder.hg2_l1)),
     c(1:length(cellorder.hg2_l1)),
     col = colors.time[src.endoderm.hg2.ext_v1.1[,cellorder.hg2_l1]$Time])
plot(c(1:length(cellorder.hg2_l3)),
     c(1:length(cellorder.hg2_l3)),
     col = colors.time[src.endoderm.hg2.ext_v1.1[,cellorder.hg2_l3]$Time])
#---------------------------
dev.off()

cellorder.hg2_hg2_FPI = cellorder.hg2_hg2
cellorder.hg2_l1_FPI = cellorder.hg2_l1
cellorder.hg2_l3_FPI = cellorder.hg2_l3

#=================================================================================
#=================================================================================
pdf("figure.v08.07/organ_development_re/hg2_heatmap_marker.re_l2.pdf",6,9)
cellorder.hg2 = c(cellorder.hg2_hg2_FPI, 
                  cellorder.hg2_l1_FPI,
                  rev(cellorder.hg2_l3_FPI))

markergene.hg2.ext_v1.1 = 
  FindAllMarkers(src.endoderm.hg2.ext_v1.1[,cellorder.hg2], assay = "RNA")
markergene.hg2.ext_v1.1 = markergene.hg2.ext_v1.1[
  markergene.hg2.ext_v1.1$avg_log2FC>0.2&markergene.hg2.ext_v1.1$p_val_adj<0.1, "gene"]
selectgene = unique(c(markergene.hg2.ext_v1.1))

selectgene = src.endoderm.hg2.ext_v1.1.filtergene_fin
#-------------------
hg2.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.hg2.ext_v1.1@assays$RNA@data[selectgene, cellorder.hg2]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.hg2.ext_v1.1$Time[cellorder.hg2], colors.time),
              MyName2Col(src.endoderm.hg2.ext_v1.1$cluster.v06.26.re[cellorder.hg2], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.hg2.ext_v1.1$cluster.extract.v1.1[cellorder.hg2], cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            # return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.hg2.ext_v1.1@meta.data[cell_hg2_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

hg2.ext_v1.1.re.rowtree.1 = as.dendrogram(hg2.ext_v1.1.re.rowtree)
src.endoderm.hg2.ext_v1.1.filtergene_fin = c(
  setdiff(
    c(
      labels(hg2.ext_v1.1.re.rowtree.1[[2]][[2]][[1]][[1]])
    ),
    c()
  ),
  setdiff(
    c(
      labels(hg2.ext_v1.1.re.rowtree.1[[2]][[2]][[2]]),
      labels(hg2.ext_v1.1.re.rowtree.1[[2]][[2]][[1]][[2]])
    ),
    c()
  ),
  setdiff(
    c(
      labels(hg2.ext_v1.1.re.rowtree.1[[2]][[1]])
    ),
      c()
  ),
  setdiff(
    c(
      labels(hg2.ext_v1.1.re.rowtree.1[[1]][[2]][[2]][[1]]),
      labels(hg2.ext_v1.1.re.rowtree.1[[1]][[2]][[1]]),
      labels(hg2.ext_v1.1.re.rowtree.1[[1]][[1]])
    ),
    c()
  )
)

names(src.endoderm.hg2.ext_v1.1.filtergene_fin) = c(
  rep(3,length(  
    setdiff(
      c(
        labels(hg2.ext_v1.1.re.rowtree.1[[2]][[2]][[1]][[1]])
      ),
      c()
    ))),
  rep(4,length(  
    setdiff(
      c(
        labels(hg2.ext_v1.1.re.rowtree.1[[2]][[2]][[2]]),
        labels(hg2.ext_v1.1.re.rowtree.1[[2]][[2]][[1]][[2]])
      ),
      c()
    ))),
  rep(9,length(  
    setdiff(
      c(
        labels(hg2.ext_v1.1.re.rowtree.1[[2]][[1]])
      ),
      c()
    ))),
  rep(7,length(  
    setdiff(
      c(
        labels(hg2.ext_v1.1.re.rowtree.1[[1]][[2]][[2]][[1]]),
        labels(hg2.ext_v1.1.re.rowtree.1[[1]][[2]][[1]]),
        labels(hg2.ext_v1.1.re.rowtree.1[[1]][[1]])
      ),
      c()
    )))
)


#--------------------------------------------------------------------------------


save(
  # Gene 
  src.endoderm.hg2.ext.v1.1.selectgene,
  src.endoderm.hg2.ext_v1.1.gene.fin,
  src.endoderm.hg2.ext_v1.1.filtergene_fin,
  # Cell type
  cellorder.hg2_hg2_FPI, 
  cellorder.hg2_l1_FPI,
  cellorder.hg2_l3_FPI,
  file= "figure.v08.07/organ_development_re/hg2_heatmap_parameter.Rdata")


