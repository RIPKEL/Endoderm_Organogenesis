#----------------
#   MG.1
#----------------
src.endoderm.mg1.ext_v1.1 = FindVariableFeatures(src.endoderm.mg1.ext_v1.1, nfeatures = 2000)
src.endoderm.mg1.ext_v1.1 = ScaleData(src.endoderm.mg1.ext_v1.1, split.by = "batch_phase",
                                      features = rownames(src.endoderm.mg1.ext_v1.1))
src.endoderm.mg1.ext_v1.1.filtergene = 
  Myfilter(as.matrix(src.endoderm.mg1.ext_v1.1@assays$RNA@data),
           gene = src.endoderm.mg1.ext_v1.1@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
#rm(src.endoderm.mg1.ext_v1.1.filtergene)
src.endoderm.mg1.ext.v1.1.selectgene = unique(c(
  select_gene_mg1, src.endoderm.mg1.ext_v1.1.filtergene))

#=============================
# Row
pdf("figure.v08.07/try.pdf")
#-----------------
mg1.ext.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.mg1.ext_v1.1@assays$RNA@data[
    src.endoderm.mg1.ext.v1.1.selectgene,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.endoderm.mg1.ext_v1.1$Time,colors.time),
      MyName2Col(src.endoderm.mg1.ext_v1.1$cluster.extract.v1.1,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.mg1.ext_v1.1$cluster.v06.26.re,cluster.endoderm.color.v5)
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
plot(mg1.ext.re.rowtree)
dev.off()

mg1.ext.re.rowtree_1 = as.dendrogram(mg1.ext.re.rowtree)
src.endoderm.mg1.ext_v1.1.gene.fin = setdiff(
  src.endoderm.mg1.ext.v1.1.selectgene,
  c(labels(mg1.ext.re.rowtree_1[[2]][[1]]),
    labels(mg1.ext.re.rowtree_1[[2]][[2]][[1]])))

#--------------------------------------------------------------------------------
src.endoderm.mg1.ext_v1.1 = RunPCA(src.endoderm.mg1.ext_v1.1,
                                   features = src.endoderm.mg1.ext_v1.1.gene.fin)
src.endoderm.mg1.ext_v1.1 = 
  seurat_mnn(seurat_name = "src.endoderm.mg1.ext_v1.1", dims = 1:30, 
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.mg1.ext_v1.1.gene.fin)
src.endoderm.mg1.ext_v1.1 = 
  seurat_int(seurat_name = "src.endoderm.mg1.ext_v1.1", dims = 1:30,
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.mg1.ext_v1.1.gene.fin)

#---------------
src.endoderm.mg1.ext_v1.1 = 
  FindNeighbors(src.endoderm.mg1.ext_v1.1, dims=1:30, 
                reduction = "pca_integrated", assay = "integrated")
src.endoderm.mg1.ext_v1.1 = 
  FindNeighbors(src.endoderm.mg1.ext_v1.1, dims=1:30, 
                reduction = "pca", assay = "RNA")
src.endoderm.mg1.ext_v1.1 = 
  FindClusters(src.endoderm.mg1.ext_v1.1, resolution = 2)
src.endoderm.mg1.ext_v1.1 = 
  FindClusters(src.endoderm.mg1.ext_v1.1, resolution = 4)

src.endoderm.mg1.ext_v1.1 = 
  FindNeighbors(src.endoderm.mg1.ext_v1.1, dims=1:2, 
                reduction = "umap", assay = "RNA")
src.endoderm.mg1.ext_v1.1 = 
  FindClusters(src.endoderm.mg1.ext_v1.1, resolution = 1)
#---------------

pdf("figure.v08.07/try.pdf")
DimPlot(src.endoderm.mg1.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.2', reduction = "umap_mnn")
DimPlot(src.endoderm.mg1.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.4', reduction = "umap_mnn")

DimPlot(src.endoderm.mg1.ext_v1.1, group.by = 'Time',
        reduction = "umap_mnn", cols = colors.time)
DimPlot(src.endoderm.mg1.ext_v1.1, group.by = 'Time', 
        reduction = "umap_integrated", cols = colors.time)
DimPlot(src.endoderm.mg1.ext_v1.1, group.by = 'batch', 
        reduction = "umap_integrated") #, cols = colors.time)

DimPlot(src.endoderm.mg1.ext_v1.1, group.by = 'cluster.extract.v1.1', 
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.mg1.ext_v1.1, group.by = 'cluster.extract.v1.1', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.mg1.ext_v1.1, group.by = 'cluster.v06.26.re', 
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.mg1.ext_v1.1, group.by = 'cluster.v06.26.re', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
dev.off()

src.endoderm.mg1.ext_v1.1@reductions$pca_integrated@cell.embeddings = 
  src.endoderm.mg1.ext_v1.1@reductions$pca_integrated@cell.embeddings[colnames(src.endoderm.mg1.ext_v1.1),]
#--------------------
src.endoderm.mg1.ext_v1.1$cluster.v06.26.re_raw = NA
src.endoderm.mg1.ext_v1.1@meta.data[
  intersect(colnames(src.endoderm.mg1.ext_v1.1),colnames(src.endoderm.mg1.ext.re)),]$cluster.v06.26.re_raw = 
  src.endoderm.mg1.ext.re@meta.data[
    intersect(colnames(src.endoderm.mg1.ext_v1.1),colnames(src.endoderm.mg1.ext.re)),]$cluster.v06.26

#---------------------
src.endoderm.mg1.ext_v1.1$cluster.v06.26.re = NA
src.endoderm.mg1.ext_v1.1@meta.data[
  intersect(colnames(src.endoderm.mg1.ext_v1.1),colnames(src.endoderm.mg1.ext.re)),]$cluster.v06.26.re = 
  src.endoderm.mg1.ext.re@meta.data[
    intersect(colnames(src.endoderm.mg1.ext_v1.1),colnames(src.endoderm.mg1.ext.re)),]$cluster.v06.26

src.endoderm.mg1.ext_v1.1@meta.data[
  src.endoderm.mg1.ext_v1.1$RNA_snn_res.2%in%
    c(11,12,14,15,10,5,6),]$cluster.v06.26.re = "MG.1"
src.endoderm.mg1.ext_v1.1@meta.data[
  src.endoderm.mg1.ext_v1.1$RNA_snn_res.2%in%
    c(1,2,7,8,9,16),]$cluster.v06.26.re = "Stomach"
src.endoderm.mg1.ext_v1.1@meta.data[
  src.endoderm.mg1.ext_v1.1$RNA_snn_res.2%in%
    c(0,3,4),]$cluster.v06.26.re = "Small.intestine.1"

src.endoderm.mg1.ext_v1.1@meta.data[
  src.endoderm.mg1.ext_v1.1$cluster.extract.v1.1%in%
    c("FG.4-MG.1/3"),]$cluster.v06.26.re = "Stomach"
src.endoderm.mg1.ext_v1.1@meta.data[
  src.endoderm.mg1.ext_v1.1$cluster.extract.v1.1%in%
    c("MG.1/2/3","AL.3-MG.1/2/3"),]$cluster.v06.26.re = "Small.intestine.1"

a = FNN::knn(
  src.endoderm.mg1.ext_v1.1@reductions$pca_integrated@cell.embeddings[
    rownames(src.endoderm.mg1.ext_v1.1@meta.data[!src.endoderm.mg1.ext_v1.1$cluster.v06.26.re%in%NA,]),],
  src.endoderm.mg1.ext_v1.1@reductions$pca_integrated@cell.embeddings[
    rownames(src.endoderm.mg1.ext_v1.1@meta.data[src.endoderm.mg1.ext_v1.1$cluster.v06.26.re%in%NA,]),],
  src.endoderm.mg1.ext_v1.1@meta.data[!src.endoderm.mg1.ext_v1.1$cluster.v06.26.re%in%NA,]$cluster.v06.26.re, k = 10)
src.endoderm.mg1.ext_v1.1@meta.data[
  src.endoderm.mg1.ext_v1.1$cluster.v06.26.re%in%NA,]$cluster.v06.26.re = as.character(a)

DimPlot(src.endoderm.mg1.ext_v1.1, group.by = 'cluster.v06.26.re', 
        reduction = "umap_mnn", cols = cluster.endoderm.color.v5)

# FDL
src.endoderm.mg1.ext_v1.1 = 
  seurat_fdl("src.endoderm.mg1.ext_v1.1", "pca","RNA")
src.endoderm.mg1.ext_v1.1 = 
  seurat_fdl("src.endoderm.mg1.ext_v1.1", "mnn","RNA")
src.endoderm.mg1.ext_v1.1 = 
  seurat_fdl("src.endoderm.mg1.ext_v1.1", "pca_integrated","integrated")


# Graph
pdf("figure.v08.07/organ_development_re/mg1_summary_graph.re_snn_PI.pdf",12,7)
# SNN
for(red in c('pca',"mnn","pca_integrated")){
  
  assay = ifelse(red=="mnn","RNA","integrated")
  
  graph.mg1.ext_v1.1 = 
    seurat_GCN(seurat_name = "src.endoderm.mg1.ext_v1.1",
               src.endoderm.mg1.ext_v1.1.gene.fin,
               reduction = red, assay =  assay)
  graph.mg1_cluster = cluster_walktrap(graph.mg1.ext_v1.1, steps = 1)
  
  for(color in c("cluster.extract.v1.1","cluster.v06.26.re",
                 "Time")){
    for(edge.color in c(NA,"lightgray")){
      for(edge.color in c(NA,"lightgray")){
        set.seed(1)
        plot.igraph(graph.mg1.ext_v1.1,
                    layout = layout_with_fr(graph.mg1.ext_v1.1),
                    edge.color = edge.color,
                    vertex.size=2.5, vertex.label=NA,
                    vertex.label.color = "black", 
                    vertex.frame.color = NA, vertex.frame.width = 0.5,
                    #vertex.color =  colors.time[src.endoderm.mg1.ext_v1.1[["Time"]][graph.mg1_cluster$names,]],
                    #vertex.color =  cluster.endoderm.color.v5[src.endoderm.mg1.ext_v1.1[["cluster.v06.26.re"]][graph.mg1_cluster$names,]],
                    vertex.color =  colors_bg[src.endoderm.mg1.ext_v1.1[[color]][graph.mg1_cluster$names,]])
      }}}
  
}
dev.off()

pdf("figure.v08.07/organ_development_re/mg1_summary_graph.re_other.pdf",12,7)
# Reduction
for(red in c("umap_integrated","umap_mnn",
             "fdl_pca",'fdl_mnn',"fdl_pca_integrated")){
  p1 = DimPlot(src.endoderm.mg1.ext_v1.1,group.by = 'cluster.v06.26.re', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p2 = DimPlot(src.endoderm.mg1.ext_v1.1,group.by = 'cluster.extract.v1.1', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p3 = DimPlot(src.endoderm.mg1.ext_v1.1,group.by = 'Time', 
               reduction = red , cols = colors.time,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  # p4 = DimPlot(src.endoderm.mg1.ext_v1.1,group.by = 'cluster.v06.26.re_hc', 
  #              reduction = red, cols = cluster.endoderm.color.v5,
  #              pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  print(p1); print(p2); print(p3);# print(p4)
}
dev.off()


pdf("figure.v08.07/organ_development_re/mg1_summary_proportion.re.pdf",12,7)
#-----------------------------
cell_type = c("MG.1","Small.intestine.1","Stomach")
meta_mg1 = src.endoderm.mg1.ext_v1.1@meta.data
meta_mg1$Time = factor(meta_mg1$Time, levels = names(colors.time))
meta_mg1$cluster.extract.v1.1 = factor(meta_mg1$cluster.extract.v1.1, 
                                       levels = c("MG.1","MG.1/2/3","AL.3-MG.1/2/3","FG.4-MG.1/3"))
meta_mg1$cluster.v06.26.re = factor(meta_mg1$cluster.v06.26.re, levels = rev(cell_type))

meta_mg1_B0 = meta_mg1
meta_mg1_B1 = meta_mg1[meta_mg1$batch%in%1&!meta_mg1$Time%in%"ss9",]
meta_mg1_B2 = meta_mg1[meta_mg1$batch%in%2,]

for(data in c("meta_mg1_B0","meta_mg1_B1","meta_mg1_B2")){
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

save(src.endoderm.mg1.ext_v1.1,
     file = "figure.v08.07/organ_development_re/src.endoderm.mg1.ext_v1.1.re.Rdata")



# Marker
#------------------------
src.endoderm.mg1.ext_v1.1 = 
  SetIdent(src.endoderm.mg1.ext_v1.1,
           value = src.endoderm.mg1.ext_v1.1$cluster.v06.26.re)
markergene.mg1.ext_v1.1 = 
  FindAllMarkers(src.endoderm.mg1.ext_v1.1, assay = "RNA")
markergene.mg1.ext_v1.1 = markergene.mg1.ext_v1.1[
  markergene.mg1.ext_v1.1$avg_log2FC>0.2&
    markergene.mg1.ext_v1.1$p_val_adj<0.1, "gene"]

#------------------------
mg1.umap.embedding = src.endoderm.mg1.ext_v1.1[["fdl_pca_integrated"]]@cell.embeddings[,1:2]
mg1.umap.embedding = src.endoderm.mg1.ext_v1.1[["fdl_mnn"]]@cell.embeddings[,1:2]
mg1.umap.embedding = src.endoderm.mg1.ext_v1.1[["umap_mnn"]]@cell.embeddings[,1:3]
mg1.umap.embedding = src.endoderm.mg1.ext_v1.1[["umap_integrated"]]@cell.embeddings[,1:3]

# Cell order :: Princurve
#---------------------------
# MG.1
cellorder.mg1_mg1 = 
  rownames(src.endoderm.mg1.ext_v1.1@meta.data[
    src.endoderm.mg1.ext_v1.1$cluster.v06.26.re%in%"MG.1",])
princurve.mg1_mg1 = princurve::principal_curve(
  x = mg1.umap.embedding[cellorder.mg1_mg1,],  smoother = "smooth.spline")
cellorder.mg1_mg1 =  names(
  princurve.mg1_mg1$lambda[cellorder.mg1_mg1][order(princurve.mg1_mg1$lambda[cellorder.mg1_mg1])])
# Stomach
cellorder.mg1_sto = 
  rownames(src.endoderm.mg1.ext_v1.1@meta.data[
    src.endoderm.mg1.ext_v1.1$cluster.v06.26.re%in%"Stomach",])
princurve.mg1_sto = princurve::principal_curve(
  x = mg1.umap.embedding[cellorder.mg1_sto,],  smoother = "smooth.spline")
cellorder.mg1_sto =  names(
  princurve.mg1_sto$lambda[cellorder.mg1_sto][order(princurve.mg1_sto$lambda[cellorder.mg1_sto])])
# Small.intestine.1
cellorder.mg1_s1 = 
  rownames(src.endoderm.mg1.ext_v1.1@meta.data[
    src.endoderm.mg1.ext_v1.1$cluster.v06.26.re%in%"Small.intestine.1",])
princurve.mg1_s1 = princurve::principal_curve(
  x = mg1.umap.embedding[cellorder.mg1_s1,],  smoother = "smooth.spline")
cellorder.mg1_s1 =  names(
  princurve.mg1_s1$lambda[cellorder.mg1_s1][order(princurve.mg1_s1$lambda[cellorder.mg1_s1])])
#---------------------------

pdf("figure.v08.07/try.pdf")
#---------------------------
plot(c(1:length(cellorder.mg1_mg1)),
     c(1:length(cellorder.mg1_mg1)),
     col = colors.time[src.endoderm.mg1.ext_v1.1[,cellorder.mg1_mg1]$Time])
plot(c(1:length(cellorder.mg1_sto)),
     c(1:length(cellorder.mg1_sto)),
     col = colors.time[src.endoderm.mg1.ext_v1.1[,cellorder.mg1_sto]$Time])
plot(c(1:length(cellorder.mg1_s1)),
     c(1:length(cellorder.mg1_s1)),
     col = colors.time[src.endoderm.mg1.ext_v1.1[,cellorder.mg1_s1]$Time])
#---------------------------
dev.off()

cellorder.mg1_mg1_FPI = cellorder.mg1_mg1
cellorder.mg1_sto_FM = cellorder.mg1_sto
cellorder.mg1_s1_FM = cellorder.mg1_s1

#=================================================================================
#=================================================================================

pdf("figure.v08.07/organ_development_re/mg1_heatmap_marker.re.pdf",6,9)
cellorder.mg1 = c(cellorder.mg1_mg1_FPI, 
                  rev(cellorder.mg1_sto_FM),
                  cellorder.mg1_s1_FM)
markergene.mg1.ext_v1.1 = 
  FindAllMarkers(src.endoderm.mg1.ext_v1.1[,cellorder.mg1], assay = "RNA")
markergene.mg1.ext_v1.1 = markergene.mg1.ext_v1.1[
  markergene.mg1.ext_v1.1$avg_log2FC>0.2&markergene.mg1.ext_v1.1$p_val_adj<0.1, "gene"]
#selectgene = unique(c(markergene.mg1.ext_v1.1))
selectgene = src.endoderm.mg1.ext_v1.1.filtergene_fin
#-------------------
mg1.ext_v1.1.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.mg1.ext_v1.1@assays$RNA@data[selectgene, cellorder.mg1]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.mg1.ext_v1.1$Time[cellorder.mg1], colors.time),
              MyName2Col(src.endoderm.mg1.ext_v1.1$cluster.v06.26.re[cellorder.mg1], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.mg1.ext_v1.1$cluster.extract.v1.1[cellorder.mg1], cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            #return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.mg1.ext_v1.1@meta.data[cell_mg1_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

mg1.ext_v1.1.re.rowtree.1 = as.dendrogram(mg1.ext_v1.1.re.rowtree)
src.endoderm.mg1.ext_v1.1.filtergene_fin = c(
  setdiff(
    c(
      labels(mg1.ext_v1.1.re.rowtree.1[[2]][[1]][[2]])
    ),
    c()
  ),
  setdiff(
    c(
      labels(mg1.ext_v1.1.re.rowtree.1[[2]][[2]])
    ),
    c()
  ),
  setdiff(
    c(
      labels(mg1.ext_v1.1.re.rowtree.1[[1]][[2]]),
      labels(mg1.ext_v1.1.re.rowtree.1[[1]][[1]][[1]][[1]]),
      labels(mg1.ext_v1.1.re.rowtree.1[[1]][[1]][[1]][[2]]),
      labels(mg1.ext_v1.1.re.rowtree.1[[1]][[1]][[2]])
    ),
    c()
  )
)

names(src.endoderm.mg1.ext_v1.1.filtergene_fin) = c(
  rep(4,length(  
    setdiff(
      c(
        labels(mg1.ext_v1.1.re.rowtree.1[[2]][[1]][[2]])
      ),
      c()
    ))),
  rep(3,length(  
    setdiff(
      c(
        labels(mg1.ext_v1.1.re.rowtree.1[[2]][[2]])
      ),
      c()
    ))),
  rep(7,length(  
    setdiff(
      c(
        labels(mg1.ext_v1.1.re.rowtree.1[[1]][[1]][[1]][[1]]),
        labels(mg1.ext_v1.1.re.rowtree.1[[1]][[2]]),
        labels(mg1.ext_v1.1.re.rowtree.1[[1]][[1]][[1]][[2]]),
        labels(mg1.ext_v1.1.re.rowtree.1[[1]][[1]][[2]])
      ),
      c()
    )))
)



save(
  # Gene 
  src.endoderm.mg1.ext.v1.1.selectgene,
  src.endoderm.mg1.ext_v1.1.gene.fin,
  src.endoderm.mg1.ext_v1.1.filtergene_fin,
  # Cell type
  cellorder.mg1_mg1_FPI, cellorder.mg1_sto_FM, cellorder.mg1_s1_FM, 
  file= "figure.v08.07/organ_development_re/mg1_heatmap_parameter.Rdata")



