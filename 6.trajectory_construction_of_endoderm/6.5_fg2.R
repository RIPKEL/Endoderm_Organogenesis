#----------------
#   FG.2
#----------------

src.endoderm.fg2.ext_v1.1 = FindVariableFeatures(src.endoderm.fg2.ext_v1.1, nfeatures = 2000, assay = "RNA")

src.endoderm.fg2.ext_v1.1.selectgene = 
  Myfilter(as.matrix(src.endoderm.fg2.ext_v1.1@assays$RNA@data),
           gene = src.endoderm.fg2.ext_v1.1@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.endoderm.fg2.ext.v1.1.selectgene = 
  unique(c(select_gene_fg2, src.endoderm.fg2.ext_v1.1.selectgene))

src.endoderm.fg2.ext_v1.1 = SetIdent(src.endoderm.fg2.ext_v1.1,
                                     value = src.endoderm.fg2.ext_v1.1$cluster.revise.re.v1.30.re.fin)
markergene.fg2.ext_v1.1 = FindAllMarkers(
  src.endoderm.fg2.ext_v1.1[,src.endoderm.fg2.ext_v1.1$Time%in%c("ss9","ss27")], assay = "RNA")
markergene.fg2.ext_v1.1 = markergene.fg2.ext_v1.1[
  markergene.fg2.ext_v1.1$avg_log2FC>0.25&
    markergene.fg2.ext_v1.1$p_val_adj<0.1, "gene"]

pdf("figure.v08.07/try.pdf")
#-----------------
fg2.ext.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.fg2.ext_v1.1@assays$RNA@data[#src.endoderm.fg2.ext.v1.1.selectgene
                                                                src.endoderm.fg2.ext_v1.1.gene.fin,
                                                                #markergene.fg2.ext_v1.1.re,
                                                                ]),
            type = "row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.fg2.ext_v1.1$Time,colors.time),
              MyName2Col(src.endoderm.fg2.ext_v1.1$cluster.extract.v1.1,cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.fg2.ext_v1.1$cluster.v06.26.re,cluster.endoderm.color.v5)
            ),
            ColSideColorsSize = 1.5,
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            #Rowv = "none",
            return.tree = "col",
            graph = T)
#-------------------
dev.off()

pdf("figure.v08.07/try.tree.pdf",100,20)
plot(fg2.ext.re.rowtree)
dev.off()

fg2.ext.re.rowtree_1 = as.dendrogram(fg2.ext.re.rowtree)
src.endoderm.fg2.ext_v1.1.gene.fin = setdiff(
  src.endoderm.fg2.ext_v1.1.selectgene,
  c(labels(fg2.ext.re.rowtree_1[[1]]),
    labels(fg2.ext.re.rowtree_1[[2]][[2]][[2]][[2]])))

fg2.ext.re.rowtree_2 = as.dendrogram(fg2.ext.re.rowtree)
fg2.ext.re.rowtree_3 = as.dendrogram(fg2.ext.re.rowtree)
fg2.ext.re.rowtree_4 = as.dendrogram(fg2.ext.re.rowtree)
fg2.ext.re.rowtree_5 = as.dendrogram(fg2.ext.re.rowtree)

src.endoderm.fg2.ext_v1.1.gene.fin = markergene.fg2.ext_v1.1
src.endoderm.fg2.ext_v1.1.gene.fin = setdiff(
  markergene.fg2.ext_v1.1,
  c(labels(fg2.ext.re.rowtree_2[[1]][[1]]),
    labels(fg2.ext.re.rowtree_2[[2]][[1]][[2]][[1]])))

fg2.ext.re.rowtree_a = as.dendrogram(fg2.ext.re.rowtree)
fg2.ext.re.rowtree_b = as.dendrogram(fg2.ext.re.rowtree)

src.endoderm.fg2.ext_v1.1.gene.fin = markergene.fg2.ext_v1.1.re
src.endoderm.fg2.ext_v1.1.gene.fin = setdiff(
  markergene.fg2.ext_v1.1,
  c(labels(fg2.ext.re.rowtree_a[[1]][[2]]),
    labels(fg2.ext.re.rowtree_a[[1]][[1]][[2]][[2]][[2]]),
    labels(fg2.ext.re.rowtree_a[[2]][[2]][[2]][[2]][[2]][[2]]),
    labels(fg2.ext.re.rowtree_b[[2]][[2]][[2]])))


fg2.ext.re.coltree = as.dendrogram(fg2.ext.re.rowtree)
src.endoderm.fg2.ext_v1.1$cluster.v06.26.re_hc = NA
src.endoderm.fg2.ext_v1.1@meta.data[labels(fg2.ext.re.coltree[[2]][[1]]),]$cluster.v06.26.re_hc = "FG.2"
src.endoderm.fg2.ext_v1.1@meta.data[c(labels(fg2.ext.re.coltree[[2]][[2]]),
                                      labels(fg2.ext.re.coltree[[1]])),]$cluster.v06.26.re_hc = "Pharynx.organ.1"



#-------------------------------------------------
src.endoderm.fg2.ext_v1.1 = 
  seurat_mnn(seurat_name = "src.endoderm.fg2.ext_v1.1",
             seurat.selectgene = markergene.fg2.ext_v1.1.re)

src.endoderm.fg2.ext_v1.1 = 
  seurat_int_phase(seurat_name = "src.endoderm.fg2.ext_v1.1",
             seurat.selectgene = markergene.fg2.ext_v1.1.re)

src.endoderm.fg2.ext_v1.1 = FindNeighbors(src.endoderm.fg2.ext_v1.1, dims=1:30)
src.endoderm.fg2.ext_v1.1 = FindClusters(src.endoderm.fg2.ext_v1.1, resolution = 2)

src.endoderm.fg2.ext_v1.1 = FindNeighbors(src.endoderm.fg2.ext_v1.1, dims=1:30, 
                                          assay = "integrated", reduction = "pca_integrated")
src.endoderm.fg2.ext_v1.1 = FindClusters(src.endoderm.fg2.ext_v1.1, resolution = 2)


DimPlot(src.endoderm.fg2.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.2', reduction = "umap_integrated")
DimPlot(src.endoderm.fg2.ext_v1.1, group.by = 'Time', 
        reduction = "umap", cols = colors.time)
DimPlot(src.endoderm.fg2.ext_v1.1, group.by = 'Time', 
        reduction = "umap_integrated", cols = colors.time)
DimPlot(src.endoderm.fg2.ext_v1.1, group.by = 'Phase', 
        reduction = "umap_integrated") #, cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.fg2.ext_v1.1, group.by = 'cluster.extract.v1.1', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.fg2.ext_v1.1, group.by = 'cluster.v06.26.re',
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.fg2.ext_v1.1, group.by = 'cluster.v06.26.re_hc',
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)
DimPlot(src.endoderm.fg2.ext_v1.1, group.by = 'cluster.v06.10', 
        reduction = "umap_integrated", cols = cluster.endoderm.color.v5)


src.endoderm.fg2.ext_v1.1$cluster.v06.26.re = NA
src.endoderm.fg2.ext_v1.1@meta.data[
  src.endoderm.fg2.ext_v1.1$RNA_snn_res.2%in%c(
    9,12,15,14,0,6,8),]$cluster.v06.26.re = "FG.2"
src.endoderm.fg2.ext_v1.1@meta.data[
  src.endoderm.fg2.ext_v1.1$RNA_snn_res.2%in%c(
    5,4,17,2,10,11,1,18,19,13,3,7,16),]$cluster.v06.26.re = "Pharynx.organ.1"
src.endoderm.fg2.ext_v1.1@reductions$umap_integrated@cell.embeddings = 
  src.endoderm.fg2.ext_v1.1@reductions$umap_integrated@cell.embeddings[colnames(src.endoderm.fg2.ext_v1.1),]
src.endoderm.fg2.ext_v1.1@meta.data[
  (src.endoderm.fg2.ext_v1.1[["umap_integrated"]]@cell.embeddings[,1]>0.5)&
    (src.endoderm.fg2.ext_v1.1[["umap_integrated"]]@cell.embeddings[,2]>0),]$cluster.v06.26.re = "FG.2"

# Correct
src.endoderm.fg2.ext_v1.1@meta.data[
  src.endoderm.fg2.ext_v1.1$Time%in%c("ss24","ss27"),]$cluster.v06.26.re_hc =  "Pharynx.organ.1"


src.endoderm.fg2.ext_v1.1 = SetIdent(src.endoderm.fg2.ext_v1.1,
                                     value = src.endoderm.fg2.ext_v1.1$cluster.v06.26.re)
markergene.fg2.ext_v1.1.re = FindAllMarkers(
  src.endoderm.fg2.ext_v1.1[,src.endoderm.fg2.ext_v1.1$Time%in%c("ss9","ss27")], assay = "RNA")
markergene.fg2.ext_v1.1.re = markergene.fg2.ext_v1.1.re[
  markergene.fg2.ext_v1.1.re$avg_log2FC>0.25&
    markergene.fg2.ext_v1.1.re$p_val_adj<0.05, "gene"]



# FDL
src.endoderm.fg2.ext_v1.1 = 
  seurat_fdl("src.endoderm.fg2.ext_v1.1", "pca","RNA")
src.endoderm.fg2.ext_v1.1 = 
  seurat_fdl("src.endoderm.fg2.ext_v1.1", "mnn","RNA")
src.endoderm.fg2.ext_v1.1 = 
  seurat_fdl("src.endoderm.fg2.ext_v1.1", "pca_integrated","integrated")

# Graph
graph.fg2.ext_v1.1 = 
  seurat_GCN(seurat_name = "src.endoderm.fg2.ext_v1.1",
             src.endoderm.fg2.ext_v1.1.gene.fin,
             reduction = "pca_integrated", assay = "integrated")
graph.fg2_cluster = cluster_walktrap(graph.fg2.ext_v1.1, steps = 1)
#----------------------------


pdf("figure.v08.07/organ_development_re/fg2_summary_graph.re_snn.pdf",12,7)
# SNN
for(color in c("cluster.extract.v1.1","cluster.v06.26.re","Time","cluster.v06.26.re_hc")){
  for(edge.color in c(NA,"lightgray")){
    for(edge.color in c(NA,"lightgray")){
      set.seed(1)
      plot.igraph(graph.fg2.ext_v1.1,
                  layout = layout_with_fr(graph.fg2.ext_v1.1),
                  edge.color = edge.color,
                  vertex.size=2.5, vertex.label=NA,
                  vertex.label.color = "black", 
                  vertex.frame.color = NA, vertex.frame.width = 0.5,
                  #vertex.color =  colors.time[src.endoderm.fg2.ext_v1.1[["Time"]][graph.fg2_cluster$names,]],
                  #vertex.color =  cluster.endoderm.color.v5[src.endoderm.fg2.ext_v1.1[["cluster.v06.26.re"]][graph.fg2_cluster$names,]],
                  vertex.color =  colors_bg[src.endoderm.fg2.ext_v1.1[[color]][graph.fg2_cluster$names,]])
    }}}
dev.off()

pdf("figure.v08.07/organ_development_re/fg2_summary_graph.re_other.pdf",12,7)
# Reduction
for(red in c("umap_integrated","umap_mnn",
             "fdl_pca",'fdl_mnn',"fdl_pca_integrated")){
  p1 = DimPlot(src.endoderm.fg2.ext_v1.1,group.by = 'cluster.v06.26.re', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p2 = DimPlot(src.endoderm.fg2.ext_v1.1,group.by = 'cluster.extract.v1.1', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p3 = DimPlot(src.endoderm.fg2.ext_v1.1,group.by = 'Time', 
               reduction = red , cols = colors.time,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  p4 = DimPlot(src.endoderm.fg2.ext_v1.1,group.by = 'cluster.v06.26.re_hc', 
               reduction = red, cols = cluster.endoderm.color.v5,
               pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
  print(p1); print(p2); print(p3); print(p4)
}
dev.off()


pdf("figure.v08.07/organ_development_re/fg2_summary_proportion.re.pdf",12,7)
#-----------------------------
cell_type = c("FG.2","Pharynx.organ.1")
meta_fg2 = src.endoderm.fg2.ext_v1.1@meta.data
meta_fg2$Time = factor(meta_fg2$Time, levels = names(colors.time))
meta_fg2$cluster.extract.v1.1 = factor(meta_fg2$cluster.extract.v1.1, levels = c("FG.2"))
meta_fg2$cluster.v06.26.re = factor(meta_fg2$cluster.v06.26.re, levels = rev(cell_type))
meta_fg2$cluster.v06.26.re_hc = factor(meta_fg2$cluster.v06.26.re_hc, levels = rev(cell_type))

meta_fg2_B0 = meta_fg2
meta_fg2_B1 = meta_fg2[meta_fg2$batch%in%1&!meta_fg2$Time%in%"ss9",]
meta_fg2_B2 = meta_fg2[meta_fg2$batch%in%2,]

meta_fg2_BO.fg1 = meta_fg2[meta_fg2$cluster.extract.v1.1%in%"FG.2",]
meta_fg2_B1.fg1 = meta_fg2[meta_fg2$batch%in%1&!meta_fg2$Time%in%"ss9"&
                             meta_fg2$cluster.extract.v1.1%in%"FG.2",]
meta_fg2_B2.fg1 = meta_fg2[meta_fg2$batch%in%2&
                             meta_fg2$cluster.extract.v1.1%in%"FG.2",]


for(data in c("meta_fg2_B0","meta_fg2_B1","meta_fg2_B2",
              "meta_fg2_BO.fg1","meta_fg2_B1.fg1","meta_fg2_B2.fg1")){
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

fg2.umap.embedding = src.endoderm.fg2.ext_v1.1[["fdl_pca_integrated"]]@cell.embeddings[,1:3]
#---------------------------
# FG.2
cellorder.fg2_fg2 = 
  rownames(src.endoderm.fg2.ext_v1.1@meta.data[
    src.endoderm.fg2.ext_v1.1$cluster.v06.26.re_hc%in%"FG.2",])
princurve.fg2_fg2 = princurve::principal_curve(
  x = fg2.umap.embedding[cellorder.fg2_fg2,],  smoother = "smooth.spline")
cellorder.fg2_fg2 =  names(
  princurve.fg2_fg2$lambda[cellorder.fg2_fg2][order(princurve.fg2_fg2$lambda[cellorder.fg2_fg2])])

# Pha1
cellorder.fg2_ph1 = 
  rownames(src.endoderm.fg2.ext_v1.1@meta.data[
    src.endoderm.fg2.ext_v1.1$cluster.v06.26.re%in%"Pharynx.organ.1",])
princurve.fg2_ph1 = princurve::principal_curve(
  x = fg2.umap.embedding[cellorder.fg2_ph1,],  smoother = "smooth.spline")
cellorder.fg2_ph1 =  names(
  princurve.fg2_ph1$lambda[cellorder.fg2_ph1][order(princurve.fg2_ph1$lambda[cellorder.fg2_ph1])])
#---------------------------

plot(c(1:length(cellorder.fg2_fg2)),
     c(1:length(cellorder.fg2_fg2)),
     col = colors.time[src.endoderm.fg2.ext_v1.1[,cellorder.fg2_fg2]$Time])
plot(c(1:length(cellorder.fg2_ph1)),
     c(1:length(cellorder.fg2_ph1)),
     col = colors.time[src.endoderm.fg2.ext_v1.1[,cellorder.fg2_ph1]$Time])

cellorder.fg2_fg2_FPI = rev(cellorder.fg2_fg2)
cellorder.fg2_ph1_FPI = rev(cellorder.fg2_ph1)


#cellorder.fg2 = c(cellorder.fg2_fg2_FPI, cellorder.fg2_ph1_FPI)
cellorder.fg2 = cellorder.fg2_col
#selectgene = markergene.fg2.ext_v1.1.re
selectgene = src.endoderm.fg2.ext_v1.1.re.markergene
pdf("figure.v08.07/organ_development_re/fg2_heatmap_marker.re.pdf",6,9)
#-------------------
fg2.ext_v1.1.rowtree =
  MyHeatmap(as.matrix(src.endoderm.fg2.ext_v1.1@assays$RNA@data[selectgene, cellorder.fg2]),
            type = "log.row.relat",
            hc.c.data.type = "row.relat",
            hc.r.data.type = "row.relat",
            c.cov.method = "s",
            r.cov.method = "s",
            c.hc.method = "ward.D",
            r.hc.method = "ward.D2",
            ColSideColors = cbind(
              MyName2Col(src.endoderm.fg2.ext_v1.1$Time[cellorder.fg2], colors.time),
              MyName2Col(src.endoderm.fg2.ext_v1.1$cluster.v06.26.re_hc[cellorder.fg2], cluster.endoderm.color.v5),
              MyName2Col(src.endoderm.fg2.ext_v1.1$cluster.extract.v1.1[cellorder.fg2], cluster.endoderm.color.v5)
            ),
            RowSideColors = t(cbind(
              MyName2Col(names(selectgene), colors.geneset)
            )),
            ColSideColorsSize = 2.5,
            Rowv = "none",
            Colv = "none",
            #return.tree = "row",
            labRow=selectgene,
            #column_split = data.frame(src.endoderm.fg1.ext.re@meta.data[cell_fg1_curve ,]$cluster.v06.26),
            #column_km = 4,
            margins = c(4,4),
            graph = T)
#-------------------
dev.off()

fg2.ext_v1.1.rowtree_1 = as.dendrogram(fg2.ext_v1.1.rowtree)
src.endoderm.fg2.ext_v1.1.re.markergene = c(
  labels(fg2.ext_v1.1.rowtree_1[[1]][[2]][[1]]),
  labels(fg2.ext_v1.1.rowtree_1[[1]][[1]]),
  labels(fg2.ext_v1.1.rowtree_1[[2]][[1]]),
  #labels(fg2.ext_v1.1.rowtree_1[[2]][[2]][[1]]),
  labels(fg2.ext_v1.1.rowtree_1[[2]][[2]][[2]][[2]][[2]][[1]]))

fg2.ext_v1.1.rowtree_2 = as.dendrogram(fg2.ext_v1.1.rowtree)
cell_Ext_fg2 = c("Gm2694","Sorbs2")
src.endoderm.fg2.ext_v1.1.re.markergene = c(
  rev(labels(fg2.ext_v1.1.rowtree_2[[2]][[2]])),
  setdiff(labels(fg2.ext_v1.1.rowtree_2[[1]]),
          cell_Ext_fg2))
names(src.endoderm.fg2.ext_v1.1.re.markergene) = c(
  rep(4, length(labels(fg2.ext_v1.1.rowtree_2[[2]][[2]]))),
  rep(7, length(labels(fg2.ext_v1.1.rowtree_2[[1]]))-
        length(cell_Ext_fg2)))

fg2.ext_v1.1.rowtree_3 = as.dendrogram(fg2.ext_v1.1.rowtree)
src.endoderm.fg2.ext_v1.1.re.markergene = c(
  rev(labels(fg2.ext_v1.1.rowtree_3[[1]][[2]][[2]])),
  labels(fg2.ext_v1.1.rowtree_3[[1]][[1]]),
  labels(fg2.ext_v1.1.rowtree_3[[1]][[2]][[1]]),
  labels(fg2.ext_v1.1.rowtree_3[[2]][[1]][[1]]),
  labels(fg2.ext_v1.1.rowtree_3[[2]][[2]][[1]]),
  labels(fg2.ext_v1.1.rowtree_3[[2]][[2]][[2]]),
  rev(labels(fg2.ext_v1.1.rowtree_3[[2]][[1]][[2]])))
names(src.endoderm.fg2.ext_v1.1.re.markergene) = c(
  rep(4, length(labels(fg2.ext_v1.1.rowtree_3[[1]]))),
  rep(7, length(labels(fg2.ext_v1.1.rowtree_3[[2]]))))


fg2.ext_v1.1.coltree_2 = as.dendrogram(fg2.ext_v1.1.rowtree)
cellorder.fg2_col = c(
  labels(fg2.ext_v1.1.coltree_2[[2]][[1]]),
  labels(fg2.ext_v1.1.coltree_2[[2]][[2]]),
  rev(labels(fg2.ext_v1.1.coltree_2[[1]])))



save(src.endoderm.fg2.ext_v1.1.gene.fin,
     markergene.fg2.ext_v1.1.re,
     src.endoderm.fg2.ext_v1.1.selectgene,
     cellorder.fg2_fg2_FPI, cellorder.fg2_ph1_FPI, cellorder.fg2_col,
     file= "figure.v08.07/organ_development_re/fg2_heatmap_parameter.Rdata")
save(src.endoderm.fg2.ext_v1.1,
     file = "figure.v08.07/organ_development_re/src.endoderm.fg2.ext_v1.1.Rdata")

























