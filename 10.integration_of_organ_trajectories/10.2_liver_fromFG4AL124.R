#----------------
#   Liver
#----------------

src.endoderm.liver.ext_v1.1$cluster.extract.v1.0.co =NA
src.endoderm.liver.ext_v1.1@meta.data[(!src.endoderm.liver.ext_v1.1$cluster.extract.v1.0.al12%in%NA)&
                                        (!src.endoderm.liver.ext_v1.1$cluster.extract.v1.0.fg4%in%NA),]$cluster.extract.v1.0.co = "Co_trace"
DimPlot(src.endoderm.liver.ext_v1.1, group.by = "cluster.extract.v1.0.co")
DimPlot(src.endoderm.liver.ext_v1.1[,!src.endoderm.liver.ext_v1.1$cluster.extract.v1.0.co%in%NA], group.by = "cluster.extract.v1.0.al12")
DimPlot(src.endoderm.liver.ext_v1.1[,!src.endoderm.liver.ext_v1.1$cluster.extract.v1.0.co%in%NA], group.by = "cluster.extract.v1.0.fg4")

cell_liver_non = 
  rownames(src.endoderm.liver.ext_v1.1@meta.data[!src.endoderm.liver.ext_v1.1$cluster.extract.v1.0.co%in%NA,])[
    src.endoderm.liver.ext_v1.1@meta.data[!src.endoderm.liver.ext_v1.1$cluster.extract.v1.0.co%in%NA,"cluster.extract.v1.0.fg3"] !=
      src.endoderm.liver.ext_v1.1@meta.data[!src.endoderm.liver.ext_v1.1$cluster.extract.v1.0.co%in%NA,"cluster.extract.v1.0.fg4"]]
# character(0)

src.endoderm.liver.ext_v1.1$cluster.extract.v1.0.merge = NA
src.endoderm.liver.ext_v1.1@meta.data[src.endoderm.liver.ext_v1.1$cluster.extract.v1.0.co%in%NA&
                                        !src.endoderm.liver.ext_v1.1$cluster.extract.v1.0.al12%in%NA,]$cluster.extract.v1.0.merge =
  src.endoderm.liver.ext_v1.1@meta.data[src.endoderm.liver.ext_v1.1$cluster.extract.v1.0.co%in%NA&
                                          !src.endoderm.liver.ext_v1.1$cluster.extract.v1.0.al12%in%NA,]$cluster.extract.v1.0.al12
src.endoderm.liver.ext_v1.1@meta.data[src.endoderm.liver.ext_v1.1$cluster.extract.v1.0.co%in%NA&
                                        !src.endoderm.liver.ext_v1.1$cluster.extract.v1.0.fg4%in%NA,]$cluster.extract.v1.0.merge =
  src.endoderm.liver.ext_v1.1@meta.data[src.endoderm.liver.ext_v1.1$cluster.extract.v1.0.co%in%NA&
                                          !src.endoderm.liver.ext_v1.1$cluster.extract.v1.0.fg4%in%NA,]$cluster.extract.v1.0.fg4
src.endoderm.liver.ext_v1.1@meta.data[!src.endoderm.liver.ext_v1.1$cluster.extract.v1.0.co%in%NA,]$cluster.extract.v1.0.merge = 'FG.4-AL.1/2'


src.endoderm.liver.ext_v1.1 = FindVariableFeatures(src.endoderm.liver.ext_v1.1, nfeatures = 2000)
src.endoderm.liver.ext_v1.1 = ScaleData(src.endoderm.liver.ext_v1.1, split.by = "batch_phase",
                                        features = rownames(src.endoderm.liver.ext_v1.1))
src.endoderm.liver.ext_v1.1.filtergene = 
  Myfilter(as.matrix(src.endoderm.liver.ext_v1.1@assays$RNA@data),
           gene = src.endoderm.liver.ext_v1.1@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.endoderm.liver.ext.v1.1.selectgene = 
  src.endoderm.liver.ext_v1.1.filtergene
rm(src.endoderm.liver.ext_v1.1.filtergene)

#=============================
# Row
pdf("figure.v08.07/try.pdf")
#-----------------
liver.ext.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.liver.ext_v1.1@assays$RNA@data[
    src.endoderm.liver.ext.v1.1.selectgene,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.endoderm.liver.ext_v1.1$Time,colors.time),
      MyName2Col(src.endoderm.liver.ext_v1.1$cluster.extract.v1.1..fg4,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.liver.ext_v1.1$cluster.extract.v1.1..al12,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.liver.ext_v1.1$cluster.extract.v1.1..al3,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.liver.ext_v1.1$cluster.v06.26.re..fg4,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.liver.ext_v1.1$cluster.v06.26.re..al12,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.liver.ext_v1.1$cluster.v06.26.re..al3,cluster.endoderm.color.v5)
    ),
    ColSideColorsSize = 4.5,
    c.hc.method = "ward.D",
    r.hc.method = "ward.D2",
    #Rowv = "none",
    return.tree = "row",
    graph = T)
#-------------------
dev.off()

pdf("figure.v08.07/try.tree.pdf",100,20)
plot(liver.ext.re.rowtree)
dev.off()

liver.ext.re.rowtree_1 = as.dendrogram(liver.ext.re.rowtree)
src.endoderm.liver.ext_v1.1.gene = setdiff(
  src.endoderm.liver.ext.v1.1.selectgene,
  c(labels(liver.ext.re.rowtree_1[[1]][[2]])))

# Extract:  V1.0
src.endoderm.liver.ext_v1.1.gene.fin = 
  rownames(src.endoderm.liver.ext.re@reductions$pca@feature.loadings)

#--------------------------------------------------------------------------------
src.endoderm.liver.ext_v1.1 = 
  seurat_mnn(seurat_name = "src.endoderm.liver.ext_v1.1", dims = 1:30, 
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.liver.ext_v1.1.gene.fin)
src.endoderm.liver.ext_v1.1 = 
  seurat_int(seurat_name = "src.endoderm.liver.ext_v1.1", dims = 1:30,
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.liver.ext_v1.1.gene.fin)
#-------------------------------------------------


# Neighbor Cluster
#-------------------
src.endoderm.liver.ext_v1.1 = 
  FindNeighbors(src.endoderm.liver.ext_v1.1, dims=1:30, 
                reduction = "pca_integrated", assay = "integrated")
src.endoderm.liver.ext_v1.1 = 
  FindNeighbors(src.endoderm.liver.ext_v1.1, dims=1:30, 
                reduction = "pca", assay = "RNA")
src.endoderm.liver.ext_v1.1 = 
  FindClusters(src.endoderm.liver.ext_v1.1, 
               resolution = 2, graph.name = "RNA_snn",)

src.endoderm.liver.ext_v1.1 = 
  FindNeighbors(src.endoderm.liver.ext_v1.1, dims=1:3, 
                reduction = "umap_mnn", assay = "RNA")
src.endoderm.liver.ext_v1.1 = 
  FindClusters(src.endoderm.liver.ext_v1.1, 
               resolution = 1.5, graph.name = "RNA_snn",)
#-------------------

pdf("figure.v08.07/try.pdf")
#-----------------------------
DimPlot(src.endoderm.liver.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.2', reduction = "umap_mnn")
DimPlot(src.endoderm.liver.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.1.5', reduction = "umap_mnn")

for(red in c("umap","umap_mnn","umap_integrated")){
  for(type in c('Time',
                "cluster.extract.v1.1..fg4",
                "cluster.extract.v1.1..al12","cluster.extract.v1.1..al3",
                "cluster.v06.26.re..fg4",
                "cluster.v06.26.re..al12","cluster.v06.26.re..al3")){
    print(
      DimPlot(src.endoderm.liver.ext_v1.1, group.by = type, pt.size = 1.2,
              reduction = red, cols = colors_bg) + p_add_leg
    )
  }
}
#-----------------------------
dev.off()

# FDL
src.endoderm.liver.ext_v1.1 = 
  seurat_fdl("src.endoderm.liver.ext_v1.1", "pca","RNA")
src.endoderm.liver.ext_v1.1 = 
  seurat_fdl("src.endoderm.liver.ext_v1.1", "mnn","RNA")
src.endoderm.liver.ext_v1.1 = 
  seurat_fdl("src.endoderm.liver.ext_v1.1", "pca_integrated","integrated")

# Merge cell type & tracing code
#------------------------------------

src.endoderm.liver.ext_v1.1$cluster.extract.v1.1..check = paste(
  src.endoderm.liver.ext_v1.1$cluster.v06.26.re..fg4,
  src.endoderm.liver.ext_v1.1$cluster.v06.26.re..al12,
  src.endoderm.liver.ext_v1.1$cluster.v06.26.re..al3,  sep="_")
table(src.endoderm.liver.ext_v1.1$cluster.extract.v1.1..check)

cell_ext_liver.rm = rownames(src.endoderm.liver.ext_v1.1@meta.data[
  src.endoderm.liver.ext_v1.1$cluster.extract.v1.1..check%in%c(
    "Liver_AL.2_Live",'Liver_Liver_EHBD',
    "Liver_Liver_Small.intestine.1",
    "NA_Small.intestine.1_Small.intestine.2",
    "Small.intestine.2_Small.intestine.2_Stomach")|
    
    src.endoderm.liver.ext_v1.1$cluster.v06.26.re..fg4%in%c(
      "FG.3/4")|
    
    src.endoderm.liver.ext_v1.1$cluster.v06.26.re..al3%in%c(
      "AL.3−MG.2"),])
#------------------------------------

src.endoderm.liver.ext_v1.1.re = 
  src.endoderm.liver.ext_v1.1[,!colnames(src.endoderm.liver.ext_v1.1)%in%cell_ext_liver.rm]

src.endoderm.liver.ext_v1.1.re$cluster.v06.26.re..merge = NA
src.endoderm.liver.ext_v1.1.re@meta.data[
  (!src.endoderm.liver.ext_v1.1.re$cluster.v06.26.re..fg4%in%NA)&
    (src.endoderm.liver.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..merge =
  src.endoderm.liver.ext_v1.1.re@meta.data[
    (!src.endoderm.liver.ext_v1.1.re$cluster.v06.26.re..fg4%in%NA)&
      (src.endoderm.liver.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..fg4
src.endoderm.liver.ext_v1.1.re@meta.data[
  (!src.endoderm.liver.ext_v1.1.re$cluster.v06.26.re..al12%in%NA)&
    (src.endoderm.liver.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..merge =
  src.endoderm.liver.ext_v1.1.re@meta.data[
    (!src.endoderm.liver.ext_v1.1.re$cluster.v06.26.re..al12%in%NA)&
      (src.endoderm.liver.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..al12
src.endoderm.liver.ext_v1.1.re@meta.data[
  (!src.endoderm.liver.ext_v1.1.re$cluster.v06.26.re..al3%in%NA)&
    (src.endoderm.liver.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..merge =
  src.endoderm.liver.ext_v1.1.re@meta.data[
    (!src.endoderm.liver.ext_v1.1.re$cluster.v06.26.re..al3%in%NA)&
      (src.endoderm.liver.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..al3


src.endoderm.liver.ext_v1.1.re$cluster.extract.v1.1..merge = NA
src.endoderm.liver.ext_v1.1.re@meta.data[
  (!src.endoderm.liver.ext_v1.1.re$cluster.extract.v1.1..fg4%in%NA)&
    (src.endoderm.liver.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..merge =
  src.endoderm.liver.ext_v1.1.re@meta.data[
    (!src.endoderm.liver.ext_v1.1.re$cluster.extract.v1.1..fg4%in%NA)&
      (src.endoderm.liver.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..fg4
src.endoderm.liver.ext_v1.1.re@meta.data[
  (!src.endoderm.liver.ext_v1.1.re$cluster.extract.v1.1..al12%in%NA)&
    (src.endoderm.liver.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..merge =
  src.endoderm.liver.ext_v1.1.re@meta.data[
    (!src.endoderm.liver.ext_v1.1.re$cluster.extract.v1.1..al12%in%NA)&
      (src.endoderm.liver.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..al12
src.endoderm.liver.ext_v1.1.re@meta.data[
  (!src.endoderm.liver.ext_v1.1.re$cluster.extract.v1.1..al3%in%NA)&
    (src.endoderm.liver.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..merge =
  src.endoderm.liver.ext_v1.1.re@meta.data[
    (!src.endoderm.liver.ext_v1.1.re$cluster.extract.v1.1..al3%in%NA)&
      (src.endoderm.liver.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..al3

#--------------------------------
src.endoderm.liver.ext_v1.1.re = 
  seurat_mnn(seurat_name = "src.endoderm.liver.ext_v1.1.re", dims = 1:30, 
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.liver.ext_v1.1.gene.fin)
src.endoderm.liver.ext_v1.1.re = 
  seurat_int(seurat_name = "src.endoderm.liver.ext_v1.1.re", dims = 1:30,
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.liver.ext_v1.1.gene.fin)

src.endoderm.liver.ext_v1.1.re@reductions$pca_integrated@cell.embeddings =
  src.endoderm.liver.ext_v1.1.re@reductions$pca_integrated@cell.embeddings[colnames(src.endoderm.liver.ext_v1.1.re),]
src.endoderm.liver.ext_v1.1.re@reductions$mnn@cell.embeddings =
  src.endoderm.liver.ext_v1.1.re@reductions$mnn@cell.embeddings[colnames(src.endoderm.liver.ext_v1.1.re),]
src.endoderm.liver.ext_v1.1.re@reductions$umap_mnn@cell.embeddings =
  src.endoderm.liver.ext_v1.1.re@reductions$umap_mnn@cell.embeddings[colnames(src.endoderm.liver.ext_v1.1.re),]
src.endoderm.liver.ext_v1.1.re@reductions$umap_integrated@cell.embeddings =
  src.endoderm.liver.ext_v1.1.re@reductions$umap_integrated@cell.embeddings[colnames(src.endoderm.liver.ext_v1.1.re),]


src.endoderm.liver.ext_v1.1.re = 
  seurat_fdl("src.endoderm.liver.ext_v1.1.re", "pca","RNA")
src.endoderm.liver.ext_v1.1.re = 
  seurat_fdl("src.endoderm.liver.ext_v1.1.re", "mnn","RNA")
src.endoderm.liver.ext_v1.1.re = 
  seurat_fdl("src.endoderm.liver.ext_v1.1.re", "pca_integrated","integrated")
#--------------------------------


seurat_GCN = function(seurat_name, seurat.selectgene,
                      reduction = "pca",assay = "RNA",dims = 1:30){
  seurat = get(seurat_name)
  seurat = FindNeighbors(seurat, dims = dims, reduction = reduction, assays = assay)
  cor.data = as.matrix(seurat@graphs[[paste(assay,"_snn",sep="")]])
  cor.data = as(cor.data,"dgCMatrix")
  graph.data  = igraph::graph.adjacency(cor.data,mode = "undirected",weighted = T)
  graph.data = igraph::simplify(graph.data,remove.multiple = T,remove.loops = T)
  set.seed(1)
  
  return(graph.data)
}
for(red in c("mnn",
             "pca_integrated")){
  
  pdf(paste("figure.v08.07/organ_development_merge_re/liver/liver_summary_graph.re_snn_",red,".pdf",sep = ""),12,7)
  
  assay = ifelse(red=="mnn","RNA","integrated")
  
  graph.liver.ext_v1.1 = 
    seurat_GCN(seurat_name = "src.endoderm.liver.ext_v1.1.re",
               src.endoderm.liver.ext_v1.1.gene.fin,
               reduction = red, assay =  assay)
  
  graph.liver_cluster = cluster_walktrap(graph.liver.ext_v1.1, steps = 1)
  
  for(color in c("cluster.extract.v1.1..merge",
                 "cluster.v06.26.re..merge",
                 "Time")){
    for(edge.color in c(NA,"lightgray")){
      for(edge.color in c(NA,"lightgray")){
        set.seed(1)
        plot.igraph(graph.liver.ext_v1.1,
                    layout = layout_with_fr(graph.liver.ext_v1.1),
                    edge.color = edge.color,
                    vertex.size=2.5, vertex.label=NA,
                    vertex.label.color = "black", 
                    vertex.frame.color = NA, vertex.frame.width = 0.5,
                    #vertex.color =  colors.time[src.endoderm.liver.ext_v1.1[["Time"]][graph.liver_cluster$names,]],
                    #vertex.color =  cluster.endoderm.color.v5[src.endoderm.liver.ext_v1.1[["cluster.v06.26.re"]][graph.liver_cluster$names,]],
                    vertex.color =  colors_bg[src.endoderm.liver.ext_v1.1.re[[color]][graph.liver_cluster$names,]])
      }}
  }
  dev.off()
  
}


pdf("figure.v08.07/organ_development_merge_re/liver/liver_summary_graph.re_other.pdf",12,7)
# Reduction
for(red in c("umap_integrated","umap_mnn",
             "fdl_pca",'fdl_mnn',"fdl_pca_integrated")){
  
  for(color in c("cluster.extract.v1.1..merge",
                 "cluster.extract.v1.1..fg4", "cluster.extract.v1.1..al12",
                 "cluster.extract.v1.1..al3",
                 "cluster.v06.26.re..merge",
                 "cluster.v06.26.re..fg4", "cluster.v06.26.re..al12",
                 "cluster.v06.26.re..al3",
                 "Time")){
    print(
      DimPlot(src.endoderm.liver.ext_v1.1.re, group.by = color, 
              reduction = red, cols = colors_bg,
              pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
    )
  }
}
dev.off() 


pdf("figure.v08.07/organ_development_merge_re/liver/liver_summary_proportion.re.pdf",12,7)
#-----------------------------
cell_type = c("FG.4","AL.1","AL.2","AL.3","FG.4-Liver","AL.1/2-Liver","AL.3-Liver","Liver")
ext_type = c("FG.4","AL.1","AL.2","AL.1/2",'AL.3',"FG.4-AL.1/2/3")

meta_liver = src.endoderm.liver.ext_v1.1.re@meta.data
meta_liver$Time = factor(meta_liver$Time, levels = names(colors.time))
meta_liver$cluster.extract.v1.1..merge = factor(meta_liver$cluster.extract.v1.1..merge, levels = rev(ext_type))
meta_liver$cluster.extract.v1.1..fg4 = factor(meta_liver$cluster.extract.v1.1..fg4, levels = rev(ext_type))
meta_liver$cluster.extract.v1.1..al12 = factor(meta_liver$cluster.extract.v1.1..al12, levels = rev(ext_type))
meta_liver$cluster.extract.v1.1..al3 = factor(meta_liver$cluster.extract.v1.1..al3, levels = rev(ext_type))

meta_liver$cluster.v06.26.re..merge = factor(meta_liver$cluster.v06.26.re..merge, levels = rev(cell_type))
meta_liver$cluster.v06.26.re..fg4 = factor(meta_liver$cluster.v06.26.re..fg4, levels = rev(cell_type))
meta_liver$cluster.v06.26.re..al12 = factor(meta_liver$cluster.v06.26.re..al12, levels = rev(cell_type))
meta_liver$cluster.v06.26.re..al3 = factor(meta_liver$cluster.v06.26.re..al3, levels = rev(cell_type))

meta_liver_B0 = meta_liver
meta_liver_B1 = meta_liver[meta_liver$batch%in%1&!meta_liver$Time%in%"ss9",]
meta_liver_B2 = meta_liver[meta_liver$batch%in%2,]

for(data in c("meta_liver_B0","meta_liver_B1","meta_liver_B2")){
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
                                 group = cluster.extract.v1.1..merge, fill = cluster.extract.v1.1..merge),
                   stat = "count", position = 'fill') + 
          xlab("Time") + ylab("Cluster proportion of
 tracing code")+ ggtitle(data))
  
  print(p + 
          geom_bar(data = get(data), 
                   mapping = aes(x = Time,
                                 group = cluster.v06.26.re..merge, fill = cluster.v06.26.re..merge),
                   stat = "count", position = 'fill') + 
          xlab("Time") + ylab("Cluster proportion of
 cell type")+ ggtitle(data))
}

#-----------------------------
dev.off()

save(src.endoderm.liver.ext_v1.1,
     src.endoderm.liver.ext_v1.1.re,
     file = "figure.v08.07/organ_development_merge_re/liver/src.endoderm.liver.ext_v1.1.re.Rdata")


metadata.liver.ext_v1.1.re = cbind(
  src.endoderm.liver.ext_v1.1.re@meta.data,
  src.endoderm.liver.ext_v1.1.re@reductions$umap_integrated@cell.embeddings[colnames(src.endoderm.liver.ext_v1.1.re),],
  src.endoderm.liver.ext_v1.1.re@reductions$fdl_pca_integrated@cell.embeddings[colnames(src.endoderm.liver.ext_v1.1.re),])
save(metadata.liver.ext_v1.1.re,
     file = "figure.v08.07/organ_development_merge_re/liver/metadata.liver.ext_v1.1.re_integrated.Rdata")

metadata.liver.ext_v1.1.re = cbind(
  src.endoderm.liver.ext_v1.1.re@meta.data,
  src.endoderm.liver.ext_v1.1.re@reductions$umap_mnn@cell.embeddings[colnames(src.endoderm.liver.ext_v1.1.re),],
  src.endoderm.liver.ext_v1.1.re@reductions$fdl_mnn@cell.embeddings[colnames(src.endoderm.liver.ext_v1.1.re),])
save(metadata.liver.ext_v1.1.re,
     file = "figure.v08.07/organ_development_merge_re/liver/metadata.liver.ext_v1.1.re_mnn.Rdata")

