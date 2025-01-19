#-------------------------------------
# stomach
#-------------------------------------
cell_mg3_pan_mg3m = rownames(src.mg3.integrated@meta.data[
  src.mg3.integrated$cluster.v06.26.re%in%"MG.3.M",])
intersect(colnames(src.endoderm.sto.ext_v1.1.re), cell_mg3_pan_mg3m)

src.endoderm.sto.ext_v1.1.refine = merge(
  src.endoderm.sto.ext_v1.1.re,
  src.mg3.integrated[,cell_mg3_pan_mg3m])
src.endoderm.sto.ext_v1.1.refine = FindVariableFeatures(src.endoderm.sto.ext_v1.1.refine, nfeatures = 2000)
src.endoderm.sto.ext_v1.1.refine = ScaleData(src.endoderm.sto.ext_v1.1.refine, split.by = "batch_phase",
                                             features = rownames(src.endoderm.sto.ext_v1.1.refine))
src.endoderm.sto.ext.v1.1.selectgene = src.endoderm.sto.ext_v1.1.gene.fin
src.endoderm.sto.ext_v1.1.refine.gene.fin = src.endoderm.sto.ext_v1.1.gene.fin

src.endoderm.sto.ext_v1.1.refine = RunPCA(src.endoderm.sto.ext_v1.1.refine, 
                                          features = src.endoderm.sto.ext_v1.1.refine.gene.fin)
src.endoderm.sto.ext_v1.1.refine = RunUMAP(src.endoderm.sto.ext_v1.1.refine, dims = 1:30, 
                                           reduction = "pca", n.components = 3)
DimPlot(src.endoderm.sto.ext_v1.1.refine, reduction = "umap", group.by = "cluster.v06.26.re..merge")

#--------------------------------------------------------------------------------
src.endoderm.sto.ext_v1.1.refine = 
  seurat_mnn(seurat_name = "src.endoderm.sto.ext_v1.1.refine",
             dims = 1:30,
             n.neighbors = 100,
             n.components = 3,
             seurat.selectgene = src.endoderm.sto.ext_v1.1.refine.gene.fin)

src.endoderm.sto.ext_v1.1.refine = 
  seurat_int(seurat_name = "src.endoderm.sto.ext_v1.1.refine", dims = 1:30,
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.sto.ext_v1.1.refine.gene.fin)

#-------------------------------------------------
src.endoderm.sto.ext_v1.1.refine = seurat_fdl(src.endoderm.sto.ext_v1.1.refine, "pca","RNA")
src.endoderm.sto.ext_v1.1.refine = seurat_fdl(src.endoderm.sto.ext_v1.1.refine, "mnn","RNA")
src.endoderm.sto.ext_v1.1.refine = seurat_fdl(src.endoderm.sto.ext_v1.1.refine, "pca_integrated","integrated")
src.endoderm.sto.ext_v1.1.refine@meta.data[cell_mg3_pan_mg3a,]$cluster.v06.26.re..merge = "MG.3.M"

pdf("figure.v08.07/organ_development_re_v240115/try.pancreas.re.pdf", 9, 7)
for(reduction in names(src.endoderm.sto.ext_v1.1.refine@reductions)){
  print(
    DimPlot(src.endoderm.sto.ext_v1.1.refine, reduction = reduction, 
            group.by = "cluster.v06.26.re..merge",
            cols = cluster.endoderm.color.v5) +
      theme_classic() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
      theme(legend.position = "none") +
      ggtitle(label = reduction)
  )
}
dev.off()

seurat = RunUMAP(src.endoderm.sto.ext_v1.1.refine, reduction = "umap_integrated",
                 dims = 1:30, n.neighbors = 100, n.components = 3)
DimPlot(seurat, reduction = "umap_integrated", c(1,2),
        group.by = "cluster.v06.26.re..merge",
        cols = cluster.endoderm.color.v5)
src.endoderm.sto.ext_v1.1.refine@reductions$umap_integrated@cell.embeddings =
  seurat@reductions$umap_integrated@cell.embeddings

save(src.endoderm.sto.ext_v1.1.refine,
     file = "figure.v08.07/organ_development_re_v240115/src.endoderm.sto.ext_v1.1.refine.Rdata")


#----------------
#     Remake
#----------------
src.endoderm.sm1.ext_v1.1 = FindVariableFeatures(src.endoderm.sm1.ext_v1.1, nfeatures = 2000)
src.endoderm.sm1.ext_v1.1 = ScaleData(src.endoderm.sm1.ext_v1.1, split.by = "batch_phase",
                                      features = rownames(src.endoderm.sm1.ext_v1.1))
src.endoderm.sm1.ext_v1.1.filtergene = 
  Myfilter(as.matrix(src.endoderm.sm1.ext_v1.1@assays$RNA@data),
           gene = src.endoderm.sm1.ext_v1.1@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.endoderm.sm1.ext.v1.1.selectgene = 
  src.endoderm.sm1.ext_v1.1.filtergene
rm(src.endoderm.sm1.ext_v1.1.filtergene)

#=============================
# Row
pdf("figure.v08.07/try.pdf")
#-----------------
sm1.ext.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.sm1.ext_v1.1@assays$RNA@data[
    src.endoderm.sm1.ext.v1.1.selectgene,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.endoderm.sm1.ext_v1.1$Time,colors.time),
      MyName2Col(src.endoderm.sm1.ext_v1.1$cluster.extract.v1.1..al3,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.sm1.ext_v1.1$cluster.extract.v1.1..mg1,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.sm1.ext_v1.1$cluster.extract.v1.1..mg2,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.sm1.ext_v1.1$cluster.extract.v1.1..mg3,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.sm1.ext_v1.1$cluster.v06.26.re..al3,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.sm1.ext_v1.1$cluster.v06.26.re..mg1,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.sm1.ext_v1.1$cluster.v06.26.re..mg2,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.sm1.ext_v1.1$cluster.v06.26.re..mg3,cluster.endoderm.color.v5)
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
plot(sm1.ext.re.rowtree)
dev.off()

sm1.ext.re.rowtree_1 = as.dendrogram(sm1.ext.re.rowtree)
src.endoderm.sm1.ext_v1.1.gene = setdiff(
  src.endoderm.sm1.ext.v1.1.selectgene,
  c(labels(sm1.ext.re.rowtree_1[[2]][[1]]),
    labels(sm1.ext.re.rowtree_1[[2]][[2]][[1]])))

# Extract:  V1.0
src.endoderm.sm1.ext_v1.1.gene.fin = 
  src.endoderm.sm1.ext.re.selectgene.fin 

#--------------------------------------------------------------------------------
src.endoderm.sm1.ext_v1.1 = 
  seurat_mnn(seurat_name = "src.endoderm.sm1.ext_v1.1", dims = 1:30, 
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.sm1.ext_v1.1.gene.fin)
src.endoderm.sm1.ext_v1.1 = 
  seurat_int(seurat_name = "src.endoderm.sm1.ext_v1.1", dims = 1:30,
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.sm1.ext_v1.1.gene.fin)
#-------------------------------------------------


# Neighbor Cluster
#-------------------
src.endoderm.sm1.ext_v1.1 = 
  FindNeighbors(src.endoderm.sm1.ext_v1.1, dims=1:30, 
                reduction = "pca_integrated", assay = "integrated")
src.endoderm.sm1.ext_v1.1 = 
  FindNeighbors(src.endoderm.sm1.ext_v1.1, dims=1:30, 
                reduction = "pca", assay = "RNA")
src.endoderm.sm1.ext_v1.1 = 
  FindClusters(src.endoderm.sm1.ext_v1.1, 
               resolution = 2, graph.name = "RNA_snn",)

src.endoderm.sm1.ext_v1.1 = 
  FindNeighbors(src.endoderm.sm1.ext_v1.1, dims=1:3, 
                reduction = "umap_mnn", assay = "RNA")
src.endoderm.sm1.ext_v1.1 = 
  FindClusters(src.endoderm.sm1.ext_v1.1, 
               resolution = 1.5, graph.name = "RNA_snn",)
#-------------------

pdf("figure.v08.07/try.pdf")
#-----------------------------
DimPlot(src.endoderm.sm1.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.2', reduction = "umap_mnn")
DimPlot(src.endoderm.sm1.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.1.5', reduction = "umap_mnn")

for(red in c("umap","umap_mnn","umap_integrated")){
  for(type in c('Time',
                "cluster.extract.v1.1..al3","cluster.extract.v1.1..mg1",
                "cluster.extract.v1.1..mg2","cluster.extract.v1.1..mg3",
                "cluster.v06.26.re..al3","cluster.v06.26.re..mg1",
                "cluster.v06.26.re..mg2","cluster.v06.26.re..mg3")){
    print(
      DimPlot(src.endoderm.sm1.ext_v1.1, group.by = type, pt.size = 1.2,
              reduction = red, cols = colors_bg) + p_add_leg
    )
  }
}
#-----------------------------
dev.off()

# FDL
src.endoderm.sm1.ext_v1.1 = 
  seurat_fdl("src.endoderm.sm1.ext_v1.1", "pca","RNA")
src.endoderm.sm1.ext_v1.1 = 
  seurat_fdl("src.endoderm.sm1.ext_v1.1", "mnn","RNA")
src.endoderm.sm1.ext_v1.1 = 
  seurat_fdl("src.endoderm.sm1.ext_v1.1", "pca_integrated","integrated")

# Merge cell type & tracing code
#------------------------------------
src.endoderm.sm1.ext_v1.1@meta.data[
  src.endoderm.sm1.ext_v1.1$cluster.v06.26.re..al3%in%"AL.3-EHBD/VP",]

src.endoderm.sm1.ext_v1.1$cluster.extract.v1.1..check = paste(
  src.endoderm.sm1.ext_v1.1$cluster.v06.26.re..al3,
  src.endoderm.sm1.ext_v1.1$cluster.v06.26.re..mg1,
  src.endoderm.sm1.ext_v1.1$cluster.v06.26.re..mg2,table(src.endoderm.sm1.ext_v1.1$cluster.extract.v1.1..check)
  src.endoderm.sm1.ext_v1.1$cluster.v06.26.re..mg3,  sep="_")


cell_ext_sm1.rm = rownames(src.endoderm.sm1.ext_v1.1@meta.data[
  
  src.endoderm.sm1.ext_v1.1$cluster.extract.v1.1..check%in%c(
    "AL.3-EHBD/VP_NA_Small.intestine.1_NA ","AL.3-Liver_NA_Small.intestine.1_",
    "EHBD_Small.intestine.1_Small.intestine.1_Small.intestine.2",
    "Liver_NA_Small.intestine.1_NA","NA_NA_MG.2_Small.intestine.1",
    "NA_NA_MG.2_Small.intestine.2 ",
    "NA_NA_Small.intestine.1_Small.intestine.2",
    "NA_Small.intestine.1_Small.intestine.1_DP",
    "SSmall.intestine.1_Small.intestine.1_Small.intestine.1_Small.intestine.2",
    "Small.intestine.1_Small.intestine.1_Small.intestine.1_Stomach","VP_NA_Small.intestine.1_NA",
    "VP_Small.intestine.1_Small.intestine.1_DP ",
    "VP_Small.intestine.1_Small.intestine.1_Stomach")|
    
    src.endoderm.sm1.ext_v1.1$cluster.v06.26.re..al3%in%c(
      "AL.3-EHBD/VP","AL.3-Liver","EHBD","VP","Liver")|
    
    src.endoderm.sm1.ext_v1.1$cluster.v06.26.re..mg3%in%c(
      "MG.3.A","DP","Small.intestine.2","Stomach"),])
#------------------------------------

src.endoderm.sm1.ext_v1.1.re = 
  src.endoderm.sm1.ext_v1.1[,!colnames(src.endoderm.sm1.ext_v1.1)%in%cell_ext_sm1.rm]

src.endoderm.sm1.ext_v1.1.re = 
  src.endoderm.sm1.ext_v1.1.re[,!(src.endoderm.sm1.ext_v1.1.re$cluster.extract.v1.1..al3%in%"FG.4-AL.1/2/3")]

src.endoderm.sm1.ext_v1.1.re$cluster.v06.26.re..merge = NA
src.endoderm.sm1.ext_v1.1.re@meta.data[
  (!src.endoderm.sm1.ext_v1.1.re$cluster.v06.26.re..al3%in%NA)&
    (src.endoderm.sm1.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..merge =
  src.endoderm.sm1.ext_v1.1.re@meta.data[
    (!src.endoderm.sm1.ext_v1.1.re$cluster.v06.26.re..al3%in%NA)&
      (src.endoderm.sm1.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..al3
src.endoderm.sm1.ext_v1.1.re@meta.data[
  (!src.endoderm.sm1.ext_v1.1.re$cluster.v06.26.re..mg1%in%NA)&
    (src.endoderm.sm1.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..merge =
  src.endoderm.sm1.ext_v1.1.re@meta.data[
    (!src.endoderm.sm1.ext_v1.1.re$cluster.v06.26.re..mg1%in%NA)&
      (src.endoderm.sm1.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..mg1
src.endoderm.sm1.ext_v1.1.re@meta.data[
  (!src.endoderm.sm1.ext_v1.1.re$cluster.v06.26.re..mg2%in%NA)&
    (src.endoderm.sm1.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..merge =
  src.endoderm.sm1.ext_v1.1.re@meta.data[
    (!src.endoderm.sm1.ext_v1.1.re$cluster.v06.26.re..mg2%in%NA)&
      (src.endoderm.sm1.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..mg2
src.endoderm.sm1.ext_v1.1.re@meta.data[
  (!src.endoderm.sm1.ext_v1.1.re$cluster.v06.26.re..mg3%in%NA)&
    (src.endoderm.sm1.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..merge =
  src.endoderm.sm1.ext_v1.1.re@meta.data[
    (!src.endoderm.sm1.ext_v1.1.re$cluster.v06.26.re..mg3%in%NA)&
      (src.endoderm.sm1.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..mg3


src.endoderm.sm1.ext_v1.1.re$cluster.extract.v1.1..merge = NA
src.endoderm.sm1.ext_v1.1.re@meta.data[
  (!src.endoderm.sm1.ext_v1.1.re$cluster.extract.v1.1..al3%in%NA)&
    (src.endoderm.sm1.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..merge =
  src.endoderm.sm1.ext_v1.1.re@meta.data[
    (!src.endoderm.sm1.ext_v1.1.re$cluster.extract.v1.1..al3%in%NA)&
      (src.endoderm.sm1.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..al3
src.endoderm.sm1.ext_v1.1.re@meta.data[
  (!src.endoderm.sm1.ext_v1.1.re$cluster.extract.v1.1..mg1%in%NA)&
    (src.endoderm.sm1.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..merge =
  src.endoderm.sm1.ext_v1.1.re@meta.data[
    (!src.endoderm.sm1.ext_v1.1.re$cluster.extract.v1.1..mg1%in%NA)&
      (src.endoderm.sm1.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..mg1
src.endoderm.sm1.ext_v1.1.re@meta.data[
  (!src.endoderm.sm1.ext_v1.1.re$cluster.extract.v1.1..mg2%in%NA)&
    (src.endoderm.sm1.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..merge =
  src.endoderm.sm1.ext_v1.1.re@meta.data[
    (!src.endoderm.sm1.ext_v1.1.re$cluster.extract.v1.1..mg2%in%NA)&
      (src.endoderm.sm1.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..mg2
src.endoderm.sm1.ext_v1.1.re@meta.data[
  (!src.endoderm.sm1.ext_v1.1.re$cluster.extract.v1.1..mg3%in%NA)&
    (src.endoderm.sm1.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..merge =
  src.endoderm.sm1.ext_v1.1.re@meta.data[
    (!src.endoderm.sm1.ext_v1.1.re$cluster.extract.v1.1..mg3%in%NA)&
      (src.endoderm.sm1.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..mg3

#--------------------------------
src.endoderm.sm1.ext_v1.1.re = 
  seurat_mnn(seurat_name = "src.endoderm.sm1.ext_v1.1.re", dims = 1:30, 
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.sm1.ext_v1.1.gene.fin)
src.endoderm.sm1.ext_v1.1.re = 
  seurat_int(seurat_name = "src.endoderm.sm1.ext_v1.1.re", dims = 1:30,
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.sm1.ext_v1.1.gene.fin)

src.endoderm.sm1.ext_v1.1.re@reductions$pca_integrated@cell.embeddings =
  src.endoderm.sm1.ext_v1.1.re@reductions$pca_integrated@cell.embeddings[colnames(src.endoderm.sm1.ext_v1.1.re),]
src.endoderm.sm1.ext_v1.1.re@reductions$mnn@cell.embeddings =
  src.endoderm.sm1.ext_v1.1.re@reductions$mnn@cell.embeddings[colnames(src.endoderm.sm1.ext_v1.1.re),]

src.endoderm.sm1.ext_v1.1.re = 
  seurat_fdl("src.endoderm.sm1.ext_v1.1.re", "pca","RNA")
src.endoderm.sm1.ext_v1.1.re = 
  seurat_fdl("src.endoderm.sm1.ext_v1.1.re", "mnn","RNA")
src.endoderm.sm1.ext_v1.1.re = 
  seurat_fdl("src.endoderm.sm1.ext_v1.1.re", "pca_integrated","integrated")
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


pdf("figure.v08.07/organ_development_merge_re/sm1/sm1_summary_graph.re_other.pdf",12,7)
# Reduction
for(red in c("umap_integrated","umap_mnn",
             "fdl_pca",'fdl_mnn',"fdl_pca_integrated")){
  
  for(color in c("cluster.extract.v1.1..merge",
                 "cluster.extract.v1.1..al3","cluster.extract.v1.1..mg1",
                 "cluster.extract.v1.1..mg2","cluster.extract.v1.1..mg3",
                 "cluster.v06.26.re..merge",
                 "cluster.v06.26.re..al3","cluster.v06.26.re..mg1",
                 "cluster.v06.26.re..mg2","cluster.v06.26.re..mg3",
                 "Time")){
    print(
      DimPlot(src.endoderm.sm1.ext_v1.1.re, group.by = color, 
              reduction = red, cols = colors_bg,
              pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
    )
  }
}
dev.off()


pdf("figure.v08.07/organ_development_merge_re/sm1/sm1_summary_proportion.re.pdf",12,7)
#-----------------------------
cell_type = c("AL.3","AL.3-Small.intestine.1","MG.1","MG.2","MG.3","MG.3.P","Small.intestine.1")
ext_type = c("AL.3","MG.1","MG.2","MG.3",'AL.3-MG.2',"MG.2/3","MG.1/2/3","AL.3-MG.1/2/3")

meta_sm1 = src.endoderm.sm1.ext_v1.1.re@meta.data
meta_sm1$Time = factor(meta_sm1$Time, levels = names(colors.time))
meta_sm1$cluster.extract.v1.1..merge = factor(meta_sm1$cluster.extract.v1.1..merge, levels = rev(ext_type))
meta_sm1$cluster.extract.v1.1..al3 = factor(meta_sm1$cluster.extract.v1.1..al3, levels = rev(ext_type))
meta_sm1$cluster.extract.v1.1..mg1 = factor(meta_sm1$cluster.extract.v1.1..mg1, levels = rev(ext_type))
meta_sm1$cluster.extract.v1.1..mg2 = factor(meta_sm1$cluster.extract.v1.1..mg2, levels = rev(ext_type))
meta_sm1$cluster.extract.v1.1..mg3 = factor(meta_sm1$cluster.extract.v1.1..mg3, levels = rev(ext_type))

meta_sm1$cluster.v06.26.re..merge = factor(meta_sm1$cluster.v06.26.re..merge, levels = rev(cell_type))
meta_sm1$cluster.v06.26.re..al3 = factor(meta_sm1$cluster.v06.26.re..al3, levels = rev(cell_type))
meta_sm1$cluster.v06.26.re..mg1 = factor(meta_sm1$cluster.v06.26.re..mg1, levels = rev(cell_type))
meta_sm1$cluster.v06.26.re..mg2 = factor(meta_sm1$cluster.v06.26.re..mg2, levels = rev(cell_type))
meta_sm1$cluster.v06.26.re..mg3 = factor(meta_sm1$cluster.v06.26.re..mg3, levels = rev(cell_type))

meta_sm1_B0 = meta_sm1
meta_sm1_B1 = meta_sm1[meta_sm1$batch%in%1&!meta_sm1$Time%in%"ss9",]
meta_sm1_B2 = meta_sm1[meta_sm1$batch%in%2,]

for(data in c("meta_sm1_B0","meta_sm1_B1","meta_sm1_B2")){
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

save(src.endoderm.sm1.ext_v1.1,
     src.endoderm.sm1.ext_v1.1.re,
     file = "figure.v08.07/organ_development_merge_re/sm1/src.endoderm.sm1.ext_v1.1.re.Rdata")


metadata.sm1.ext_v1.1.re = cbind(
  src.endoderm.sm1.ext_v1.1.re@meta.data,
  src.endoderm.sm1.ext_v1.1.re@reductions$umap_mnn@cell.embeddings[colnames(src.endoderm.sm1.ext_v1.1.re),],
  src.endoderm.sm1.ext_v1.1.re@reductions$fdl_mnn@cell.embeddings[colnames(src.endoderm.sm1.ext_v1.1.re),])
save(metadata.sm1.ext_v1.1.re,
     file = "figure.v08.07/organ_development_merge_re/sm1/metadata.sm1.ext_v1.1.re.Rdata")

#----------------
#     Remake
#----------------
src.endoderm.sm2.ext_v1.1 = FindVariableFeatures(src.endoderm.sm2.ext_v1.1, nfeatures = 2000)
src.endoderm.sm2.ext_v1.1 = ScaleData(src.endoderm.sm2.ext_v1.1, split.by = "batch_phase",
                                      features = rownames(src.endoderm.sm2.ext_v1.1))
src.endoderm.sm2.ext_v1.1.filtergene = 
  Myfilter(as.matrix(src.endoderm.sm2.ext_v1.1@assays$RNA@data),
           gene = src.endoderm.sm2.ext_v1.1@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.endoderm.sm2.ext.v1.1.selectgene = 
  src.endoderm.sm2.ext_v1.1.filtergene
rm(src.endoderm.sm2.ext_v1.1.filtergene)

#=============================
# Row
pdf("figure.v08.07/try.pdf")
#-----------------
sm2.ext.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.sm2.ext_v1.1@assays$RNA@data[
    src.endoderm.sm2.ext.v1.1.selectgene,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.endoderm.sm2.ext_v1.1$Time,colors.time),
      MyName2Col(src.endoderm.sm2.ext_v1.1$cluster.extract.v1.1..hg1,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.sm2.ext_v1.1$cluster.extract.v1.1..mg2,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.sm2.ext_v1.1$cluster.extract.v1.1..mg3,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.sm2.ext_v1.1$cluster.v06.26.re..hg1,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.sm2.ext_v1.1$cluster.v06.26.re..mg2,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.sm2.ext_v1.1$cluster.v06.26.re..mg3,cluster.endoderm.color.v5)
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
plot(sm2.ext.re.rowtree)
dev.off()

sm2.ext.re.rowtree_1 = as.dendrogram(sm2.ext.re.rowtree)
src.endoderm.sm2.ext_v1.1.gene = setdiff(
  src.endoderm.sm2.ext.v1.1.selectgene,
  c(labels(sm2.ext.re.rowtree_1[[1]][[1]][[1]][[2]][[1]]),
    labels(sm2.ext.re.rowtree_1[[2]][[1]]),
    labels(sm2.ext.re.rowtree_1[[2]][[2]][[1]])))

# Extract:  V1.0
src.endoderm.sm2.ext_v1.1.gene.fin = 
  src.endoderm.sm2.ext.selectgene.fin 



#--------------------------------------------------------------------------------
src.endoderm.sm2.ext_v1.1 = 
  seurat_mnn(seurat_name = "src.endoderm.sm2.ext_v1.1", dims = 1:30, 
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.sm2.ext_v1.1.gene.fin)
src.endoderm.sm2.ext_v1.1 = 
  seurat_int(seurat_name = "src.endoderm.sm2.ext_v1.1", dims = 1:30,
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.sm2.ext_v1.1.gene.fin)
#-------------------------------------------------


# Neighbor Cluster
#-------------------
src.endoderm.sm2.ext_v1.1 = 
  FindNeighbors(src.endoderm.sm2.ext_v1.1, dims=1:30, 
                reduction = "pca_integrated", assay = "integrated")
src.endoderm.sm2.ext_v1.1 = 
  FindNeighbors(src.endoderm.sm2.ext_v1.1, dims=1:30, 
                reduction = "pca", assay = "RNA")
src.endoderm.sm2.ext_v1.1 = 
  FindClusters(src.endoderm.sm2.ext_v1.1, 
               resolution = 2, graph.name = "RNA_snn",)

src.endoderm.sm2.ext_v1.1 = 
  FindNeighbors(src.endoderm.sm2.ext_v1.1, dims=1:3, 
                reduction = "umap_mnn", assay = "RNA")
src.endoderm.sm2.ext_v1.1 = 
  FindClusters(src.endoderm.sm2.ext_v1.1, 
               resolution = 1.5, graph.name = "RNA_snn",)
#-------------------

pdf("figure.v08.07/try.pdf")
#-----------------------------
DimPlot(src.endoderm.sm2.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.2', reduction = "umap_mnn")
DimPlot(src.endoderm.sm2.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.1.5', reduction = "umap_mnn")

for(red in c("umap","umap_mnn","umap_integrated")){
  for(type in c('Time',
                "cluster.extract.v1.1..hg1",
                "cluster.extract.v1.1..mg2","cluster.extract.v1.1..mg3",
                "cluster.v06.26.re..hg1",
                "cluster.v06.26.re..mg2","cluster.v06.26.re..mg3")){
    print(
      DimPlot(src.endoderm.sm2.ext_v1.1, group.by = type, pt.size = 1.2,
              reduction = red, cols = colors_bg) + p_add_leg
    )
  }
}
#-----------------------------
dev.off()

# FDL
src.endoderm.sm2.ext_v1.1 = 
  seurat_fdl("src.endoderm.sm2.ext_v1.1", "pca","RNA")
src.endoderm.sm2.ext_v1.1 = 
  seurat_fdl("src.endoderm.sm2.ext_v1.1", "mnn","RNA")
src.endoderm.sm2.ext_v1.1 = 
  seurat_fdl("src.endoderm.sm2.ext_v1.1", "pca_integrated","integrated")

# Merge cell type & tracing code
#------------------------------------

src.endoderm.sm2.ext_v1.1$cluster.extract.v1.1..check = paste(
  src.endoderm.sm2.ext_v1.1$cluster.v06.26.re..hg1,
  src.endoderm.sm2.ext_v1.1$cluster.v06.26.re..mg2,
  src.endoderm.sm2.ext_v1.1$cluster.v06.26.re..mg3,  sep="_")
table(src.endoderm.sm2.ext_v1.1$cluster.extract.v1.1..check)

cell_ext_sm2.rm = rownames(src.endoderm.sm2.ext_v1.1@meta.data[
  src.endoderm.sm2.ext_v1.1$cluster.extract.v1.1..check%in%c(
    "Large.intestine.1_Small.intestine.2_Small.intestine.2",
    'Large.intestine.2_Small.intestine.2_Small.intestine.2',
    "NA_MG.2_Small.intestine.1","NA_MG.2_Small.intestine.2",
    "NA_Small.intestine.1_Small.intestine.2",
    "Small.intestine.2_Small.intestine.2_Stomach")|
    
    src.endoderm.sm2.ext_v1.1$cluster.v06.26.re..hg1%in%c(
      "Large.intestine.1","Large.intestine.2")|
    
    src.endoderm.sm2.ext_v1.1$cluster.v06.26.re..hg1%in%c(
      "Small.intestine.1")|
    
    src.endoderm.sm2.ext_v1.1$cluster.v06.26.re..mg3%in%c(
      "Small.intestine.1","Stomach"),])
#------------------------------------

src.endoderm.sm2.ext_v1.1.re = 
  src.endoderm.sm2.ext_v1.1[,!colnames(src.endoderm.sm2.ext_v1.1)%in%cell_ext_sm2.rm]
src.endoderm.sm2.ext_v1.1.re = 
  src.endoderm.sm2.ext_v1.1.re[,!src.endoderm.sm2.ext_v1.1.re$cluster.extract.v1.1..merge%in%"FG.4-MG.1/3"]
table(src.endoderm.sm2.ext_v1.1.re$cluster.extract.v1.1..merge)

src.endoderm.sm2.ext_v1.1.re$cluster.v06.26.re..merge = NA
src.endoderm.sm2.ext_v1.1.re@meta.data[
  (!src.endoderm.sm2.ext_v1.1.re$cluster.v06.26.re..hg1%in%NA)&
    (src.endoderm.sm2.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..merge =
  src.endoderm.sm2.ext_v1.1.re@meta.data[
    (!src.endoderm.sm2.ext_v1.1.re$cluster.v06.26.re..hg1%in%NA)&
      (src.endoderm.sm2.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..hg1
src.endoderm.sm2.ext_v1.1.re@meta.data[
  (!src.endoderm.sm2.ext_v1.1.re$cluster.v06.26.re..mg2%in%NA)&
    (src.endoderm.sm2.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..merge =
  src.endoderm.sm2.ext_v1.1.re@meta.data[
    (!src.endoderm.sm2.ext_v1.1.re$cluster.v06.26.re..mg2%in%NA)&
      (src.endoderm.sm2.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..mg2
src.endoderm.sm2.ext_v1.1.re@meta.data[
  (!src.endoderm.sm2.ext_v1.1.re$cluster.v06.26.re..mg3%in%NA)&
    (src.endoderm.sm2.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..merge =
  src.endoderm.sm2.ext_v1.1.re@meta.data[
    (!src.endoderm.sm2.ext_v1.1.re$cluster.v06.26.re..mg3%in%NA)&
      (src.endoderm.sm2.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..mg3


src.endoderm.sm2.ext_v1.1.re$cluster.extract.v1.1..merge = NA
src.endoderm.sm2.ext_v1.1.re@meta.data[
  (!src.endoderm.sm2.ext_v1.1.re$cluster.extract.v1.1..hg1%in%NA)&
    (src.endoderm.sm2.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..merge =
  src.endoderm.sm2.ext_v1.1.re@meta.data[
    (!src.endoderm.sm2.ext_v1.1.re$cluster.extract.v1.1..hg1%in%NA)&
      (src.endoderm.sm2.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..hg1
src.endoderm.sm2.ext_v1.1.re@meta.data[
  (!src.endoderm.sm2.ext_v1.1.re$cluster.extract.v1.1..mg2%in%NA)&
    (src.endoderm.sm2.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..merge =
  src.endoderm.sm2.ext_v1.1.re@meta.data[
    (!src.endoderm.sm2.ext_v1.1.re$cluster.extract.v1.1..mg2%in%NA)&
      (src.endoderm.sm2.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..mg2
src.endoderm.sm2.ext_v1.1.re@meta.data[
  (!src.endoderm.sm2.ext_v1.1.re$cluster.extract.v1.1..mg3%in%NA)&
    (src.endoderm.sm2.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..merge =
  src.endoderm.sm2.ext_v1.1.re@meta.data[
    (!src.endoderm.sm2.ext_v1.1.re$cluster.extract.v1.1..mg3%in%NA)&
      (src.endoderm.sm2.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..mg3

#--------------------------------
src.endoderm.sm2.ext_v1.1.re = 
  seurat_mnn(seurat_name = "src.endoderm.sm2.ext_v1.1.re", dims = 1:30, 
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.sm2.ext_v1.1.gene.fin)
src.endoderm.sm2.ext_v1.1.re = 
  seurat_int(seurat_name = "src.endoderm.sm2.ext_v1.1.re", dims = 1:30,
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.sm2.ext_v1.1.gene.fin)

src.endoderm.sm2.ext_v1.1.re@reductions$pca_integrated@cell.embeddings =
  src.endoderm.sm2.ext_v1.1.re@reductions$pca_integrated@cell.embeddings[colnames(src.endoderm.sm2.ext_v1.1.re),]
src.endoderm.sm2.ext_v1.1.re@reductions$mnn@cell.embeddings =
  src.endoderm.sm2.ext_v1.1.re@reductions$mnn@cell.embeddings[colnames(src.endoderm.sm2.ext_v1.1.re),]
src.endoderm.sm2.ext_v1.1.re@reductions$umap_mnn@cell.embeddings =
  src.endoderm.sm2.ext_v1.1.re@reductions$umap_mnn@cell.embeddings[colnames(src.endoderm.sm2.ext_v1.1.re),]
src.endoderm.sm2.ext_v1.1.re@reductions$umap_pca_integrated@cell.embeddings =
  src.endoderm.sm2.ext_v1.1.re@reductions$umap_pca_integrated@cell.embeddings[colnames(src.endoderm.sm2.ext_v1.1.re),]


src.endoderm.sm2.ext_v1.1.re = 
  seurat_fdl("src.endoderm.sm2.ext_v1.1.re", "pca","RNA")
src.endoderm.sm2.ext_v1.1.re = 
  seurat_fdl("src.endoderm.sm2.ext_v1.1.re", "mnn","RNA")
src.endoderm.sm2.ext_v1.1.re = 
  seurat_fdl("src.endoderm.sm2.ext_v1.1.re", "pca_integrated","integrated")
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
  
  pdf(paste("figure.v08.07/organ_development_merge_re/sm2/sm2_summary_graph.re_snn_",red,".pdf",sep = ""),12,7)
  
  assay = ifelse(red=="mnn","RNA","integrated")
  
  graph.sm2.ext_v1.1 = 
    seurat_GCN(seurat_name = "src.endoderm.sm2.ext_v1.1.re",
               src.endoderm.sm2.ext_v1.1.gene.fin,
               reduction = red, assay =  assay)
  
  graph.sm2_cluster = cluster_walktrap(graph.sm2.ext_v1.1, steps = 1)
  
  for(color in c("cluster.extract.v1.1..merge",
                 "cluster.v06.26.re..merge",
                 "Time")){
    for(edge.color in c(NA,"lightgray")){
      for(edge.color in c(NA,"lightgray")){
        set.seed(1)
        plot.igraph(graph.sm2.ext_v1.1,
                    layout = layout_with_fr(graph.sm2.ext_v1.1),
                    edge.color = edge.color,
                    vertex.size=2.5, vertex.label=NA,
                    vertex.label.color = "black", 
                    vertex.frame.color = NA, vertex.frame.width = 0.5,
                    #vertex.color =  colors.time[src.endoderm.sm2.ext_v1.1[["Time"]][graph.sm2_cluster$names,]],
                    #vertex.color =  cluster.endoderm.color.v5[src.endoderm.sm2.ext_v1.1[["cluster.v06.26.re"]][graph.sm2_cluster$names,]],
                    vertex.color =  colors_bg[src.endoderm.sm2.ext_v1.1.re[[color]][graph.sm2_cluster$names,]])
      }}
  }
  dev.off()
  
}


pdf("figure.v08.07/organ_development_merge_re/sm2/sm2_summary_graph.re_other.pdf",12,7)
# Reduction
for(red in c("umap_integrated","umap_mnn",
             "fdl_pca",'fdl_mnn',"fdl_pca_integrated")){
  
  for(color in c("cluster.extract.v1.1..merge",
                 "cluster.extract.v1.1..hg1", "cluster.extract.v1.1..mg2",
                 "cluster.extract.v1.1..mg3",
                 "cluster.v06.26.re..merge",
                 "cluster.v06.26.re..hg1", "cluster.v06.26.re..mg2",
                 "cluster.v06.26.re..mg3",
                 "Time")){
    print(
      DimPlot(src.endoderm.sm2.ext_v1.1.re, group.by = color, 
              reduction = red, cols = colors_bg,
              pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
    )
  }
}
dev.off() 

save(src.endoderm.sm2.ext_v1.1,
     src.endoderm.sm2.ext_v1.1.re,
     file = "figure.v08.07/organ_development_merge_re/sm2/src.endoderm.sm2.ext_v1.1.re.Rdata")


metadata.sm2.ext_v1.1.re = cbind(
  src.endoderm.sm2.ext_v1.1.re@meta.data,
  src.endoderm.sm2.ext_v1.1.re@reductions$umap_integrated@cell.embeddings[colnames(src.endoderm.sm2.ext_v1.1.re),],
  src.endoderm.sm2.ext_v1.1.re@reductions$fdl_pca_integrated@cell.embeddings[colnames(src.endoderm.sm2.ext_v1.1.re),])
save(metadata.sm2.ext_v1.1.re,
     file = "figure.v08.07/organ_development_merge_re/sm2/metadata.sm2.ext_v1.1.re.Rdata")

#----------------
#     Remake
#----------------
src.endoderm.lar1.ext_v1.1 = FindVariableFeatures(src.endoderm.lar1.ext_v1.1, nfeatures = 2000)
src.endoderm.lar1.ext_v1.1 = ScaleData(src.endoderm.lar1.ext_v1.1, split.by = "batch_phase",
                                       features = rownames(src.endoderm.lar1.ext_v1.1))
src.endoderm.lar1.ext_v1.1.filtergene = 
  Myfilter(as.matrix(src.endoderm.lar1.ext_v1.1@assays$RNA@data),
           gene = src.endoderm.lar1.ext_v1.1@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.endoderm.lar1.ext.v1.1.selectgene = 
  src.endoderm.lar1.ext_v1.1.filtergene
rm(src.endoderm.lar1.ext_v1.1.filtergene)

#=============================
# Row
pdf("figure.v08.07/try.pdf")
#-----------------
lar1.ext.re.rowtree =
  MyHeatmap(as.matrix(src.endoderm.lar1.ext_v1.1@assays$RNA@data[
    src.endoderm.lar1.ext.v1.1.selectgene,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    hc.r.data.type = "row.relat",
    c.cov.method = "s",
    r.cov.method = "s",
    ColSideColors = cbind(
      MyName2Col(src.endoderm.lar1.ext_v1.1$Time,colors.time),
      MyName2Col(src.endoderm.lar1.ext_v1.1$cluster.extract.v1.1..hg1,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.lar1.ext_v1.1$cluster.extract.v1.1..hg2,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.lar1.ext_v1.1$cluster.v06.26.re..hg1,cluster.endoderm.color.v5),
      MyName2Col(src.endoderm.lar1.ext_v1.1$cluster.v06.26.re..hg2,cluster.endoderm.color.v5)
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
plot(lar1.ext.re.rowtree)
dev.off()

lar1.ext.re.rowtree_1 = as.dendrogram(lar1.ext.re.rowtree)
src.endoderm.lar1.ext_v1.1.gene = setdiff(
  src.endoderm.lar1.ext.v1.1.selectgene,
  c(labels(lar1.ext.re.rowtree_1[[1]][[2]]),
    labels(lar1.ext.re.rowtree_1[[2]][[1]])))

# Extract:  V1.0
src.endoderm.lar1.ext_v1.1.gene.fin = 
  src.endoderm.lar1.ext.selectgene.fin 


#--------------------------------------------------------------------------------
src.endoderm.lar1.ext_v1.1 = 
  seurat_mnn(seurat_name = "src.endoderm.lar1.ext_v1.1", dims = 1:30, 
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.lar1.ext_v1.1.gene.fin)
src.endoderm.lar1.ext_v1.1 = 
  seurat_int(seurat_name = "src.endoderm.lar1.ext_v1.1", dims = 1:30,
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.lar1.ext_v1.1.gene.fin)
#-------------------------------------------------


# Neighbor Cluster
#-------------------
src.endoderm.lar1.ext_v1.1 = 
  FindNeighbors(src.endoderm.lar1.ext_v1.1, dims=1:30, 
                reduction = "pca_integrated", assay = "integrated")
src.endoderm.lar1.ext_v1.1 = 
  FindNeighbors(src.endoderm.lar1.ext_v1.1, dims=1:30, 
                reduction = "pca", assay = "RNA")
src.endoderm.lar1.ext_v1.1 = 
  FindClusters(src.endoderm.lar1.ext_v1.1, 
               resolution = 2, graph.name = "RNA_snn",)

src.endoderm.lar1.ext_v1.1 = 
  FindNeighbors(src.endoderm.lar1.ext_v1.1, dims=1:3, 
                reduction = "umap_mnn", assay = "RNA")
src.endoderm.lar1.ext_v1.1 = 
  FindClusters(src.endoderm.lar1.ext_v1.1, 
               resolution = 1.5, graph.name = "RNA_snn",)
#-------------------

pdf("figure.v08.07/try.pdf")
#-----------------------------
DimPlot(src.endoderm.lar1.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.2', reduction = "umap_mnn")
DimPlot(src.endoderm.lar1.ext_v1.1, label = T, label.size = 4, cols = colors.num,
        group.by = 'RNA_snn_res.1.5', reduction = "umap_mnn")

for(red in c("umap","umap_mnn","umap_integrated")){
  for(type in c('Time',
                "cluster.extract.v1.1..hg1", "cluster.extract.v1.1..hg2",
                "cluster.v06.26.re..hg1", "cluster.v06.26.re..hg2")){
    print(
      DimPlot(src.endoderm.lar1.ext_v1.1, group.by = type, pt.size = 1.2,
              reduction = red, cols = colors_bg) + p_add_leg
    )
  }
}
#-----------------------------
dev.off()

# FDL
src.endoderm.lar1.ext_v1.1 = 
  seurat_fdl("src.endoderm.lar1.ext_v1.1", "pca","RNA")
src.endoderm.lar1.ext_v1.1 = 
  seurat_fdl("src.endoderm.lar1.ext_v1.1", "mnn","RNA")
src.endoderm.lar1.ext_v1.1 = 
  seurat_fdl("src.endoderm.lar1.ext_v1.1", "pca_integrated","integrated")

# Merge cell type & tracing code
#------------------------------------

src.endoderm.lar1.ext_v1.1$cluster.extract.v1.1..check = paste(
  src.endoderm.lar1.ext_v1.1$cluster.v06.26.re..hg1,
  src.endoderm.lar1.ext_v1.1$cluster.v06.26.re..hg2,  sep="_")
table(src.endoderm.lar1.ext_v1.1$cluster.extract.v1.1..check)

cell_ext_lar1.rm = rownames(src.endoderm.lar1.ext_v1.1@meta.data[
  src.endoderm.lar1.ext_v1.1$cluster.extract.v1.1..hg1%in%c(
    "MG.2/3−HG.1"),])
#------------------------------------

src.endoderm.lar1.ext_v1.1.re = 
  src.endoderm.lar1.ext_v1.1[,!colnames(src.endoderm.lar1.ext_v1.1)%in%cell_ext_lar1.rm]

src.endoderm.lar1.ext_v1.1.re$cluster.v06.26.re..merge = NA
src.endoderm.lar1.ext_v1.1.re@meta.data[
  (!src.endoderm.lar1.ext_v1.1.re$cluster.v06.26.re..hg1%in%NA)&
    (src.endoderm.lar1.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..merge =
  src.endoderm.lar1.ext_v1.1.re@meta.data[
    (!src.endoderm.lar1.ext_v1.1.re$cluster.v06.26.re..hg1%in%NA)&
      (src.endoderm.lar1.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..hg1
src.endoderm.lar1.ext_v1.1.re@meta.data[
  (!src.endoderm.lar1.ext_v1.1.re$cluster.v06.26.re..hg2%in%NA)&
    (src.endoderm.lar1.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..merge =
  src.endoderm.lar1.ext_v1.1.re@meta.data[
    (!src.endoderm.lar1.ext_v1.1.re$cluster.v06.26.re..hg2%in%NA)&
      (src.endoderm.lar1.ext_v1.1.re$cluster.v06.26.re..merge%in%NA),]$cluster.v06.26.re..hg2

src.endoderm.lar1.ext_v1.1.re$cluster.extract.v1.1..merge = NA
src.endoderm.lar1.ext_v1.1.re@meta.data[
  (!src.endoderm.lar1.ext_v1.1.re$cluster.extract.v1.1..hg1%in%NA)&
    (src.endoderm.lar1.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..merge =
  src.endoderm.lar1.ext_v1.1.re@meta.data[
    (!src.endoderm.lar1.ext_v1.1.re$cluster.extract.v1.1..hg1%in%NA)&
      (src.endoderm.lar1.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..hg1
src.endoderm.lar1.ext_v1.1.re@meta.data[
  (!src.endoderm.lar1.ext_v1.1.re$cluster.extract.v1.1..hg2%in%NA)&
    (src.endoderm.lar1.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..merge =
  src.endoderm.lar1.ext_v1.1.re@meta.data[
    (!src.endoderm.lar1.ext_v1.1.re$cluster.extract.v1.1..hg2%in%NA)&
      (src.endoderm.lar1.ext_v1.1.re$cluster.extract.v1.1..merge%in%NA),]$cluster.extract.v1.1..hg2

#--------------------------------
src.endoderm.lar1.ext_v1.1.re = 
  seurat_mnn(seurat_name = "src.endoderm.lar1.ext_v1.1.re", dims = 1:30, 
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.lar1.ext_v1.1.gene.fin)
src.endoderm.lar1.ext_v1.1.re = 
  seurat_int(seurat_name = "src.endoderm.lar1.ext_v1.1.re", dims = 1:30,
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.lar1.ext_v1.1.gene.fin)

src.endoderm.lar1.ext_v1.1.re@reductions$pca_integrated@cell.embeddings =
  src.endoderm.lar1.ext_v1.1.re@reductions$pca_integrated@cell.embeddings[colnames(src.endoderm.lar1.ext_v1.1.re),]
src.endoderm.lar1.ext_v1.1.re@reductions$mnn@cell.embeddings =
  src.endoderm.lar1.ext_v1.1.re@reductions$mnn@cell.embeddings[colnames(src.endoderm.lar1.ext_v1.1.re),]
src.endoderm.lar1.ext_v1.1.re@reductions$umap_mnn@cell.embeddings =
  src.endoderm.lar1.ext_v1.1.re@reductions$umap_mnn@cell.embeddings[colnames(src.endoderm.lar1.ext_v1.1.re),]
src.endoderm.lar1.ext_v1.1.re@reductions$umap_integrated@cell.embeddings =
  src.endoderm.lar1.ext_v1.1.re@reductions$umap_integrated@cell.embeddings[colnames(src.endoderm.lar1.ext_v1.1.re),]


src.endoderm.lar1.ext_v1.1.re = 
  seurat_fdl("src.endoderm.lar1.ext_v1.1.re", "pca","RNA")
src.endoderm.lar1.ext_v1.1.re = 
  seurat_fdl("src.endoderm.lar1.ext_v1.1.re", "mnn","RNA")
src.endoderm.lar1.ext_v1.1.re = 
  seurat_fdl("src.endoderm.lar1.ext_v1.1.re", "pca_integrated","integrated")
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
  
  pdf(paste("figure.v08.07/organ_development_merge_re/lar1/lar1_summary_graph.re_snn_",red,".pdf",sep = ""),12,7)
  
  assay = ifelse(red=="mnn","RNA","integrated")
  
  graph.lar1.ext_v1.1 = 
    seurat_GCN(seurat_name = "src.endoderm.lar1.ext_v1.1.re",
               src.endoderm.lar1.ext_v1.1.gene.fin,
               reduction = red, assay =  assay)
  
  graph.lar1_cluster = cluster_walktrap(graph.lar1.ext_v1.1, steps = 1)
  
  for(color in c("cluster.extract.v1.1..merge",
                 "cluster.v06.26.re..merge",
                 "Time")){
    for(edge.color in c(NA,"lightgray")){
      for(edge.color in c(NA,"lightgray")){
        set.seed(1)
        plot.igraph(graph.lar1.ext_v1.1,
                    layout = layout_with_fr(graph.lar1.ext_v1.1),
                    edge.color = edge.color,
                    vertex.size=2.5, vertex.label=NA,
                    vertex.label.color = "black", 
                    vertex.frame.color = NA, vertex.frame.width = 0.5,
                    #vertex.color =  colors.time[src.endoderm.lar1.ext_v1.1[["Time"]][graph.lar1_cluster$names,]],
                    #vertex.color =  cluster.endoderm.color.v5[src.endoderm.lar1.ext_v1.1[["cluster.v06.26.re"]][graph.lar1_cluster$names,]],
                    vertex.color =  colors_bg[src.endoderm.lar1.ext_v1.1.re[[color]][graph.lar1_cluster$names,]])
      }}
  }
  dev.off()
  
}


pdf("figure.v08.07/organ_development_merge_re/lar1/lar1_summary_graph.re_other.pdf",12,7)
# Reduction
for(red in c("umap_integrated","umap_mnn",
             "fdl_pca",'fdl_mnn',"fdl_pca_integrated")){
  
  for(color in c("cluster.extract.v1.1..merge",
                 "cluster.extract.v1.1..hg1", "cluster.extract.v1.1..hg2",
                 "cluster.v06.26.re..merge",
                 "cluster.v06.26.re..hg1", "cluster.v06.26.re..hg2",
                 "Time")){
    print(
      DimPlot(src.endoderm.lar1.ext_v1.1.re, group.by = color, 
              reduction = red, cols = colors_bg,
              pt.size = 1.15,  na.value = "darkgray")+ p_add + theme(aspect.ratio = 1)
    )
  }
}
dev.off() 


pdf("figure.v08.07/organ_development_merge_re/lar1/lar1_summary_proportion.re.pdf",12,7)
#-----------------------------
cell_type = c("HG.1","HG.2","Large.intestine.1")
ext_type = c("HG.1","HG.2")

meta_lar1 = src.endoderm.lar1.ext_v1.1.re@meta.data
meta_lar1$Time = factor(meta_lar1$Time, levels = names(colors.time))
meta_lar1$cluster.extract.v1.1..merge = factor(meta_lar1$cluster.extract.v1.1..merge, levels = rev(ext_type))
meta_lar1$cluster.extract.v1.1..hg1 = factor(meta_lar1$cluster.extract.v1.1..hg1, levels = rev(ext_type))
meta_lar1$cluster.extract.v1.1..hg2 = factor(meta_lar1$cluster.extract.v1.1..hg2, levels = rev(ext_type))

meta_lar1$cluster.v06.26.re..merge = factor(meta_lar1$cluster.v06.26.re..merge, levels = rev(cell_type))
meta_lar1$cluster.v06.26.re..hg1 = factor(meta_lar1$cluster.v06.26.re..hg1, levels = rev(cell_type))
meta_lar1$cluster.v06.26.re..hg2 = factor(meta_lar1$cluster.v06.26.re..hg2, levels = rev(cell_type))

meta_lar1_B0 = meta_lar1
meta_lar1_B1 = meta_lar1[meta_lar1$batch%in%1&!meta_lar1$Time%in%"ss9",]
meta_lar1_B2 = meta_lar1[meta_lar1$batch%in%2,]

for(data in c("meta_lar1_B0","meta_lar1_B1","meta_lar1_B2")){
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

save(src.endoderm.lar1.ext_v1.1,
     src.endoderm.lar1.ext_v1.1.re,
     file = "figure.v08.07/organ_development_merge_re/lar1/src.endoderm.lar1.ext_v1.1.re.Rdata")


metadata.lar1.ext_v1.1.re = cbind(
  src.endoderm.lar1.ext_v1.1.re@meta.data,
  src.endoderm.lar1.ext_v1.1.re@reductions$umap_integrated@cell.embeddings[colnames(src.endoderm.lar1.ext_v1.1.re),],
  src.endoderm.lar1.ext_v1.1.re@reductions$fdl_pca_integrated@cell.embeddings[colnames(src.endoderm.lar1.ext_v1.1.re),])
save(metadata.lar1.ext_v1.1.re,
     file = "figure.v08.07/organ_development_merge_re/lar1/metadata.lar1.ext_v1.1.re_integrated.Rdata")

metadata.lar1.ext_v1.1.re = cbind(
  src.endoderm.lar1.ext_v1.1.re@meta.data,
  src.endoderm.lar1.ext_v1.1.re@reductions$umap_mnn@cell.embeddings[colnames(src.endoderm.lar1.ext_v1.1.re),],
  src.endoderm.lar1.ext_v1.1.re@reductions$fdl_mnn@cell.embeddings[colnames(src.endoderm.lar1.ext_v1.1.re),])
save(metadata.lar1.ext_v1.1.re,
     file = "figure.v08.07/organ_development_merge_re/lar1/metadata.lar1.ext_v1.1.re_mnn.Rdata")



