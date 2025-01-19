#-------------------------------------
# Pancreas
#-------------------------------------
cell_mg3_pan_mg3a = rownames(src.mg3.integrated@meta.data[
  src.mg3.integrated$cluster.v06.26.re%in%"MG.3.A",])
intersect(colnames(src.endoderm.pan.ext_v1.1.re), cell_mg3_pan_mg3a)

src.endoderm.pan.ext_v1.1.refine = merge(
  src.endoderm.pan.ext_v1.1.re,
  src.mg3.integrated[,cell_mg3_pan_mg3a])
src.endoderm.pan.ext_v1.1.refine = FindVariableFeatures(src.endoderm.pan.ext_v1.1.refine, nfeatures = 2000)
src.endoderm.pan.ext_v1.1.refine = ScaleData(src.endoderm.pan.ext_v1.1.refine, split.by = "batch_phase",
                                             features = rownames(src.endoderm.pan.ext_v1.1.refine))
src.endoderm.pan.ext.v1.1.selectgene = src.endoderm.pan.ext_v1.1.gene.fin
src.endoderm.pan.ext_v1.1.refine.gene.fin = src.endoderm.pan.ext_v1.1.gene.fin

src.endoderm.pan.ext_v1.1.refine = RunPCA(src.endoderm.pan.ext_v1.1.refine, 
                                          features = src.endoderm.pan.ext_v1.1.refine.gene.fin)
src.endoderm.pan.ext_v1.1.refine = RunUMAP(src.endoderm.pan.ext_v1.1.refine, dims = 1:30, 
                                           reduction = "pca", n.components = 3)
DimPlot(src.endoderm.pan.ext_v1.1.refine, reduction = "umap", group.by = "cluster.v06.26.re..merge")

#--------------------------------------------------------------------------------
src.endoderm.pan.ext_v1.1.refine = 
  seurat_mnn(seurat_name = "src.endoderm.pan.ext_v1.1.refine",
             dims = 1:30,
             n.neighbors = 100,
             n.components = 3,
             seurat.selectgene = src.endoderm.pan.ext_v1.1.refine.gene.fin)

src.endoderm.pan.ext_v1.1.refine = 
  seurat_int(seurat_name = "src.endoderm.pan.ext_v1.1.refine", dims = 1:30,
             n.neighbors = 100,n.components = 3,
             seurat.selectgene = src.endoderm.pan.ext_v1.1.refine.gene.fin)

#-------------------------------------------------
src.endoderm.pan.ext_v1.1.refine = seurat_fdl(src.endoderm.pan.ext_v1.1.refine, "pca","RNA")
src.endoderm.pan.ext_v1.1.refine = seurat_fdl(src.endoderm.pan.ext_v1.1.refine, "mnn","RNA")
src.endoderm.pan.ext_v1.1.refine = seurat_fdl(src.endoderm.pan.ext_v1.1.refine, "pca_integrated","integrated")
src.endoderm.pan.ext_v1.1.refine@meta.data[cell_mg3_pan_mg3a,]$cluster.v06.26.re..merge = "MG.3.A"

pdf("figure.v08.07/organ_development_re_v240115/try.pancreas.re.pdf", 9, 7)
for(reduction in names(src.endoderm.pan.ext_v1.1.refine@reductions)){
  print(
    DimPlot(src.endoderm.pan.ext_v1.1.refine, reduction = reduction, 
            group.by = "cluster.v06.26.re..merge",
            cols = cluster.endoderm.color.v5) +
      theme_classic() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
      theme(legend.position = "none") +
      ggtitle(label = reduction)
  )
}
dev.off()

seurat = RunUMAP(src.endoderm.pan.ext_v1.1.refine, reduction = "pca_integrated",
                 dims = 1:30, n.neighbors = 100, #min.dist = 0.05,
                 n.components = 3 #, learning.rate = 0.01
                 )
DimPlot(seurat, reduction = "umap", c(1,2),
        group.by = "cluster.v06.26.re..merge",
        cols = cluster.endoderm.color.v5)
src.endoderm.pan.ext_v1.1.refine@reductions$umap_integrated@cell.embeddings =
  seurat@reductions$umap_integrated@cell.embeddings

save(src.endoderm.pan.ext_v1.1.refine,
     file = "figure.v08.07/organ_development_re_v240115/src.endoderm.pan.ext_v1.1.refine.Rdata")



