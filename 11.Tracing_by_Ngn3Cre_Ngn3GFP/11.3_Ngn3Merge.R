#-----------------------------------------------
#>> 11.3 Merge Ngn3-Cre and Ngn3-GFP
#------------------------------------------------


src.sm3.update.rawdata_Ngn3 = merge(
  src.sm3.update.rawdata_Ngn3Cre,
  src.sm3.update.rawdata_Ngn3GFP)

src.sm3.update.rawdata_Ngn3 = FindVariableFeatures(src.sm3.update.rawdata_Ngn3, nfeatures = 2000)
src.sm3.update.rawdata_Ngn3 = ScaleData(src.sm3.update.rawdata_Ngn3, features = rownames(src.sm3.update.rawdata_Ngn3))

src.sm3.update.rawdata_Ngn3.selectgene = 
  Myfilter(as.matrix(src.sm3.update.rawdata_Ngn3@assays$RNA@data),
           gene = src.sm3.update.rawdata_Ngn3@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.sm3.update.rawdata_Ngn3.selectgene.re = src.sm3.update.rawdata_Ngn3.selectgene
rm(src.sm3.update.rawdata_Ngn3.selectgene)
src.sm3.update.rawdata_Ngn3.selectgene = src.sm3.update.rawdata_Ngn3.selectgene.re

src.sm3.update.rawdata_Ngn3 = RunPCA(src.sm3.update.rawdata_Ngn3, 
                                     features = src.sm3.update.rawdata_Ngn3.selectgene)
src.sm3.update.rawdata_Ngn3 = RunUMAP(src.sm3.update.rawdata_Ngn3, dims = 1:30)
src.sm3.update.rawdata_Ngn3 = RunTSNE(src.sm3.update.rawdata_Ngn3, dims = 1:30)

src.sm3.update.rawdata_Ngn3$lineage = "Ngn3Cre"
src.sm3.update.rawdata_Ngn3@meta.data[
  rownames(src.sm3.update.rawdata_Ngn3GFP@meta.data[
    src.sm3.update.rawdata_Ngn3GFP$index%in%c("HM","HZ"),]),]$lineage = 
  paste("Ngn3GFP", 
        src.sm3.update.rawdata_Ngn3GFP@meta.data[
          src.sm3.update.rawdata_Ngn3GFP$index%in%c("HM","HZ"),]$index, sep="")

src.sm3.update.rawdata_Ngn3@meta.data[
  rownames(src.sm3.update.rawdata_Ngn3GFP@meta.data[
    src.sm3.update.rawdata_Ngn3GFP$Treatment%in%c("HM","HZ"),]),]$lineage = 
  paste("Ngn3GFP", 
        src.sm3.update.rawdata_Ngn3GFP@meta.data[
          src.sm3.update.rawdata_Ngn3GFP$Treatment%in%c("HM","HZ"),]$Treatment, sep="")


pdf("figure.v08.07/Ngn3_process/qc.ngn3.pdf",9,7)
data = src.sm3.update.rawdata_Ngn3@meta.data
ggplot() +
  geom_histogram(data[data$nFeature_RNA<=2500,], 
                 mapping = aes(x = nFeature_RNA),
                 binwidth = 5, colour = "darkgray") +
  geom_histogram(data[data$nFeature_RNA>2500,], 
                 mapping = aes(x = nFeature_RNA),
                 binwidth = 5, colour = "black") +
  geom_line(mapping = aes(x=c(6000,6000), y=c(0,6.5)), 
            linetype = "dashed", size = 1.5)+
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks =  element_line(linewidth=1.2),
        axis.line.x = element_line(linewidth=1.2),
        axis.line.y = element_line(linewidth=1.2),
        aspect.ratio=1/3)

DimPlot(src.sm3.update.rawdata_Ngn3, group.by = "cluster", reduction = "tsne",
        cols = qc.color2, pt.size = 1.25)+
  theme_classic() + theme(aspect.ratio = 1) + p_add

DimPlot(src.sm3.update.rawdata_Ngn3, group.by = "lineage", reduction = "tsne",
        cols = c(colors.type.lineage, color.lineage), 
        pt.size = 1.25)+
  theme_classic() + theme(aspect.ratio = 1) + p_add
dev.off()



#===============================================================================
src.sm3.endoderm_Ngn3Cre = FindVariableFeatures(src.sm3.endoderm_Ngn3Cre, nfeatures = 2000)
src.sm3.endoderm_Ngn3Cre = ScaleData(src.sm3.endoderm_Ngn3Cre, features = rownames(src.sm3.endoderm_Ngn3Cre))

src.sm3.endoderm_Ngn3Cre.selectgene = 
  Myfilter(as.matrix(src.sm3.endoderm_Ngn3Cre@assays$RNA@data),
           gene = src.sm3.endoderm_Ngn3Cre@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.sm3.endoderm_Ngn3Cre.selectgene.re = src.sm3.endoderm_Ngn3Cre.selectgene
rm(src.sm3.endoderm_Ngn3Cre.selectgene)
src.sm3.endoderm_Ngn3Cre.selectgene = src.sm3.endoderm_Ngn3Cre.selectgene.re
src.sm3.endoderm_Ngn3Cre = RunPCA(src.sm3.endoderm_Ngn3Cre, 
                                  features = src.sm3.endoderm_Ngn3Cre.selectgene)
src.sm3.endoderm_Ngn3Cre = RunUMAP(src.sm3.endoderm_Ngn3Cre, dims = 1:30)
src.sm3.endoderm_Ngn3Cre = RunTSNE(src.sm3.endoderm_Ngn3Cre, dims = 1:30)


src.sm3.endoderm_Ngn3GFP = FindVariableFeatures(src.sm3.endoderm_Ngn3GFP, nfeatures = 2000)
src.sm3.endoderm_Ngn3GFP = ScaleData(src.sm3.endoderm_Ngn3GFP, 
                                     features = rownames(src.sm3.endoderm_Ngn3GFP), split.by = "SeqDate")
src.sm3.endoderm_Ngn3GFP.selectgene = 
  Myfilter(as.matrix(src.sm3.endoderm_Ngn3GFP@assays$RNA@data),
           gene = src.sm3.endoderm_Ngn3GFP@assays$RNA@var.features,
           pearson.threshold = 0.15,partner.threshold = 5,
           bottom.dispersion.interval = 0.1)
src.sm3.endoderm_Ngn3GFP.selectgene.re = src.sm3.endoderm_Ngn3GFP.selectgene
rm(src.sm3.endoderm_Ngn3GFP.selectgene)
src.sm3.endoderm_Ngn3GFP.selectgene = src.sm3.endoderm_Ngn3GFP.selectgene.re
src.sm3.endoderm_Ngn3GFP = RunPCA(src.sm3.endoderm_Ngn3GFP, 
                                  features = src.sm3.endoderm_Ngn3GFP.selectgene)
src.sm3.endoderm_Ngn3GFP = RunUMAP(src.sm3.endoderm_Ngn3GFP, dims = 1:30)
src.sm3.endoderm_Ngn3GFP = RunTSNE(src.sm3.endoderm_Ngn3GFP, dims = 1:30)

DimPlot(src.sm3.endoderm_Ngn3Cre, reduction = "umap",
        group.by = "cluster.predict.umap_int.ext.v1.1",
        cols = cluster.endoderm.color.v5)
DimPlot(src.sm3.endoderm_Ngn3GFP, reduction = "umap",
        group.by = "cluster.predict.umap_int.ext.v1.1",
        cols = cluster.endoderm.color.v5)


src.sm3.endoderm_Ngn3 = merge(src.sm3.endoderm_Ngn3GFP, src.sm3.endoderm_Ngn3Cre) 
table(src.sm3.endoderm_Ngn3$Time)



src.sm3.endoderm_Ngn3Cre$cluster.predict.umap_int.ext.v1.1 = 
  src.sm3.endoderm_Ngn3$cluster.predict.umap_int.ext.v1.1[colnames(src.sm3.endoderm_Ngn3Cre)]
src.sm3.endoderm_Ngn3GFP$cluster.predict.umap_int.ext.v1.1 = 
  src.sm3.endoderm_Ngn3$cluster.predict.umap_int.ext.v1.1[colnames(src.sm3.endoderm_Ngn3GFP)]

src.sm3.endoderm_Ngn3Cre$lineage = src.sm3.update.rawdata_Ngn3$lineage[colnames(src.sm3.endoderm_Ngn3Cre)]
src.sm3.endoderm_Ngn3GFP$lineage = src.sm3.update.rawdata_Ngn3$lineage[colnames(src.sm3.endoderm_Ngn3GFP)]



src.15ss.integrated.merge.selectgene = rownames(src.15ss.integrated@reductions$pca@feature.loadings)
src.18ss.integrated.merge.selectgene = rownames(src.18ss.integrated@reductions$pca@feature.loadings)
src.21ss.integrated.merge.selectgene = rownames(src.21ss.integrated@reductions$pca@feature.loadings)
src.24ss.integrated.merge.selectgene = rownames(src.24ss.integrated@reductions$pca@feature.loadings)
src.27ss.integrated.merge.selectgene = rownames(src.27ss.integrated@reductions$pca@feature.loadings)

# 15ss
#--------------------------------------------------------------------------------
src.15ss.integrated.merge_20231223 = merge(src.15ss.integrated, 
                                           c(src.sm3.endoderm_Ngn3[,src.sm3.endoderm_Ngn3$Time%in%"15ss"],
                                             src.sm3.merge[,src.sm3.merge$Time%in%"15ss"]))
src.15ss.integrated.merge_20231223 = src.15ss.integrated.merge_20231223[,src.15ss.integrated.merge_20231223$nFeature_RNA>2500]

anchor.15ss.integrated.merge_20231223 = FindIntegrationAnchors(c(src.15ss.integrated.merge_20231223[,src.15ss.integrated.merge_20231223$Time%in%"ss15"&
                                                                                                      src.15ss.integrated.merge_20231223$batch%in%1],
                                                                 src.15ss.integrated.merge_20231223[,src.15ss.integrated.merge_20231223$Time%in%"ss15"&
                                                                                                      src.15ss.integrated.merge_20231223$batch%in%2],
                                                                 src.15ss.integrated.merge_20231223[,!src.15ss.integrated.merge_20231223$Time%in%"ss15"]),
                                                               reduction = "cca",
                                                               anchor.features = src.15ss.integrated.merge.selectgene,
                                                               dims = 1:30)

src.15ss.integrated.merge_20231223.re = IntegrateData(anchorset = anchor.15ss.integrated.merge_20231223 , dims = 1:30)
DefaultAssay(src.15ss.integrated.merge_20231223.re) = "integrated"
src.15ss.integrated.merge_20231223.re = ScaleData(src.15ss.integrated.merge_20231223.re, features = rownames(src.15ss.integrated.merge_20231223.re))
src.15ss.integrated.merge_20231223.re = RunPCA(src.15ss.integrated.merge_20231223.re, features = rownames(src.15ss.integrated@reductions$pca@feature.loadings))
src.15ss.integrated.merge_20231223.re = FindNeighbors(src.15ss.integrated.merge_20231223.re, dims = 1:30)
src.15ss.integrated.merge_20231223.re = FindClusters(src.15ss.integrated.merge_20231223.re, resolution = 2)
src.15ss.integrated.merge_20231223.re = FindClusters(src.15ss.integrated.merge_20231223.re, resolution = 3)
src.15ss.integrated.merge_20231223.re = FindClusters(src.15ss.integrated.merge_20231223.re, resolution = 4)
src.15ss.integrated.merge_20231223.re = RunUMAP(src.15ss.integrated.merge_20231223.re, dims = 1:20, n.neighbors = 100)

src.15ss.integrated.merge_20231223@reductions$umap_integrated = src.15ss.integrated.merge_20231223.re@reductions$umap
src.15ss.integrated.merge_20231223@reductions$pca_integrated = src.15ss.integrated.merge_20231223.re@reductions$pca
src.15ss.integrated.merge_20231223@reductions$umap_integrated@assay.used = "RNA"
src.15ss.integrated.merge_20231223@reductions$pca_integrated@assay.used = "RNA"
rm(src.15ss.integrated.merge_20231223.re)

src.15ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings = cbind(
  -src.15ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings[,1],
  -src.15ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings[,2])
colnames(src.15ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings) = c("UMAP_1","UMAP_2")
DimPlot(src.15ss.integrated.merge_20231223, group.by = "lineage", reduction = "umap_integrated",
        cols = color.lineage, na.value = "#eeeeee")

a = FNN::knn(src.15ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings[
  rownames(src.15ss.integrated.merge_20231223@meta.data[src.15ss.integrated.merge_20231223$Time%in%"ss15"&
                                                          !src.15ss.integrated.merge_20231223$cluster.revise.re.v1.30.re%in%NA,]),],
  src.15ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings[
    rownames(src.15ss.integrated.merge_20231223@meta.data[!src.15ss.integrated.merge_20231223$Time%in%"ss15",]),],
  src.15ss.integrated.merge_20231223@meta.data[src.15ss.integrated.merge_20231223$Time%in%"ss15"&
                                                 !src.15ss.integrated.merge_20231223$cluster.revise.re.v1.30.re%in%NA,]$cluster.v06.26.re..merge, k = 10)

src.15ss.integrated.merge_20231223$cluster.predict.umap_int.ext.v1.1 = src.15ss.integrated.merge_20231223$cluster.v06.26.re..merge
src.15ss.integrated.merge_20231223@meta.data[!src.15ss.integrated.merge_20231223$Time%in%"ss15",]$cluster.predict.umap_int.ext.v1.1 = as.character(a)
src.sm3.endoderm_Ngn3@meta.data[src.sm3.endoderm_Ngn3$Time%in%"15ss",]$cluster.predict.umap_int.ext.v1.1 = 
  src.15ss.integrated.merge_20231223@meta.data[rownames(
    src.sm3.endoderm_Ngn3@meta.data[src.sm3.endoderm_Ngn3$Time%in%"15ss",]),]$cluster.predict.umap_int.ext.v1.1
#--------------------------------------------------------------------------------

# 18ss
#--------------------------------------------------------------------------------
src.18ss.integrated.merge_20231223 = merge(src.18ss.integrated, 
                                           c(src.sm3.endoderm_Ngn3[,src.sm3.endoderm_Ngn3$Time%in%"18ss"],
                                             src.sm3.merge[,src.sm3.merge$Time%in%"18ss"]))
src.18ss.integrated.merge_20231223 = src.18ss.integrated.merge_20231223[,src.18ss.integrated.merge_20231223$nFeature_RNA>2500]

anchor.18ss.integrated.merge_20231223 = FindIntegrationAnchors(c(src.18ss.integrated.merge_20231223[,src.18ss.integrated.merge_20231223$Time%in%"ss18"&
                                                                                                      src.18ss.integrated.merge_20231223$batch%in%1],
                                                                 src.18ss.integrated.merge_20231223[,src.18ss.integrated.merge_20231223$Time%in%"ss18"&
                                                                                                      src.18ss.integrated.merge_20231223$batch%in%2],
                                                                 src.18ss.integrated.merge_20231223[,!src.18ss.integrated.merge_20231223$Time%in%"ss18"]),
                                                               reduction = "cca",
                                                               anchor.features = src.18ss.integrated.merge.selectgene,
                                                               dims = 1:30)

src.18ss.integrated.merge_20231223.re = IntegrateData(anchorset = anchor.18ss.integrated.merge_20231223 , dims = 1:30)
DefaultAssay(src.18ss.integrated.merge_20231223.re) = "integrated"
src.18ss.integrated.merge_20231223.re = ScaleData(src.18ss.integrated.merge_20231223.re, features = rownames(src.18ss.integrated.merge_20231223.re))
src.18ss.integrated.merge_20231223.re = RunPCA(src.18ss.integrated.merge_20231223.re, features = rownames(src.18ss.integrated@reductions$pca@feature.loadings))
src.18ss.integrated.merge_20231223.re = FindNeighbors(src.18ss.integrated.merge_20231223.re, dims = 1:30)
src.18ss.integrated.merge_20231223.re = FindClusters(src.18ss.integrated.merge_20231223.re, resolution = 2)
src.18ss.integrated.merge_20231223.re = FindClusters(src.18ss.integrated.merge_20231223.re, resolution = 3)
src.18ss.integrated.merge_20231223.re = FindClusters(src.18ss.integrated.merge_20231223.re, resolution = 4)
src.18ss.integrated.merge_20231223.re = RunUMAP(src.18ss.integrated.merge_20231223.re, dims = 1:20, n.neighbors = 100)

src.18ss.integrated.merge_20231223@reductions$umap_integrated = src.18ss.integrated.merge_20231223.re@reductions$umap
src.18ss.integrated.merge_20231223@reductions$pca_integrated = src.18ss.integrated.merge_20231223.re@reductions$pca
src.18ss.integrated.merge_20231223@reductions$umap_integrated@assay.used = "RNA"
src.18ss.integrated.merge_20231223@reductions$pca_integrated@assay.used = "RNA"
rm(src.18ss.integrated.merge_20231223.re)

src.18ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings = cbind(
  -src.18ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings[,1],
  -src.18ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings[,2])
colnames(src.18ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings) = c("UMAP_1","UMAP_2")
DimPlot(src.18ss.integrated.merge_20231223, group.by = "lineage", reduction = "umap_integrated",
        cols = color.lineage, na.value = "#eeeeee")

a = FNN::knn(src.18ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings[
  rownames(src.18ss.integrated.merge_20231223@meta.data[src.18ss.integrated.merge_20231223$Time%in%"ss18"&
                                                          !src.18ss.integrated.merge_20231223$cluster.revise.re.v1.30.re%in%NA,]),],
  src.18ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings[
    rownames(src.18ss.integrated.merge_20231223@meta.data[!src.18ss.integrated.merge_20231223$Time%in%"ss18",]),],
  src.18ss.integrated.merge_20231223@meta.data[src.18ss.integrated.merge_20231223$Time%in%"ss18"&
                                                 !src.18ss.integrated.merge_20231223$cluster.revise.re.v1.30.re%in%NA,]$cluster.revise.re.v1.30.re, k = 10)

src.18ss.integrated.merge_20231223$cluster.predict.umap_int.ext.v1.1 = src.18ss.integrated.merge_20231223$cluster.revise.re.v1.30.re
src.18ss.integrated.merge_20231223@meta.data[!src.18ss.integrated.merge_20231223$Time%in%"ss18",]$cluster.predict.umap_int.ext.v1.1 = as.character(a)
src.sm3.endoderm_Ngn3@meta.data[src.sm3.endoderm_Ngn3$Time%in%"18ss",]$cluster.predict.umap_int.ext.v1.1 = 
  src.18ss.integrated.merge_20231223@meta.data[rownames(
    src.sm3.endoderm_Ngn3@meta.data[src.sm3.endoderm_Ngn3$Time%in%"18ss",]),]$cluster.predict.umap_int.ext.v1.1
#--------------------------------------------------------------------------------

# 21ss
#--------------------------------------------------------------------------------
src.21ss.integrated.merge_20231223 = merge(src.21ss.integrated, 
                                           c(src.sm3.endoderm_Ngn3[,src.sm3.endoderm_Ngn3$Time%in%"21ss"],
                                             src.sm3.merge[,src.sm3.merge$Time%in%"21ss"]))
src.21ss.integrated.merge_20231223 = src.21ss.integrated.merge_20231223[,src.21ss.integrated.merge_20231223$nFeature_RNA>2500]

anchor.21ss.integrated.merge_20231223 = FindIntegrationAnchors(c(src.21ss.integrated.merge_20231223[,src.21ss.integrated.merge_20231223$Time%in%"ss21"&
                                                                                                      src.21ss.integrated.merge_20231223$batch%in%1],
                                                                 src.21ss.integrated.merge_20231223[,src.21ss.integrated.merge_20231223$Time%in%"ss21"&
                                                                                                      src.21ss.integrated.merge_20231223$batch%in%2],
                                                                 src.21ss.integrated.merge_20231223[,!src.21ss.integrated.merge_20231223$Time%in%"ss21"]),
                                                               reduction = "cca",
                                                               anchor.features = src.21ss.integrated.merge.selectgene,
                                                               dims = 1:30)

src.21ss.integrated.merge_20231223.re = IntegrateData(anchorset = anchor.21ss.integrated.merge_20231223 , dims = 1:30)
DefaultAssay(src.21ss.integrated.merge_20231223.re) = "integrated"
src.21ss.integrated.merge_20231223.re = ScaleData(src.21ss.integrated.merge_20231223.re, features = rownames(src.21ss.integrated.merge_20231223.re))
src.21ss.integrated.merge_20231223.re = RunPCA(src.21ss.integrated.merge_20231223.re, features = rownames(src.21ss.integrated@reductions$pca@feature.loadings))
src.21ss.integrated.merge_20231223.re = FindNeighbors(src.21ss.integrated.merge_20231223.re, dims = 1:30)
src.21ss.integrated.merge_20231223.re = FindClusters(src.21ss.integrated.merge_20231223.re, resolution = 2)
src.21ss.integrated.merge_20231223.re = FindClusters(src.21ss.integrated.merge_20231223.re, resolution = 3)
src.21ss.integrated.merge_20231223.re = FindClusters(src.21ss.integrated.merge_20231223.re, resolution = 4)
src.21ss.integrated.merge_20231223.re = RunUMAP(src.21ss.integrated.merge_20231223.re, dims = 1:20, n.neighbors = 100)

src.21ss.integrated.merge_20231223@reductions$umap_integrated = src.21ss.integrated.merge_20231223.re@reductions$umap
src.21ss.integrated.merge_20231223@reductions$pca_integrated = src.21ss.integrated.merge_20231223.re@reductions$pca
src.21ss.integrated.merge_20231223@reductions$umap_integrated@assay.used = "RNA"
src.21ss.integrated.merge_20231223@reductions$pca_integrated@assay.used = "RNA"
rm(src.21ss.integrated.merge_20231223.re)

src.21ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings = cbind(
  -src.21ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings[,1],
  -src.21ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings[,2])
colnames(src.21ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings) = c("UMAP_1","UMAP_2")
DimPlot(src.21ss.integrated.merge_20231223, group.by = "lineage", reduction = "umap_integrated",
        cols = color.lineage, na.value = "#eeeeee")

a = FNN::knn(src.21ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings[
  rownames(src.21ss.integrated.merge_20231223@meta.data[src.21ss.integrated.merge_20231223$Time%in%"ss21"&
                                                          !src.21ss.integrated.merge_20231223$cluster.revise.re.v1.30.re%in%NA,]),],
  src.21ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings[
    rownames(src.21ss.integrated.merge_20231223@meta.data[!src.21ss.integrated.merge_20231223$Time%in%"ss21",]),],
  src.21ss.integrated.merge_20231223@meta.data[src.21ss.integrated.merge_20231223$Time%in%"ss21"&
                                                 !src.21ss.integrated.merge_20231223$cluster.revise.re.v1.30.re%in%NA,]$cluster.revise.re.v1.30.re, k = 10)

src.21ss.integrated.merge_20231223$cluster.predict.umap_int.ext.v1.1 = src.21ss.integrated.merge_20231223$cluster.revise.re.v1.30.re
src.21ss.integrated.merge_20231223@meta.data[!src.21ss.integrated.merge_20231223$Time%in%"ss21",]$cluster.predict.umap_int.ext.v1.1 = as.character(a)
src.sm3.endoderm_Ngn3@meta.data[src.sm3.endoderm_Ngn3$Time%in%"21ss",]$cluster.predict.umap_int.ext.v1.1 = 
  src.21ss.integrated.merge_20231223@meta.data[rownames(
    src.sm3.endoderm_Ngn3@meta.data[src.sm3.endoderm_Ngn3$Time%in%"21ss",]),]$cluster.predict.umap_int.ext.v1.1
#--------------------------------------------------------------------------------

# 24ss
#--------------------------------------------------------------------------------
src.24ss.integrated.merge_20231223 = merge(src.24ss.integrated, 
                                           c(src.sm3.endoderm_Ngn3[,src.sm3.endoderm_Ngn3$Time%in%"24ss"],
                                             src.sm3.merge[,src.sm3.merge$Time%in%"24ss"]))
src.24ss.integrated.merge_20231223 = src.24ss.integrated.merge_20231223[,src.24ss.integrated.merge_20231223$nFeature_RNA>2500]

anchor.24ss.integrated.merge_20231223 = FindIntegrationAnchors(c(src.24ss.integrated.merge_20231223[,src.24ss.integrated.merge_20231223$Time%in%"ss24"&
                                                                                                      src.24ss.integrated.merge_20231223$batch%in%1],
                                                                 src.24ss.integrated.merge_20231223[,src.24ss.integrated.merge_20231223$Time%in%"ss24"&
                                                                                                      src.24ss.integrated.merge_20231223$batch%in%2],
                                                                 src.24ss.integrated.merge_20231223[,!src.24ss.integrated.merge_20231223$Time%in%"ss24"]),
                                                               reduction = "cca",
                                                               anchor.features = src.24ss.integrated.merge.selectgene,
                                                               dims = 1:30)

src.24ss.integrated.merge_20231223.re = IntegrateData(anchorset = anchor.24ss.integrated.merge_20231223 , dims = 1:30)
DefaultAssay(src.24ss.integrated.merge_20231223.re) = "integrated"
src.24ss.integrated.merge_20231223.re = ScaleData(src.24ss.integrated.merge_20231223.re, features = rownames(src.24ss.integrated.merge_20231223.re))
src.24ss.integrated.merge_20231223.re = RunPCA(src.24ss.integrated.merge_20231223.re, features = rownames(src.24ss.integrated@reductions$pca@feature.loadings))
src.24ss.integrated.merge_20231223.re = FindNeighbors(src.24ss.integrated.merge_20231223.re, dims = 1:30)
src.24ss.integrated.merge_20231223.re = FindClusters(src.24ss.integrated.merge_20231223.re, resolution = 2)
src.24ss.integrated.merge_20231223.re = FindClusters(src.24ss.integrated.merge_20231223.re, resolution = 3)
src.24ss.integrated.merge_20231223.re = FindClusters(src.24ss.integrated.merge_20231223.re, resolution = 4)
src.24ss.integrated.merge_20231223.re = RunUMAP(src.24ss.integrated.merge_20231223.re, dims = 1:20, n.neighbors = 100)

src.24ss.integrated.merge_20231223@reductions$umap_integrated = src.24ss.integrated.merge_20231223.re@reductions$umap
src.24ss.integrated.merge_20231223@reductions$pca_integrated = src.24ss.integrated.merge_20231223.re@reductions$pca
src.24ss.integrated.merge_20231223@reductions$umap_integrated@assay.used = "RNA"
src.24ss.integrated.merge_20231223@reductions$pca_integrated@assay.used = "RNA"
rm(src.24ss.integrated.merge_20231223.re)

src.24ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings = cbind(
  -src.24ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings[,1],
  -src.24ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings[,2])
colnames(src.24ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings) = c("UMAP_1","UMAP_2")
DimPlot(src.24ss.integrated.merge_20231223, group.by = "lineage", reduction = "umap_integrated",
        cols = color.lineage, na.value = "#eeeeee")

a = FNN::knn(src.24ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings[
  rownames(src.24ss.integrated.merge_20231223@meta.data[src.24ss.integrated.merge_20231223$Time%in%"ss24"&
                                                          !src.24ss.integrated.merge_20231223$cluster.revise.re.v1.30.re%in%NA,]),],
  src.24ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings[
    rownames(src.24ss.integrated.merge_20231223@meta.data[!src.24ss.integrated.merge_20231223$Time%in%"ss24",]),],
  src.24ss.integrated.merge_20231223@meta.data[src.24ss.integrated.merge_20231223$Time%in%"ss24"&
                                                 !src.24ss.integrated.merge_20231223$cluster.revise.re.v1.30.re%in%NA,]$cluster.revise.re.v1.30.re, k = 10)

src.24ss.integrated.merge_20231223$cluster.predict.umap_int.ext.v1.1 = src.24ss.integrated.merge_20231223$cluster.revise.re.v1.30.re
src.24ss.integrated.merge_20231223@meta.data[!src.24ss.integrated.merge_20231223$Time%in%"ss24",]$cluster.predict.umap_int.ext.v1.1 = as.character(a)
src.sm3.endoderm_Ngn3@meta.data[src.sm3.endoderm_Ngn3$Time%in%"24ss",]$cluster.predict.umap_int.ext.v1.1 = 
  src.24ss.integrated.merge_20231223@meta.data[rownames(
    src.sm3.endoderm_Ngn3@meta.data[src.sm3.endoderm_Ngn3$Time%in%"24ss",]),]$cluster.predict.umap_int.ext.v1.1
#--------------------------------------------------------------------------------

# 27ss
#--------------------------------------------------------------------------------
src.27ss.integrated.merge_20231223 = merge(src.27ss.integrated, 
                                           c(src.sm3.endoderm_Ngn3[,src.sm3.endoderm_Ngn3$Time%in%"27ss"],
                                             src.sm3.merge[,src.sm3.merge$Time%in%"27ss"]))
src.27ss.integrated.merge_20231223 = src.27ss.integrated.merge_20231223[,src.27ss.integrated.merge_20231223$nFeature_RNA>2500]

anchor.27ss.integrated.merge_20231223 = FindIntegrationAnchors(c(src.27ss.integrated.merge_20231223[,src.27ss.integrated.merge_20231223$Time%in%"ss27"&
                                                                                                      src.27ss.integrated.merge_20231223$batch%in%1],
                                                                 src.27ss.integrated.merge_20231223[,src.27ss.integrated.merge_20231223$Time%in%"ss27"&
                                                                                                      src.27ss.integrated.merge_20231223$batch%in%2],
                                                                 src.27ss.integrated.merge_20231223[,!src.27ss.integrated.merge_20231223$Time%in%"ss27"]),
                                                               reduction = "cca",
                                                               anchor.features = src.27ss.integrated.merge.selectgene,
                                                               dims = 1:30)

src.27ss.integrated.merge_20231223.re = IntegrateData(anchorset = anchor.27ss.integrated.merge_20231223 , dims = 1:30)
DefaultAssay(src.27ss.integrated.merge_20231223.re) = "integrated"
src.27ss.integrated.merge_20231223.re = ScaleData(src.27ss.integrated.merge_20231223.re, features = rownames(src.27ss.integrated.merge_20231223.re))
src.27ss.integrated.merge_20231223.re = RunPCA(src.27ss.integrated.merge_20231223.re, features = rownames(src.27ss.integrated@reductions$pca@feature.loadings))
src.27ss.integrated.merge_20231223.re = FindNeighbors(src.27ss.integrated.merge_20231223.re, dims = 1:30)
src.27ss.integrated.merge_20231223.re = FindClusters(src.27ss.integrated.merge_20231223.re, resolution = 2)
src.27ss.integrated.merge_20231223.re = FindClusters(src.27ss.integrated.merge_20231223.re, resolution = 3)
src.27ss.integrated.merge_20231223.re = FindClusters(src.27ss.integrated.merge_20231223.re, resolution = 4)
src.27ss.integrated.merge_20231223.re = RunUMAP(src.27ss.integrated.merge_20231223.re, dims = 1:20, n.neighbors = 100)

src.27ss.integrated.merge_20231223@reductions$umap_integrated = src.27ss.integrated.merge_20231223.re@reductions$umap
src.27ss.integrated.merge_20231223@reductions$pca_integrated = src.27ss.integrated.merge_20231223.re@reductions$pca
src.27ss.integrated.merge_20231223@reductions$umap_integrated@assay.used = "RNA"
src.27ss.integrated.merge_20231223@reductions$pca_integrated@assay.used = "RNA"
rm(src.27ss.integrated.merge_20231223.re)

src.27ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings = cbind(
  -src.27ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings[,1],
  -src.27ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings[,2])
colnames(src.27ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings) = c("UMAP_1","UMAP_2")
DimPlot(src.27ss.integrated.merge_20231223, group.by = "lineage", reduction = "umap_integrated",
        cols = color.lineage, na.value = "#eeeeee")

a = FNN::knn(src.27ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings[
  rownames(src.27ss.integrated.merge_20231223@meta.data[src.27ss.integrated.merge_20231223$Time%in%"ss27"&
                                                          !src.27ss.integrated.merge_20231223$cluster.revise.re.v1.30.re%in%NA,]),],
  src.27ss.integrated.merge_20231223@reductions$umap_integrated@cell.embeddings[
    rownames(src.27ss.integrated.merge_20231223@meta.data[!src.27ss.integrated.merge_20231223$Time%in%"ss27",]),],
  src.27ss.integrated.merge_20231223@meta.data[src.27ss.integrated.merge_20231223$Time%in%"ss27"&
                                                 !src.27ss.integrated.merge_20231223$cluster.revise.re.v1.30.re%in%NA,]$cluster.v06.26.re..merge, k = 10)

src.27ss.integrated.merge_20231223$cluster.predict.umap_int.ext.v1.1 = src.27ss.integrated.merge_20231223$cluster.v06.26.re..merge
src.27ss.integrated.merge_20231223@meta.data[!src.27ss.integrated.merge_20231223$Time%in%"ss27",]$cluster.predict.umap_int.ext.v1.1 = as.character(a)
src.sm3.endoderm_Ngn3@meta.data[src.sm3.endoderm_Ngn3$Time%in%"27ss",]$cluster.predict.umap_int.ext.v1.1 = 
  src.27ss.integrated.merge_20231223@meta.data[rownames(
    src.sm3.endoderm_Ngn3@meta.data[src.sm3.endoderm_Ngn3$Time%in%"27ss",]),]$cluster.predict.umap_int.ext.v1.1
#--------------------------------------------------------------------------------

rm(src.15ss.integrated.merge_20231223,anchor.15ss.integrated.merge_20231223,
   src.18ss.integrated.merge_20231223,anchor.18ss.integrated.merge_20231223,
   src.21ss.integrated.merge_20231223,anchor.21ss.integrated.merge_20231223,
   src.24ss.integrated.merge_20231223,anchor.24ss.integrated.merge_20231223,
   src.27ss.integrated.merge_20231223,anchor.27ss.integrated.merge_20231223,
   src.27ss.integrated.merge_20231223,anchor.27ss.integrated.merge_20231223)



src.sm3.endoderm_Ngn3Cre = FindNeighbors(src.sm3.endoderm_Ngn3Cre, dims = 1:30)
src.sm3.endoderm_Ngn3Cre = FindClusters(src.sm3.endoderm_Ngn3Cre, resolution = 1.75)
DimPlot(src.sm3.endoderm_Ngn3Cre, group.by = "RNA_snn_res.1.5", 
        pt.size = 2, label = T, label.size = 6) + 
  DimPlot(src.sm3.endoderm_Ngn3Cre, group.by = "Time", 
          pt.size = 2, cols = colors.time.2) +
  FeaturePlot(src.sm3.endoderm_Ngn3Cre, features = c("Gcg","Pdx1","Neurog3","Sim1","Mnx1"))
table(src.sm3.endoderm_Ngn3Cre$cluster.refine)

src.sm3.endoderm_Ngn3Cre$cluster.refine = NA
src.sm3.endoderm_Ngn3Cre@meta.data[src.sm3.endoderm_Ngn3Cre$RNA_snn_res.1.5%in%c(2,3,4),]$cluster.refine = "EP.1"
src.sm3.endoderm_Ngn3Cre@meta.data[src.sm3.endoderm_Ngn3Cre$RNA_snn_res.1.5%in%c(0,5,6),]$cluster.refine = "EP.2"
src.sm3.endoderm_Ngn3Cre@meta.data[src.sm3.endoderm_Ngn3Cre$RNA_snn_res.1.5%in%c(1),]$cluster.refine = "DP"


pdf("figure.v08.07/Ngn3_process/Ngn3Cre_summary.pdf",9,7)
for(seurat in c("src.sm3.endoderm_Ngn3Cre")){
  seurat = get(seurat)
  print(DimPlot(seurat, group.by = "cluster.refine",
                reduction = "umap",
                cols = c(colors.type.lineage, color.lineage,
                         cluster.endoderm.color.v5), 
                pt.size = 6)+
          theme_classic() + theme(aspect.ratio = 1) + p_add) 
  
  print(DimPlot(seurat, group.by = "lineage",
                reduction = "umap",
                cols = c(colors.type.lineage, color.lineage,
                         cluster.endoderm.color.v5), 
                pt.size = 6)+
          theme_classic() + theme(aspect.ratio = 1) + p_add) 
  
  print(DimPlot(seurat, group.by = "Time",
                reduction = "umap",
                cols = c(colors.type.lineage, color.lineage,
                         cluster.endoderm.color.v5,colors.time.2), 
                pt.size = 6)+
          theme_classic() + theme(aspect.ratio = 1) + p_add) 
}
dev.off()



src.sm3.endoderm_Ngn3GFP = FindNeighbors(src.sm3.endoderm_Ngn3GFP, dims = 1:30)
src.sm3.endoderm_Ngn3GFP = FindClusters(src.sm3.endoderm_Ngn3GFP, resolution = 3.5)
DimPlot(src.sm3.endoderm_Ngn3GFP, group.by = "RNA_snn_res.3", 
        pt.size = 2, label = T, label.size = 6) + 
  DimPlot(src.sm3.endoderm_Ngn3GFP, group.by = "SeqDate", 
          pt.size = 2) 
FeaturePlot(src.sm3.endoderm_Ngn3GFP,
            features = c("Gcg","Pdx1","Neurog3", "Sim1","Mnx1",
                         "Osr2","Pitx2","Cdx2"))
table(src.sm3.endoderm_Ngn3GFP$cluster.predict.umap_int.ext.v1.1)
FeaturePlot(src.21ss.integrated, reduction = 'umap', features = "Neurog3")
table(src.sm3.endoderm_Ngn3GFP$lineage)

src.sm3.endoderm_Ngn3GFP$cluster.refine = NA
src.sm3.endoderm_Ngn3GFP@meta.data[src.sm3.endoderm_Ngn3GFP$RNA_snn_res.3%in%c(3,11,10,7,8,2,14),]$cluster.refine = "Small.intestine.1"
src.sm3.endoderm_Ngn3GFP@meta.data[src.sm3.endoderm_Ngn3GFP$RNA_snn_res.3%in%c(5),]$cluster.refine = "Small.intestine.2"
src.sm3.endoderm_Ngn3GFP@meta.data[src.sm3.endoderm_Ngn3GFP$RNA_snn_res.3%in%c(0)|
                                     src.sm3.endoderm_Ngn3GFP$RNA_snn_res.3.5%in%c(4),]$cluster.refine = "MG.3.M"
src.sm3.endoderm_Ngn3GFP@meta.data[src.sm3.endoderm_Ngn3GFP$RNA_snn_res.3%in%c(9),]$cluster.refine = "EP.1"
src.sm3.endoderm_Ngn3GFP@meta.data[src.sm3.endoderm_Ngn3GFP$RNA_snn_res.3%in%c(6),]$cluster.refine = "EP.2"
src.sm3.endoderm_Ngn3GFP@meta.data[src.sm3.endoderm_Ngn3GFP$RNA_snn_res.3%in%c(1,4,13,12),]$cluster.refine = "DP"

DimPlot(src.sm3.endoderm_Ngn3GFP, group.by = "cluster.refine",
        cols = cluster.endoderm.color.v5,
        pt.size = 2, label = F, label.size = 6) +
  theme_classic() + theme(aspect.ratio = 1) + p_add



pdf("figure.v08.07/Ngn3_process/Ngn3GFP_summary.pdf",9,7)
for(seurat in c("src.sm3.endoderm_Ngn3GFP")){
  seurat = get(seurat)
  print(DimPlot(seurat, group.by = "cluster.refine",
                reduction = "umap",
                cols = c(colors.type.lineage, color.lineage,
                         cluster.endoderm.color.v5), 
                pt.size = 4)+
          theme_classic() + theme(aspect.ratio = 1) + p_add) 
  
  print(DimPlot(seurat, group.by = "lineage",
                reduction = "umap",
                cols = c(colors.type.lineage, color.lineage,
                         cluster.endoderm.color.v5), 
                pt.size = 4)+
          theme_classic() + theme(aspect.ratio = 1) + p_add) 
  
  print(DimPlot(seurat, group.by = "Time",
                reduction = "umap",
                cols = c(colors.type.lineage, color.lineage,
                         cluster.endoderm.color.v5,colors.time.2), 
                pt.size = 4)+
          theme_classic() + theme(aspect.ratio = 1) + p_add) 
}
dev.off()








