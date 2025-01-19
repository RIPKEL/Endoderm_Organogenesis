#-------------------------------
#>> 11.2 Ngn3-GFP
#--------------------------------

#>> after pre-process for Smartseq-3

src.sm3.update.rawdata_Ngn3GFP$index = "HM"
src.sm3.update.rawdata_Ngn3GFP@meta.data[
  src.sm3.update.rawdata_Ngn3GFP$Mouse%in%c("Neurog3GFP_HZ","Ngn3GFP_HMHZ","Ngn3GFP_HZ"),]$index = "HZ"

mat_sm3_raw = 
  src.sm3.update.rawdata_Ngn3GFP@meta.data[,c("index","nFeature_RNA","cluster")]
mat_sm3_raw$index = factor(mat_sm3_raw$index,
                           levels = c("HZ","HM"))
mat_sm3_raw$cluster = factor(mat_sm3_raw$cluster,
                             levels = names(qc.color))


pdf("figure.v08.07/Ngn3GFP/ngn3gfp_qc.pdf",6,6)
ggplot() +
  geom_hline(yintercept = 6000, linetype = "dashed", colour = "black")+
  geom_violin(data = mat_sm3_raw, 
              mapping = aes(x = index, y = nFeature_RNA, 
                            group = index, fill = index))+
  scale_fill_manual(values = colors.type)+
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=90, size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        plot.title = element_blank(),
        aspect.ratio=1/0.3) +
  xlab("Type") + ylab("The number of\n detected gene")
ggplot() +
  #geom_vline(xintercept = 2.5, linetype = "dashed", colour = "black")+
  geom_bar(data = mat_sm3_raw, 
           mapping = aes(x = index,
                         group = cluster, fill = cluster),
           stat = "count", position = 'fill')+
  # geom_text(mapping = aes(
  #   x = rep(1.0,12), y = c(1:12),
  #   label = paste("(",sapply(table(mat_10x_raw$index), as.character)[1:12],")",sep="")))+
  annotate("text",  y = rep(1.1,2), x = c(1:2), angle=90,
           label = paste("(",sapply(table(mat_sm3_raw$index), as.character)[c(1,2)],")",sep="")) +
  scale_fill_manual(values = qc.color)+
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=90, size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        plot.title = element_blank(),
        aspect.ratio=1/0.3) +
  xlab("Type") + ylab("Cluster Proportion of QC")
dev.off()


src.sm3.Ngn3GFP = 
  src.sm3.update.rawdata_Ngn3GFP[,src.sm3.update.rawdata_Ngn3GFP$cluster%in%"Endoderm"]
for(seurat in c("src.sm3.Ngn3GFP")){
  
  assign(seurat,
         FindVariableFeatures(get(seurat), selection.method = "vst", nfeatures = 2000))
  VariableFeaturePlot(get(seurat))
  assign(seurat,
         ScaleData(get(seurat),features = VariableFeatures(get(seurat))))
  # assign(seurat,
  #        MyDoubletFinder(get(seurat),round(ncol(get(seurat)) * 0.1,0)))
  assign(seurat,
         RunPCA(get(seurat), features = VariableFeatures(get(seurat))))
  assign(seurat,
         RunTSNE(get(seurat), dims = 1:30,check_duplicates = F,
                 perplexity=round((30+ncol(get(seurat))/100)),
                 max_iter = round((ncol(get(seurat))/12))))
  assign(seurat,
         RunUMAP(get(seurat), dims = 1:30,
                 n.neighbors = 30))
  
  assign(seurat,
         FindNeighbors(get(seurat), dims = 1:30))
  assign(seurat,
         FindClusters(get(seurat), resolution = 1.2))
  assign(seurat,
         FindClusters(get(seurat), resolution = 3))
  assign(seurat,
         FindClusters(get(seurat), resolution = 1.5))
}


##-------------------------------------------------------------------------------
#                  Annotation 
##-------------------------------------------------------------------------------
src.sm3.Ngn3GFP_add = merge(src.sm3.Ngn3GFP, src.sm3.Mnx1)
src.sm3.update.9ss = src.sm3.Ngn3GFP_add[,src.sm3.Ngn3GFP_add$Time%in%"9ss"]
src.sm3.update.12ss = src.sm3.Ngn3GFP_add[,src.sm3.Ngn3GFP_add$Time%in%"12ss"]
src.sm3.update.15ss = src.sm3.Ngn3GFP_add[,src.sm3.Ngn3GFP_add$Time%in%"15ss"]
src.sm3.update.18ss = src.sm3.Ngn3GFP_add[,src.sm3.Ngn3GFP_add$Time%in%"18ss"]
src.sm3.update.21ss = src.sm3.Ngn3GFP_add[,src.sm3.Ngn3GFP_add$Time%in%"21ss"]
src.sm3.update.24ss = src.sm3.Ngn3GFP_add[,src.sm3.Ngn3GFP_add$Time%in%"24ss"]
src.sm3.update.27ss = src.sm3.Ngn3GFP_add[,src.sm3.Ngn3GFP_add$Time%in%"27ss"]

########################
## 1.9ss 
########################
src.9ss.integrated.merge_230909 = merge(src.9ss.integrated, src.sm3.update.9ss)
src.9ss.integrated.merge_230909 = ScaleData(src.9ss.integrated.merge_230909, 
                                            features = rownames(src.9ss.integrated.merge_230909))
src.9ss.integrated.merge_230909 = RunPCA(src.9ss.integrated.merge_230909, features = src.9ss.integrated.selectgene)
src.9ss.integrated.merge_230909 = RunUMAP(src.9ss.integrated.merge_230909, dims = 1:20, n.neighbors = 100)

#-------- MNN corrected in RNA level--------
src.9ss.integrated.merge_230909.selectgene = src.9ss.integrated.selectgene
MNN.res = 
  mnnCorrect(as.matrix(src.9ss.integrated.merge_230909@assays$RNA@data[src.9ss.integrated.merge_230909.selectgene, src.9ss.integrated.merge_230909$Time%in%"ss9"]),
             as.matrix(src.9ss.integrated.merge_230909@assays$RNA@data[src.9ss.integrated.merge_230909.selectgene, src.9ss.integrated.merge_230909$Time%in%"9ss"]),
             k = 5,cos.norm.out=F)
src.9ss.integrated.merge_230909@assays$mnnRNA = CreateAssayObject(data = MNN.res@assays@data$corrected)
src.9ss.integrated.merge_230909@assays$mnnRNA@key = "mnn_"
src.9ss.integrated.merge_230909 = ScaleData(src.9ss.integrated.merge_230909, 
                                            rownames(src.9ss.integrated.merge_230909@assays$mnnRNA),
                                            assay = "mnnRNA")
#----------------

anchor.9ss.integrated = 
  FindTransferAnchors(reference = src.9ss.integrated.merge_230909[,colnames(src.9ss.integrated)],
                      query = src.9ss.integrated.merge_230909[,colnames(src.sm3.update.9ss)],
                      reference.assay = "mnnRNA",
                      query.assay = "mnnRNA",scale = T,
                      features = src.9ss.integrated.selectgene)


#############
# pca
pca.transfer.9ss.integrated = TransferData(anchor.9ss.integrated,
                                           t(src.9ss.integrated@reductions$pca@cell.embeddings))
src.9ss.integrated.merge_230909@reductions$pca_fta = src.9ss.integrated.merge_230909@reductions$pca
src.9ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.9ss.integrated),] = src.9ss.integrated@reductions$pca@cell.embeddings
src.9ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.sm3.update.9ss),] = as.matrix(t(pca.transfer.9ss.integrated@data))

# tsne
tsne.transfer.9ss.integrated = TransferData(anchor.9ss.integrated,
                                            t(src.9ss.integrated@reductions$tsne@cell.embeddings))
src.9ss.integrated.merge_230909@reductions$tsne_fta = src.9ss.integrated.merge_230909@reductions$umap
src.9ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.9ss.integrated),] = src.9ss.integrated@reductions$tsne@cell.embeddings
src.9ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.sm3.update.9ss),] = as.matrix(t(tsne.transfer.9ss.integrated@data))

# umap
umap.transfer.9ss.integrated = TransferData(anchor.9ss.integrated,
                                            t(src.9ss.integrated@reductions$umap@cell.embeddings))
src.9ss.integrated.merge_230909@reductions$umap_fta = src.9ss.integrated.merge_230909@reductions$umap
src.9ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.9ss.integrated),] = src.9ss.integrated@reductions$umap@cell.embeddings
src.9ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.sm3.update.9ss),] = as.matrix(t(umap.transfer.9ss.integrated@data))
#############

#################
# predict.v06.26.re
#################
# pca
pca.transfer.9ss.integrated = TransferData(anchor.9ss.integrated,
                                           t(src.9ss.integrated@reductions$pca@cell.embeddings))
src.9ss.integrated.merge_230909@reductions$pca_fta = src.9ss.integrated.merge_230909@reductions$pca
src.9ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.9ss.integrated),] = src.9ss.integrated@reductions$pca@cell.embeddings
src.9ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.sm3.update.9ss),] = as.matrix(t(pca.transfer.9ss.integrated@data))

src.9ss.integrated.merge_230909@meta.data$cluster.predict.pca.v06.26.re =
  as.character(FNN::knn(src.9ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.9ss.integrated),],
                        src.9ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings,
                        as.character(src.9ss.integrated.merge_230909@meta.data[colnames(src.9ss.integrated),]$cluster.v06.26.re..merge), k = 10))
src.9ss.integrated.merge_230909@meta.data[colnames(src.9ss.integrated),]$cluster.predict.pca.v06.26.re =
  src.9ss.integrated.merge_230909@meta.data[colnames(src.9ss.integrated),]$cluster.predict.pca.v06.26.re  

src.sm3.update.9ss$cluster.9ss.predict.pca.v06.26.re = 
  src.9ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.9ss),]$cluster.predict.pca.v06.26.re
src.9ss.integrated.merge_230909$cluster.9ss.predict.pca.v06.26.re = NA
src.9ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.9ss),]$cluster.predict.pca.v06.26.re = 
  src.sm3.update.9ss$cluster.9ss.predict.pca.v06.26.re

# tsne
tsne.transfer.9ss.integrated = TransferData(anchor.9ss.integrated,
                                            t(src.9ss.integrated@reductions$tsne@cell.embeddings))
src.9ss.integrated.merge_230909@reductions$tsne_fta = src.9ss.integrated.merge_230909@reductions$tsne
src.9ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.9ss.integrated),] = src.9ss.integrated@reductions$tsne@cell.embeddings
src.9ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.sm3.update.9ss),] = as.matrix(t(tsne.transfer.9ss.integrated@data))

src.9ss.integrated.merge_230909@meta.data$cluster.predict.tsne.v06.26.re =
  as.character(FNN::knn(src.9ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.9ss.integrated),],
                        src.9ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings,
                        as.character(src.9ss.integrated.merge_230909@meta.data[colnames(src.9ss.integrated),]$cluster.v06.26.re..merge), k = 10))
src.9ss.integrated.merge_230909@meta.data[colnames(src.9ss.integrated),]$cluster.predict.tsne.v06.26.re =
  src.9ss.integrated.merge_230909@meta.data[colnames(src.9ss.integrated),]$cluster.predict.tsne.v06.26.re  

src.sm3.update.9ss$cluster.9ss.predict.tsne.v06.26.re = 
  src.9ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.9ss),]$cluster.predict.tsne.v06.26.re
src.9ss.integrated.merge_230909$cluster.9ss.predict.tsne.v06.26.re = NA
src.9ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.9ss),]$cluster.predict.tsne.v06.26.re = 
  src.sm3.update.9ss$cluster.9ss.predict.tsne.v06.26.re

# umap
umap.transfer.9ss.integrated = TransferData(anchor.9ss.integrated,
                                            t(src.9ss.integrated@reductions$umap@cell.embeddings))
src.9ss.integrated.merge_230909@reductions$umap_fta = src.9ss.integrated.merge_230909@reductions$umap
src.9ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.9ss.integrated),] = src.9ss.integrated@reductions$umap@cell.embeddings
src.9ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.sm3.update.9ss),] = as.matrix(t(umap.transfer.9ss.integrated@data))

src.9ss.integrated.merge_230909@meta.data$cluster.v06.26.re..merge =
  as.character(FNN::knn(src.9ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.9ss.integrated),],
                        src.9ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings,
                        as.character(src.9ss.integrated.merge_230909@meta.data[colnames(src.9ss.integrated),]$cluster.v06.26.re..merge), k = 10))
src.9ss.integrated.merge_230909@meta.data[colnames(src.9ss.integrated),]$cluster.v06.26.re..merge =
  src.9ss.integrated.merge_230909@meta.data[colnames(src.9ss.integrated),]$cluster.v06.26.re..merge  

src.sm3.update.9ss$cluster.9ss.predict.umap.v06.26.re = 
  src.9ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.9ss),]$cluster.v06.26.re..merge

src.9ss.integrated.merge_230909$cluster.9ss.predict.umap.v06.26.re = NA
src.9ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.9ss),]$cluster.9ss.predict.umap.v06.26.re = 
  src.sm3.update.9ss$cluster.9ss.predict.umap.v06.26.re
##################


########################
## 2.12ss 
########################
src.12ss.integrated.merge_230909 = merge(src.12ss.integrated, src.sm3.update.12ss)
src.12ss.integrated.merge_230909 = ScaleData(src.12ss.integrated.merge_230909, features = rownames(src.12ss.integrated.merge_230909))
src.12ss.integrated.merge_230909 = RunPCA(src.12ss.integrated.merge_230909, features = src.12ss.integrated.selectgene)
src.12ss.integrated.merge_230909 = RunUMAP(src.12ss.integrated.merge_230909, dims = 1:20, n.neighbors = 100)

#-------- MNN corrected in RNA level--------
src.12ss.integrated.merge_230909.selectgene = src.12ss.integrated.selectgene
MNN.res = 
  mnnCorrect(as.matrix(src.12ss.integrated.merge_230909@assays$RNA@data[src.12ss.integrated.merge_230909.selectgene, src.12ss.integrated.merge_230909$Time%in%"ss12"]),
             as.matrix(src.12ss.integrated.merge_230909@assays$RNA@data[src.12ss.integrated.merge_230909.selectgene, src.12ss.integrated.merge_230909$Time%in%"12ss"]),
             k = 5,cos.norm.out=F)
src.12ss.integrated.merge_230909@assays$mnnRNA = CreateAssayObject(data = MNN.res@assays@data$corrected)
src.12ss.integrated.merge_230909@assays$mnnRNA@key = "mnn_"
src.12ss.integrated.merge_230909 = ScaleData(src.12ss.integrated.merge_230909, 
                                             rownames(src.12ss.integrated.merge_230909@assays$mnnRNA),
                                             assay = "mnnRNA")
#----------------

anchor.12ss.integrated = 
  FindTransferAnchors(reference = src.12ss.integrated.merge_230909[,colnames(src.12ss.integrated)],
                      query = src.12ss.integrated.merge_230909[,colnames(src.sm3.update.12ss)],
                      reference.assay = "mnnRNA",
                      query.assay = "mnnRNA",scale = T,
                      features = src.12ss.integrated.selectgene)

#############
# pca
pca.transfer.12ss.integrated = TransferData(anchor.12ss.integrated,
                                            t(src.12ss.integrated@reductions$pca@cell.embeddings))
src.12ss.integrated.merge_230909@reductions$pca_fta = src.12ss.integrated.merge_230909@reductions$pca
src.12ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.12ss.integrated),] = src.12ss.integrated@reductions$pca@cell.embeddings
src.12ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.sm3.update.12ss),] = as.matrix(t(pca.transfer.12ss.integrated@data))

# tsne
tsne.transfer.12ss.integrated = TransferData(anchor.12ss.integrated,
                                             t(src.12ss.integrated@reductions$tsne@cell.embeddings))
src.12ss.integrated.merge_230909@reductions$tsne_fta = src.12ss.integrated.merge_230909@reductions$umap
src.12ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.12ss.integrated),] = src.12ss.integrated@reductions$tsne@cell.embeddings
src.12ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.sm3.update.12ss),] = as.matrix(t(tsne.transfer.12ss.integrated@data))

# umap
umap.transfer.12ss.integrated = TransferData(anchor.12ss.integrated,
                                             t(src.12ss.integrated@reductions$umap@cell.embeddings))
src.12ss.integrated.merge_230909@reductions$umap_fta = src.12ss.integrated.merge_230909@reductions$umap
src.12ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.12ss.integrated),] = src.12ss.integrated@reductions$umap@cell.embeddings
src.12ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.sm3.update.12ss),] = as.matrix(t(umap.transfer.12ss.integrated@data))
#############

#################
# predict.v06.26.re
#################
# pca
pca.transfer.12ss.integrated = TransferData(anchor.12ss.integrated,
                                            t(src.12ss.integrated@reductions$pca@cell.embeddings))
src.12ss.integrated.merge_230909@reductions$pca_fta = src.12ss.integrated.merge_230909@reductions$pca
src.12ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.12ss.integrated),] = src.12ss.integrated@reductions$pca@cell.embeddings
src.12ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.sm3.update.12ss),] = as.matrix(t(pca.transfer.12ss.integrated@data))

src.12ss.integrated.merge_230909@meta.data$cluster.predict.pca.v06.26.re =
  as.character(FNN::knn(src.12ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.12ss.integrated),],
                        src.12ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings,
                        as.character(src.12ss.integrated.merge_230909@meta.data[colnames(src.12ss.integrated),]$cluster.v06.26.re..merge), k = 10))
src.12ss.integrated.merge_230909@meta.data[colnames(src.12ss.integrated),]$cluster.predict.pca.v06.26.re =
  src.12ss.integrated.merge_230909@meta.data[colnames(src.12ss.integrated),]$cluster.predict.pca.v06.26.re  

src.sm3.update.12ss$cluster.12ss.predict.pca.v06.26.re = 
  src.12ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.12ss),]$cluster.predict.pca.v06.26.re
src.12ss.integrated.merge_230909$cluster.12ss.predict.pca.v06.26.re = NA
src.12ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.12ss),]$cluster.predict.pca.v06.26.re = 
  src.sm3.update.12ss$cluster.12ss.predict.pca.v06.26.re

# tsne
tsne.transfer.12ss.integrated = TransferData(anchor.12ss.integrated,
                                             t(src.12ss.integrated@reductions$tsne@cell.embeddings))
src.12ss.integrated.merge_230909@reductions$tsne_fta = src.12ss.integrated.merge_230909@reductions$tsne
src.12ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.12ss.integrated),] = src.12ss.integrated@reductions$tsne@cell.embeddings
src.12ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.sm3.update.12ss),] = as.matrix(t(tsne.transfer.12ss.integrated@data))

src.12ss.integrated.merge_230909@meta.data$cluster.predict.tsne.v06.26.re =
  as.character(FNN::knn(src.12ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.12ss.integrated),],
                        src.12ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings,
                        as.character(src.12ss.integrated.merge_230909@meta.data[colnames(src.12ss.integrated),]$cluster.v06.26.re..merge), k = 10))
src.12ss.integrated.merge_230909@meta.data[colnames(src.12ss.integrated),]$cluster.predict.tsne.v06.26.re =
  src.12ss.integrated.merge_230909@meta.data[colnames(src.12ss.integrated),]$cluster.predict.tsne.v06.26.re  

src.sm3.update.12ss$cluster.12ss.predict.tsne.v06.26.re = 
  src.12ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.12ss),]$cluster.predict.tsne.v06.26.re
src.12ss.integrated.merge_230909$cluster.12ss.predict.tsne.v06.26.re = NA
src.12ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.12ss),]$cluster.predict.tsne.v06.26.re = 
  src.sm3.update.12ss$cluster.12ss.predict.tsne.v06.26.re

# umap
umap.transfer.12ss.integrated = TransferData(anchor.12ss.integrated,
                                             t(src.12ss.integrated@reductions$umap@cell.embeddings))
src.12ss.integrated.merge_230909@reductions$umap_fta = src.12ss.integrated.merge_230909@reductions$umap
src.12ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.12ss.integrated),] = src.12ss.integrated@reductions$umap@cell.embeddings
src.12ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.sm3.update.12ss),] = as.matrix(t(umap.transfer.12ss.integrated@data))

src.12ss.integrated.merge_230909@meta.data$cluster.v06.26.re..merge =
  as.character(FNN::knn(src.12ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.12ss.integrated),],
                        src.12ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings,
                        as.character(src.12ss.integrated.merge_230909@meta.data[colnames(src.12ss.integrated),]$cluster.v06.26.re..merge), k = 10))
src.12ss.integrated.merge_230909@meta.data[colnames(src.12ss.integrated),]$cluster.v06.26.re..merge =
  src.12ss.integrated.merge_230909@meta.data[colnames(src.12ss.integrated),]$cluster.v06.26.re..merge  

src.sm3.update.12ss$cluster.12ss.predict.umap.v06.26.re = 
  src.12ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.12ss),]$cluster.v06.26.re..merge

src.12ss.integrated.merge_230909$cluster.12ss.predict.umap.v06.26.re = NA
src.12ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.12ss),]$cluster.12ss.predict.umap.v06.26.re = 
  src.sm3.update.12ss$cluster.12ss.predict.umap.v06.26.re
##################


########################
## 3.15ss 
########################
src.15ss.integrated.merge_230909 = merge(src.15ss.integrated, src.sm3.update.15ss)
src.15ss.integrated.merge_230909 = ScaleData(src.15ss.integrated.merge_230909, features = rownames(src.15ss.integrated.merge_230909))
src.15ss.integrated.merge_230909 = RunPCA(src.15ss.integrated.merge_230909, features = src.15ss.integrated.selectgene)
src.15ss.integrated.merge_230909 = RunUMAP(src.15ss.integrated.merge_230909, dims = 1:20, n.neighbors = 100)

#-------- MNN corrected in RNA level--------
src.15ss.integrated.merge_230909.selectgene = src.15ss.integrated.selectgene
MNN.res = 
  mnnCorrect(as.matrix(src.15ss.integrated.merge_230909@assays$RNA@data[src.15ss.integrated.merge_230909.selectgene, src.15ss.integrated.merge_230909$Time%in%"ss15"]),
             as.matrix(src.15ss.integrated.merge_230909@assays$RNA@data[src.15ss.integrated.merge_230909.selectgene, src.15ss.integrated.merge_230909$Time%in%"15ss"]),
             k = 5,cos.norm.out=F)
src.15ss.integrated.merge_230909@assays$mnnRNA = CreateAssayObject(data = MNN.res@assays@data$corrected)
src.15ss.integrated.merge_230909@assays$mnnRNA@key = "mnn_"
src.15ss.integrated.merge_230909 = ScaleData(src.15ss.integrated.merge_230909, 
                                             rownames(src.15ss.integrated.merge_230909@assays$mnnRNA),
                                             assay = "mnnRNA")
#----------------

anchor.15ss.integrated = 
  FindTransferAnchors(reference = src.15ss.integrated.merge_230909[,colnames(src.15ss.integrated)],
                      query = src.15ss.integrated.merge_230909[,colnames(src.sm3.update.15ss)],
                      reference.assay = "mnnRNA",
                      query.assay = "mnnRNA",scale = T,
                      features = src.15ss.integrated.selectgene)

#############
# pca
pca.transfer.15ss.integrated = TransferData(anchor.15ss.integrated,
                                            t(src.15ss.integrated@reductions$pca@cell.embeddings))
src.15ss.integrated.merge_230909@reductions$pca_fta = src.15ss.integrated.merge_230909@reductions$pca
src.15ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.15ss.integrated),] = src.15ss.integrated@reductions$pca@cell.embeddings
src.15ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.sm3.update.15ss),] = as.matrix(t(pca.transfer.15ss.integrated@data))

# tsne
tsne.transfer.15ss.integrated = TransferData(anchor.15ss.integrated,
                                             t(src.15ss.integrated@reductions$tsne@cell.embeddings))
src.15ss.integrated.merge_230909@reductions$tsne_fta = src.15ss.integrated.merge_230909@reductions$umap
src.15ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.15ss.integrated),] = src.15ss.integrated@reductions$tsne@cell.embeddings
src.15ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.sm3.update.15ss),] = as.matrix(t(tsne.transfer.15ss.integrated@data))

# umap
umap.transfer.15ss.integrated = TransferData(anchor.15ss.integrated,
                                             t(src.15ss.integrated@reductions$umap@cell.embeddings))
src.15ss.integrated.merge_230909@reductions$umap_fta = src.15ss.integrated.merge_230909@reductions$umap
src.15ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.15ss.integrated),] = src.15ss.integrated@reductions$umap@cell.embeddings
src.15ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.sm3.update.15ss),] = as.matrix(t(umap.transfer.15ss.integrated@data))
#############

#################
# predict.v06.26.re
#################
# pca
pca.transfer.15ss.integrated = TransferData(anchor.15ss.integrated,
                                            t(src.15ss.integrated@reductions$pca@cell.embeddings))
src.15ss.integrated.merge_230909@reductions$pca_fta = src.15ss.integrated.merge_230909@reductions$pca
src.15ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.15ss.integrated),] = src.15ss.integrated@reductions$pca@cell.embeddings
src.15ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.sm3.update.15ss),] = as.matrix(t(pca.transfer.15ss.integrated@data))

src.15ss.integrated.merge_230909@meta.data$cluster.predict.pca.v06.26.re =
  as.character(FNN::knn(src.15ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.15ss.integrated),],
                        src.15ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings,
                        as.character(src.15ss.integrated.merge_230909@meta.data[colnames(src.15ss.integrated),]$cluster.v06.26.re..merge), k = 10))
src.15ss.integrated.merge_230909@meta.data[colnames(src.15ss.integrated),]$cluster.predict.pca.v06.26.re =
  src.15ss.integrated.merge_230909@meta.data[colnames(src.15ss.integrated),]$cluster.predict.pca.v06.26.re  

src.sm3.update.15ss$cluster.15ss.predict.pca.v06.26.re = 
  src.15ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.15ss),]$cluster.predict.pca.v06.26.re
src.15ss.integrated.merge_230909$cluster.15ss.predict.pca.v06.26.re = NA
src.15ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.15ss),]$cluster.predict.pca.v06.26.re = 
  src.sm3.update.15ss$cluster.15ss.predict.pca.v06.26.re

# tsne
tsne.transfer.15ss.integrated = TransferData(anchor.15ss.integrated,
                                             t(src.15ss.integrated@reductions$tsne@cell.embeddings))
src.15ss.integrated.merge_230909@reductions$tsne_fta = src.15ss.integrated.merge_230909@reductions$tsne
src.15ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.15ss.integrated),] = src.15ss.integrated@reductions$tsne@cell.embeddings
src.15ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.sm3.update.15ss),] = as.matrix(t(tsne.transfer.15ss.integrated@data))

src.15ss.integrated.merge_230909@meta.data$cluster.predict.tsne.v06.26.re =
  as.character(FNN::knn(src.15ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.15ss.integrated),],
                        src.15ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings,
                        as.character(src.15ss.integrated.merge_230909@meta.data[colnames(src.15ss.integrated),]$cluster.v06.26.re..merge), k = 10))
src.15ss.integrated.merge_230909@meta.data[colnames(src.15ss.integrated),]$cluster.predict.tsne.v06.26.re =
  src.15ss.integrated.merge_230909@meta.data[colnames(src.15ss.integrated),]$cluster.predict.tsne.v06.26.re  

src.sm3.update.15ss$cluster.15ss.predict.tsne.v06.26.re = 
  src.15ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.15ss),]$cluster.predict.tsne.v06.26.re
src.15ss.integrated.merge_230909$cluster.15ss.predict.tsne.v06.26.re = NA
src.15ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.15ss),]$cluster.predict.tsne.v06.26.re = 
  src.sm3.update.15ss$cluster.15ss.predict.tsne.v06.26.re

# umap
umap.transfer.15ss.integrated = TransferData(anchor.15ss.integrated,
                                             t(src.15ss.integrated@reductions$umap@cell.embeddings))
src.15ss.integrated.merge_230909@reductions$umap_fta = src.15ss.integrated.merge_230909@reductions$umap
src.15ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.15ss.integrated),] = src.15ss.integrated@reductions$umap@cell.embeddings
src.15ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.sm3.update.15ss),] = as.matrix(t(umap.transfer.15ss.integrated@data))

src.15ss.integrated.merge_230909@meta.data$cluster.v06.26.re..merge =
  as.character(FNN::knn(src.15ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.15ss.integrated),],
                        src.15ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings,
                        as.character(src.15ss.integrated.merge_230909@meta.data[colnames(src.15ss.integrated),]$cluster.v06.26.re..merge), k = 10))
src.15ss.integrated.merge_230909@meta.data[colnames(src.15ss.integrated),]$cluster.v06.26.re..merge =
  src.15ss.integrated.merge_230909@meta.data[colnames(src.15ss.integrated),]$cluster.v06.26.re..merge  

src.sm3.update.15ss$cluster.15ss.predict.umap.v06.26.re = 
  src.15ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.15ss),]$cluster.v06.26.re..merge

src.15ss.integrated.merge_230909$cluster.15ss.predict.umap.v06.26.re = NA
src.15ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.15ss),]$cluster.15ss.predict.umap.v06.26.re = 
  src.sm3.update.15ss$cluster.15ss.predict.umap.v06.26.re
##################


########################
## 4.18ss 
########################
src.18ss.integrated.merge_230909 = merge(src.18ss.integrated, src.sm3.update.18ss)
src.18ss.integrated.merge_230909 = ScaleData(src.18ss.integrated.merge_230909, features = rownames(src.18ss.integrated.merge_230909))
src.18ss.integrated.merge_230909 = RunPCA(src.18ss.integrated.merge_230909, features = src.18ss.integrated.selectgene)
src.18ss.integrated.merge_230909 = RunUMAP(src.18ss.integrated.merge_230909, dims = 1:20, n.neighbors = 100)

#-------- MNN corrected in RNA level--------
src.18ss.integrated.merge_230909.selectgene = src.18ss.integrated.selectgene
MNN.res = 
  mnnCorrect(as.matrix(src.18ss.integrated.merge_230909@assays$RNA@data[src.18ss.integrated.merge_230909.selectgene, src.18ss.integrated.merge_230909$Time%in%"ss18"]),
             as.matrix(src.18ss.integrated.merge_230909@assays$RNA@data[src.18ss.integrated.merge_230909.selectgene, src.18ss.integrated.merge_230909$Time%in%"18ss"]),
             k = 5,cos.norm.out=F)
src.18ss.integrated.merge_230909@assays$mnnRNA = CreateAssayObject(data = MNN.res@assays@data$corrected)
src.18ss.integrated.merge_230909@assays$mnnRNA@key = "mnn_"
src.18ss.integrated.merge_230909 = ScaleData(src.18ss.integrated.merge_230909, 
                                             rownames(src.18ss.integrated.merge_230909@assays$mnnRNA),
                                             assay = "mnnRNA")
#----------------

anchor.18ss.integrated = 
  FindTransferAnchors(reference = src.18ss.integrated.merge_230909[,colnames(src.18ss.integrated)],
                      query = src.18ss.integrated.merge_230909[,colnames(src.sm3.update.18ss)],
                      reference.assay = "mnnRNA",
                      query.assay = "mnnRNA",scale = T,
                      features = src.18ss.integrated.selectgene)
#############
# pca
pca.transfer.18ss.integrated = TransferData(anchor.18ss.integrated,
                                            t(src.18ss.integrated@reductions$pca@cell.embeddings))
src.18ss.integrated.merge_230909@reductions$pca_fta = src.18ss.integrated.merge_230909@reductions$pca
src.18ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.18ss.integrated),] = src.18ss.integrated@reductions$pca@cell.embeddings
src.18ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.sm3.update.18ss),] = as.matrix(t(pca.transfer.18ss.integrated@data))

# tsne
tsne.transfer.18ss.integrated = TransferData(anchor.18ss.integrated,
                                             t(src.18ss.integrated@reductions$tsne@cell.embeddings))
src.18ss.integrated.merge_230909@reductions$tsne_fta = src.18ss.integrated.merge_230909@reductions$umap
src.18ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.18ss.integrated),] = src.18ss.integrated@reductions$tsne@cell.embeddings
src.18ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.sm3.update.18ss),] = as.matrix(t(tsne.transfer.18ss.integrated@data))

# umap
umap.transfer.18ss.integrated = TransferData(anchor.18ss.integrated,
                                             t(src.18ss.integrated@reductions$umap@cell.embeddings))
src.18ss.integrated.merge_230909@reductions$umap_fta = src.18ss.integrated.merge_230909@reductions$umap
src.18ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.18ss.integrated),] = src.18ss.integrated@reductions$umap@cell.embeddings
src.18ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.sm3.update.18ss),] = as.matrix(t(umap.transfer.18ss.integrated@data))
#############

#################
# predict.v06.26.re
#################
# pca
pca.transfer.18ss.integrated = TransferData(anchor.18ss.integrated,
                                            t(src.18ss.integrated@reductions$pca@cell.embeddings))
src.18ss.integrated.merge_230909@reductions$pca_fta = src.18ss.integrated.merge_230909@reductions$pca
src.18ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.18ss.integrated),] = src.18ss.integrated@reductions$pca@cell.embeddings
src.18ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.sm3.update.18ss),] = as.matrix(t(pca.transfer.18ss.integrated@data))

src.18ss.integrated.merge_230909@meta.data$cluster.predict.pca.v06.26.re =
  as.character(FNN::knn(src.18ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.18ss.integrated),],
                        src.18ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings,
                        as.character(src.18ss.integrated.merge_230909@meta.data[colnames(src.18ss.integrated),]$cluster.v06.26.re..merge), k = 10))
src.18ss.integrated.merge_230909@meta.data[colnames(src.18ss.integrated),]$cluster.predict.pca.v06.26.re =
  src.18ss.integrated.merge_230909@meta.data[colnames(src.18ss.integrated),]$cluster.predict.pca.v06.26.re  

src.sm3.update.18ss$cluster.18ss.predict.pca.v06.26.re = 
  src.18ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.18ss),]$cluster.predict.pca.v06.26.re
src.18ss.integrated.merge_230909$cluster.18ss.predict.pca.v06.26.re = NA
src.18ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.18ss),]$cluster.predict.pca.v06.26.re = 
  src.sm3.update.18ss$cluster.18ss.predict.pca.v06.26.re

# tsne
tsne.transfer.18ss.integrated = TransferData(anchor.18ss.integrated,
                                             t(src.18ss.integrated@reductions$tsne@cell.embeddings))
src.18ss.integrated.merge_230909@reductions$tsne_fta = src.18ss.integrated.merge_230909@reductions$tsne
src.18ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.18ss.integrated),] = src.18ss.integrated@reductions$tsne@cell.embeddings
src.18ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.sm3.update.18ss),] = as.matrix(t(tsne.transfer.18ss.integrated@data))

src.18ss.integrated.merge_230909@meta.data$cluster.predict.tsne.v06.26.re =
  as.character(FNN::knn(src.18ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.18ss.integrated),],
                        src.18ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings,
                        as.character(src.18ss.integrated.merge_230909@meta.data[colnames(src.18ss.integrated),]$cluster.v06.26.re..merge), k = 10))
src.18ss.integrated.merge_230909@meta.data[colnames(src.18ss.integrated),]$cluster.predict.tsne.v06.26.re =
  src.18ss.integrated.merge_230909@meta.data[colnames(src.18ss.integrated),]$cluster.predict.tsne.v06.26.re  

src.sm3.update.18ss$cluster.18ss.predict.tsne.v06.26.re = 
  src.18ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.18ss),]$cluster.predict.tsne.v06.26.re
src.18ss.integrated.merge_230909$cluster.18ss.predict.tsne.v06.26.re = NA
src.18ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.18ss),]$cluster.predict.tsne.v06.26.re = 
  src.sm3.update.18ss$cluster.18ss.predict.tsne.v06.26.re

# umap
umap.transfer.18ss.integrated = TransferData(anchor.18ss.integrated,
                                             t(src.18ss.integrated@reductions$umap@cell.embeddings))
src.18ss.integrated.merge_230909@reductions$umap_fta = src.18ss.integrated.merge_230909@reductions$umap
src.18ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.18ss.integrated),] = src.18ss.integrated@reductions$umap@cell.embeddings
src.18ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.sm3.update.18ss),] = as.matrix(t(umap.transfer.18ss.integrated@data))

src.18ss.integrated.merge_230909@meta.data$cluster.v06.26.re..merge =
  as.character(FNN::knn(src.18ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.18ss.integrated),],
                        src.18ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings,
                        as.character(src.18ss.integrated.merge_230909@meta.data[colnames(src.18ss.integrated),]$cluster.v06.26.re..merge), k = 10))
src.18ss.integrated.merge_230909@meta.data[colnames(src.18ss.integrated),]$cluster.v06.26.re..merge =
  src.18ss.integrated.merge_230909@meta.data[colnames(src.18ss.integrated),]$cluster.v06.26.re..merge  

src.sm3.update.18ss$cluster.18ss.predict.umap.v06.26.re = 
  src.18ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.18ss),]$cluster.v06.26.re..merge

src.18ss.integrated.merge_230909$cluster.18ss.predict.umap.v06.26.re = NA
src.18ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.18ss),]$cluster.18ss.predict.umap.v06.26.re = 
  src.sm3.update.18ss$cluster.18ss.predict.umap.v06.26.re
##################


########################
## 5.21ss 
########################
src.21ss.integrated.merge_230909 = merge(src.21ss.integrated, src.sm3.update.21ss)
src.21ss.integrated.merge_230909 = ScaleData(src.21ss.integrated.merge_230909, features = rownames(src.21ss.integrated.merge_230909))
src.21ss.integrated.merge_230909 = RunPCA(src.21ss.integrated.merge_230909, features = src.21ss.integrated.selectgene)
src.21ss.integrated.merge_230909 = RunUMAP(src.21ss.integrated.merge_230909, dims = 1:20, n.neighbors = 100)

#-------- MNN corrected in RNA level--------
src.21ss.integrated.merge_230909.selectgene = src.21ss.integrated.selectgene
MNN.res = 
  mnnCorrect(as.matrix(src.21ss.integrated.merge_230909@assays$RNA@data[src.21ss.integrated.merge_230909.selectgene, src.21ss.integrated.merge_230909$Time%in%"ss21"]),
             as.matrix(src.21ss.integrated.merge_230909@assays$RNA@data[src.21ss.integrated.merge_230909.selectgene, src.21ss.integrated.merge_230909$Time%in%"21ss"]),
             k = 5,cos.norm.out=F)
src.21ss.integrated.merge_230909@assays$mnnRNA = CreateAssayObject(data = MNN.res@assays@data$corrected)
src.21ss.integrated.merge_230909@assays$mnnRNA@key = "mnn_"
src.21ss.integrated.merge_230909 = ScaleData(src.21ss.integrated.merge_230909, 
                                             rownames(src.21ss.integrated.merge_230909@assays$mnnRNA),
                                             assay = "mnnRNA")
#----------------

anchor.21ss.integrated = 
  FindTransferAnchors(reference = src.21ss.integrated.merge_230909[,colnames(src.21ss.integrated)],
                      query = src.21ss.integrated.merge_230909[,colnames(src.sm3.update.21ss)],
                      reference.assay = "mnnRNA",
                      query.assay = "mnnRNA",scale = T,
                      features = src.21ss.integrated.selectgene)

#############
# pca
pca.transfer.21ss.integrated = TransferData(anchor.21ss.integrated,
                                            t(src.21ss.integrated@reductions$pca@cell.embeddings))
src.21ss.integrated.merge_230909@reductions$pca_fta = src.21ss.integrated.merge_230909@reductions$pca
src.21ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.21ss.integrated),] = src.21ss.integrated@reductions$pca@cell.embeddings
src.21ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.sm3.update.21ss),] = as.matrix(t(pca.transfer.21ss.integrated@data))

# tsne
tsne.transfer.21ss.integrated = TransferData(anchor.21ss.integrated,
                                             t(src.21ss.integrated@reductions$tsne@cell.embeddings))
src.21ss.integrated.merge_230909@reductions$tsne_fta = src.21ss.integrated.merge_230909@reductions$umap
src.21ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.21ss.integrated),] = src.21ss.integrated@reductions$tsne@cell.embeddings
src.21ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.sm3.update.21ss),] = as.matrix(t(tsne.transfer.21ss.integrated@data))

# umap
umap.transfer.21ss.integrated = TransferData(anchor.21ss.integrated,
                                             t(src.21ss.integrated@reductions$umap@cell.embeddings))
src.21ss.integrated.merge_230909@reductions$umap_fta = src.21ss.integrated.merge_230909@reductions$umap
src.21ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.21ss.integrated),] = src.21ss.integrated@reductions$umap@cell.embeddings
src.21ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.sm3.update.21ss),] = as.matrix(t(umap.transfer.21ss.integrated@data))
#############

#################
# predict.v06.26.re
#################
# pca
pca.transfer.21ss.integrated = TransferData(anchor.21ss.integrated,
                                            t(src.21ss.integrated@reductions$pca@cell.embeddings))
src.21ss.integrated.merge_230909@reductions$pca_fta = src.21ss.integrated.merge_230909@reductions$pca
src.21ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.21ss.integrated),] = src.21ss.integrated@reductions$pca@cell.embeddings
src.21ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.sm3.update.21ss),] = as.matrix(t(pca.transfer.21ss.integrated@data))

src.21ss.integrated.merge_230909@meta.data$cluster.predict.pca.v06.26.re =
  as.character(FNN::knn(src.21ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.21ss.integrated),],
                        src.21ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings,
                        as.character(src.21ss.integrated.merge_230909@meta.data[colnames(src.21ss.integrated),]$cluster.v06.26.re..merge), k = 10))
src.21ss.integrated.merge_230909@meta.data[colnames(src.21ss.integrated),]$cluster.predict.pca.v06.26.re =
  src.21ss.integrated.merge_230909@meta.data[colnames(src.21ss.integrated),]$cluster.predict.pca.v06.26.re  

src.sm3.update.21ss$cluster.21ss.predict.pca.v06.26.re = 
  src.21ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.21ss),]$cluster.predict.pca.v06.26.re
src.21ss.integrated.merge_230909$cluster.21ss.predict.pca.v06.26.re = NA
src.21ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.21ss),]$cluster.predict.pca.v06.26.re = 
  src.sm3.update.21ss$cluster.21ss.predict.pca.v06.26.re

# tsne
tsne.transfer.21ss.integrated = TransferData(anchor.21ss.integrated,
                                             t(src.21ss.integrated@reductions$tsne@cell.embeddings))
src.21ss.integrated.merge_230909@reductions$tsne_fta = src.21ss.integrated.merge_230909@reductions$tsne
src.21ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.21ss.integrated),] = src.21ss.integrated@reductions$tsne@cell.embeddings
src.21ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.sm3.update.21ss),] = as.matrix(t(tsne.transfer.21ss.integrated@data))

src.21ss.integrated.merge_230909@meta.data$cluster.predict.tsne.v06.26.re =
  as.character(FNN::knn(src.21ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.21ss.integrated),],
                        src.21ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings,
                        as.character(src.21ss.integrated.merge_230909@meta.data[colnames(src.21ss.integrated),]$cluster.v06.26.re..merge), k = 10))
src.21ss.integrated.merge_230909@meta.data[colnames(src.21ss.integrated),]$cluster.predict.tsne.v06.26.re =
  src.21ss.integrated.merge_230909@meta.data[colnames(src.21ss.integrated),]$cluster.predict.tsne.v06.26.re  

src.sm3.update.21ss$cluster.21ss.predict.tsne.v06.26.re = 
  src.21ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.21ss),]$cluster.predict.tsne.v06.26.re
src.21ss.integrated.merge_230909$cluster.21ss.predict.tsne.v06.26.re = NA
src.21ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.21ss),]$cluster.predict.tsne.v06.26.re = 
  src.sm3.update.21ss$cluster.21ss.predict.tsne.v06.26.re

# umap
umap.transfer.21ss.integrated = TransferData(anchor.21ss.integrated,
                                             t(src.21ss.integrated@reductions$umap@cell.embeddings))
src.21ss.integrated.merge_230909@reductions$umap_fta = src.21ss.integrated.merge_230909@reductions$umap
src.21ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.21ss.integrated),] = src.21ss.integrated@reductions$umap@cell.embeddings
src.21ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.sm3.update.21ss),] = as.matrix(t(umap.transfer.21ss.integrated@data))

src.21ss.integrated.merge_230909@meta.data$cluster.v06.26.re..merge =
  as.character(FNN::knn(src.21ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.21ss.integrated),],
                        src.21ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings,
                        as.character(src.21ss.integrated.merge_230909@meta.data[colnames(src.21ss.integrated),]$cluster.v06.26.re..merge), k = 10))
src.21ss.integrated.merge_230909@meta.data[colnames(src.21ss.integrated),]$cluster.v06.26.re..merge =
  src.21ss.integrated.merge_230909@meta.data[colnames(src.21ss.integrated),]$cluster.v06.26.re..merge  

src.sm3.update.21ss$cluster.21ss.predict.umap.v06.26.re = 
  src.21ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.21ss),]$cluster.v06.26.re..merge

src.21ss.integrated.merge_230909$cluster.21ss.predict.umap.v06.26.re = NA
src.21ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.21ss),]$cluster.21ss.predict.umap.v06.26.re = 
  src.sm3.update.21ss$cluster.21ss.predict.umap.v06.26.re
##################


########################
## 6.24ss 
########################
src.24ss.integrated.merge_230909 = merge(src.24ss.integrated, src.sm3.update.24ss)
src.24ss.integrated.merge_230909 = ScaleData(src.24ss.integrated.merge_230909, features = rownames(src.24ss.integrated.merge_230909))
src.24ss.integrated.merge_230909 = RunPCA(src.24ss.integrated.merge_230909, features = src.24ss.integrated.selectgene)
src.24ss.integrated.merge_230909 = RunUMAP(src.24ss.integrated.merge_230909, dims = 1:20, n.neighbors = 100)

#-------- MNN corrected in RNA level--------
src.24ss.integrated.merge_230909.selectgene = src.24ss.integrated.selectgene
MNN.res = 
  mnnCorrect(as.matrix(src.24ss.integrated.merge_230909@assays$RNA@data[src.24ss.integrated.merge_230909.selectgene, src.24ss.integrated.merge_230909$Time%in%"ss24"]),
             as.matrix(src.24ss.integrated.merge_230909@assays$RNA@data[src.24ss.integrated.merge_230909.selectgene, src.24ss.integrated.merge_230909$Time%in%"24ss"]),
             k = 5,cos.norm.out=F)
src.24ss.integrated.merge_230909@assays$mnnRNA = CreateAssayObject(data = MNN.res@assays@data$corrected)
src.24ss.integrated.merge_230909@assays$mnnRNA@key = "mnn_"
src.24ss.integrated.merge_230909 = ScaleData(src.24ss.integrated.merge_230909, 
                                             rownames(src.24ss.integrated.merge_230909@assays$mnnRNA),
                                             assay = "mnnRNA")
#----------------

anchor.24ss.integrated = 
  FindTransferAnchors(reference = src.24ss.integrated.merge_230909[,colnames(src.24ss.integrated)],
                      query = src.24ss.integrated.merge_230909[,colnames(src.sm3.update.24ss)],
                      reference.assay = "mnnRNA",
                      query.assay = "mnnRNA",scale = T,
                      features = src.24ss.integrated.selectgene)

#############
# pca
pca.transfer.24ss.integrated = TransferData(anchor.24ss.integrated,
                                            t(src.24ss.integrated@reductions$pca@cell.embeddings))
src.24ss.integrated.merge_230909@reductions$pca_fta = src.24ss.integrated.merge_230909@reductions$pca
src.24ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.24ss.integrated),] = src.24ss.integrated@reductions$pca@cell.embeddings
src.24ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.sm3.update.24ss),] = as.matrix(t(pca.transfer.24ss.integrated@data))

# tsne
tsne.transfer.24ss.integrated = TransferData(anchor.24ss.integrated,
                                             t(src.24ss.integrated@reductions$tsne@cell.embeddings))
src.24ss.integrated.merge_230909@reductions$tsne_fta = src.24ss.integrated.merge_230909@reductions$umap
src.24ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.24ss.integrated),] = src.24ss.integrated@reductions$tsne@cell.embeddings
src.24ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.sm3.update.24ss),] = as.matrix(t(tsne.transfer.24ss.integrated@data))

# umap
umap.transfer.24ss.integrated = TransferData(anchor.24ss.integrated,
                                             t(src.24ss.integrated@reductions$umap@cell.embeddings))
src.24ss.integrated.merge_230909@reductions$umap_fta = src.24ss.integrated.merge_230909@reductions$umap
src.24ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.24ss.integrated),] = src.24ss.integrated@reductions$umap@cell.embeddings
src.24ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.sm3.update.24ss),] = as.matrix(t(umap.transfer.24ss.integrated@data))
#############

#################
# predict.v06.26.re
#################
# pca
pca.transfer.24ss.integrated = TransferData(anchor.24ss.integrated,
                                            t(src.24ss.integrated@reductions$pca@cell.embeddings))
src.24ss.integrated.merge_230909@reductions$pca_fta = src.24ss.integrated.merge_230909@reductions$pca
src.24ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.24ss.integrated),] = src.24ss.integrated@reductions$pca@cell.embeddings
src.24ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.sm3.update.24ss),] = as.matrix(t(pca.transfer.24ss.integrated@data))

src.24ss.integrated.merge_230909@meta.data$cluster.predict.pca.v06.26.re =
  as.character(FNN::knn(src.24ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.24ss.integrated),],
                        src.24ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings,
                        as.character(src.24ss.integrated.merge_230909@meta.data[colnames(src.24ss.integrated),]$cluster.v06.26.re..merge), k = 10))
src.24ss.integrated.merge_230909@meta.data[colnames(src.24ss.integrated),]$cluster.predict.pca.v06.26.re =
  src.24ss.integrated.merge_230909@meta.data[colnames(src.24ss.integrated),]$cluster.predict.pca.v06.26.re  

src.sm3.update.24ss$cluster.24ss.predict.pca.v06.26.re = 
  src.24ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.24ss),]$cluster.predict.pca.v06.26.re
src.24ss.integrated.merge_230909$cluster.24ss.predict.pca.v06.26.re = NA
src.24ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.24ss),]$cluster.predict.pca.v06.26.re = 
  src.sm3.update.24ss$cluster.24ss.predict.pca.v06.26.re

# tsne
tsne.transfer.24ss.integrated = TransferData(anchor.24ss.integrated,
                                             t(src.24ss.integrated@reductions$tsne@cell.embeddings))
src.24ss.integrated.merge_230909@reductions$tsne_fta = src.24ss.integrated.merge_230909@reductions$tsne
src.24ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.24ss.integrated),] = src.24ss.integrated@reductions$tsne@cell.embeddings
src.24ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.sm3.update.24ss),] = as.matrix(t(tsne.transfer.24ss.integrated@data))

src.24ss.integrated.merge_230909@meta.data$cluster.predict.tsne.v06.26.re =
  as.character(FNN::knn(src.24ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.24ss.integrated),],
                        src.24ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings,
                        as.character(src.24ss.integrated.merge_230909@meta.data[colnames(src.24ss.integrated),]$cluster.v06.26.re..merge), k = 10))
src.24ss.integrated.merge_230909@meta.data[colnames(src.24ss.integrated),]$cluster.predict.tsne.v06.26.re =
  src.24ss.integrated.merge_230909@meta.data[colnames(src.24ss.integrated),]$cluster.predict.tsne.v06.26.re  

src.sm3.update.24ss$cluster.24ss.predict.tsne.v06.26.re = 
  src.24ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.24ss),]$cluster.predict.tsne.v06.26.re
src.24ss.integrated.merge_230909$cluster.24ss.predict.tsne.v06.26.re = NA
src.24ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.24ss),]$cluster.predict.tsne.v06.26.re = 
  src.sm3.update.24ss$cluster.24ss.predict.tsne.v06.26.re

# umap
umap.transfer.24ss.integrated = TransferData(anchor.24ss.integrated,
                                             t(src.24ss.integrated@reductions$umap@cell.embeddings))
src.24ss.integrated.merge_230909@reductions$umap_fta = src.24ss.integrated.merge_230909@reductions$umap
src.24ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.24ss.integrated),] = src.24ss.integrated@reductions$umap@cell.embeddings
src.24ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.sm3.update.24ss),] = as.matrix(t(umap.transfer.24ss.integrated@data))

src.24ss.integrated.merge_230909@meta.data$cluster.v06.26.re..merge =
  as.character(FNN::knn(src.24ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.24ss.integrated),],
                        src.24ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings,
                        as.character(src.24ss.integrated.merge_230909@meta.data[colnames(src.24ss.integrated),]$cluster.v06.26.re..merge), k = 10))
src.24ss.integrated.merge_230909@meta.data[colnames(src.24ss.integrated),]$cluster.v06.26.re..merge =
  src.24ss.integrated.merge_230909@meta.data[colnames(src.24ss.integrated),]$cluster.v06.26.re..merge  

src.sm3.update.24ss$cluster.24ss.predict.umap.v06.26.re = 
  src.24ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.24ss),]$cluster.v06.26.re..merge

src.24ss.integrated.merge_230909$cluster.24ss.predict.umap.v06.26.re = NA
src.24ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.24ss),]$cluster.24ss.predict.umap.v06.26.re = 
  src.sm3.update.24ss$cluster.24ss.predict.umap.v06.26.re
##################


########################
## 7.27ss 
########################
src.27ss.integrated.merge_230909 = merge(src.27ss.integrated, src.sm3.update.27ss)
src.27ss.integrated.merge_230909 = ScaleData(src.27ss.integrated.merge_230909, features = rownames(src.27ss.integrated.merge_230909))
src.27ss.integrated.merge_230909 = RunPCA(src.27ss.integrated.merge_230909, features = src.27ss.integrated.selectgene)
src.27ss.integrated.merge_230909 = RunUMAP(src.27ss.integrated.merge_230909, dims = 1:20, n.neighbors = 100)

#-------- MNN corrected in RNA level--------
src.27ss.integrated.merge_230909.selectgene = src.27ss.integrated.selectgene
MNN.res = 
  mnnCorrect(as.matrix(src.27ss.integrated.merge_230909@assays$RNA@data[src.27ss.integrated.merge_230909.selectgene, src.27ss.integrated.merge_230909$Time%in%"ss27"]),
             as.matrix(src.27ss.integrated.merge_230909@assays$RNA@data[src.27ss.integrated.merge_230909.selectgene, src.27ss.integrated.merge_230909$Time%in%"27ss"]),
             k = 5,cos.norm.out=F)
src.27ss.integrated.merge_230909@assays$mnnRNA = CreateAssayObject(data = MNN.res@assays@data$corrected)
src.27ss.integrated.merge_230909@assays$mnnRNA@key = "mnn_"
src.27ss.integrated.merge_230909 = ScaleData(src.27ss.integrated.merge_230909, 
                                             rownames(src.27ss.integrated.merge_230909@assays$mnnRNA),
                                             assay = "mnnRNA")
#----------------

anchor.27ss.integrated = 
  FindTransferAnchors(reference = src.27ss.integrated.merge_230909[,colnames(src.27ss.integrated)],
                      query = src.27ss.integrated.merge_230909[,colnames(src.sm3.update.27ss)],
                      reference.assay = "mnnRNA",
                      query.assay = "mnnRNA",scale = T,
                      features = src.27ss.integrated.selectgene)

#############
# pca
pca.transfer.27ss.integrated = TransferData(anchor.27ss.integrated,
                                            t(src.27ss.integrated@reductions$pca@cell.embeddings))
src.27ss.integrated.merge_230909@reductions$pca_fta = src.27ss.integrated.merge_230909@reductions$pca
src.27ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.27ss.integrated),] = src.27ss.integrated@reductions$pca@cell.embeddings
src.27ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.sm3.update.27ss),] = as.matrix(t(pca.transfer.27ss.integrated@data))

# tsne
tsne.transfer.27ss.integrated = TransferData(anchor.27ss.integrated,
                                             t(src.27ss.integrated@reductions$tsne@cell.embeddings))
src.27ss.integrated.merge_230909@reductions$tsne_fta = src.27ss.integrated.merge_230909@reductions$umap
src.27ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.27ss.integrated),] = src.27ss.integrated@reductions$tsne@cell.embeddings
src.27ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.sm3.update.27ss),] = as.matrix(t(tsne.transfer.27ss.integrated@data))

# umap
umap.transfer.27ss.integrated = TransferData(anchor.27ss.integrated,
                                             t(src.27ss.integrated@reductions$umap@cell.embeddings))
src.27ss.integrated.merge_230909@reductions$umap_fta = src.27ss.integrated.merge_230909@reductions$umap
src.27ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.27ss.integrated),] = src.27ss.integrated@reductions$umap@cell.embeddings
src.27ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.sm3.update.27ss),] = as.matrix(t(umap.transfer.27ss.integrated@data))
#############

#################
# predict.v06.26.re
#################
# pca
pca.transfer.27ss.integrated = TransferData(anchor.27ss.integrated,
                                            t(src.27ss.integrated@reductions$pca@cell.embeddings))
src.27ss.integrated.merge_230909@reductions$pca_fta = src.27ss.integrated.merge_230909@reductions$pca
src.27ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.27ss.integrated),] = src.27ss.integrated@reductions$pca@cell.embeddings
src.27ss.integrated.merge_230909@reductions$pca_fta@cell.embeddings[colnames(src.sm3.update.27ss),] = as.matrix(t(pca.transfer.27ss.integrated@data))

src.27ss.integrated.merge_230909@meta.data$cluster.predict.pca.v06.26.re =
  as.character(FNN::knn(src.27ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.27ss.integrated),],
                        src.27ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings,
                        as.character(src.27ss.integrated.merge_230909@meta.data[colnames(src.27ss.integrated),]$cluster.v06.26.re..merge), k = 10))
src.27ss.integrated.merge_230909@meta.data[colnames(src.27ss.integrated),]$cluster.predict.pca.v06.26.re =
  src.27ss.integrated.merge_230909@meta.data[colnames(src.27ss.integrated),]$cluster.predict.pca.v06.26.re  

src.sm3.update.27ss$cluster.27ss.predict.pca.v06.26.re = 
  src.27ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.27ss),]$cluster.predict.pca.v06.26.re
src.27ss.integrated.merge_230909$cluster.27ss.predict.pca.v06.26.re = NA
src.27ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.27ss),]$cluster.predict.pca.v06.26.re = 
  src.sm3.update.27ss$cluster.27ss.predict.pca.v06.26.re

# tsne
tsne.transfer.27ss.integrated = TransferData(anchor.27ss.integrated,
                                             t(src.27ss.integrated@reductions$tsne@cell.embeddings))
src.27ss.integrated.merge_230909@reductions$tsne_fta = src.27ss.integrated.merge_230909@reductions$tsne
src.27ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.27ss.integrated),] = src.27ss.integrated@reductions$tsne@cell.embeddings
src.27ss.integrated.merge_230909@reductions$tsne_fta@cell.embeddings[colnames(src.sm3.update.27ss),] = as.matrix(t(tsne.transfer.27ss.integrated@data))

src.27ss.integrated.merge_230909@meta.data$cluster.predict.tsne.v06.26.re =
  as.character(FNN::knn(src.27ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.27ss.integrated),],
                        src.27ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings,
                        as.character(src.27ss.integrated.merge_230909@meta.data[colnames(src.27ss.integrated),]$cluster.v06.26.re..merge), k = 10))
src.27ss.integrated.merge_230909@meta.data[colnames(src.27ss.integrated),]$cluster.predict.tsne.v06.26.re =
  src.27ss.integrated.merge_230909@meta.data[colnames(src.27ss.integrated),]$cluster.predict.tsne.v06.26.re  

src.sm3.update.27ss$cluster.27ss.predict.tsne.v06.26.re = 
  src.27ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.27ss),]$cluster.predict.tsne.v06.26.re
src.27ss.integrated.merge_230909$cluster.27ss.predict.tsne.v06.26.re = NA
src.27ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.27ss),]$cluster.predict.tsne.v06.26.re = 
  src.sm3.update.27ss$cluster.27ss.predict.tsne.v06.26.re

# umap
umap.transfer.27ss.integrated = TransferData(anchor.27ss.integrated,
                                             t(src.27ss.integrated@reductions$umap@cell.embeddings))
src.27ss.integrated.merge_230909@reductions$umap_fta = src.27ss.integrated.merge_230909@reductions$umap
src.27ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.27ss.integrated),] = src.27ss.integrated@reductions$umap@cell.embeddings
src.27ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.sm3.update.27ss),] = as.matrix(t(umap.transfer.27ss.integrated@data))

src.27ss.integrated.merge_230909@meta.data$cluster.v06.26.re..merge =
  as.character(FNN::knn(src.27ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings[colnames(src.27ss.integrated),],
                        src.27ss.integrated.merge_230909@reductions$umap_fta@cell.embeddings,
                        as.character(src.27ss.integrated.merge_230909@meta.data[colnames(src.27ss.integrated),]$cluster.v06.26.re..merge), k = 10))
src.27ss.integrated.merge_230909@meta.data[colnames(src.27ss.integrated),]$cluster.v06.26.re..merge =
  src.27ss.integrated.merge_230909@meta.data[colnames(src.27ss.integrated),]$cluster.v06.26.re..merge  

src.sm3.update.27ss$cluster.27ss.predict.umap.v06.26.re = 
  src.27ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.27ss),]$cluster.v06.26.re..merge

src.27ss.integrated.merge_230909$cluster.27ss.predict.umap.v06.26.re = NA
src.27ss.integrated.merge_230909@meta.data[colnames(src.sm3.update.27ss),]$cluster.27ss.predict.umap.v06.26.re = 
  src.sm3.update.27ss$cluster.27ss.predict.umap.v06.26.re
#===============================================================================



##################
# v06.26.re
##########
src.sm3.update.9ss$cluster.9ss.fin.pca.v06.26.re  = src.sm3.update.9ss$cluster.9ss.predict.pca.v06.26.re 
src.sm3.update.12ss$cluster.12ss.fin.pca.v06.26.re  = src.sm3.update.12ss$cluster.12ss.predict.pca.v06.26.re 
src.sm3.update.15ss$cluster.15ss.fin.pca.v06.26.re  = src.sm3.update.15ss$cluster.15ss.predict.pca.v06.26.re 
src.sm3.update.18ss$cluster.18ss.fin.pca.v06.26.re  = src.sm3.update.18ss$cluster.18ss.predict.pca.v06.26.re 
src.sm3.update.21ss$cluster.21ss.fin.pca.v06.26.re  = src.sm3.update.21ss$cluster.21ss.predict.pca.v06.26.re 
src.sm3.update.24ss$cluster.24ss.fin.pca.v06.26.re  = src.sm3.update.24ss$cluster.24ss.predict.pca.v06.26.re 
src.sm3.update.27ss$cluster.27ss.fin.pca.v06.26.re  = src.sm3.update.27ss$cluster.27ss.predict.pca.v06.26.re 

src.sm3.update.9ss$cluster.9ss.fin.tsne.v06.26.re  = src.sm3.update.9ss$cluster.9ss.predict.tsne.v06.26.re 
src.sm3.update.12ss$cluster.12ss.fin.tsne.v06.26.re  = src.sm3.update.12ss$cluster.12ss.predict.tsne.v06.26.re 
src.sm3.update.15ss$cluster.15ss.fin.tsne.v06.26.re  = src.sm3.update.15ss$cluster.15ss.predict.tsne.v06.26.re 
src.sm3.update.18ss$cluster.18ss.fin.tsne.v06.26.re  = src.sm3.update.18ss$cluster.18ss.predict.tsne.v06.26.re 
src.sm3.update.21ss$cluster.21ss.fin.tsne.v06.26.re  = src.sm3.update.21ss$cluster.21ss.predict.tsne.v06.26.re 
src.sm3.update.24ss$cluster.24ss.fin.tsne.v06.26.re  = src.sm3.update.24ss$cluster.24ss.predict.tsne.v06.26.re 
src.sm3.update.27ss$cluster.27ss.fin.tsne.v06.26.re  = src.sm3.update.27ss$cluster.27ss.predict.tsne.v06.26.re 

src.sm3.update.9ss$cluster.9ss.fin.umap.v06.26.re  = src.sm3.update.9ss$cluster.9ss.predict.umap.v06.26.re 
src.sm3.update.12ss$cluster.12ss.fin.umap.v06.26.re  = src.sm3.update.12ss$cluster.12ss.predict.umap.v06.26.re 
src.sm3.update.15ss$cluster.15ss.fin.umap.v06.26.re  = src.sm3.update.15ss$cluster.15ss.predict.umap.v06.26.re 
src.sm3.update.18ss$cluster.18ss.fin.umap.v06.26.re  = src.sm3.update.18ss$cluster.18ss.predict.umap.v06.26.re 
src.sm3.update.21ss$cluster.21ss.fin.umap.v06.26.re  = src.sm3.update.21ss$cluster.21ss.predict.umap.v06.26.re 
src.sm3.update.24ss$cluster.24ss.fin.umap.v06.26.re  = src.sm3.update.24ss$cluster.24ss.predict.umap.v06.26.re 
src.sm3.update.27ss$cluster.27ss.fin.umap.v06.26.re  = src.sm3.update.27ss$cluster.27ss.predict.umap.v06.26.re 


src.sm3.Ngn3GFP$cluster.12ss.fin.umap.v06.26.re = "na"
src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"12ss",]$cluster.12ss.fin.umap.v06.26.re = 
  src.sm3.update.12ss@meta.data[
    rownames(src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"12ss",]),]$cluster.12ss.fin.umap.v06.26.re
src.sm3.Ngn3GFP$cluster.15ss.fin.umap.v06.26.re = "na"
src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"15ss",]$cluster.15ss.fin.umap.v06.26.re = 
  src.sm3.update.15ss@meta.data[
    rownames(src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"15ss",]),]$cluster.15ss.fin.umap.v06.26.re
src.sm3.Ngn3GFP$cluster.18ss.fin.umap.v06.26.re = "na"
src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"18ss",]$cluster.18ss.fin.umap.v06.26.re = 
  src.sm3.update.18ss@meta.data[
    rownames(src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"18ss",]),]$cluster.18ss.fin.umap.v06.26.re
src.sm3.Ngn3GFP$cluster.21ss.fin.umap.v06.26.re = "na"
src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"21ss",]$cluster.21ss.fin.umap.v06.26.re = 
  src.sm3.update.21ss@meta.data[
    rownames(src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"21ss",]),]$cluster.21ss.fin.umap.v06.26.re
src.sm3.Ngn3GFP$cluster.24ss.fin.umap.v06.26.re = "na"
src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"24ss",]$cluster.24ss.fin.umap.v06.26.re = 
  src.sm3.update.24ss@meta.data[
    rownames(src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"24ss",]),]$cluster.24ss.fin.umap.v06.26.re
src.sm3.Ngn3GFP$cluster.27ss.fin.umap.v06.26.re = "na"
src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"27ss",]$cluster.27ss.fin.umap.v06.26.re = 
  src.sm3.update.27ss@meta.data[
    rownames(src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"27ss",]),]$cluster.27ss.fin.umap.v06.26.re



src.sm3.Ngn3GFP$cluster.predict.fin.umap.v06.26.re = "na"
src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"12ss",]$cluster.predict.fin.umap.v06.26.re = 
  src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"12ss",]$cluster.12ss.fin.umap.v06.26.re
src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"15ss",]$cluster.predict.fin.umap.v06.26.re = 
  src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"15ss",]$cluster.15ss.fin.umap.v06.26.re
src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"18ss",]$cluster.predict.fin.umap.v06.26.re = 
  src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"18ss",]$cluster.18ss.fin.umap.v06.26.re
src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"21ss",]$cluster.predict.fin.umap.v06.26.re = 
  src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"21ss",]$cluster.21ss.fin.umap.v06.26.re
src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"24ss",]$cluster.predict.fin.umap.v06.26.re = 
  src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"24ss",]$cluster.24ss.fin.umap.v06.26.re
src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"27ss",]$cluster.predict.fin.umap.v06.26.re = 
  src.sm3.Ngn3GFP@meta.data[src.sm3.Ngn3GFP$Time%in%"27ss",]$cluster.27ss.fin.umap.v06.26.re

src.sm3.Ngn3GFP$lineage = "Ngn3GFP"

pdf("figure.v08.07/Ngn3GFP/ngn3gfp_summary.pdf",9,5)
DimPlot(src.sm3.Ngn3GFP, group.by = "cluster.predict.fin.umap.v06.26.re",
        cols = cluster.endoderm.color.v5, pt.size = 1.2) + p_add_leg 
DimPlot(src.sm3.Ngn3GFP, group.by = "Time",
        cols = colors.time.2, pt.size = 1.2) + p_add_leg 
DimPlot(src.sm3.Ngn3GFP, group.by = "index",
        cols = colors.type, pt.size = 1.2) + p_add_leg 
dev.off()


meta_ngn3gfp = src.sm3.Ngn3GFP@meta.data
meta_ngn3gfp$Time = factor(x = meta_ngn3gfp$Time,
                           levels = c('15ss',"18ss","21ss","24ss","27ss"))
meta_ngn3gfp$cluster.predict.fin.umap.v06.26.re = 
  factor(x = meta_ngn3gfp$cluster.predict.fin.umap.v06.26.re,
         levels = c("MG.1","MG.2","MG.3.A","MG.3.M","MG.3.P",
                    "Stomach","DP","EP.1","EP.2","Small.intestine.1","Small.intestine.2",
                    "Esophagus","Pharynx.organ.2","Thyroid",
                    "AL.3-Liver","Extrahepatic biliary tract","Large.intestine.2"))


pdf("figure.v08.07/Ngn3GFP/ngn3gfp_index_cluster.pdf",6,6)
for(i in c("HZ","HM")){
  meta_ngn3gfp_re = meta_ngn3gfp[meta_ngn3gfp$index%in%i,]
  print(
    ggplot() +
      #geom_vline(xintercept = 2.5, linetype = "dashed", colour = "black")+
      geom_bar(data = meta_ngn3gfp_re, 
               mapping = aes(x = Time,
                             group = cluster.predict.fin.umap.v06.26.re, 
                             fill = cluster.predict.fin.umap.v06.26.re),
               stat = "count", position = 'fill')+
      scale_fill_manual(values = cluster.endoderm.color.v5)+
      theme_classic() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle=90, size = 13),
            axis.text.y = element_text(size = 13),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            legend.position = "none",
            legend.title = element_blank(),
            legend.text = element_text(size = 10),
            plot.title = element_blank(),
            aspect.ratio=1/0.3) +
      xlab("Time") + ylab("Cell type")
  )
}
dev.off()

#------------------------
src.sm3.Ngn3GFP.re = src.sm3.Ngn3GFP[,!src.sm3.Ngn3GFP$cluster.predict.fin.umap.v06.26.re%in%c(
  "Thyroid","AL.3-Liver","Extrahepatic biliary tract","Large.intestine.2")]

src.sm3.Ngn3GFP.re = ScaleData(src.sm3.Ngn3GFP.re,
                               split.by = "index")
src.sm3.Ngn3GFP.re = RunPCA(src.sm3.Ngn3GFP.re,
                            features = src.endoderm.mg3.ext.re.selectgene.fin)
src.sm3.Ngn3GFP.re = RunUMAP(src.sm3.Ngn3GFP.re, dims = 1:20)

pdf("figure.v08.07/Ngn3GFP/ngn3gfp.re_summary.pdf",9,5)
DimPlot(src.sm3.Ngn3GFP.re, group.by = "cluster.predict.fin.umap.v06.26.re",
        cols = cluster.endoderm.color.v5, pt.size = 1.2) + p_add_leg 
DimPlot(src.sm3.Ngn3GFP.re, group.by = "Time",
        cols = colors.time.2, pt.size = 1.2) + p_add_leg 
DimPlot(src.sm3.Ngn3GFP.re, group.by = "index",
        cols = colors.type, pt.size = 1.2) + p_add_leg 
dev.off()


meta_ngn3gfp = src.sm3.Ngn3GFP.re@meta.data
meta_ngn3gfp$Time = factor(x = meta_ngn3gfp$Time,
                           levels = c('15ss',"18ss","21ss","24ss","27ss"))
meta_ngn3gfp$cluster.predict.fin.umap.v06.26.re = 
  factor(x = meta_ngn3gfp$cluster.predict.fin.umap.v06.26.re,
         levels = c("MG.1","MG.2","MG.3.A","MG.3.M","MG.3.P",
                    "Stomach","DP","EP.1","EP.2","Small.intestine.1","Small.intestine.2",
                    "Esophagus","Pharynx.organ.2"))


pdf("figure.v08.07/Ngn3GFP/ngn3gfp.re_index_cluster.pdf",6,6)
for(i in c("HZ","HM")){
  meta_ngn3gfp_re = meta_ngn3gfp[meta_ngn3gfp$index%in%i&
                                   !meta_ngn3gfp$cluster.predict.fin.umap.v06.26.re%in%NA,]
  print(
    ggplot() +
      #geom_vline(xintercept = 2.5, linetype = "dashed", colour = "black")+
      geom_bar(data = meta_ngn3gfp_re, 
               mapping = aes(x = Time,
                             group = cluster.predict.fin.umap.v06.26.re, 
                             fill = cluster.predict.fin.umap.v06.26.re),
               stat = "count", position = 'fill')+
      scale_fill_manual(values = cluster.endoderm.color.v5)+
      theme_classic() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle=90, size = 13),
            axis.text.y = element_text(size = 13),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            legend.position = "none",
            legend.title = element_blank(),
            legend.text = element_text(size = 10),
            plot.title = element_blank(),
            aspect.ratio=1/0.3) +
      xlab("Time") + ylab("Cell type")
  )
}
dev.off()
#------------------------

# MG.3 merge Ngn3
#===============================================================================
src.endoderm.mg3.merge_v1.1_Ngn3GFP = merge(
  src.sm3.Ngn3GFP[,src.sm3.Ngn3GFP$cluster.predict.fin.umap.v06.26.re%in%c(
    "DP",'EP.1','EP.2',"MG.3.A",'MG.3.M',"MG.3.P","Stomach","Small.intestine.1",'Small.intestine.2'
  )],
  src.endoderm.mg3.ext_v1.1)

src.sm3.endoderm_Ngn3GFP =  src.sm3.Ngn3GFP[,src.sm3.Ngn3GFP$cluster.predict.fin.umap.v06.26.re%in%c(
  "DP",'EP.1','EP.2',"MG.3.A",'MG.3.M',"MG.3.P","Stomach","Small.intestine.1",'Small.intestine.2'
)]


src.endoderm.mg3.merge_v1.1_Ngn3GFP$type = 1
src.endoderm.mg3.merge_v1.1_Ngn3GFP@meta.data[src.endoderm.mg3.merge_v1.1_Ngn3GFP$lineage%in%NA,]$type = 2
src.endoderm.mg3.merge_v1.1_Ngn3GFP$batch_type = 
  paste(src.endoderm.mg3.merge_v1.1_Ngn3GFP$batch,
        src.endoderm.mg3.merge_v1.1_Ngn3GFP$type, sep = "_")

src.endoderm.mg3.merge_v1.1_Ngn3GFP = ScaleData(src.endoderm.mg3.merge_v1.1_Ngn3GFP, 
                                                features = rownames(src.endoderm.mg3.merge_v1.1_Ngn3GFP),
                                                split.by = "batch_type")

src.endoderm.mg3.merge_v1.1_Ngn3GFP.selectgene = src.endoderm.mg3.ext_v1.1.gene.fin
src.endoderm.mg3.merge_v1.1_Ngn3GFP = RunPCA(src.endoderm.mg3.merge_v1.1_Ngn3GFP,
                                             features = src.endoderm.mg3.ext_v1.1.gene.fin)
src.endoderm.mg3.merge_v1.1_Ngn3GFP = RunUMAP(src.endoderm.mg3.merge_v1.1_Ngn3GFP,
                                              dims = 1:20,n.neighbors = 50)

src.endoderm.mg3.merge_v1.1_Ngn3GFP@reductions$umap@cell.embeddings = 
  cbind(src.endoderm.mg3.merge_v1.1_Ngn3GFP@reductions$umap@cell.embeddings[,2],
        src.endoderm.mg3.merge_v1.1_Ngn3GFP@reductions$umap@cell.embeddings[,1])
colnames(src.endoderm.mg3.merge_v1.1_Ngn3GFP@reductions$umap@cell.embeddings) = c("UMAP_1","UMAP_2")

# MNN
#---------------
MNN.res = 
  mnnCorrect(as.matrix(src.endoderm.mg3.merge_v1.1_Ngn3GFP@assays$RNA@data[src.endoderm.mg3.merge_v1.1_Ngn3GFP.selectgene, 
                                                                           src.endoderm.mg3.merge_v1.1_Ngn3GFP$batch_type%in%"1_1"]),
             as.matrix(src.endoderm.mg3.merge_v1.1_Ngn3GFP@assays$RNA@data[src.endoderm.mg3.merge_v1.1_Ngn3GFP.selectgene, 
                                                                           src.endoderm.mg3.merge_v1.1_Ngn3GFP$batch_type%in%"2_1"]),
             as.matrix(src.endoderm.mg3.merge_v1.1_Ngn3GFP@assays$RNA@data[src.endoderm.mg3.merge_v1.1_Ngn3GFP.selectgene, 
                                                                           src.endoderm.mg3.merge_v1.1_Ngn3GFP$batch_type%in%"1_2"]),
             as.matrix(src.endoderm.mg3.merge_v1.1_Ngn3GFP@assays$RNA@data[src.endoderm.mg3.merge_v1.1_Ngn3GFP.selectgene, 
                                                                           src.endoderm.mg3.merge_v1.1_Ngn3GFP$batch_type%in%"2_2"]),
             as.matrix(src.endoderm.mg3.merge_v1.1_Ngn3GFP@assays$RNA@data[src.endoderm.mg3.merge_v1.1_Ngn3GFP.selectgene, 
                                                                           src.endoderm.mg3.merge_v1.1_Ngn3GFP$batch_type%in%"3_1"]),
             k = 5,cos.norm.out=F)
src.endoderm.mg3.merge_v1.1_Ngn3GFP@assays$mnnRNA = CreateAssayObject(data = MNN.res@assays@data$corrected)
src.endoderm.mg3.merge_v1.1_Ngn3GFP@assays$mnnRNA@key = "mnn_"
src.endoderm.mg3.merge_v1.1_Ngn3GFP = ScaleData(src.endoderm.mg3.merge_v1.1_Ngn3GFP, 
                                                rownames(src.endoderm.mg3.merge_v1.1_Ngn3GFP@assays$mnnRNA),
                                                assay = "mnnRNA")

src.endoderm.mg3.merge_v1.1_Ngn3GFP@meta.data = src.endoderm.mg3.merge_v1.1_Ngn3GFP@meta.data[!src.endoderm.mg3.merge_v1.1_Ngn3GFP$Time%in%NA,]

anchor.mg3.merge_Ngn3GFP = 
  FindTransferAnchors(reference = src.endoderm.mg3.merge_v1.1_Ngn3GFP[,src.endoderm.mg3.merge_v1.1_Ngn3GFP$type%in%2],
                      query = src.endoderm.mg3.merge_v1.1_Ngn3GFP[,src.endoderm.mg3.merge_v1.1_Ngn3GFP$type%in%1],
                      reference.assay = "mnnRNA",
                      query.assay = "mnnRNA",scale = T,
                      features = src.endoderm.mg3.ext_v1.1.gene.fin)
umap.transfer.mg3.merge_Ngn3GFP = TransferData(anchor.mg3.merge_Ngn3GFP,
                                               t(src.endoderm.mg3.ext_v1.1@reductions$umap_mnn@cell.embeddings))
# 
# src.endoderm.mg3.ext_v1.1@reductions$umap_mnn@cell.embeddings = 
#   src.endoderm.mg3.ext_v1.1@reductions$umap_mnn@cell.embeddings[colnames(src.endoderm.mg3.ext_v1.1),]

src.endoderm.mg3.merge_v1.1_Ngn3GFP@reductions$umap_fta = src.endoderm.mg3.merge_v1.1_Ngn3GFP@reductions$umap
src.endoderm.mg3.merge_v1.1_Ngn3GFP@reductions$umap_fta@cell.embeddings[colnames(src.endoderm.mg3.ext_v1.1),] = 
  src.endoderm.mg3.ext_v1.1@reductions$umap_mnn@cell.embeddings
src.endoderm.mg3.merge_v1.1_Ngn3GFP@reductions$umap_fta@cell.embeddings[colnames(src.sm3.endoderm_Ngn3GFP),] = 
  as.matrix(t(umap.transfer.mg3.merge_Ngn3GFP@data))
src.endoderm.mg3.merge_v1.1_Ngn3GFP@reductions$umap_fta@key = "Coord_"
colnames(src.endoderm.mg3.merge_v1.1_Ngn3GFP@reductions$umap_fta@cell.embeddings) = c("Coord_1","Coord_2")
DimPlot(src.endoderm.mg3.merge_v1.1_Ngn3GFP, reduction = "umap_fta",
        group.by = "lineage", cols = color.lineage)

a = FNN::knn(
  src.endoderm.mg3.merge_v1.1_Ngn3GFP@reductions$umap_fta@cell.embeddings[
    rownames(src.endoderm.mg3.merge_v1.1_Ngn3GFP@meta.data[src.endoderm.mg3.merge_v1.1_Ngn3GFP$type%in%2,]),],
  src.endoderm.mg3.merge_v1.1_Ngn3GFP@reductions$umap_fta@cell.embeddings[
    rownames(src.endoderm.mg3.merge_v1.1_Ngn3GFP@meta.data[src.endoderm.mg3.merge_v1.1_Ngn3GFP$type%in%1,]),],
  src.endoderm.mg3.merge_v1.1_Ngn3GFP@meta.data[src.endoderm.mg3.merge_v1.1_Ngn3GFP$type%in%2,]$cluster.v06.26.re, k = 10)
src.endoderm.mg3.merge_v1.1_Ngn3GFP$cluster.predict.umap_int.ext.v1.1 = src.endoderm.mg3.merge_v1.1_Ngn3GFP$cluster.v06.26.re
src.endoderm.mg3.merge_v1.1_Ngn3GFP@meta.data[src.endoderm.mg3.merge_v1.1_Ngn3GFP$type%in%1,]$cluster.predict.umap_int.ext.v1.1 = as.character(a)
DimPlot(src.endoderm.mg3.merge_v1.1_Ngn3GFP, group.by = "cluster.predict.umap_int.ext.v1.1")


p = ggplot()+
  geom_point(data = 
               cbind(src.endoderm.mg3.ext_v1.1@meta.data,
                     src.endoderm.mg3.ext_v1.1@reductions$umap_mnn@cell.embeddings),
             mapping = aes(x=UMAP_1, y=UMAP_2), colour="#eeeeee") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_classic2() + p_add
p0 = p + geom_point(data = 
                      cbind(src.endoderm.mg3.merge_v1.1_Ngn3GFP@meta.data[colnames(src.sm3.endoderm_Ngn3GFP),],
                            src.endoderm.mg3.merge_v1.1_Ngn3GFP@reductions$umap_fta@cell.embeddings[colnames(src.sm3.endoderm_Ngn3GFP),]),
                    mapping = aes(x=Coord_1, y=Coord_2, color = index)) +
  scale_color_manual(values = colors.type)
p1 = p + geom_point(data = 
                      cbind(src.endoderm.mg3.merge_v1.1_Ngn3GFP@meta.data[colnames(src.sm3.endoderm_Ngn3GFP),],
                            src.endoderm.mg3.merge_v1.1_Ngn3GFP@reductions$umap_fta@cell.embeddings[colnames(src.sm3.endoderm_Ngn3GFP),]),
                    mapping = aes(x=Coord_1, y=Coord_2, color = Time)) +
  scale_color_manual(values = colors.time.2)
p2 = p + geom_point(data = 
                      cbind(src.endoderm.mg3.merge_v1.1_Ngn3GFP@meta.data[colnames(src.sm3.endoderm_Ngn3GFP),],
                            src.endoderm.mg3.merge_v1.1_Ngn3GFP@reductions$umap_fta@cell.embeddings[colnames(src.sm3.endoderm_Ngn3GFP),]),
                    mapping = aes(x=Coord_1, y=Coord_2, color = cluster.predict.umap_int.ext.v1.1)) +
  scale_color_manual(values = cluster.endoderm.color.v5)

mat = src.endoderm.mg3.merge_v1.1_Ngn3GFP@meta.data[colnames(src.sm3.endoderm_Ngn3GFP),]
mat$Time = factor(mat$Time,
                  levels = c("27ss","24ss","21ss","18ss"))
time_list = table(src.endoderm.mg3.merge_v1.1_Ngn3GFP@meta.data[colnames(src.sm3.endoderm_Ngn3GFP),"Time"])


pdf("figure.v08.07/ngn3gfp_fta.pdf",9,7)
p0;p1;p2;p3;
dev.off()










