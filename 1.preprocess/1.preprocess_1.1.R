#===============================================================================
#>> 1.Pre-processing 
#>
#>  1.1 Pre-processing of 10x-v3
#===============================================================================

setwd("~/Bioinformatic/project_20220204/endoderm 10x/cluster/project_20220327")
source("~/bin/MyFunction.R")
source("~/bin/MyFunction.Seurat3.R")
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(rgl)
library(scran)
setwd("D:/project180311/endoderm 10x/cluster/")
gi=read.csv("~/genome/mm10/mm10.gene.inf.merge.10X.v1.csv",stringsAsFactors = F)
gi$gene=gi$Symbol2
rownames(gi)=gi$Symbol2

marker.sym <- c(
  #-------------------------------------------
  "Tbx6", # neuralmesoderm
  "Nkx1-2", # neural tube 
  "T", # notochord
  "Noto", # notochord
  "Lmx1a", # neural tube
  "Tdh", # yolk sac endoderm
  "Lgals2", # yolk sac endoderm
  "Cubn", # yolk sac endoderm
  "Pla2g12b", # yolk sac endoderm
  "Apoc2", # yolk sac endoderm
  "Foxf1", # mesoderm
  "Snai2", # mesoderm
  "Hand1", # mesoderm
  "Twist2", # mesoderm
  "Msgn1", #presomitic mesoderm
  "Wnt8a", #presomitic mesoderm
  "Sox7", # vessel
  "Gata1", # blood
  "Epcam", # endoderm
  "Otx2", # foregut
  "Prrx2", # foregut
  "Pyy", # foregut
  "Trh", # foregut
  "Tbx1", # foregut
  "Nepn", # midgut
  "Cyp26a1", # hindgut
  "Hoxb6", # hindgut
  "Cdx1", # hindgut
  "Six1", # FA
  "Pax1", # FA
  "Pax9", # FA
  "Hhex", # FL ventral
  "Foxa1","Foxa2","Foxa3","Tbx3","Cdx2",
  "Afp","Pdx1","Mnx1","Neurod1","Gcg","Ins1",
  "Tfap2a", # neural crest progenitor
  "Foxi2", # Cranial Ectoderm?
  "Tfap2c", # neural ectoderm
  "Dlx5", # surface ectoderm
  "Meox1", # Somite
  "Cxcr4", # Endoderm
  "Gjb3", # migrated neural crest
  "Trap1a", # yolk sac boundary
  "Bhlha9", # Limb ?
  "Dppa3", #germline
  # "Myl3", #Parietal Endoderm myocyte
  "Six1","Has2","Sox2","Hhex","Syt6","Nepn","Hoxb9","T","Wnt5b","Hoxa9",
  "Nkx2-1","Mgll","Nkx2-3","Pitx2",
  "Pax9","Fgf8","Jag1","Fgf1",
  "Igf1","Chrd","Hhip",
  "Osr1","Hoxb1",
  "Afp","Apom", #Liver
  "Ghrl","Pdx1","Ly6h", #Pancreas
  "Pitx2","Fabp1","Upk1a","Adcy8","Tspan8", #Intestine
  "Mnx1","Lgals2",
  "Cdx2","Hoxa9","Hoxb9","Hoxd10","Fzd10",
  "Pax9", "Has2", "Hhex", "Nepn", "Cdx2", 
  "Wnt5b", "Smoc2", "Foxg1", "Cdc42ep3", "Prrx2", 
  "Asb4", "Nkx2-1", "Nkx2-3", "Fgf8", "Tbx1", 
  "Cldn11", "Ly6h", "Dab2", "Afp", 
  "Hoxb1", "Upk1a", "Mnx1", "2610528A11Rik"
  #-------------------------------------------
  )

#==========================================================
#>>>> (1.1.1) Load Rawdata
#==========================================================
src.12ss <- Read10X(data.dir = "../rawdata/ss12/")
src.15ss <- Read10X(data.dir = "../rawdata/ss15/")
src.18ss <- Read10X(data.dir = "../rawdata/ss18/")
src.21ss <- Read10X(data.dir = "../rawdata/ss21/")
src.24ss <- Read10X(data.dir = "../rawdata/ss24/")
src.27ss <- Read10X(data.dir = "../rawdata/ss27/")
src.EpcamE8.5 <- Read10X(data.dir = "../rawdata/EpcamE8.5/")
src.Epcam12SS <- Read10X(data.dir = "../rawdata/Epcamss12/")
src.Epcam15SS <- Read10X(data.dir = "../rawdata/Epcamss15/")
src.Epcam18SS <- Read10X(data.dir = "../rawdata/Epcamss18/")
src.Epcam21SS <- Read10X(data.dir = "../rawdata/Epcamss21/")
src.Epcam24SS <- Read10X(data.dir = "../rawdata/Epcamss24/")
src.EpcamE9.5 <- Read10X(data.dir = "../rawdata/EpcamE9.5/")

src.12ss=CreateSeuratObject(src.12ss)
src.15ss=CreateSeuratObject(src.15ss)
src.18ss=CreateSeuratObject(src.18ss)
src.21ss=CreateSeuratObject(src.21ss)
src.24ss=CreateSeuratObject(src.24ss)
src.27ss=CreateSeuratObject(src.27ss)
src.EpcamE8.5=CreateSeuratObject(src.EpcamE8.5)
src.Epcam12SS <- CreateSeuratObject(src.Epcam12SS)
src.Epcam15SS <- CreateSeuratObject(src.Epcam15SS)
src.Epcam18SS <- CreateSeuratObject(src.Epcam18SS)
src.Epcam21SS <- CreateSeuratObject(src.Epcam21SS)
src.Epcam24SS <- CreateSeuratObject(src.Epcam24SS)
src.EpcamE9.5=CreateSeuratObject(src.EpcamE9.5)

src.12ss$Time="ss12"
src.15ss$Time="ss15"
src.18ss$Time="ss18"
src.21ss$Time="ss21"
src.24ss$Time="ss24"
src.27ss$Time="ss27"
src.EpcamE8.5$Time="ss9"
src.Epcam12SS$Time="ss12"
src.Epcam15SS$Time="ss15"
src.Epcam18SS$Time="ss18"
src.Epcam21SS$Time="ss21"
src.Epcam24SS$Time="ss24"
src.EpcamE9.5$Time="ss27"

src.12ss$Time2="/"
src.15ss$Time2="/"
src.18ss$Time2="/"
src.21ss$Time2="/"
src.24ss$Time2="/"
src.27ss$Time2="/"
src.EpcamE8.5$Time2="E8.5"
src.Epcam12SS$Time2="/"
src.Epcam15SS$Time2="/"
src.Epcam18SS$Time2="/"
src.Epcam21SS$Time2="/"
src.Epcam24SS$Time2="/"
src.EpcamE9.5$Time2="E9.5"

src.12ss$injecttime=NA
src.15ss$injecttime=NA
src.18ss$injecttime=NA
src.21ss$injecttime=NA
src.24ss$injecttime=NA
src.27ss$injecttime=NA
src.EpcamE8.5$injecttime=NA
src.Epcam12SS$injecttime=NA
src.Epcam15SS$injecttime=NA
src.Epcam18SS$injecttime=NA
src.Epcam21SS$injecttime=NA
src.Epcam24SS$injecttime=NA
src.EpcamE9.5$injecttime=NA

src.12ss$Sort="UnderSomite"
src.15ss$Sort="UnderSomite"
src.18ss$Sort="UnderSomite"
src.21ss$Sort="UnderSomite"
src.24ss$Sort="UnderSomite"
src.27ss$Sort="UnderSomite"
src.EpcamE8.5$Sort="Whole"
src.Epcam12SS$Sort="Whole"
src.Epcam15SS$Sort="Whole"
src.Epcam18SS$Sort="Whole"
src.Epcam21SS$Sort="Whole"
src.Epcam24SS$Sort="Whole"
src.EpcamE9.5$Sort="Whole"

src.12ss$Batch=1
src.15ss$Batch=1
src.18ss$Batch=1
src.21ss$Batch=1
src.24ss$Batch=1
src.27ss$Batch=1
src.EpcamE8.5$Batch=2
src.Epcam12SS$Batch=2
src.Epcam15SS$Batch=2
src.Epcam18SS$Batch=2
src.Epcam21SS$Batch=2
src.Epcam24SS$Batch=2
src.EpcamE9.5$Batch=2

src.12ss$Cre=NA
src.15ss$Cre=NA
src.18ss$Cre=NA
src.21ss$Cre=NA
src.24ss$Cre=NA
src.27ss$Cre=NA
src.EpcamE8.5$Cre=NA
src.Epcam12SS$Cre=NA
src.Epcam15SS$Cre=NA
src.Epcam18SS$Cre=NA
src.Epcam21SS$Cre=NA
src.Epcam24SS$Cre=NA
src.EpcamE9.5$Cre=NA
#==========================================================


#==========================================================
#>>>> (1.1.2) QC
#==========================================================
src.list.name.10x = c("src.12ss",
                      "src.15ss",
                      "src.18ss",
                      "src.21ss",
                      "src.24ss",
                      "src.27ss",
                      "src.EpcamE8.5",
                      "src.Epcam12SS",
                      "src.Epcam15SS",
                      "src.Epcam18SS",
                      "src.Epcam21SS",
                      "src.Epcam24SS",
                      "src.EpcamE9.5")


for (seurat in c(src.list.name.10x)) {
  assign(seurat,
         NormalizeData(get(seurat),
                       normalization.method = "LogNormalize",
                       scale.factor = 10000))
  mt.gene=grep("^mt-",rownames(get(seurat)),value = T)
  mt.gene.exp=colMeans(get(seurat)@assays$RNA@data[mt.gene,])
  assign(seurat,
         AddMetaData(get(seurat),mt.gene.exp,"Mito.gene"))
  assign(seurat,
         FindVariableFeatures(get(seurat), selection.method = "vst", nfeatures = 2000))
  VariableFeaturePlot(get(seurat))
  assign(seurat,
         ScaleData(get(seurat),features = VariableFeatures(get(seurat))))
  assign(seurat,
         MyDoubletFinder(get(seurat),round(ncol(get(seurat)) * 0.1,0)))
  assign(seurat,
         RunPCA(get(seurat), features = VariableFeatures(get(seurat))))
  assign(seurat,
         RunTSNE(get(seurat), dims = 1:30,check_duplicates = F,
                 perplexity=round((30+ncol(get(seurat))/100)),
                 max_iter = round((ncol(get(seurat))/12))))
  assign(seurat,
         FindNeighbors(get(seurat), dims = 1:30))
  assign(seurat,
         FindClusters(get(seurat), resolution = 1.2))
  assign(seurat,
         FindClusters(get(seurat), resolution = 3))
  assign(seurat,
         FindClusters(get(seurat), resolution = 1.5))
  pdf(paste(seurat,".cluster.pdf",sep = ""),8,7)
  if (seurat%in%c("src.12ss",
                  "src.15ss",
                  "src.18ss",
                  "src.21ss",
                  "src.24ss",
                  "src.27ss")) {
    assign(seurat,
           SetIdent(get(seurat), value = get(seurat)$RNA_snn_res.1.2))
  }
  if (seurat%in%c("src.EpcamE8.5",
                  "src.Epcam12SS",
                  "src.Epcam15SS",
                  "src.Epcam18SS",
                  "src.Epcam21SS",
                  "src.Epcam24SS",
                  "src.EpcamE9.5",
                  "src.E9.0EpcamE10.5")) {
    assign(seurat,
           SetIdent(get(seurat), value = get(seurat)$RNA_snn_res.1.5))
  }
  print(TSNEPlot(get(seurat),label = T)+
          theme(axis.ticks =  element_blank(),axis.text = element_blank()) +
          guides(colour = guide_legend(title = "",override.aes = list(size = 8)))+
          theme_bw() +
          scale_color_manual(values = c(brewer.pal(9,"Set1"),
                                        brewer.pal(12,"Set3"),
                                        brewer.pal(8,"Set2"),
                                        brewer.pal(9,"Set1"),
                                        brewer.pal(12,"Set3"),
                                        brewer.pal(8,"Set2"))) +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                aspect.ratio=1))
  assign(seurat,
         SetIdent(get(seurat), value = get(seurat)$pANNPredictions))
  print(TSNEPlot(get(seurat))+
          theme(axis.ticks =  element_blank(),axis.text = element_blank()) +
          guides(colour = guide_legend(title = "",override.aes = list(size = 8)))+
          theme_bw() +
          scale_color_manual(values = c(brewer.pal(9,"Set1"))) +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                aspect.ratio=1))
  dev.off()
  pdf(paste(seurat,".marker.pdf",sep = ""),8,7)
  for (gene in c(marker.sym,"nFeature_RNA")) {
    print(FeaturePlot(get(seurat),gene,cols = c("#C7E8CC","#FF0000"))+
            theme(axis.ticks =  element_blank(),axis.text = element_blank()) +
            guides(colour = guide_colorbar(title = "",barwidth = 3,barheight = 15))+
            theme_bw() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  aspect.ratio=1))
  }
  dev.off()
}


src.12ss@meta.data$cluster=as.character(src.12ss@meta.data$RNA_snn_res.1.2)
src.12ss@meta.data[src.12ss@meta.data$cluster%in%c(23),]$cluster="Erythrocyte"
src.12ss@meta.data[src.12ss@meta.data$cluster%in%c(16),]$cluster="Endothelial"
src.12ss@meta.data[src.12ss@meta.data$cluster%in%c(4,17,11,13,19,6),]$cluster="Mesoderm"
src.12ss@meta.data[src.12ss@meta.data$cluster%in%c(21),]$cluster="Ectoderm"
src.12ss@meta.data[src.12ss@meta.data$cluster%in%c(1,10,12,8,20,2,5,9,7),]$cluster="Endoderm"
src.12ss@meta.data[src.12ss@meta.data$cluster%in%c(15,3),]$cluster="Yolk sac"
src.12ss@meta.data[src.12ss@meta.data$cluster%in%c(0,14,18),]$cluster="Low Quality"
src.12ss@meta.data[src.12ss@meta.data$cluster%in%c(24),]$cluster="Primordial Germ Cell"
src.12ss@meta.data[src.12ss@meta.data$cluster%in%c(22)|
                     src.12ss@meta.data$pANNPredictions=="Doublet",]$cluster="Doublet"

src.15ss@meta.data$cluster=as.character(src.15ss@meta.data$RNA_snn_res.1.2)
src.15ss@meta.data[src.15ss@meta.data$cluster%in%c(25),]$cluster="Erythrocyte"
src.15ss@meta.data[src.15ss@meta.data$cluster%in%c(19),]$cluster="Endothelial"
src.15ss@meta.data[src.15ss@meta.data$cluster%in%c(2,18,17,14,15,21),]$cluster="Mesoderm"
src.15ss@meta.data[src.15ss@meta.data$cluster%in%c(22,23,26),]$cluster="Ectoderm"
src.15ss@meta.data[src.15ss@meta.data$cluster%in%c(3,6,20,1,7,4,12,9,10,11,0),]$cluster="Endoderm"
src.15ss@meta.data[src.15ss@meta.data$cluster%in%c(16,13),]$cluster="Yolk sac"
src.15ss@meta.data[src.15ss@meta.data$cluster%in%c(5,8),]$cluster="Low Quality"
src.15ss@meta.data[src.15ss@meta.data$cluster%in%c(24),]$cluster="Primordial Germ Cell"
src.15ss@meta.data[src.15ss@meta.data$pANNPredictions=="Doublet",]$cluster="Doublet"

pdf("src.18ss.cluster.pdf",8,7)
src.18ss=SetIdent(src.18ss, value = src.18ss$RNA_snn_res.1.5)
print(TSNEPlot(src.18ss,label = T)+
        theme(axis.ticks =  element_blank(),axis.text = element_blank()) +
        guides(colour = guide_legend(title = "",override.aes = list(size = 8)))+
        theme_bw() +
        scale_color_manual(values = c(brewer.pal(9,"Set1"),
                                      brewer.pal(12,"Set3"),
                                      brewer.pal(8,"Set2"),
                                      brewer.pal(9,"Set1"),
                                      brewer.pal(12,"Set3"),
                                      brewer.pal(8,"Set2"))) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              aspect.ratio=1))
src.18ss=SetIdent(src.18ss, value = src.18ss$pANNPredictions)
print(TSNEPlot(src.18ss)+
        theme(axis.ticks =  element_blank(),axis.text = element_blank()) +
        guides(colour = guide_legend(title = "",override.aes = list(size = 8)))+
        theme_bw() +
        scale_color_manual(values = c(brewer.pal(9,"Set1"))) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              aspect.ratio=1))
dev.off()
src.18ss@meta.data$cluster=as.character(src.18ss@meta.data$RNA_snn_res.1.5)
src.18ss@meta.data[src.18ss@meta.data$cluster%in%c(26),]$cluster="Erythrocyte"
src.18ss@meta.data[src.18ss@meta.data$cluster%in%c(20),]$cluster="Endothelial"
src.18ss@meta.data[src.18ss@meta.data$cluster%in%c(8,10,16,17,18,22),]$cluster="Mesoderm"
src.18ss@meta.data[src.18ss@meta.data$cluster%in%c(19),]$cluster="Ectoderm"
src.18ss@meta.data[src.18ss@meta.data$cluster%in%c(12,9,4,11,1,21,27,6,0,15,7,24,25,5,2),]$cluster="Endoderm"
src.18ss@meta.data[src.18ss@meta.data$cluster%in%c(14),]$cluster="Yolk sac"
src.18ss@meta.data[src.18ss@meta.data$cluster%in%c(3,28),]$cluster="Low Quality"
src.18ss@meta.data[src.18ss@meta.data$cluster%in%c(23),]$cluster="Primordial Germ Cell"
src.18ss@meta.data[(src.18ss@meta.data$cluster%in%c(13)|
                      src.18ss@meta.data$pANNPredictions=="Doublet")&
                     !src.18ss@meta.data$RNA_snn_res.1.5%in%c(2),]$cluster="Doublet"

src.21ss@meta.data$cluster=as.character(src.21ss@meta.data$RNA_snn_res.1.2)
src.21ss@meta.data[src.21ss@meta.data$cluster%in%c(25),]$cluster="Erythrocyte"
src.21ss@meta.data[src.21ss@meta.data$cluster%in%c(19),]$cluster="Endothelial"
src.21ss@meta.data[src.21ss@meta.data$cluster%in%c(21,5,11,15,24,9),]$cluster="Mesoderm"
src.21ss@meta.data[src.21ss@meta.data$cluster%in%c(27,16,22),]$cluster="Ectoderm"
src.21ss@meta.data[src.21ss@meta.data$cluster%in%c(6,8,20,4,12,14,0,7,18,26,1,10,3,13),]$cluster="Endoderm"
src.21ss@meta.data[src.21ss@meta.data$cluster%in%c(17),]$cluster="Yolk sac"
src.21ss@meta.data[src.21ss@meta.data$cluster%in%c(2),]$cluster="Low Quality"
src.21ss@meta.data[src.21ss@meta.data$cluster%in%c(23),]$cluster="Primordial Germ Cell"
src.21ss@meta.data[(src.21ss@meta.data$pANNPredictions=="Doublet")&
                     !src.21ss@meta.data$RNA_snn_res.1.2%in%c(20),]$cluster="Doublet"

src.24ss$Dbcluster=dbscan::dbscan(src.24ss@reductions$tsne@cell.embeddings,eps = 0.5)$cluster
pdf("src.24ss.cluster.pdf",8,7)
src.24ss=SetIdent(src.24ss, value = src.24ss$RNA_snn_res.1.5)
print(TSNEPlot(src.24ss,label = T)+
        theme(axis.ticks =  element_blank(),axis.text = element_blank()) +
        guides(colour = guide_legend(title = "",override.aes = list(size = 8)))+
        theme_bw() +
        scale_color_manual(values = c(brewer.pal(9,"Set1"),
                                      brewer.pal(12,"Set3"),
                                      brewer.pal(8,"Set2"),
                                      brewer.pal(9,"Set1"),
                                      brewer.pal(12,"Set3"),
                                      brewer.pal(8,"Set2"))) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              aspect.ratio=1))
src.24ss=SetIdent(src.24ss, value = src.24ss$pANNPredictions)
print(TSNEPlot(src.24ss)+
        theme(axis.ticks =  element_blank(),axis.text = element_blank()) +
        guides(colour = guide_legend(title = "",override.aes = list(size = 8)))+
        theme_bw() +
        scale_color_manual(values = c(brewer.pal(9,"Set1"))) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              aspect.ratio=1))
src.24ss=SetIdent(src.24ss, value = src.24ss$Dbcluster)
print(TSNEPlot(src.24ss,label = T)+
        theme(axis.ticks =  element_blank(),axis.text = element_blank()) +
        guides(colour = guide_legend(title = "",override.aes = list(size = 8)))+
        theme_bw() +
        scale_color_manual(values = rep(c(brewer.pal(9,"Set1"),
                                          brewer.pal(12,"Set3"),
                                          brewer.pal(8,"Set2")),
                                        10)) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              aspect.ratio=1))
dev.off()
src.24ss@meta.data$cluster=as.character(src.24ss@meta.data$RNA_snn_res.1.5)
src.24ss@meta.data[src.24ss@meta.data$cluster%in%c(31),]$cluster="Erythrocyte"
src.24ss@meta.data[src.24ss@meta.data$cluster%in%c(21),]$cluster="Endothelial"
src.24ss@meta.data[src.24ss@meta.data$cluster%in%c(29,12,4,11),]$cluster="Mesoderm"
src.24ss@meta.data[src.24ss@meta.data$cluster%in%c(16,30,26),]$cluster="Ectoderm"
src.24ss@meta.data[src.24ss@meta.data$cluster%in%c(18,19,1,20,8,3,17,22,2,9,13,14,5,6,7,25,28),]$cluster="Endoderm"
src.24ss@meta.data[src.24ss@meta.data$cluster%in%c(23),]$cluster="Yolk sac"
src.24ss@meta.data[src.24ss@meta.data$Dbcluster%in%c(76,58,54)|
                     src.24ss@meta.data$cluster%in%c(0,10,15,27),]$cluster="Low Quality"
src.24ss@meta.data[src.24ss@meta.data$cluster%in%c(24),]$cluster="Primordial Germ Cell"
src.24ss@meta.data[(src.24ss@meta.data$cluster%in%c(32)|
                      src.24ss@meta.data$pANNPredictions=="Doublet")&
                     !src.24ss@meta.data$RNA_snn_res.1.5%in%c(25,2),]$cluster="Doublet"

src.27ss=FindClusters(src.27ss, resolution = 2.4)
src.27ss$Dbcluster=dbscan::dbscan(src.27ss@reductions$tsne@cell.embeddings,eps = 0.5)$cluster
pdf("src.27ss.cluster.pdf",8,7)
src.27ss=SetIdent(src.27ss, value = src.27ss$RNA_snn_res.2.4)
print(TSNEPlot(src.27ss,label = T)+
        theme(axis.ticks =  element_blank(),axis.text = element_blank()) +
        guides(colour = guide_legend(title = "",override.aes = list(size = 8)))+
        theme_bw() +
        scale_color_manual(values = c(brewer.pal(9,"Set1"),
                                      brewer.pal(12,"Set3"),
                                      brewer.pal(8,"Set2"),
                                      brewer.pal(9,"Set1"),
                                      brewer.pal(12,"Set3"),
                                      brewer.pal(8,"Set2"))) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              aspect.ratio=1))
src.27ss=SetIdent(src.27ss, value = src.27ss$Dbcluster)
print(TSNEPlot(src.27ss,label = T)+
        theme(axis.ticks =  element_blank(),axis.text = element_blank()) +
        guides(colour = guide_legend(title = "",override.aes = list(size = 8)))+
        theme_bw() +
        scale_color_manual(values = c(brewer.pal(9,"Set1"),
                                      brewer.pal(12,"Set3"),
                                      brewer.pal(8,"Set2"),
                                      brewer.pal(9,"Set1"),
                                      brewer.pal(12,"Set3"),
                                      brewer.pal(8,"Set2"),
                                      brewer.pal(9,"Set1"),
                                      brewer.pal(12,"Set3"),
                                      brewer.pal(8,"Set2"),
                                      brewer.pal(9,"Set1"),
                                      brewer.pal(12,"Set3"),
                                      brewer.pal(8,"Set2"))) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              aspect.ratio=1))
src.27ss=SetIdent(src.27ss, value = src.27ss$pANNPredictions)
print(TSNEPlot(src.27ss)+
        theme(axis.ticks =  element_blank(),axis.text = element_blank()) +
        guides(colour = guide_legend(title = "",override.aes = list(size = 8)))+
        theme_bw() +
        scale_color_manual(values = c(brewer.pal(9,"Set1"))) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              aspect.ratio=1))
dev.off()
src.27ss@meta.data$cluster=as.character(src.27ss@meta.data$RNA_snn_res.2.4)
src.27ss@meta.data[src.27ss@meta.data$cluster%in%c(35),]$cluster="Erythrocyte"
src.27ss@meta.data[src.27ss@meta.data$cluster%in%c(28),]$cluster="Endothelial"
src.27ss@meta.data[src.27ss@meta.data$cluster%in%c(31,17,14,8,3),]$cluster="Mesoderm"
src.27ss@meta.data[src.27ss@meta.data$cluster%in%c(30,11,26,25),]$cluster="Ectoderm"
src.27ss@meta.data[src.27ss@meta.data$cluster%in%c(21,5,9,16,20,10,13,4,27,12,19,15,2,6,29,23,34,1,18,24),]$cluster="Endoderm"
src.27ss@meta.data[src.27ss@meta.data$Dbcluster%in%c(47),]$cluster="Yolk sac"
src.27ss@meta.data[src.27ss@meta.data$cluster%in%c(0,7,22,32)|
                     src.27ss@meta.data$Dbcluster%in%c(89),]$cluster="Low Quality"
src.27ss@meta.data[src.27ss@meta.data$cluster%in%c(33),]$cluster="Primordial Germ Cell"
src.27ss@meta.data[src.27ss@meta.data$pANNPredictions=="Doublet"&
                     !src.27ss@meta.data$RNA_snn_res.2.4%in%c(34)&
                     !src.27ss@meta.data$RNA_snn_res.2.4%in%c(24),]$cluster="Doublet"

pdf("src.EpcamE8.5.cluster.pdf",8,7)
src.EpcamE8.5=SetIdent(src.EpcamE8.5, value = src.EpcamE8.5$RNA_snn_res.1.5)
print(TSNEPlot(src.EpcamE8.5,label = T)+
        theme(axis.ticks =  element_blank(),axis.text = element_blank()) +
        guides(colour = guide_legend(title = "",override.aes = list(size = 8)))+
        theme_bw() +
        scale_color_manual(values = c(brewer.pal(9,"Set1"),
                                      brewer.pal(12,"Set3"),
                                      brewer.pal(8,"Set2"),
                                      brewer.pal(9,"Set1"),
                                      brewer.pal(12,"Set3"),
                                      brewer.pal(8,"Set2"))) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              aspect.ratio=1))
src.EpcamE8.5=SetIdent(src.EpcamE8.5, value = src.EpcamE8.5$RNA_snn_res.3)
print(TSNEPlot(src.EpcamE8.5,label = T)+
        theme(axis.ticks =  element_blank(),axis.text = element_blank()) +
        guides(colour = guide_legend(title = "",override.aes = list(size = 8)))+
        theme_bw() +
        scale_color_manual(values = c(brewer.pal(9,"Set1"),
                                      brewer.pal(12,"Set3"),
                                      brewer.pal(8,"Set2"),
                                      brewer.pal(9,"Set1"),
                                      brewer.pal(12,"Set3"),
                                      brewer.pal(8,"Set2"))) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              aspect.ratio=1))
src.EpcamE8.5=SetIdent(src.EpcamE8.5, value = src.EpcamE8.5$pANNPredictions)
print(TSNEPlot(src.EpcamE8.5)+
        theme(axis.ticks =  element_blank(),axis.text = element_blank()) +
        guides(colour = guide_legend(title = "",override.aes = list(size = 8)))+
        theme_bw() +
        scale_color_manual(values = c(brewer.pal(9,"Set1"))) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              aspect.ratio=1))
dev.off()
src.EpcamE8.5@meta.data$cluster=as.character(src.EpcamE8.5@meta.data$RNA_snn_res.1.5)
src.EpcamE8.5@meta.data[src.EpcamE8.5@meta.data$cluster%in%c(11,16,13,8)|
                          src.EpcamE8.5$RNA_snn_res.3%in%c(41),]$cluster="Low Quality"
src.EpcamE8.5@meta.data[src.EpcamE8.5@meta.data$cluster%in%c(19,26,23,22),]$cluster="Yolk sac"
src.EpcamE8.5@meta.data[src.EpcamE8.5@meta.data$cluster%in%c(25),]$cluster="Ectoderm"
src.EpcamE8.5@meta.data[src.EpcamE8.5@meta.data$cluster%in%c(18),]$cluster="Primordial Germ Cell"
src.EpcamE8.5@meta.data[src.EpcamE8.5@meta.data$cluster%in%c(0,1,2,3,4,5,6,7,9,10,12,14,15,17,20,21,24,27),]$cluster="Endoderm"

src.Epcam12SS@meta.data$cluster=as.character(src.Epcam12SS@meta.data$RNA_snn_res.1.5)
src.Epcam12SS@meta.data[src.Epcam12SS@meta.data$cluster%in%c(15,10,24,19,4,6,14,23,3,2,0,11,18,20,17,1,21,8),]$cluster="Endoderm"
src.Epcam12SS@meta.data[src.Epcam12SS@meta.data$cluster%in%c(12,27),]$cluster="Ectoderm"
src.Epcam12SS@meta.data[src.Epcam12SS@meta.data$cluster%in%c(16,29,22,26)|
                          src.Epcam12SS@assays$RNA@data["Trap1a",]>0,]$cluster="Yolk sac"
src.Epcam12SS@meta.data[src.Epcam12SS@meta.data$cluster%in%c(28),]$cluster="Mesoderm"
src.Epcam12SS@meta.data[src.Epcam12SS@meta.data$cluster%in%c(30),]$cluster="Notochord"
src.Epcam12SS@meta.data[src.Epcam12SS@meta.data$cluster%in%c(9,7,13,5),]$cluster="Low Quality"
src.Epcam12SS@meta.data[src.Epcam12SS@meta.data$cluster%in%c(25),]$cluster="Primordial Germ Cell"

src.Epcam15SS@meta.data$cluster=as.character(src.Epcam15SS@meta.data$RNA_snn_res.1.5)
src.Epcam15SS@meta.data[src.Epcam15SS@meta.data$cluster%in%c(1,3,10,12,16),]$cluster="Low Quality"
src.Epcam15SS@meta.data[src.Epcam15SS@meta.data$cluster%in%c(30,20,15),]$cluster="Yolk sac"
src.Epcam15SS@meta.data[src.Epcam15SS@meta.data$cluster%in%c(4,5,21,14,17,26,22,6,0,24,7,19,29,18,2,8),]$cluster="Endoderm"
src.Epcam15SS@meta.data[src.Epcam15SS@meta.data$cluster%in%c(25),]$cluster="Primordial Germ Cell"
src.Epcam15SS@meta.data[src.Epcam15SS@meta.data$cluster%in%c(23,28,9,11,31,27,13),]$cluster="Ectoderm"

src.Epcam18SS@meta.data$cluster=as.character(src.Epcam18SS@meta.data$RNA_snn_res.1.5)
src.Epcam18SS@meta.data[src.Epcam18SS@meta.data$cluster%in%c(15,21,0,9,3,29),]$cluster="Low Quality"
src.Epcam18SS@meta.data[src.Epcam18SS@meta.data$cluster%in%c(26,22),]$cluster="Yolk sac"
src.Epcam18SS@meta.data[src.Epcam18SS@meta.data$cluster%in%c(1,17,28,4,7,27),]$cluster="Ectoderm"
src.Epcam18SS@meta.data[src.Epcam18SS@meta.data$cluster%in%c(30),]$cluster="Primordial Germ Cell"
src.Epcam18SS@meta.data[src.Epcam18SS@meta.data$cluster%in%c(2,24,8,25,20,13,10,11,6,23,5,12,18,19,16,31,14),]$cluster="Endoderm"
src.Epcam18SS@meta.data[src.Epcam18SS@meta.data$cluster%in%c(32),]$cluster="Doublet"

pdf("src.Epcam21SS.cluster.pdf",8,7)
src.Epcam21SS=SetIdent(src.Epcam21SS, value = src.Epcam21SS$RNA_snn_res.3)
print(TSNEPlot(src.Epcam21SS,label = T)+
        theme(axis.ticks =  element_blank(),axis.text = element_blank()) +
        guides(colour = guide_legend(title = "",override.aes = list(size = 8)))+
        theme_bw() +
        scale_color_manual(values = c(brewer.pal(9,"Set1"),
                                      brewer.pal(12,"Set3"),
                                      brewer.pal(8,"Set2"),
                                      brewer.pal(9,"Set1"),
                                      brewer.pal(12,"Set3"),
                                      brewer.pal(8,"Set2"))) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              aspect.ratio=1))
src.Epcam21SS=SetIdent(src.Epcam21SS, value = src.Epcam21SS$pANNPredictions)
print(TSNEPlot(src.Epcam21SS)+
        theme(axis.ticks =  element_blank(),axis.text = element_blank()) +
        guides(colour = guide_legend(title = "",override.aes = list(size = 8)))+
        theme_bw() +
        scale_color_manual(values = c(brewer.pal(9,"Set1"))) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              aspect.ratio=1))
dev.off()
src.Epcam21SS@meta.data$cluster=as.character(src.Epcam21SS@meta.data$RNA_snn_res.3)
src.Epcam21SS@meta.data[src.Epcam21SS@meta.data$cluster%in%c(0,1,12,3,28,11,35,37,2,10),]$cluster="Low Quality"
src.Epcam21SS@meta.data[src.Epcam21SS@meta.data$cluster%in%c(32,38),]$cluster="Yolk sac"
src.Epcam21SS@meta.data[src.Epcam21SS@meta.data$cluster%in%c(4,6,20,19,36,27,7,5,33,30),]$cluster="Ectoderm"
src.Epcam21SS@meta.data[src.Epcam21SS@meta.data$cluster%in%c(17,21,18,25,22,24,31,34,15,14,23,13,9,8,16,26,39,29),]$cluster="Endoderm"

src.Epcam24SS@meta.data$cluster=as.character(src.Epcam24SS@meta.data$RNA_snn_res.1.5)
src.Epcam24SS@meta.data[src.Epcam24SS@meta.data$cluster%in%c(2,4,22,1,23),]$cluster="Low Quality"
src.Epcam24SS@meta.data[src.Epcam24SS@meta.data$cluster%in%c()|
                          src.Epcam24SS@assays$RNA@data["Trap1a",]>0,]$cluster="Yolk sac"
src.Epcam24SS@meta.data[src.Epcam24SS@meta.data$cluster%in%c(16,5,14,12,8,29,0,9,25),]$cluster="Ectoderm"
src.Epcam24SS@meta.data[src.Epcam24SS@meta.data$cluster%in%c(30),]$cluster="Primordial Germ Cell"
src.Epcam24SS@meta.data[src.Epcam24SS@meta.data$cluster%in%c(24,17,11,10,21,20,3,28,15,19,27,13,6,26,7,31,18),]$cluster="Endoderm"

src.EpcamE9.5@meta.data$cluster=as.character(src.EpcamE9.5@meta.data$RNA_snn_res.1.5)
src.EpcamE9.5@meta.data[src.EpcamE9.5@meta.data$cluster%in%c(13,0),]$cluster="Low Quality"
src.EpcamE9.5@meta.data[src.EpcamE9.5@meta.data$cluster%in%c(20,2,12,6,11,21,25),]$cluster="Ectoderm"
src.EpcamE9.5@meta.data[src.EpcamE9.5@meta.data$cluster%in%c(26),]$cluster="Primordial Germ Cell"
src.EpcamE9.5@meta.data[src.EpcamE9.5@meta.data$cluster%in%c(3),]$cluster="Mesoderm"
src.EpcamE9.5@meta.data[src.EpcamE9.5@meta.data$cluster%in%c(1,4,5,7,8,9,10,14,15,17,18,19,22,23,24,16),]$cluster="Endoderm"


#-->> Add-Lowa quality standard: less than 2k5 detected genes
for (seurat in c(src.list.name.10x)) {
  seurat.temp = get(seurat)
  seurat.temp@assays[seurat.temp$nFeature_RNA<2500,]$cluster="Low Quality"
  assign(seurat, seurat.temp)
}

#-->> merge-to-plot
qc.color2=c(
  c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072"),
  c("#80B1D3","#FDB462","#B3DE69","#FCCDE5"),
  c("#BC80BD","#CCEBC5","#FFED6F",'#9370db'),
  c("#00ad1c","#D9D9D9","#D9D9Da","black"))
names(qc.color2)=c(
  "Endoderm","Mesoderm","Ectoderm","Erythrocyte",
  "Endothelial","Yolk sac","Primordial Germ Cell","Erythroid cell",
  "Neuromesoderm","Paraxial Mesoderm","Notochord","Somite",
  "Doublet","Low Quality","Low quality","Contamination")

src.merge.rawdata.no9SS = merge(
  src.12ss,  # src.EpcamE8.5, 
  c(src.Epcam12SS,
    src.15ss, src.Epcam15SS,
    src.18ss, src.Epcam18SS,
    src.21ss, src.Epcam21SS,
    src.24ss, src.Epcam24SS,
    src.27ss, src.Epcam27SS),
  add.cell.ids = c("ss12","ss12","ss15","ss15",
                   "ss18","ss18","ss21","ss21",
                   "ss24","ss24","ss27","ss27"))

for (seurat in c("src.merge.rawdata.no9SS")) {
  assign(seurat,
         NormalizeData(get(seurat),
                       normalization.method = "LogNormalize",
                       scale.factor = 10000))
  mt.gene=grep("^mt-",rownames(get(seurat)),value = T)
  mt.gene.exp=colMeans(get(seurat)@assays$RNA@data[mt.gene,])
  assign(seurat,
         AddMetaData(get(seurat),mt.gene.exp,"Mito.gene"))
  assign(seurat,
         FindVariableFeatures(get(seurat), selection.method = "vst", nfeatures = 2000))
  VariableFeaturePlot(get(seurat))
  assign(seurat,
         ScaleData(get(seurat),features = VariableFeatures(get(seurat))))
  assign(seurat,
         MyDoubletFinder(get(seurat),round(ncol(get(seurat)) * 0.1,0)))
  assign(seurat,
         RunPCA(get(seurat), features = VariableFeatures(get(seurat))))
  assign(seurat,
         RunTSNE(get(seurat), dims = 1:30,check_duplicates = F,
                 perplexity=round((30+ncol(get(seurat))/100)),
                 max_iter = round((ncol(get(seurat))/12))))
}

pdf(paste("qc.cluster.src.merge.rawdata.no9SS.pdf",sep = ""),8,7)
for(i.src in c("src.EpcamE8.5", "src.merge.rawdata.no9SS")){
  print(DimPlot(get(i.src), group.by = "cluster", pt.size = 1.25)+
          scale_color_manual(values = qc.color2) +
          theme(axis.ticks =  element_blank(),axis.text = element_blank()) +
          # guides(colour = guide_colorbar(title = "",barwidth = 3,barheight = 15))+
          theme_bw() +
          theme(legend.position = "none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                aspect.ratio=1)) 
}
dev.off()



#==========================================================


#==========================================================
#>>>> (1.1.3) endoderm
#==========================================================
src.list.name.10x = c(
  "src.12ss", "src.15ss", "src.18ss",
  "src.21ss", "src.24ss", "src.27ss",
  "src.EpcamE8.5", "src.Epcam12SS",
  "src.Epcam15SS", "src.Epcam18SS",
  "src.Epcam21SS", "src.Epcam24SS",
  "src.EpcamE9.5")

for (seurat in c(src.list.name.10x)) {
  assign(paste(seurat,".endoderm",sep = ""),
         get(seurat)[,get(seurat)$cluster=="Endoderm"])
}

# rm(src.12ss,src.12ss.marker.all,
#    src.15ss,src.15ss.marker.all,
#    src.18ss,src.18ss.marker.all,
#    src.21ss,src.21ss.marker.all,
#    src.24ss,src.24ss.marker.all,
#    src.27ss,src.27ss.marker.all,
#    src.EpcamE8.5,src.EpcamE8.5.marker.all,
#    src.Epcam12SS,src.Epcam12SS.marker.all,
#    src.Epcam15SS,src.Epcam15SS.marker.all,
#    src.Epcam18SS,src.Epcam18SS.marker.all,
#    src.Epcam21SS,src.Epcam21SS.marker.all,
#    src.Epcam24SS,src.Epcam24SS.marker.all,
#    src.EpcamE9.5,src.EpcamE9.5.marker.all,olddata)

src.list.name.10x.endoderm = c(
  "src.12ss.endoderm",
  "src.15ss.endoderm",
  "src.18ss.endoderm",
  "src.21ss.endoderm",
  "src.24ss.endoderm",
  "src.27ss.endoderm",
  "src.EpcamE8.5.endoderm",
  "src.Epcam12SS.endoderm",
  "src.Epcam15SS.endoderm",
  "src.Epcam18SS.endoderm",
  "src.Epcam21SS.endoderm",
  "src.Epcam24SS.endoderm",
  "src.EpcamE9.5.endoderm")

for (seurat in c(src.list.name.10x.endoderm)) {
  assign(seurat,
         FindVariableFeatures(get(seurat), selection.method = "vst", nfeatures = 3000))
  VariableFeaturePlot(get(seurat))
  assign(seurat,
         ScaleData(get(seurat),features = VariableFeatures(get(seurat))))
  assign(paste(seurat,".selectgene",sep = ""),
         Myfilter(as.matrix(get(seurat)@assays$RNA@data),
                  gene = get(seurat)@assays$RNA@var.features,
                  pearson.threshold = 0.2,partner.threshold = 5,
                  bottom.dispersion.interval = 0.1))
  assign(paste(seurat,".row.tree",sep = ""),
         MyHeatmap(as.matrix(get(seurat)@assays$RNA@data[get(paste(seurat,".selectgene",sep = "")),]),
                   type = "row.relat",
                   hc.c.data.type = "row.relat",
                   hc.r.data.type = "row.relat",
                   c.cov.method = "s",
                   r.cov.method = "s",
                   c.hc.method = "ward.D",
                   r.hc.method = "ward.D2",
                   return.tree = "row",
                   graph = F))
  assign(paste(seurat,".row.tree",sep = ""),
         as.dendrogram(get(paste(seurat,".row.tree",sep = ""))))
  plot(get(paste(seurat,".row.tree",sep = "")))
}


src.12ss.endoderm.selectgene=setdiff(src.12ss.endoderm.selectgene,
                                     c(labels(src.12ss.endoderm.row.tree[[2]][[1]]),
                                       labels(src.12ss.endoderm.row.tree[[2]][[2]][[2]][[1]][[1]]),
                                       labels(src.12ss.endoderm.row.tree[[2]][[2]][[2]][[2]][[2]][[1]])))
src.15ss.endoderm.selectgene=setdiff(src.15ss.endoderm.selectgene,
                                     c(labels(src.15ss.endoderm.row.tree[[2]][[1]]),
                                       labels(src.15ss.endoderm.row.tree[[2]][[2]][[2]][[2]][[2]])))
src.18ss.endoderm.selectgene=setdiff(src.18ss.endoderm.selectgene,
                                     c(labels(src.18ss.endoderm.row.tree[[2]][[2]][[1]]),
                                       labels(src.18ss.endoderm.row.tree[[2]][[2]][[2]][[2]][[2]])))
src.21ss.endoderm.selectgene=setdiff(src.21ss.endoderm.selectgene,
                                     c(labels(src.21ss.endoderm.row.tree[[2]][[2]][[1]]),
                                       labels(src.21ss.endoderm.row.tree[[2]][[2]][[2]][[2]][[1]])))
src.24ss.endoderm.selectgene=setdiff(src.24ss.endoderm.selectgene,
                                     c(labels(src.24ss.endoderm.row.tree[[2]][[2]][[1]]),
                                       labels(src.24ss.endoderm.row.tree[[2]][[2]][[2]][[1]])))
src.27ss.endoderm.selectgene=setdiff(src.27ss.endoderm.selectgene,
                                     c(labels(src.27ss.endoderm.row.tree[[2]][[1]][[1]]),
                                       labels(src.27ss.endoderm.row.tree[[2]][[1]][[2]][[1]])))

src.EpcamE8.5.endoderm.selectgene=setdiff(src.EpcamE8.5.endoderm.selectgene,
                                          c(labels(src.EpcamE8.5.endoderm.row.tree[[2]][[1]]),
                                            labels(src.EpcamE8.5.endoderm.row.tree[[2]][[2]][[1]][[2]])))
src.Epcam12SS.endoderm.selectgene=setdiff(src.Epcam12SS.endoderm.selectgene,
                                          c(labels(src.Epcam12SS.endoderm.row.tree[[2]][[1]]),
                                            labels(src.Epcam12SS.endoderm.row.tree[[2]][[2]][[1]][[2]]),
                                            labels(src.Epcam12SS.endoderm.row.tree[[2]][[2]][[2]][[2]][[2]][[1]])))
src.Epcam15SS.endoderm.selectgene=setdiff(src.Epcam15SS.endoderm.selectgene,
                                          c(labels(src.Epcam15SS.endoderm.row.tree[[2]][[1]]),
                                            labels(src.Epcam15SS.endoderm.row.tree[[2]][[2]][[2]][[1]])))
src.Epcam18SS.endoderm.selectgene=setdiff(src.Epcam18SS.endoderm.selectgene,
                                          c(labels(src.Epcam18SS.endoderm.row.tree[[2]][[1]]),
                                            labels(src.Epcam18SS.endoderm.row.tree[[2]][[2]][[2]][[1]])))
src.Epcam21SS.endoderm.selectgene=setdiff(src.Epcam21SS.endoderm.selectgene,
                                          c(labels(src.Epcam21SS.endoderm.row.tree[[2]][[1]]),
                                            labels(src.Epcam21SS.endoderm.row.tree[[2]][[2]][[2]][[1]])))
src.Epcam24SS.endoderm.selectgene=setdiff(src.Epcam24SS.endoderm.selectgene,
                                          c(labels(src.Epcam24SS.endoderm.row.tree[[2]][[2]][[1]]),
                                            labels(src.Epcam24SS.endoderm.row.tree[[2]][[2]][[2]][[1]])))
src.EpcamE9.5.endoderm.selectgene=setdiff(src.EpcamE9.5.endoderm.selectgene,
                                          c(labels(src.EpcamE9.5.endoderm.row.tree[[1]][[1]]),
                                            labels(src.EpcamE9.5.endoderm.row.tree[[1]][[2]])))


for (seurat in c(src.list.name.10x.endoderm)) {
  assign(seurat,
         RunPCA(get(seurat), features = get(paste(seurat,".selectgene",sep = ""))))
  assign(seurat,
         RunTSNE(get(seurat), dims = 1:30,
                 perplexity= round((30+ncol(get(seurat))/100)),
                 max_iter = round((ncol(get(seurat))/12)),
                 dim.embed = 2,seed.use = 10))
  assign(seurat,
         RunUMAP(get(seurat), dims = 1:30,n.neighbors = 200))
  assign(seurat,
         FindNeighbors(get(seurat), dims = 1:30))
  assign(seurat,
         FindClusters(get(seurat), resolution = 1.5))
  assign(seurat,
         FindClusters(get(seurat), resolution = 3.5))
}

src.24ss.endoderm=RunTSNE(src.24ss.endoderm, dims = 1:20,
                          perplexity= round((30+ncol(src.24ss.endoderm)/100)),
                          dim.embed = 2,seed.use = 10)
src.27ss.endoderm=RunTSNE(src.27ss.endoderm, dims = 1:20,
                          perplexity= round((30+ncol(src.27ss.endoderm)/100)),
                          dim.embed = 2,seed.use = 10)
src.Epcam21SS.endoderm=RunTSNE(src.Epcam21SS.endoderm, dims = 1:30,
                               perplexity= round((30+ncol(src.Epcam21SS.endoderm)/100)))
src.Epcam24SS.endoderm=RunTSNE(src.Epcam24SS.endoderm, dims = 1:30,
                               perplexity= round((30+ncol(src.Epcam24SS.endoderm)/100)))
src.EpcamE9.5.endoderm=RunTSNE(src.EpcamE9.5.endoderm, dims = 1:30,
                               perplexity= round((30+ncol(src.EpcamE9.5.endoderm)/100)))

for(i.just.rotated in c("i.just.rotated")){
  
  src.12ss.endoderm@reductions$tsne@cell.embeddings=
    cbind(-src.12ss.endoderm@reductions$tsne@cell.embeddings[,2],
          src.12ss.endoderm@reductions$tsne@cell.embeddings[,1])
  colnames(src.12ss.endoderm@reductions$tsne@cell.embeddings)=c("tSNE_1","tSNE_2")
  src.12ss.endoderm@reductions$umap@cell.embeddings=
    cbind(src.12ss.endoderm@reductions$umap@cell.embeddings[,1],
          -src.12ss.endoderm@reductions$umap@cell.embeddings[,2])
  colnames(src.12ss.endoderm@reductions$umap@cell.embeddings)=c("UMAP_1","UMAP_2")
  
  src.15ss.endoderm@reductions$tsne@cell.embeddings=
    cbind(src.15ss.endoderm@reductions$tsne@cell.embeddings[,2],
          src.15ss.endoderm@reductions$tsne@cell.embeddings[,1])
  colnames(src.15ss.endoderm@reductions$tsne@cell.embeddings)=c("tSNE_1","tSNE_2")
  src.15ss.endoderm@reductions$umap@cell.embeddings=
    cbind(-src.15ss.endoderm@reductions$umap@cell.embeddings[,1],
          -src.15ss.endoderm@reductions$umap@cell.embeddings[,2])
  colnames(src.15ss.endoderm@reductions$umap@cell.embeddings)=c("UMAP_1","UMAP_2")
  
  src.18ss.endoderm@reductions$tsne@cell.embeddings=
    cbind(src.18ss.endoderm@reductions$tsne@cell.embeddings[,1],
          src.18ss.endoderm@reductions$tsne@cell.embeddings[,2])
  colnames(src.18ss.endoderm@reductions$tsne@cell.embeddings)=c("tSNE_1","tSNE_2")
  src.18ss.endoderm@reductions$umap@cell.embeddings=
    cbind(src.18ss.endoderm@reductions$umap@cell.embeddings[,1],
          src.18ss.endoderm@reductions$umap@cell.embeddings[,2])
  colnames(src.18ss.endoderm@reductions$umap@cell.embeddings)=c("UMAP_1","UMAP_2")
  
  src.21ss.endoderm@reductions$tsne@cell.embeddings=
    cbind(src.21ss.endoderm@reductions$tsne@cell.embeddings[,2],
          -src.21ss.endoderm@reductions$tsne@cell.embeddings[,1])
  colnames(src.21ss.endoderm@reductions$tsne@cell.embeddings)=c("tSNE_1","tSNE_2")
  src.21ss.endoderm@reductions$umap@cell.embeddings=
    cbind(src.21ss.endoderm@reductions$umap@cell.embeddings[,1],
          -src.21ss.endoderm@reductions$umap@cell.embeddings[,2])
  colnames(src.21ss.endoderm@reductions$umap@cell.embeddings)=c("UMAP_1","UMAP_2")
  
  src.24ss.endoderm@reductions$tsne@cell.embeddings=
    cbind(src.24ss.endoderm@reductions$tsne@cell.embeddings[,2],
          -src.24ss.endoderm@reductions$tsne@cell.embeddings[,1])
  colnames(src.24ss.endoderm@reductions$tsne@cell.embeddings)=c("tSNE_1","tSNE_2")
  src.24ss.endoderm@reductions$umap@cell.embeddings=
    cbind(-src.24ss.endoderm@reductions$umap@cell.embeddings[,1],
          src.24ss.endoderm@reductions$umap@cell.embeddings[,2])
  colnames(src.24ss.endoderm@reductions$umap@cell.embeddings)=c("UMAP_1","UMAP_2")
  
  src.27ss.endoderm@reductions$tsne@cell.embeddings=
    cbind(src.27ss.endoderm@reductions$tsne@cell.embeddings[,1],
          src.27ss.endoderm@reductions$tsne@cell.embeddings[,2])
  colnames(src.27ss.endoderm@reductions$tsne@cell.embeddings)=c("tSNE_1","tSNE_2")
  src.27ss.endoderm@reductions$umap@cell.embeddings=
    cbind(src.27ss.endoderm@reductions$umap@cell.embeddings[,1],
          -src.27ss.endoderm@reductions$umap@cell.embeddings[,2])
  colnames(src.27ss.endoderm@reductions$umap@cell.embeddings)=c("UMAP_1","UMAP_2")
  
  src.EpcamE8.5.endoderm@reductions$umap@cell.embeddings=
    cbind(src.EpcamE8.5.endoderm@reductions$umap@cell.embeddings[,2],
          src.EpcamE8.5.endoderm@reductions$umap@cell.embeddings[,1])
  colnames(src.EpcamE8.5.endoderm@reductions$umap@cell.embeddings)=c("UMAP_1","UMAP_2")
  
  src.Epcam12SS.endoderm@reductions$tsne@cell.embeddings=
    cbind(-src.Epcam12SS.endoderm@reductions$tsne@cell.embeddings[,2],
          src.Epcam12SS.endoderm@reductions$tsne@cell.embeddings[,1])
  colnames(src.Epcam12SS.endoderm@reductions$tsne@cell.embeddings)=c("tSNE_1","tSNE_2")
  src.Epcam12SS.endoderm@reductions$umap@cell.embeddings=
    cbind(-src.Epcam12SS.endoderm@reductions$umap@cell.embeddings[,2],
          src.Epcam12SS.endoderm@reductions$umap@cell.embeddings[,1])
  colnames(src.Epcam12SS.endoderm@reductions$umap@cell.embeddings)=c("UMAP_1","UMAP_2")
  
  src.Epcam15SS.endoderm@reductions$tsne@cell.embeddings=
    cbind(src.Epcam15SS.endoderm@reductions$tsne@cell.embeddings[,1],
          src.Epcam15SS.endoderm@reductions$tsne@cell.embeddings[,2])
  colnames(src.Epcam15SS.endoderm@reductions$tsne@cell.embeddings)=c("tSNE_1","tSNE_2")
  src.Epcam15SS.endoderm@reductions$umap@cell.embeddings=
    cbind(-src.Epcam15SS.endoderm@reductions$umap@cell.embeddings[,2],
          -src.Epcam15SS.endoderm@reductions$umap@cell.embeddings[,1])
  colnames(src.Epcam15SS.endoderm@reductions$umap@cell.embeddings)=c("UMAP_1","UMAP_2")
  
  src.Epcam18SS.endoderm@reductions$tsne@cell.embeddings=
    cbind(-src.Epcam18SS.endoderm@reductions$tsne@cell.embeddings[,2],
          src.Epcam18SS.endoderm@reductions$tsne@cell.embeddings[,1])
  colnames(src.Epcam18SS.endoderm@reductions$tsne@cell.embeddings)=c("tSNE_1","tSNE_2")
  src.Epcam18SS.endoderm@reductions$umap@cell.embeddings=
    cbind(src.Epcam18SS.endoderm@reductions$umap@cell.embeddings[,1],
          -src.Epcam18SS.endoderm@reductions$umap@cell.embeddings[,2])
  colnames(src.Epcam18SS.endoderm@reductions$umap@cell.embeddings)=c("UMAP_1","UMAP_2")
  
  src.Epcam21SS.endoderm@reductions$tsne@cell.embeddings=
    cbind(-src.Epcam21SS.endoderm@reductions$tsne@cell.embeddings[,1],
          -src.Epcam21SS.endoderm@reductions$tsne@cell.embeddings[,2])
  colnames(src.Epcam21SS.endoderm@reductions$tsne@cell.embeddings)=c("tSNE_1","tSNE_2")
  src.Epcam21SS.endoderm@reductions$umap@cell.embeddings=
    cbind(-src.Epcam21SS.endoderm@reductions$umap@cell.embeddings[,1],
          -src.Epcam21SS.endoderm@reductions$umap@cell.embeddings[,2])
  colnames(src.Epcam21SS.endoderm@reductions$umap@cell.embeddings)=c("UMAP_1","UMAP_2")
  
  src.Epcam24SS.endoderm@reductions$tsne@cell.embeddings=
    cbind(-src.Epcam24SS.endoderm@reductions$tsne@cell.embeddings[,2],
          -src.Epcam24SS.endoderm@reductions$tsne@cell.embeddings[,1])
  colnames(src.Epcam24SS.endoderm@reductions$tsne@cell.embeddings)=c("tSNE_1","tSNE_2")
  src.Epcam24SS.endoderm@reductions$umap@cell.embeddings=
    cbind(-src.Epcam24SS.endoderm@reductions$umap@cell.embeddings[,2],
          -src.Epcam24SS.endoderm@reductions$umap@cell.embeddings[,1])
  colnames(src.Epcam24SS.endoderm@reductions$umap@cell.embeddings)=c("UMAP_1","UMAP_2")
  
  src.EpcamE9.5.endoderm@reductions$tsne@cell.embeddings=
    cbind(src.EpcamE9.5.endoderm@reductions$tsne@cell.embeddings[,1],
          -src.EpcamE9.5.endoderm@reductions$tsne@cell.embeddings[,2])
  colnames(src.EpcamE9.5.endoderm@reductions$tsne@cell.embeddings)=c("tSNE_1","tSNE_2")
  src.EpcamE9.5.endoderm@reductions$umap@cell.embeddings=
    cbind(src.EpcamE9.5.endoderm@reductions$umap@cell.embeddings[,2],
          src.EpcamE9.5.endoderm@reductions$umap@cell.embeddings[,1])
  colnames(src.EpcamE9.5.endoderm@reductions$umap@cell.embeddings)=c("UMAP_1","UMAP_2")
  
}
#==========================================================


#==========================================================
#>>>> (1.1.4) endoderm-integration 
#==========================================================
src.12ss.integrated.selectgene=union(src.12ss.endoderm.selectgene,src.Epcam12SS.endoderm.selectgene)
src.15ss.integrated.selectgene=union(src.15ss.endoderm.selectgene,src.Epcam15SS.endoderm.selectgene)
src.18ss.integrated.selectgene=union(src.18ss.endoderm.selectgene,src.Epcam18SS.endoderm.selectgene)
src.21ss.integrated.selectgene=union(src.21ss.endoderm.selectgene,src.Epcam21SS.endoderm.selectgene)
src.24ss.integrated.selectgene=union(src.24ss.endoderm.selectgene,src.Epcam24SS.endoderm.selectgene)
src.27ss.integrated.selectgene=union(src.27ss.endoderm.selectgene,src.EpcamE9.5.endoderm.selectgene)

src.12ss.integrated=merge(src.12ss.endoderm, src.Epcam12SS.endoderm)
src.15ss.integrated=merge(src.15ss.endoderm, src.Epcam15SS.endoderm)
src.18ss.integrated=merge(src.18ss.endoderm, src.Epcam18SS.endoderm)
src.21ss.integrated=merge(src.21ss.endoderm, src.Epcam21SS.endoderm)
src.24ss.integrated=merge(src.24ss.endoderm, src.Epcam24SS.endoderm)
src.27ss.integrated=merge(src.27ss.endoderm, src.EpcamE9.5.endoderm)

for (seurat in c("src.12ss.integrated",
                 "src.15ss.integrated",
                 "src.18ss.integrated",
                 "src.21ss.integrated",
                 "src.24ss.integrated",
                 "src.27ss.integrated")) {
  assign(seurat,
         ScaleData(get(seurat),features = get(paste(seurat,".selectgene",sep = ""))))
  assign(seurat,
         RunPCA(get(seurat), features = get(paste(seurat,".selectgene",sep = ""))))
  assign(paste("mnnPCA.",seurat,sep = ""),
         fastMNN(get(seurat)@reductions$pca@cell.embeddings[!is.na(get(seurat)$Sort),1:30],
                 get(seurat)@reductions$pca@cell.embeddings[is.na(get(seurat)$Sort),1:30],
                 pc.input=TRUE))
  assign(paste("mnnPCA.",seurat,sep = ""),
         get(paste("mnnPCA.",seurat,sep = ""))$corrected[colnames(get(seurat)),])
  eval(parse(text = paste(seurat,"@reductions$mnnpca=CreateDimReducObject(mnnPCA.",seurat,",key = 'PC')",sep = "")))
  assign(seurat,
         RunTSNE(get(seurat), reduction = "mnnpca", dims = 1:30,check_duplicates = F,
                 perplexity=round((30+ncol(get(seurat))/100)),
                 dim.embed = 2,seed.use = 10))
  assign(seurat,
         RunUMAP(get(seurat), reduction = "mnnpca", 
                 dims = 1:30, n.neighbors=100, n.components=2))
}

save(src.12ss.integrated, src.12ss.integrated.selectgene, mnnPCA.src.12ss.integrated,
     src.15ss.integrated, src.15ss.integrated.selectgene, mnnPCA.src.15ss.integrated,
     src.18ss.integrated, src.18ss.integrated.selectgene, mnnPCA.src.18ss.integrated,
     src.21ss.integrated, src.21ss.integrated.selectgene, mnnPCA.src.21ss.integrated,
     src.24ss.integrated, src.24ss.integrated.selectgene, mnnPCA.src.24ss.integrated,
     src.27ss.integrated, src.27ss.integrated.selectgene, mnnPCA.src.27ss.integrated,
     file = "integrated.RData")


#---->>>>>>   merge from 9ss to 27ss

src.endoderm=merge(src.EpcamE8.5.endoderm,
                   c(src.12ss.integrated[,!is.na(src.12ss.integrated$Sort)],
                     src.12ss.integrated[,is.na(src.12ss.integrated$Sort)],
                     src.15ss.integrated[,!is.na(src.15ss.integrated$Sort)],
                     src.15ss.integrated[,is.na(src.15ss.integrated$Sort)],
                     src.18ss.integrated[,!is.na(src.18ss.integrated$Sort)],
                     src.18ss.integrated[,is.na(src.18ss.integrated$Sort)],
                     src.21ss.integrated[,!is.na(src.21ss.integrated$Sort)],
                     src.21ss.integrated[,is.na(src.21ss.integrated$Sort)],
                     src.24ss.integrated[,!is.na(src.24ss.integrated$Sort)],
                     src.24ss.integrated[,is.na(src.24ss.integrated$Sort)],
                     src.27ss.integrated[,!is.na(src.27ss.integrated$Sort)],
                     src.27ss.integrated[,is.na(src.27ss.integrated$Sort)]),
                   add.cell.ids = c("ss9","ss12","ss12","ss15","ss15","ss18","ss18","ss21","ss21","ss24","ss24","ss27","ss27"))

src.endoderm=FindVariableFeatures(src.endoderm, selection.method = "vst", nfeatures = 3000)
VariableFeaturePlot(src.endoderm)
src.endoderm.selectgene=
  Myfilter(as.matrix(src.endoderm@assays$RNA@data),
           gene = src.endoderm@assays$RNA@var.features,
           pearson.threshold = 0.2,partner.threshold = 5,
           bottom.dispersion.interval = 0.05)
r.cov <- WGCNA::cor(t(as.matrix(src.endoderm@assays$RNA@data[src.endoderm.selectgene,] /
                                  rowMax(as.matrix(src.endoderm@assays$RNA@data[src.endoderm.selectgene,])))),
                    method = "s")
r.hc <- hclust(as.dist(abs(1 - r.cov)),
               method = "ward.D2")
src.endoderm.rowtree=as.dendrogram(r.hc)
plot(src.endoderm.rowtree)
src.endoderm.selectgene=setdiff(src.endoderm.selectgene,
                                c(labels(src.endoderm.rowtree[[2]][[1]]),
                                  labels(src.endoderm.rowtree[[2]][[2]][[2]][[2]])))
for (seurat in c("src.12ss.integrated",
                 "src.15ss.integrated",
                 "src.18ss.integrated",
                 "src.21ss.integrated",
                 "src.24ss.integrated",
                 "src.27ss.integrated")) {
  assign(paste("mnnlogtpm.",seurat,sep = ""),
         mnnCorrect(as.matrix(get(seurat)@assays$RNA@data[src.endoderm.selectgene,!is.na(get(seurat)$Sort)]),
                    as.matrix(get(seurat)@assays$RNA@data[src.endoderm.selectgene,is.na(get(seurat)$Sort)]),
                    k = 5,cos.norm.out=F))
  assign(paste("mnnlogtpm.",seurat,sep = ""),
         do.call("cbind", get(paste("mnnlogtpm.",seurat,sep = ""))[[1]]))
}
mnnlogtpm=cbind(as.matrix(src.EpcamE8.5.endoderm@assays$RNA@data[src.endoderm.selectgene,]),
                mnnlogtpm.src.12ss.integrated,
                mnnlogtpm.src.15ss.integrated,
                mnnlogtpm.src.18ss.integrated,
                mnnlogtpm.src.21ss.integrated,
                mnnlogtpm.src.24ss.integrated,
                mnnlogtpm.src.27ss.integrated)
colnames(mnnlogtpm)=colnames(src.endoderm)
src.endoderm@assays$mnnRNA=CreateAssayObject(data = mnnlogtpm)
src.endoderm@assays$mnnRNA@var.features=src.endoderm.selectgene
DefaultAssay(src.endoderm)="mnnRNA"
src.endoderm=ScaleData(src.endoderm, features = src.endoderm.selectgene)
src.endoderm=RunPCA(src.endoderm, features = src.endoderm.selectgene)
src.endoderm=FindNeighbors(src.endoderm,dims = 1:30)
fdl.endoderm=layout_with_fr(graph.adjacency(src.endoderm@graphs$mnnRNA_nn),
                            dim = 3)
src.endoderm=RunUMAP(src.endoderm,dims = 1:30,n.neighbors = 200,n.components = 3)
src.endoderm@reductions$umap@cell.embeddings[,1]=-src.endoderm@reductions$umap@cell.embeddings[,1]
src.endoderm@reductions$umap@cell.embeddings[,2]=-src.endoderm@reductions$umap@cell.embeddings[,2]
src.endoderm@reductions$umap@cell.embeddings[,3]=-src.endoderm@reductions$umap@cell.embeddings[,3]

save(src.endoderm, file = "src.endoderm.Rdata")

#==========================================================


















