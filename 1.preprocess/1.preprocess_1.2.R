#===============================================================================
#>> 1.Pre-processing of scRNA-seq data
#>  1.2 Pre-processing of Smart-seq3
#===============================================================================
source("~/bin/MyFunction.R")
source("~/bin/MyFunction.Seurat3.R")
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(rgl)
library(scran)
library(batchelor)
library(patchwork)
library(SeuratWrappers)

gi=read.csv("~/genome/mm10/mm10.gene.inf.merge.10X.v1.csv",stringsAsFactors = F)
gi$gene=gi$Symbol
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
#>>>> (1.2.1) Load RawData
#==========================================================
#> Load the latest batch of updated data. for example: 
load("mm10_20230530_ypl.RData") 

data.all=data.all
sample.inf=sample.inf

data.name=as.data.frame(data.all)
data.name$EnsemblGeneID=rownames(data.name)
data.name=left_join(data.name,gi,by="EnsemblGeneID")
data.name.select=data.name[!data.name$Symbol2%in%NA,]
rownames(data.name.select)=data.name.select$Symbol2
data.name.select=data.name.select[1:length(colnames(data.all))]

src=CreateSeuratObject(data.name.select)
src=AddMetaData(src,sample.inf)
src$Samplename=sample.inf$SampleName
src$Time=sample.inf$Time
src$Contributor=sample.inf$Contributor
src$SeqData=sample.inf$SeqDate
src$Mouse=sample.inf$Mouse
factor(src$Contributor)

#-- Select our data
src.update.pre = subset(src,cells=rownames(src@meta.data[src$Contributor%in%c("ZQQ",'LKR'),])) 
rm(src)

update_count = src.update.pre@assays$RNA@counts
update_metadata = src.update.pre@meta.data
update_count_tpm = 
  MyCalTpm(read.count.data = as.matrix(update_count),
           gene.length = as.numeric(gi[rownames(update_count),]$GeneLength))
#rownames(MyGeneExp(update_count_tpm,100000,10))
#update_count_tpm.log2 = log2(update_count_tpm + 1)
update_count_tpm.ln = log10(update_count_tpm + 1)

cell_exclude_update = 
  colnames(update_count_tpm.ln)[update_count_tpm.ln[1,]=="NaN"]
update_count = update_count[,!colnames(update_count)%in%cell_exclude_update]
update_count_tpm = update_count_tpm[,!colnames(update_count_tpm)%in%cell_exclude_update]
#update_count_tpm.log2 = update_count_tpm.log2[,!colnames(update_count_tpm.log2)%in%cell_exclude_update]
update_count_tpm.ln = update_count_tpm.ln[,!colnames(update_count_tpm.ln)%in%cell_exclude_update]
update_metadata = update_metadata[!rownames(update_metadata)%in%cell_exclude_update,]

src.sm3.update.rawdata = CreateSeuratObject(update_count_tpm)
src.sm3.update.rawdata[["RNA"]]@data = 
  as(as.matrix(update_count_tpm.ln), "dgCMatrix")
rm(update_count,update_count_tpm,update_count_tpm.ln,src.update.pre)
DefaultAssay(src.sm3.update.rawdata) = "RNA"
src.sm3.update.rawdata = 
  FindVariableFeatures(
    src.sm3.update.rawdata, 
    selection.method = "vst", nfeatures = 2000)# vst
src.sm3.update.rawdata = AddMetaData(
  src.sm3.update.rawdata,
  metadata = update_metadata)

src.sm3.update.rawdata = src.sm3.update.rawdata[,src.sm3.update.rawdata$SeqDate%in%c("20230530")]                                          
#==========================================================


#==========================================================
#>>>> (1.2.2) QC
#==========================================================
#>>> provided a QC example for each batch.
#>  
#> Mnx1-CreER-GFP;Rosa-RFP 
src.sm3.Mnx1GFP.rawdata.230530 = src.sm3.update.rawdata[,src.sm3.update.rawdata$Mouse%in%"Mnx1CreERGFP;RossaRFP"]
src.sm3.Mnx1GFP.rawdata = merge(src.sm3.Mnx1GFP.rawdata,src.sm3.Mnx1GFP.rawdata.230530)
src.sm3.Mnx1GFP.rawdata$lineage = "Mnx1"

#> Ngn3-CreER-Ai6  
src.sm3.Ngn3Ai6.rawdata.230530 = src.sm3.update.rawdata[,src.sm3.update.rawdata$Mouse%in%"Neurog3CreAi6"]
src.sm3.Ngn3Ai6.rawdata = merge(src.sm3.Ngn3Ai6.rawdata,src.sm3.Ngn3Ai6.rawdata.230530)
src.sm3.Ngn3Ai6.rawdata$lineage = "Neurog3"
unique(src.sm3.Ngn3Ai6.rawdata$SeqDate)

#> Ngn3-GFP-HM&HZ
src.sm3.Ngn3GFP.rawdata.230530 = src.sm3.update.rawdata[,src.sm3.update.rawdata$Mouse%in%c("Neurog3GFP_HM","Neurog3GFP_HZ")]
src.sm3.Ngn3GFP.rawdata = merge(src.sm3.Ngn3GFP.rawdata, src.sm3.Ngn3GFP.rawdata.230530)
src.sm3.Ngn3GFP.rawdata$lineage = "Neurog3"
unique(src.sm3.Ngn3GFP.rawdata$SeqDate)

#> Pdx1GFP  
src.sm3.Pdx1.rawdata.230530 = src.sm3.update.rawdata[,src.sm3.update.rawdata$Mouse%in%"Pdx1GFP"]
src.sm3.Pdx1.rawdata = merge(src.sm3.Pdx1.rawdata,src.sm3.Pdx1.rawdata.230530)
src.sm3.Pdx1.rawdata$lineage = "Pdx1"


#>>> re-QC based on merge-data
for(seurat in c("src.sm3.Mnx1GFP.rawdata",
                "src.sm3.Ngn3Ai6.rawdata",
                "src.sm3.Ngn3GFP.rawdata",
                "src.sm3.Pdx1.rawdata")){
  
  mt.gene=grep("^mt-",rownames(get(seurat)),value = T)
  mt.gene.exp=colMeans(get(seurat)@assays$RNA@data[mt.gene,])
  assign(seurat,
         AddMetaData(get(seurat),mt.gene.exp,"Mito.gene"))
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


#> Mnx1-CreER-GFP;Rosa-RFP 
DimPlot(src.sm3.Mnx1GFP.rawdata,group.by = "RNA_snn_res.3",label = T)
DimPlot(src.sm3.Mnx1GFP.rawdata,group.by = "cluster")
DimPlot(src.sm3.Mnx1GFP.rawdata,group.by = "SeqDate")
FeaturePlot(src.sm3.Mnx1GFP.rawdata,
            features = c("Epcam","Trap1a",
                         "Foxi2","Dlx5","Tfap2a",
                         "Noto", "Nkx1-2","Snai1","Pou5f1"))

src.sm3.Mnx1GFP.rawdata$cluster="Endoderm"
src.sm3.Mnx1GFP.rawdata@meta.data[src.sm3.Mnx1GFP.rawdata@assays$RNA@data["Trap1a",]>1.5,]$cluster="Yolk sac"
src.sm3.Mnx1GFP.rawdata@meta.data[src.sm3.Mnx1GFP.rawdata$nFeature_RNA<6000,]$cluster='Low quality'
src.sm3.Mnx1GFP.rawdata@meta.data[src.sm3.Mnx1GFP.rawdata$RNA_snn_res.3%in%c(17),]$cluster='Mesoderm'
src.sm3.Mnx1GFP.rawdata@meta.data[src.sm3.Mnx1GFP.rawdata$RNA_snn_res.3%in%c(0),]$cluster='Ectoderm'

#> Ngn3-CreER-Ai6  
DimPlot(src.sm3.Ngn3Ai6.rawdata,group.by = "RNA_snn_res.1.2",label = T)
DimPlot(src.sm3.Ngn3Ai6.rawdata,group.by = "cluster")
#DimPlot(src.sm3.Ngn3Ai6.rawdata,group.by = "Time")
FeaturePlot(src.sm3.Ngn3Ai6.rawdata,
            features = c("Epcam","Trap1a",
                         "Foxi2","Dlx5","Tfap2a",
                         "Noto", "Nkx1-2","Snai1","Pou5f1"))

src.sm3.Ngn3Ai6.rawdata$cluster="Endoderm"
src.sm3.Ngn3Ai6.rawdata@meta.data[src.sm3.Ngn3Ai6.rawdata@assays$RNA@data["Trap1a",]>1.5,]$cluster="Yolk sac"
src.sm3.Ngn3Ai6.rawdata@meta.data[src.sm3.Ngn3Ai6.rawdata$nFeature_RNA<6000,]$cluster='Low quality'
src.sm3.Ngn3Ai6.rawdata@meta.data[src.sm3.Ngn3Ai6.rawdata$RNA_snn_res.1.2%in%c(0,2),]$cluster='Ectoderm'

#> Ngn3-GFP-HM&HZ
DimPlot(src.sm3.Ngn3GFP.rawdata,group.by = "RNA_snn_res.3",label = T)
DimPlot(src.sm3.Ngn3GFP.rawdata,group.by = "cluster")
FeaturePlot(src.sm3.Ngn3GFP.rawdata,
            features = c("Epcam","Trap1a",
                         "Foxi2","Dlx5","Tfap2a",
                         "Noto", "Nkx1-2","Snai1","Pou5f1"))

src.sm3.Ngn3GFP.rawdata$cluster="Endoderm"
src.sm3.Ngn3GFP.rawdata@meta.data[src.sm3.Ngn3GFP.rawdata@assays$RNA@data["Trap1a",]>1.5,]$cluster="Yolk sac"
src.sm3.Ngn3GFP.rawdata@meta.data[src.sm3.Ngn3GFP.rawdata$nFeature_RNA<6000,]$cluster='Low quality'
src.sm3.Ngn3GFP.rawdata@meta.data[src.sm3.Ngn3GFP.rawdata$RNA_snn_res.3%in%c(5,12),]$cluster='Ectoderm'
src.sm3.Ngn3GFP.rawdata@meta.data[src.sm3.Ngn3GFP.rawdata$RNA_snn_res.3%in%c(3),]$cluster='Yolk sac'
src.sm3.Ngn3GFP.rawdata@meta.data[src.sm3.Ngn3GFP.rawdata$RNA_snn_res.3%in%c(5,12),]$cluster='Ectoderm'

#> Pdx1GFP  
DimPlot(src.sm3.Pdx1.rawdata,group.by = "RNA_snn_res.3",label = T)
DimPlot(src.sm3.Pdx1.rawdata,group.by = "cluster")
DimPlot(src.sm3.Pdx1.rawdata,group.by = "SeqDate")
FeaturePlot(src.sm3.Pdx1.rawdata,
            features = c("Noto","Epcam","T","Nkx1-2","Gata1","Snai1"))

src.sm3.Pdx1.rawdata$cluster="Endoderm"
src.sm3.Pdx1.rawdata@meta.data[src.sm3.Pdx1.rawdata@assays$RNA@data["Trap1a",]>1.5,]$cluster="Yolk sac"
src.sm3.Pdx1.rawdata@meta.data[src.sm3.Pdx1.rawdata$nFeature_RNA<6000,]$cluster='Low quality'
src.sm3.Pdx1.rawdata@meta.data[src.sm3.Pdx1.rawdata$RNA_snn_res.3%in%c(2,4),]$cluster='Mesoderm'
src.sm3.Pdx1.rawdata@meta.data[src.sm3.Pdx1.rawdata$RNA_snn_res.3%in%c(21),]$cluster='Erythrocyte'
src.sm3.Pdx1.rawdata@meta.data[src.sm3.Pdx1.rawdata$RNA_snn_res.3%in%c(4),]$cluster='Notochord'
src.sm3.Pdx1.rawdata@meta.data[src.sm3.Pdx1.rawdata@assays$RNA@data["Pou5f1",]>0.5,]$cluster='Primordial Germ Cell'

save(src.sm3.Mnx1GFP.rawdata,
     src.sm3.Ngn3Ai6.rawdata,
     src.sm3.Ngn3GFP.rawdata,
     src.sm3.Pdx1.rawdata,
     file = "Update.20230530_20230602.Rdata")


#>>> update merged rawdata: 9ss & non-9ss
src.sm3.rawdata.temp = merge(src.sm3.Mnx1GFP.rawdata,
                             c(src.sm3.Ngn3Ai6.rawdata,
                               src.sm3.Ngn3GFP.rawdata,
                               src.sm3.Pdx1.rawdata))
src.sm3.rawdata.temp = src.sm3.rawdata.temp[,src.sm3.rawdata.temp$SeqDate%in%"20230530"]
src.sm3.rawdata.merge = merge(
  src.sm3.rawdata.merge[,setdiff(colnames(src.sm3.rawdata.merge), colnames(src.sm3.rawdata.temp))],
  src.sm3.rawdata.temp)

src.sm3.rawdata.merge.9ss = src.sm3.rawdata.merge[,src.sm3.rawdata.merge$Time%in%"9ss"]
src.sm3.rawdata.merge.no9ss = src.sm3.rawdata.merge[,!src.sm3.rawdata.merge$Time%in%"9ss"]

for(seurat in c("src.sm3.rawdata.merge.9ss",
                "src.sm3.rawdata.merge.no9ss")){
  
  mt.gene=grep("^mt-",rownames(get(seurat)),value = T)
  mt.gene.exp=colMeans(get(seurat)@assays$RNA@data[mt.gene,])
  assign(seurat,
         AddMetaData(get(seurat),mt.gene.exp,"Mito.gene"))
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


pdf("update.qc.summary.tracing.pdf",9,7)
for (seurat in c("src.sm3.rawdata.merge.9ss",
                 "src.sm3.rawdata.merge.no9ss")){
  print(
    DimPlot(get(seurat), reduction = "umap",
            group.by = "cluster", cols = qc.color,pt.size = 2.5)+
      theme(axis.ticks =  element_blank(),axis.text = element_blank()) +
      guides(colour = guide_legend(title = "",override.aes = list(size = 4)))+
      theme_bw() +
      #scale_color_manual(values = c(colors.seq)) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            aspect.ratio=1)+
      theme(text=element_text(size=16,face = 'bold'))
  )
}
dev.off()
#==========================================================


#==========================================================
#>>>> (1.2.3) Endoderm
#==========================================================
src.sm3.Pax9.endoderm = src.sm3.Pax9.rawdata[,src.sm3.Pax9.rawdata$cluster%in%"Endoderm"]
src.sm3.Sox2.endoderm = src.sm3.Sox2.rawdata[,src.sm3.Sox2.rawdata$cluster%in%"Endoderm"]
src.sm3.Nepn.endoderm = src.sm3.Nepn.rawdata[,src.sm3.Nepn.rawdata$cluster%in%"Endoderm"]
src.sm3.Wnt5b.endoderm = src.sm3.Wnt5b.rawdata[,src.sm3.Wnt5b.rawdata$cluster%in%"Endoderm"]
src.sm3.Nkx2_3.endoderm = src.sm3.Nkx2_3.rawdata[,src.sm3.Nkx2_3.rawdata$cluster%in%"Endoderm"]
src.sm3.Hhex.endoderm = src.sm3.Hhex.rawdata[,src.sm3.Hhex.rawdata$cluster%in%"Endoderm"]
src.sm3.Pdx1.endoderm = src.sm3.Pdx1.rawdata[,src.sm3.Pdx1.rawdata$cluster%in%"Endoderm"]
src.sm3.Mnx1GFP.endoderm = src.sm3.Mnx1GFP.rawdata[,src.sm3.Mnx1GFP.rawdata$cluster%in%"Endoderm"]
src.sm3.Mnx1.endoderm = src.sm3.Mnx1.rawdata[,src.sm3.Mnx1.rawdata$cluster%in%"Endoderm"]

src.sm3.Ngn3Ai6.endoderm = src.sm3.Ngn3Ai6.rawdata[,src.sm3.Ngn3Ai6.rawdata$cluster%in%"Endoderm"]
src.sm3.Ngn3GFP.endoderm = src.sm3.Ngn3GFP.rawdata[,src.sm3.Ngn3GFP.rawdata$cluster%in%"Endoderm"]

src.sm3.merge = merge(src.sm3.Pax9.endoderm, c(
  src.sm3.Sox2.endoderm, src.sm3.Nepn.endoderm, src.sm3.Wnt5b.endoderm,
  src.sm3.Nkx2_3.endoderm, src.sm3.Hhex.endoderm, src.sm3.Pdx1.endoderm,
  src.sm3.Mnx1GFP.endoderm, src.sm3.Mnx1.endoderm))

src.sm3.endoderm_Ngn3 = merge(src.sm3.Ngn3Ai6.endoderm, src.sm3.Ngn3GFP.endoderm)

for(seurat in c("src.sm3.Pax9.endoderm",
                "src.sm3.Sox2.endoderm",
                "src.sm3.Nepn.endoderm",
                "src.sm3.Wnt5b.endoderm",
                "src.sm3.Nkx2_3.endoderm",
                "src.sm3.Hhex.endoderm",
                "src.sm3.Pdx1.endoderm",
                "src.sm3.Mnx1GFP.endoderm",
                "src.sm3.Mnx1.endoderm",
                
                "src.sm3.Ngn3Ai6.endoderm",
                "src.sm3.Ngn3GFP.endoderm",
                
                "src.sm3.merge",
                "src.sm3.endoderm_Ngn3")){
  
  mt.gene=grep("^mt-",rownames(get(seurat)),value = T)
  mt.gene.exp=colMeans(get(seurat)@assays$RNA@data[mt.gene,])
  assign(seurat,
         AddMetaData(get(seurat),mt.gene.exp,"Mito.gene"))
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

save(src.sm3.Pax9.endoderm,
     src.sm3.Sox2.endoderm,
     src.sm3.Nepn.endoderm,
     src.sm3.Wnt5b.endoderm,
     src.sm3.Nkx2_3.endoderm,
     src.sm3.Hhex.endoderm,
     src.sm3.Pdx1.endoderm,
     src.sm3.Mnx1GFP.endoderm,
     src.sm3.Mnx1.endoderm,
     file = "src.sm3.tracing.summary.Rdata")

save(src.sm3.Ngn3Ai6.endoderm,
     src.sm3.Ngn3GFP.endoderm,
     file = "src.sm3.tracing.Ngn3.summary.Rdata")

save(src.sm3.merge, file = "src.sm3.merge.Rdata")
save(src.sm3.endoderm_Ngn3, file = "src.sm3.endoderm_Ngn3.Rdata")

#==========================================================


















