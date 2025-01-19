#===============================================================================
#>> 5.Simulation of LAI
#===============================================================================

#==========================================================
#>>>> (5.0) Set color-list
#==========================================================
color.pancreatic.type = c(c("#cb8d71","#8e9658",'#dc8bb1',"#a05f38","#a05f39","#c93832"),
                          c("#417ab2","#87548b","#76ad99","#dd803e","#62a05a"))
names(color.pancreatic.type) = c("EP1","EP2","EP3","EP4",'EP',"alpha/PP Pro.",
                                 "epsilon","beta","delta","alpha","PP")

color.pancreatic.time = c(c("#cd603e","#5c6ba1","#4c905d","#e09c4d","#d7b1c9","#4799bb","#798689",'#765f4a'))
names(color.pancreatic.time) = c('E13.5',"E14.5","E15.5","E16.5","E17.5","E18.5","P0","P3")

color.pancreatic.mouse = c(c("#5079b9","#dc7c3b","#966432","#966433","#d05a5f"),
                           c("#bac265","#8db9be","#669d56","#d2719c","#beaece"))
names(color.pancreatic.mouse) = c("Pdx1-Cre;Rosa-RFP","Ngn3-GFP","Ngn3-CreER-1D;Rosa-RFP","Ngn3-CreER-2D;Rosa-RFP","Ins1-RFP",
                                  "Gcg-Cre;Rosa-RFP","Gcg-GFP","Sst-BFP","Ppy-mNeptune","Ghrl-CFP")

color.pancreatic.merge = c(color.pancreatic.type, color.pancreatic.time, color.pancreatic.mouse)


color.pancreatic.type.plot = c("EP","alpha/PP Pro.","epsilon","beta","delta","alpha","PP")
color.pancreatic.mouse.plot = c("Pdx1-Cre;Rosa-RFP","Ngn3-CreER;Rosa-RFP","Ngn3-GFP",
                                "Ins1-RFP","Gcg-Cre;Rosa-RFP","Gcg-GFP",
                                "Sst-BFP","Ppy-mNeptune","Ghrl-CFP")

list.type.pancreatic = c("EP","alpha/PP Pro.","beta","alpha","delta","PP","epsilon")
matirx.type.pancreatic = cbind.data.frame(
  c(1,1,1,0,0,0,0,0,0), # EP
  c(1,1,1,0,0,0,0,0,0), # alpha/PP Pro.
  c(1,1,0,1,0,0,0,0,0), # beta
  c(1,1,0,0,1,1,0,0,0), # alpha
  c(1,1,0,0,0,0,1,0,0), # delta
  c(1,1,0,0,0,0,0,1,0),  # PP
  c(1,1,0,0,0,0,0,0,1) # epsilon
)
colnames(matirx.type.pancreatic) = c("EP","alpha/PP Pro.",'epsilon',"beta","delta","alpha","PP")
rownames(matirx.type.pancreatic) = c("Pdx1-Cre;Rosa-RFP","Ngn3-GFP","Ngn3-CreER;Rosa-RFP",
                                     "Ins1-RFP","Gcg-Cre;Rosa-RFP","Gcg-GFP",
                                     "Sst-BFP","Ppy-mNeptune","Ghrl-CFP")
matirx.type.pancreatic.plot = as.data.frame(melt(matirx.type.pancreatic))
matirx.type.pancreatic.plot$lineage = rep(c(1:9),7)
matirx.type.pancreatic.plot$lineage.value = matirx.type.pancreatic.plot$lineage * matirx.type.pancreatic.plot$value

pdf('Simulation/color.tile.lineage.pancreatic.pdf',5,5)
ggplot() + 
  geom_tile(data = matirx.type.pancreatic.plot,
            mapping = aes(x=lineage, y=rev(variable), fill=as.character(lineage.value)),
            width = 0.80, height= 0.8) +
  scale_fill_manual(values = c("#eeeeee", as.character(color.pancreatic.mouse[color.pancreatic.mouse.plot]))) +
  theme_void()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks =  element_blank(),
        axis.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        aspect.ratio = 1.1/1) +
  theme(legend.position = "none")
dev.off()
#==========================================================


#================================================================
#>>>> (5.1) Load data: pancreatic endocrine cell differentation
#================================================================
#> GEO: GSE139627, GSE115931, GSE87375
#> 
#> Cell selection Strategie was guided by: Figure 1a (Cell research, 2021, PMID: 33692492)
#> 
#> ref: Yu XX, Qiu WL, Yang L, Wang YC, He MY, Wang D, Zhang Y, Li LC, Zhang J, Wang Y, Xu CR. \
#> Sequential progenitor states mark the generation of pancreatic endocrine lineages in mice and humans. \
#> Cell Res. 2021 Aug;31(8):886-903. (DOI: 10.1038/s41422-021-00486-w) \
#> 
#> 
#> src.seu.Yu.all.10x: store 10x-dataset 
#> src.seu.data.sm2: store Smart-seq2-dataset 

src.seu.Yu.all.10x <- readRDS("Simulation/data_cell.res_pan/seu.Yu.all.10x.rds")
src.seu.data.sm2 <- readRDS("Simulation/data_cell.res_pan/seurat.data.sm2.rds")

fdl.src.seu.Yu.all.10x = layout_with_fr(graph_from_adjacency_matrix(src.seu.Yu.all.10x@graphs$RNA_nn), dim = 3)
fdl.src.seu.data.sm2 = layout_with_fr(graph_from_adjacency_matrix(src.seu.data.sm2@graphs$RNA_nn), dim = 3)

save(fdl.src.seu.Yu.all.10x, file = "Simulation/fdl.src.seu.Yu.all.10x.Rdata")
save(fdl.src.seu.data.sm2, file = "Simulation/fdl.src.seu.data.sm2.Rdata")

#-- FDL
for(i.src in c("src.seu.Yu.all.10x", "src.seu.data.sm2")){
  src = get(i.src)
  fdl.embedding = as.matrix(get(paste("fdl.", i.src, sep = "")))
  rownames(fdl.embedding) = colnames(src)
  colnames(fdl.embedding) = paste("FDL_", c(1:3), sep = "")
  
  src@reductions$fdl = CreateDimReducObject(embeddings = fdl.embedding, key = "FDL_")
  rm(fdl.embedding)
  
  assign(i.src, src)
}

save(src.seu.Yu.all.10x, file = "Simulation/src.seu.Yu.all.10x.Rdata")
save(src.seu.data.sm2, file = 'Simulation/src.seu.data.sm2.Rdata')


#-- FDL used Sm2-genes
src.seu.data.sm2.filtergene = 
  Myfilter(as.matrix(src.seu.data.sm2@assays$RNA@data),
           gene = src.seu.data.sm2@assays$RNA@var.features,
           bottom.dispersion.interval = 0.1,
           pearson.threshold = 0.15, partner.threshold = 5)

for(i.plot in c("tree.seu.data.sm2.filtergene")){
  src = src.seu.data.sm2
  gene.src = src.seu.data.sm2.filtergene
  pdf("Simulation/try.tree.pdf",9,7)
  for(i.tree in c("row")){
    src.rowtree =
      MyHeatmap(as.matrix(src@assays$RNA@data[gene.src,]),
                type = "row.relat",
                hc.c.data.type = "row.relat",
                hc.r.data.type = "row.relat",
                c.cov.method = "s",
                r.cov.method = "s",
                ColSideColors = cbind(
                  MyName2Col(src$Time, color.pancreatic.merge),
                  MyName2Col(src$CellType, color.pancreatic.merge)),
                ColSideColorsSize = 4,
                c.hc.method = "ward.D",
                r.hc.method = "ward.D2",
                #Rowv = "none",
                return.tree = "row",
                graph = T) 
  }
  dev.off()  
  assign(i.plot, as.dendrogram(src.rowtree))
  
  pdf("Simulation/try.tree.pdf",100,20)
  plot(as.dendrogram(src.rowtree))
  dev.off()
}

gene.src.seu.data.sm2.rowtree = setdiff(
  src.seu.data.sm2.filtergene,
  c(labels(tree.seu.data.sm2.filtergene[[1]][[2]][[1]])))

src.seu.Yu.all.10x.pca = RunPCA(src.seu.Yu.all.10x, features = intersect(gene.src.seu.data.sm2.rowtree, rownames(src.seu.Yu.all.10x)))
src.seu.Yu.all.10x.pca = FindNeighbors(src.seu.Yu.all.10x.pca, k.param = 30, assay = 'RNA')
fdl.src.seu.Yu.all.10x.pca = layout_with_fr(graph_from_adjacency_matrix(src.seu.Yu.all.10x.pca@graphs$RNA_nn), dim = 3)
save(fdl.src.seu.Yu.all.10x.pca, file = "Simulation/fdl.src.seu.Yu.all.10x.pca.Rdata")
rm(src.seu.Yu.all.10x.pca);gc()

#-- FDL-rotated (selected)
#--------------------------------------------------
#> Set src.seu.Yu.all.10x as example and
#> running this program locally
#> 
library("rgl")
data = cbind(src.seu.Yu.all.10x@meta.data,
             src.seu.Yu.all.10x@reductions$fdl@cell.embeddings)
colnames(data) = gsub("fdl","Coord",colnames(data))

open3d() 
view = par3d(family="arial", cex=20, font=1)
par3d(userMatrix = view)

plot3d(data[,c("Coord_1","Coord_2","Coord_3")],
       col=cluster.endoderm.color.v5[data$cluster.v06.26.re..correct..un],
       lwd = 1, type="p", size=5, axes= F,
       xlab = "", ylab = "", zlab = "",
       ticktype = "detailed")

axes3d(col = "black", lwd = 1, family = "serif", font = 1, marklen = 0,marklen.rel = F)
axes3d(c("x","y","z")) 
grid3d(c("x","y","z"),n = 8) 
close3d()

view_data = as.matrix(dput(par3d("userMatrix")))
view.fdl.src.seu.Yu.all.10x = view_data
src.seu.Yu.all.10x@reductions$fdl.rotated = src.seu.Yu.all.10x@reductions$fdl
data.temp =  view.fdl.src.seu.Yu.all.10x[1:3,1:3] %*% 
  t(as.matrix(src.seu.Yu.all.10x@reductions$fdl@cell.embeddings)); data.temp = t(data.temp)
colnames(data.temp) = paste("UMAP_", c(1:3), sep = "")
src.seu.Yu.all.10x@reductions$fdl.rotated@cell.embeddings = data.temp

fdl.rotated.src.seu.Yu.all.10x = src.seu.Yu.all.10x@reductions$fdl.rotated

# load("Simulation/fdl.rotated.src.seu.data.sm2.summary.Rdata")
# load("Simulation/fdl.rotated.src.seu.Yu.all.10x.summary.Rdata")

src.seu.Yu.all.10x@reductions$fdl.rotated = fdl.rotated.src.seu.Yu.all.10x
src.seu.data.sm2@reductions$fdl.rotated = fdl.rotated.src.seu.data.sm2
#--------------------------------------------------

pdf("Simulation/cell.res.data.pdf", 9, 7)
for(i.type in c("Time","CellTypeEnd","CellType")){
  print(DimPlot(src.seu.Yu.all.10x, group.by = i.type,  
                reduction = "fdl.rotated", dims = c(1,2),
                cols = color.pancreatic.merge)+
          theme_classic()+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.title.x = element_text(color="black", size=20, face="bold"),
                axis.title.y = element_text(color="black", size=20, face="bold"),
                axis.text.x = element_text(size = 15),
                axis.text.y = element_text(size = 15),
                axis.line.x = element_line(linetype=1, color="black", size=1.5),
                axis.line.y = element_line(linetype=1, color="black", size=1.5),
                axis.ticks.x = element_line(linetype=1, color="black", size=1.5),
                axis.ticks.y = element_line(linetype=1, color="black", size=1.5),
                axis.ticks.length  = unit(0.2, "cm"),
                aspect.ratio=1))
}
for(i.type in c("Time","CellType","Source","Source4")){
  print(DimPlot(src.seu.data.sm2, group.by = i.type, 
                reduction = "fdl.rotated", dims = c(1,2),
                cols = color.pancreatic.merge)+
          theme_classic()+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.title.x = element_text(color="black", size=20, face="bold"),
                axis.title.y = element_text(color="black", size=20, face="bold"),
                axis.text.x = element_text(size = 15),
                axis.text.y = element_text(size = 15),
                axis.line.x = element_line(linetype=1, color="black", size=1.5),
                axis.line.y = element_line(linetype=1, color="black", size=1.5),
                axis.ticks.x = element_line(linetype=1, color="black", size=1.5),
                axis.ticks.y = element_line(linetype=1, color="black", size=1.5),
                axis.ticks.length  = unit(0.2, "cm"),
                aspect.ratio=1))
}
dev.off()

#-- Add UMAP for src.seu.Yu.all.10x
src.seu.Yu.all.10x = RunUMAP(src.seu.Yu.all.10x, dims = 1:30,
                             reduction = "pca", n.neighbors = 100)
for(i.type in c("Time","CellTypeEnd","CellType")){
  print(DimPlot(src.seu.Yu.all.10x, group.by = i.type,  
                reduction = "umap", dims = c(1,2),
                cols = color.pancreatic.merge)+
          theme_classic()+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.title.x = element_text(color="black", size=20, face="bold"),
                axis.title.y = element_text(color="black", size=20, face="bold"),
                axis.text.x = element_text(size = 15),
                axis.text.y = element_text(size = 15),
                axis.line.x = element_line(linetype=1, color="black", size=1.5),
                axis.line.y = element_line(linetype=1, color="black", size=1.5),
                axis.ticks.x = element_line(linetype=1, color="black", size=1.5),
                axis.ticks.y = element_line(linetype=1, color="black", size=1.5),
                axis.ticks.length  = unit(0.2, "cm"),
                aspect.ratio=1))
}
#===============================================================================


#================================================================
#>>>> (5.2) Integration and LAI calculation
#================================================================

#-- Integration of 10x with Sm2 (whole datasets)
#--------------------------------------------------------------------------------
src.seu.data.sm2 = UpdateSeuratObject(src.seu.data.sm2)
src.seu.Yu.all.10x = UpdateSeuratObject(src.seu.Yu.all.10x)

src.seu.data.sm2.red = src.seu.data.sm2
src.seu.data.sm2.red = ScaleData(src.seu.data.sm2.red, features = rownames(src.seu.data.sm2.red))
src.seu.data.sm2.red@assays$RNA@counts = src.seu.data.sm2.red@assays$RNA@counts[rownames(src.seu.Yu.all.10x),]
src.seu.data.sm2.red@assays$RNA@data = src.seu.data.sm2.red@assays$RNA@data[rownames(src.seu.Yu.all.10x),]
src.seu.data.sm2.red@assays$RNA@scale.data = src.seu.data.sm2.red@assays$RNA@scale.data[rownames(src.seu.Yu.all.10x),]

src.cellres.pan.integrated.merge = merge(src.seu.Yu.all.10x, src.seu.data.sm2.red)
src.cellres.pan.integrated.merge = src.cellres.pan.integrated.merge[,src.cellres.pan.integrated.merge$Time%in%c("E13.5","E14.5","E15.5","E16.5","E17.5","E18.5")]

src.cellres.pan.integrated.merge.selectgene = intersect(
  unique(c(rownames(src.seu.data.sm2@reductions$pca@feature.loadings))),
  rownames(src.seu.Yu.all.10x))

src.cellres.pan.integrated.merge$Tech = ifelse(
  colnames(src.cellres.pan.integrated.merge)%in%colnames(src.seu.data.sm2), "Sm2", "10x")

MNN.res = mnnCorrect(
  as.matrix(src.cellres.pan.integrated.merge@assays$RNA@data[
    src.cellres.pan.integrated.merge.selectgene, 
    src.cellres.pan.integrated.merge$Tech%in%"10x"]),
  as.matrix(src.cellres.pan.integrated.merge@assays$RNA@data[
    src.cellres.pan.integrated.merge.selectgene,
    src.cellres.pan.integrated.merge$Tech%in%"Sm2"]),
  k = 5,cos.norm.out=F)

src.cellres.pan.integrated.merge@assays$mnnRNA = CreateAssayObject(data = MNN.res@assays@data$corrected)
src.cellres.pan.integrated.merge@assays$mnnRNA@key = "mnn_"
src.cellres.pan.integrated.merge = ScaleData(src.cellres.pan.integrated.merge, rownames(src.cellres.pan.integrated.merge@assays$mnnRNA), assay = "mnnRNA")


transfer.cellres.pan.integrated.merge = 
  FindTransferAnchors(reference = src.cellres.pan.integrated.merge[, intersect(colnames(src.cellres.pan.integrated.merge), colnames(src.seu.Yu.all.10x))],
                      query = src.cellres.pan.integrated.merge[, intersect(colnames(src.cellres.pan.integrated.merge), colnames(src.seu.data.sm2))],
                      reference.assay = "mnnRNA",query.assay = "mnnRNA", scale = T,
                      features = rownames(src.cellres.pan.integrated.merge@assays$mnnRNA@data))
fdl.transfer.cellres.pan.integrated.merge = TransferData(
  transfer.cellres.pan.integrated.merge,
  t(src.seu.Yu.all.10x@reductions$fdl.rotated@cell.embeddings[
    intersect(colnames(src.cellres.pan.integrated.merge), colnames(src.seu.Yu.all.10x)),]))

embedding.pan.integrated.merge.fdl.rotated = rbind(
  src.seu.Yu.all.10x@reductions$fdl.rotated@cell.embeddings[
    intersect(colnames(src.cellres.pan.integrated.merge), colnames(src.seu.Yu.all.10x)),],
  t(fdl.transfer.cellres.pan.integrated.merge@data))

src.cellres.pan.integrated.merge@reductions$fdl.rotated.fta = 
  CreateDimReducObject(embeddings = as.matrix(embedding.pan.integrated.merge.fdl.rotated),
                       key = "FDL_")
#--------------------------------------------------------------------------------

#-- Set function for Integration
batch_process_integration_index.pan = function(seurat, time){
  
  timess = "Sm2"
  sstime = "10x"
  seurat$Tech = ifelse(colnames(seurat)%in%colnames(src.seu.Yu.all.10x), sstime, timess)
  seurat.selectgene = get(paste("src.cellres.pan.", time, ".merge.selectgene",sep=""))
  seurat = ScaleData(seurat, split.by = "Tech")
  seurat = RunPCA(seurat, features = seurat.selectgene)
  seurat = RunUMAP(seurat, dims = 1:30,n.neighbors = 100)
  
  seurat@reductions$umap@cell.embeddings = 
    cbind(seurat@reductions$umap@cell.embeddings[,2],
          seurat@reductions$umap@cell.embeddings[,1])
  colnames(seurat@reductions$umap@cell.embeddings) = c("UMAP_1","UMAP_2")
  
  # MNN on Counts
  #---------------
  MNN.res = 
    mnnCorrect(as.matrix(seurat@assays$RNA@data[seurat.selectgene, seurat$Tech%in%sstime]),
               as.matrix(seurat@assays$RNA@data[seurat.selectgene, seurat$Tech%in%timess]),
               k = 5,cos.norm.out=F)
  seurat@assays$mnnRNA = CreateAssayObject(data = MNN.res@assays@data$corrected)
  seurat@assays$mnnRNA@key = "mnn_"
  seurat = ScaleData(seurat, 
                     rownames(seurat@assays$mnnRNA),
                     assay = "mnnRNA")
  
  seurat@meta.data = seurat@meta.data[!seurat$Tech%in%NA,]
  #---------------
  
  # MNN on PCA
  #---------------
  fMNN.res =  fastMNN(as.matrix(seurat@assays$RNA@data[seurat.selectgene, seurat$Tech%in%sstime]),
                      as.matrix(seurat@assays$RNA@data[seurat.selectgene, seurat$Tech%in%timess]),
                      
                      d=50,
                      subset.row = seurat.selectgene)
  seurat@reductions$mnn = seurat@reductions$umap
  seurat@reductions$mnn@cell.embeddings = fMNN.res@assays@data$reconstructed@seed@components
  seurat@reductions$mnn@feature.loadings = fMNN.res@assays@data$reconstructed@seed@rotation
  colnames(seurat@reductions$mnn@cell.embeddings) = paste("MNN",1:50,sep="_")
  colnames(seurat@reductions$mnn@feature.loadings) = paste("MNN",1:50,sep="_")
  seurat@reductions$mnn@key = "MNN_"
  #---------------
  
  # Finde integrated on RNA
  #--------------------------
  anchor.timess.integrated.merge = FindIntegrationAnchors(c(seurat[,seurat$Tech%in%timess],
                                                            seurat[,seurat$Tech%in%sstime]),
                                                          reduction = "cca",
                                                          assay = c("RNA","RNA"),
                                                          anchor.features = seurat.selectgene,
                                                          dims = 1:30)
  
  seurat.re = IntegrateData(anchorset = anchor.timess.integrated.merge , dims = 1:30)
  DefaultAssay(seurat.re) = "integrated"
  seurat.re = ScaleData(seurat.re, features = rownames(seurat.re), split.by = "Tech")
  seurat.re = RunPCA(seurat.re, features = seurat.selectgene)
  seurat.re = RunUMAP(seurat.re, dims = 1:30, n.neighbors = 100)
  
  seurat@reductions$umap_integrated = seurat.re@reductions$umap
  seurat@reductions$pca_integrated = seurat.re@reductions$pca
  
  seurat@reductions$umap_integrated@assay.used = "RNA"
  seurat@reductions$pca_integrated@assay.used = "RNA"
  
  rm(seurat.re,  anchor.timess.integrated.merge)
  #---------------
  
  # Finde integrated on mnnRNA
  #--------------------------
  anchor.timess.integrated.merge = FindIntegrationAnchors(c(seurat[,seurat$Tech%in%timess],
                                                            seurat[,seurat$Tech%in%sstime]),
                                                          reduction = "cca",
                                                          assay = c("mnnRNA","mnnRNA"),
                                                          anchor.features = seurat.selectgene,
                                                          dims = 1:30)
  
  seurat.re = IntegrateData(anchorset = anchor.timess.integrated.merge , dims = 1:30)
  DefaultAssay(seurat.re) = "integrated"
  seurat.re = ScaleData(seurat.re, features = rownames(seurat.re), split.by = "Tech")
  seurat.re = RunPCA(seurat.re, features = seurat.selectgene)
  seurat.re = RunUMAP(seurat.re, dims = 1:30, n.neighbors = 100)
  
  seurat@reductions$mnn_umap_integrated = seurat.re@reductions$umap
  seurat@reductions$mnn_pca_integrated = seurat.re@reductions$pca
  seurat@reductions$mnn_umap_integrated@assay.used = "mnnRNA"
  seurat@reductions$mnn_pca_integrated@assay.used = "mnnRNA"
  
  rm(seurat.re,  anchor.timess.integrated.merge)
  #---------------
  
  return(seurat)
  
}

#-- Integration of 10x with Sm2 at E15.5, E16.5 and E17.5
#--------------------------------------------------------------------------------
#>>>> select E15.5 E16.5, E17.5 for various lineages
src.cellres.pan.E15.5.merge = src.cellres.pan.integrated.merge[,src.cellres.pan.integrated.merge$Time%in%"E15.5"]
src.cellres.pan.E16.5.merge = src.cellres.pan.integrated.merge[,src.cellres.pan.integrated.merge$Time%in%"E16.5"]
src.cellres.pan.E17.5.merge = src.cellres.pan.integrated.merge[,src.cellres.pan.integrated.merge$Time%in%"E17.5"]

src.cellres.pan.E15.5.merge.selectgene = src.cellres.pan.integrated.merge.selectgene
src.cellres.pan.E16.5.merge.selectgene = src.cellres.pan.integrated.merge.selectgene
src.cellres.pan.E17.5.merge.selectgene = src.cellres.pan.integrated.merge.selectgene

src.cellres.pan.E15.5.merge = batch_process_integration_index.pan(seurat = src.cellres.pan.E15.5.merge, time = "E15.5")
src.cellres.pan.E16.5.merge = batch_process_integration_index.pan(seurat = src.cellres.pan.E16.5.merge, time = "E16.5")
src.cellres.pan.E17.5.merge = batch_process_integration_index.pan(seurat = src.cellres.pan.E17.5.merge, time = "E17.5")

src.cellres.pan.E15.5.merge.re = RunUMAP(src.cellres.pan.E15.5.merge, reduction = "mnn", dims = 1:30,n.neighbors = 100)
src.cellres.pan.E16.5.merge.re = RunUMAP(src.cellres.pan.E16.5.merge, reduction = "mnn", dims = 1:30,n.neighbors = 100)
src.cellres.pan.E17.5.merge.re = RunUMAP(src.cellres.pan.E17.5.merge, reduction = "mnn", dims = 1:30,n.neighbors = 100)
src.cellres.pan.E15.5.merge@reductions$mnn_umap = src.cellres.pan.E15.5.merge.re@reductions$umap
src.cellres.pan.E16.5.merge@reductions$mnn_umap = src.cellres.pan.E16.5.merge.re@reductions$umap
src.cellres.pan.E17.5.merge@reductions$mnn_umap = src.cellres.pan.E17.5.merge.re@reductions$umap
rm(src.cellres.pan.E15.5.merge.re,
   src.cellres.pan.E16.5.merge.re,
   src.cellres.pan.E17.5.merge.re)

pdf("Simulation/cell.res.pan.merge.pdf", 9, 7)
for(i.src in c("src.cellres.pan.E15.5.merge",
               "src.cellres.pan.E16.5.merge",
               "src.cellres.pan.E17.5.merge")){
  src = get(i.src)
  for(i.reduction in c("mnn_umap")){
    for(i.type in c("Time","CellTypeEnd","CellType","Source4")){
      print(DimPlot(src, group.by = i.type,  
                    reduction = i.reduction, dims = c(1,2),
                    cols = color.pancreatic.merge)+
              theme_classic()+
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.title.x = element_text(color="black", size=20, face="bold"),
                    axis.title.y = element_text(color="black", size=20, face="bold"),
                    axis.text.x = element_text(size = 15),
                    axis.text.y = element_text(size = 15),
                    axis.line.x = element_line(linetype=1, color="black", size=1.5),
                    axis.line.y = element_line(linetype=1, color="black", size=1.5),
                    axis.ticks.x = element_line(linetype=1, color="black", size=1.5),
                    axis.ticks.y = element_line(linetype=1, color="black", size=1.5),
                    axis.ticks.length  = unit(0.2, "cm"),
                    aspect.ratio=1))
    }
  }
}
dev.off()
#--------------------------------------------------------------------------------


#-- Set function for LAI calculation
network_extract_construct = function(seurat, reduction.list = c('mnn',"pca_integrated")){
  
  for (reduction in reduction.list) {
    seurat = FindNeighbors(seurat, reduction = reduction, dims = 1:30)
    graph_name = ifelse( reduction %in% c("pca_integrated"), "RNA", "RNA")
    cor.data = as.matrix(seurat@graphs[[ paste(graph_name,"_snn",sep="") ]])
    cor.data = as(cor.data, "dgCMatrix")
    assign(paste("cor.data",reduction,sep="_"), cor.data)
  }
  
  matrix = list(
    mnn = cor.data_mnn,
    pca_integrated = cor.data_pca_integrated,
    mnn_pca_integrated = cor.data_mnn_pca_integrated)
  
  return(matrix)
}


calculate_tracing_index = function(seurat, time, 
                                   type = c("Net_MinPath"), reduction = "mnn",
                                   lineage.list = c("Pdx1","Neurog3","Neruog3GFP",
                                                    "Ins1","Gcg","GcgGFP","Sst","Ppy","Ghrl")){
  
  seurat.merge = seurat # integrated.merge
  
  #-- load distance
  if(type == "Net_MinPath"){
    
    if(time %in% c("E15.5")){
      distance_index1 = get(paste("distance.pan.",time,".merge.index.1", sep=""))
      distance_index2 = get(paste("distance.pan.",time,".merge.index.2", sep=""))
      distance_index3 = get(paste("distance.pan.",time,".merge.index.3", sep=""))
      distance_index4 = get(paste("distance.pan.",time,".merge.index.4", sep=""))
    }else if(time %in% c("E16.5")){
      distance_index1 = get(paste("distance.pan.",time,".merge.index.1", sep=""))
      distance_index2 = get(paste("distance.pan.",time,".merge.index.2", sep=""))
      distance_index3 = get(paste("distance.pan.",time,".merge.index.3", sep=""))
      distance_index4 = get(paste("distance.pan.",time,".merge.index.4", sep=""))
      distance_index5 = get(paste("distance.pan.",time,".merge.index.5", sep=""))
      distance_index6 = get(paste("distance.pan.",time,".merge.index.6", sep=""))
      distance_index7 = get(paste("distance.pan.",time,".merge.index.7", sep=""))
      distance_index8 = get(paste("distance.pan.",time,".merge.index.8", sep=""))
      distance_index9 = get(paste("distance.pan.",time,".merge.index.9", sep=""))
    }else if(time %in% c("E17.5")){
      distance_index1 = get(paste("distance.pan.",time,".merge.index.1", sep=""))
      distance_index2 = get(paste("distance.pan.",time,".merge.index.2", sep=""))
      distance_index3 = get(paste("distance.pan.",time,".merge.index.3", sep=""))
      distance_index4 = get(paste("distance.pan.",time,".merge.index.4", sep=""))
      distance_index5 = get(paste("distance.pan.",time,".merge.index.5", sep=""))
    }else{
      return("Change Time !")
    }
    
  }else{
    return("Change type !")
  }
  
  #-- load seurat
  if(time %in% c("E15.5")){
    seurat_index1 = get(paste("src.cellres.pan.",time,".merge.index.1", sep=""))
    seurat_index2 = get(paste("src.cellres.pan.",time,".merge.index.2", sep=""))
    seurat_index3 = get(paste("src.cellres.pan.",time,".merge.index.3", sep=""))
    seurat_index4 = get(paste("src.cellres.pan.",time,".merge.index.4", sep=""))
  }else if(time %in% c("E16.5")){
    seurat_index1 = get(paste("src.cellres.pan.",time,".merge.index.1", sep=""))
    seurat_index2 = get(paste("src.cellres.pan.",time,".merge.index.2", sep=""))
    seurat_index3 = get(paste("src.cellres.pan.",time,".merge.index.3", sep=""))
    seurat_index4 = get(paste("src.cellres.pan.",time,".merge.index.4", sep=""))
    seurat_index5 = get(paste("src.cellres.pan.",time,".merge.index.5", sep=""))
    seurat_index6 = get(paste("src.cellres.pan.",time,".merge.index.6", sep=""))
    seurat_index7 = get(paste("src.cellres.pan.",time,".merge.index.7", sep=""))
    seurat_index8 = get(paste("src.cellres.pan.",time,".merge.index.8", sep=""))
    seurat_index9 = get(paste("src.cellres.pan.",time,".merge.index.9", sep=""))
  }else if(time %in% c("E17.5")){
    seurat_index1 = get(paste("src.cellres.pan.",time,".merge.index.1", sep=""))
    seurat_index2 = get(paste("src.cellres.pan.",time,".merge.index.2", sep=""))
    seurat_index3 = get(paste("src.cellres.pan.",time,".merge.index.3", sep=""))
    seurat_index4 = get(paste("src.cellres.pan.",time,".merge.index.4", sep=""))
    seurat_index5 = get(paste("src.cellres.pan.",time,".merge.index.5", sep=""))
  }else{
    return("Change Time !")
  }
  
  #-- Progress note (1)
  index_type = paste("tracing.index_",type,"_",reduction, sep="")
  print(paste(index_type, '\n'))
  
  #-- Set seurat_index_list & lineage_list 
  lineage_list = lineage.list
  
  if(time %in% c("E15.5")){
    seurat_index_list = c("seurat_index1","seurat_index2","seurat_index3","seurat_index4")
  }else if(time %in% c("E16.5")){
    seurat_index_list = c("seurat_index1","seurat_index2","seurat_index3",
                          "seurat_index4","seurat_index5","seurat_index6",
                          "seurat_index7","seurat_index8","seurat_index9")
  }else if(time %in% c("E17.5")){
    seurat_index_list = c("seurat_index1","seurat_index2","seurat_index3",
                          "seurat_index4","seurat_index5")
  }else{}
  
  
  if(type == "Net_MinPath"){
    for(seurat_temp_index in seurat_index_list){
      seurat_temp = get(seurat_temp_index)
      nwn_data = get(paste("distance", gsub("seurat","",seurat_temp_index), sep=""))
      nwn_data = nwn_data[[reduction]]
      for(lineage in lineage_list){
        name_list = rownames(seurat_temp@meta.data[seurat_temp$lineage%in%lineage,])
        seurat_temp[[paste(index_type, lineage, sep="_")]] = 0
        for(cell in name_list){
          if(cell %in% c(black_list_summary)){
            next()
          }
          cell_list = names(nwn_data[,cell][nwn_data[,cell]==1])
          seurat_temp@meta.data[cell_list, paste(index_type, lineage, sep="_")] = 
            seurat_temp@meta.data[cell_list, paste(index_type, lineage, sep="_")] + 1
        }
      }
      assign(seurat_temp_index, seurat_temp)
    }
  }else{}
  
  for(lineage in lineage_list){
    seurat.merge[[paste(index_type, lineage, sep="_")]] = 0
    for(seurat_temp_index in seurat_index_list){
      seurat_temp = get(seurat_temp_index)
      seurat.merge@meta.data[colnames(seurat_temp), paste(index_type, lineage, sep="_")] =
        seurat_temp@meta.data[colnames(seurat_temp), paste(index_type, lineage, sep="_")] 
    }
  }
  
  return(seurat.merge)
}

list.lineage.affinity.index_mnn = 
  paste("tracing.index_Net_MinPath_mnn_",
        c("Pdx1","Neurog3","Neruog3GFP",
          "Ins1","Gcg","GcgGFP","Sst",
          "Ppy","Ghrl"), sep = "")

list.lineage.affinity.index_pca_integrated = 
  paste("tracing.index_Net_MinPath_pca_integrated_",
        c("Pdx1","Neurog3","Neruog3GFP",
          "Ins1","Gcg","GcgGFP","Sst",
          "Ppy","Ghrl"), sep = "")


#-- LAI calculation
#--------------------------------------------------------------------------------
#> set replicate as 3
df.lineage.affinity.index = list()
for(i.replicate in c(1:3)){
  
  #-- split-refer: E15.5
  for(i.batch.merge in c("src.cellres.pan.E15.5.merge")){
    #-- E15.5
    src.cellres.pan.E15.5.merge$source_batch = ifelse(src.cellres.pan.E15.5.merge$Tech%in%"10x", "refer", 'query')
    
    #>> random set-number
    src.cellres.pan.E15.5.merge$index_refer_raw = NA
    src.cellres.pan.E15.5.merge@meta.data[
      src.cellres.pan.E15.5.merge$source_batch%in%"refer",]$index_refer_raw = 
      sample(x = c(1: table(src.cellres.pan.E15.5.merge$source_batch)["refer"]),
             size = table(src.cellres.pan.E15.5.merge$source_batch)["refer"])
    
    
    # sample by batch size
    query_E15.5_num = table(src.cellres.pan.E15.5.merge$source_batch)["query"]
    refer_E15.5_num = table(src.cellres.pan.E15.5.merge$source_batch)["refer"]
    size_E15.5 = floor(refer_E15.5_num / (query_E15.5_num * 1.5))
    print(size_E15.5)
    refer_E15.5_num_1 = floor(refer_E15.5_num/4)
    refer_E15.5_num_2 = floor(refer_E15.5_num*2/4)
    refer_E15.5_num_3 = floor(refer_E15.5_num*3/4)
    
    src.cellres.pan.E15.5.merge$index_refer = NA
    src.cellres.pan.E15.5.merge@meta.data[
      src.cellres.pan.E15.5.merge$index_refer_raw%in%c(1:refer_E15.5_num_1),]$index_refer = 1
    src.cellres.pan.E15.5.merge@meta.data[
      src.cellres.pan.E15.5.merge$index_refer_raw%in%c((refer_E15.5_num_1+1):(refer_E15.5_num_2)),]$index_refer = 2
    src.cellres.pan.E15.5.merge@meta.data[
      src.cellres.pan.E15.5.merge$index_refer_raw%in%c((refer_E15.5_num_2+1):(refer_E15.5_num_3)),]$index_refer = 3
    src.cellres.pan.E15.5.merge@meta.data[
      src.cellres.pan.E15.5.merge$index_refer_raw%in%c((refer_E15.5_num_3+1):(refer_E15.5_num)),]$index_refer = 4
    table(src.cellres.pan.E15.5.merge$index_refer)
    
    for(i.src in c("src.cellres.pan.E15.5.merge")){
      src = get(i.src)
      for(i.index in c(1,2,3,4)){
        src_index = merge(
          src[,colnames(src)%in%rownames(src@meta.data[src$index_refer%in%i.index,])],
          src[,colnames(src)%in%rownames(src@meta.data[src$source_batch%in%"query",])])
        assign(paste(i.src, ".index.", i.index, sep = ""), src_index)
      }
      rm(src, src_index)
    }
    
    src.cellres.pan.E15.5.merge.selectgene = src.cellres.pan.integrated.merge.selectgene
    
    #---- integrated for batch
    src.cellres.pan.E15.5.merge.index.1 = batch_process_integration_index.pan(src.cellres.pan.E15.5.merge.index.1, "E15.5")
    src.cellres.pan.E15.5.merge.index.2 = batch_process_integration_index.pan(src.cellres.pan.E15.5.merge.index.2, "E15.5")
    src.cellres.pan.E15.5.merge.index.3 = batch_process_integration_index.pan(src.cellres.pan.E15.5.merge.index.3, "E15.5")
    src.cellres.pan.E15.5.merge.index.4 = batch_process_integration_index.pan(src.cellres.pan.E15.5.merge.index.4, "E15.5")
  }
  
  #-- split-refer: E16.5
  for(i.batch.merge in c("src.cellres.pan.E16.5.merge")){
    #-- E16.5
    src.cellres.pan.E16.5.merge$source_batch = ifelse(src.cellres.pan.E16.5.merge$Tech%in%"10x", "refer", 'query')
    src.cellres.pan.E16.5.merge$index_refer_raw = NA
    src.cellres.pan.E16.5.merge@meta.data[
      src.cellres.pan.E16.5.merge$source_batch%in%"refer",]$index_refer_raw = 
      sample(x = c(1: table(src.cellres.pan.E16.5.merge$source_batch)["refer"]),
             size = table(src.cellres.pan.E16.5.merge$source_batch)["refer"])
    
    # sample by batch size
    query_E16.5_num = table(src.cellres.pan.E16.5.merge$source_batch)["query"]
    refer_E16.5_num = table(src.cellres.pan.E16.5.merge$source_batch)["refer"]
    size_E16.5 = floor(refer_E16.5_num / (query_E16.5_num * 1.5))
    print(size_E16.5)
    refer_E16.5_num_1 = floor(refer_E16.5_num/9)
    refer_E16.5_num_2 = floor(refer_E16.5_num*2/9)
    refer_E16.5_num_3 = floor(refer_E16.5_num*3/9)
    refer_E16.5_num_4 = floor(refer_E16.5_num*4/9)
    refer_E16.5_num_5 = floor(refer_E16.5_num*5/9)
    refer_E16.5_num_6 = floor(refer_E16.5_num*6/9)
    refer_E16.5_num_7 = floor(refer_E16.5_num*7/9)
    refer_E16.5_num_8 = floor(refer_E16.5_num*8/9)
    
    
    src.cellres.pan.E16.5.merge$index_refer = NA
    src.cellres.pan.E16.5.merge@meta.data[src.cellres.pan.E16.5.merge$index_refer_raw%in%c(1:refer_E16.5_num_1),]$index_refer = 1
    src.cellres.pan.E16.5.merge@meta.data[src.cellres.pan.E16.5.merge$index_refer_raw%in%c((refer_E16.5_num_1+1):(refer_E16.5_num_2)),]$index_refer = 2
    src.cellres.pan.E16.5.merge@meta.data[src.cellres.pan.E16.5.merge$index_refer_raw%in%c((refer_E16.5_num_2+1):(refer_E16.5_num_3)),]$index_refer = 3
    src.cellres.pan.E16.5.merge@meta.data[src.cellres.pan.E16.5.merge$index_refer_raw%in%c((refer_E16.5_num_3+1):(refer_E16.5_num_4)),]$index_refer = 4
    src.cellres.pan.E16.5.merge@meta.data[src.cellres.pan.E16.5.merge$index_refer_raw%in%c((refer_E16.5_num_4+1):(refer_E16.5_num_5)),]$index_refer = 5
    src.cellres.pan.E16.5.merge@meta.data[src.cellres.pan.E16.5.merge$index_refer_raw%in%c((refer_E16.5_num_5+1):(refer_E16.5_num_6)),]$index_refer = 6
    src.cellres.pan.E16.5.merge@meta.data[src.cellres.pan.E16.5.merge$index_refer_raw%in%c((refer_E16.5_num_6+1):(refer_E16.5_num_7)),]$index_refer = 7
    src.cellres.pan.E16.5.merge@meta.data[src.cellres.pan.E16.5.merge$index_refer_raw%in%c((refer_E16.5_num_7+1):(refer_E16.5_num_8)),]$index_refer = 8
    src.cellres.pan.E16.5.merge@meta.data[src.cellres.pan.E16.5.merge$index_refer_raw%in%c((refer_E16.5_num_8+1):(refer_E16.5_num)),]$index_refer = 9
    
    table(src.cellres.pan.E16.5.merge$index_refer)
    
    for(i.src in c("src.cellres.pan.E16.5.merge")){
      src = get(i.src)
      for(i.index in c(1,2,3,4,5,6,7,8,9)){
        src_index = merge(
          src[,colnames(src)%in%rownames(src@meta.data[src$index_refer%in%i.index,])],
          src[,colnames(src)%in%rownames(src@meta.data[src$source_batch%in%"query",])])
        assign(paste(i.src, ".index.", i.index, sep = ""), src_index)
      }
      rm(src, src_index)
    }
    
    src.cellres.pan.E16.5.merge.selectgene = src.cellres.pan.integrated.merge.selectgene
    
    #---- integrated for batch
    src.cellres.pan.E16.5.merge.index.1 = batch_process_integration_index.pan(src.cellres.pan.E16.5.merge.index.1, "E16.5")
    src.cellres.pan.E16.5.merge.index.2 = batch_process_integration_index.pan(src.cellres.pan.E16.5.merge.index.2, "E16.5")
    src.cellres.pan.E16.5.merge.index.3 = batch_process_integration_index.pan(src.cellres.pan.E16.5.merge.index.3, "E16.5")
    src.cellres.pan.E16.5.merge.index.4 = batch_process_integration_index.pan(src.cellres.pan.E16.5.merge.index.4, "E16.5")
    src.cellres.pan.E16.5.merge.index.5 = batch_process_integration_index.pan(src.cellres.pan.E16.5.merge.index.5, "E16.5")
    src.cellres.pan.E16.5.merge.index.6 = batch_process_integration_index.pan(src.cellres.pan.E16.5.merge.index.6, "E16.5")
    src.cellres.pan.E16.5.merge.index.7 = batch_process_integration_index.pan(src.cellres.pan.E16.5.merge.index.7, "E16.5")
    src.cellres.pan.E16.5.merge.index.8 = batch_process_integration_index.pan(src.cellres.pan.E16.5.merge.index.8, "E16.5")
    src.cellres.pan.E16.5.merge.index.9 = batch_process_integration_index.pan(src.cellres.pan.E16.5.merge.index.9, "E16.5")
  }
  
  #-- split-refer: E17.5
  for(i.batch.merge in c("src.cellres.pan.E17.5.merge")){
    #-- E17.5
    src.cellres.pan.E17.5.merge$source_batch = ifelse(src.cellres.pan.E17.5.merge$Tech%in%"10x", "refer", 'query')
    src.cellres.pan.E17.5.merge$index_refer_raw = NA
    src.cellres.pan.E17.5.merge@meta.data[
      src.cellres.pan.E17.5.merge$source_batch%in%"refer",]$index_refer_raw = 
      sample(x = c(1: table(src.cellres.pan.E17.5.merge$source_batch)["refer"]),
             size = table(src.cellres.pan.E17.5.merge$source_batch)["refer"])
    
    # sample by batch size
    query_E17.5_num = table(src.cellres.pan.E17.5.merge$source_batch)["query"]
    refer_E17.5_num = table(src.cellres.pan.E17.5.merge$source_batch)["refer"]
    size_E17.5 = floor(refer_E17.5_num / (query_E17.5_num * 1.5))
    print(size_E17.5)
    refer_E17.5_num_1 = floor(refer_E17.5_num/5)
    refer_E17.5_num_2 = floor(refer_E17.5_num*2/5)
    refer_E17.5_num_3 = floor(refer_E17.5_num*3/5)
    refer_E17.5_num_4 = floor(refer_E17.5_num*4/5)
    
    src.cellres.pan.E17.5.merge$index_refer = NA
    src.cellres.pan.E17.5.merge@meta.data[
      src.cellres.pan.E17.5.merge$index_refer_raw%in%c(1:refer_E17.5_num_1),]$index_refer = 1
    src.cellres.pan.E17.5.merge@meta.data[
      src.cellres.pan.E17.5.merge$index_refer_raw%in%c((refer_E17.5_num_1+1):(refer_E17.5_num_2)),]$index_refer = 2
    src.cellres.pan.E17.5.merge@meta.data[
      src.cellres.pan.E17.5.merge$index_refer_raw%in%c((refer_E17.5_num_2+1):(refer_E17.5_num_3)),]$index_refer = 3
    src.cellres.pan.E17.5.merge@meta.data[
      src.cellres.pan.E17.5.merge$index_refer_raw%in%c((refer_E17.5_num_3+1):(refer_E17.5_num_4)),]$index_refer = 4
    src.cellres.pan.E17.5.merge@meta.data[
      src.cellres.pan.E17.5.merge$index_refer_raw%in%c((refer_E17.5_num_4+1):(refer_E17.5_num)),]$index_refer = 5
    table(src.cellres.pan.E17.5.merge$index_refer)
    
    for(i.src in c("src.cellres.pan.E17.5.merge")){
      src = get(i.src)
      for(i.index in c(1,2,3,4,5)){
        src_index = merge(
          src[,colnames(src)%in%rownames(src@meta.data[src$index_refer%in%i.index,])],
          src[,colnames(src)%in%rownames(src@meta.data[src$source_batch%in%"query",])])
        assign(paste(i.src, ".index.", i.index, sep = ""), src_index)
      }
      rm(src, src_index)
    }
    
    src.cellres.pan.E17.5.merge.selectgene = src.cellres.pan.integrated.merge.selectgene
    
    #---- integrated for batch
    src.cellres.pan.E17.5.merge.index.1 = batch_process_integration_index.pan(src.cellres.pan.E17.5.merge.index.1, "E17.5")
    src.cellres.pan.E17.5.merge.index.2 = batch_process_integration_index.pan(src.cellres.pan.E17.5.merge.index.2, "E17.5")
    src.cellres.pan.E17.5.merge.index.3 = batch_process_integration_index.pan(src.cellres.pan.E17.5.merge.index.3, "E17.5")
    src.cellres.pan.E17.5.merge.index.4 = batch_process_integration_index.pan(src.cellres.pan.E17.5.merge.index.4, "E17.5")
    src.cellres.pan.E17.5.merge.index.5 = batch_process_integration_index.pan(src.cellres.pan.E17.5.merge.index.5, "E17.5")
  }
  
  for(i.network in c("src.cellres.pan.E15.5.merge")){
    network.pan.E15.5.merge.index.1 = network_extract_construct(src.cellres.pan.E15.5.merge.index.1)
    network.pan.E15.5.merge.index.2 = network_extract_construct(src.cellres.pan.E15.5.merge.index.2)
    network.pan.E15.5.merge.index.3 = network_extract_construct(src.cellres.pan.E15.5.merge.index.3)
    network.pan.E15.5.merge.index.4 = network_extract_construct(src.cellres.pan.E15.5.merge.index.4)
    
    distance.pan.E15.5.merge.index.1 = distance_NW_network_index(src.cellres.pan.E15.5.merge, network.pan.E15.5.merge.index.1)
    distance.pan.E15.5.merge.index.2 = distance_NW_network_index(src.cellres.pan.E15.5.merge, network.pan.E15.5.merge.index.2)
    distance.pan.E15.5.merge.index.3 = distance_NW_network_index(src.cellres.pan.E15.5.merge, network.pan.E15.5.merge.index.3)
    distance.pan.E15.5.merge.index.4 = distance_NW_network_index(src.cellres.pan.E15.5.merge, network.pan.E15.5.merge.index.4)
    
  }
  
  for(i.network in c("src.cellres.pan.E16.5.merge")){
    network.pan.E16.5.merge.index.1 = network_extract_construct(src.cellres.pan.E16.5.merge.index.1)
    network.pan.E16.5.merge.index.2 = network_extract_construct(src.cellres.pan.E16.5.merge.index.2)
    network.pan.E16.5.merge.index.3 = network_extract_construct(src.cellres.pan.E16.5.merge.index.3)
    network.pan.E16.5.merge.index.4 = network_extract_construct(src.cellres.pan.E16.5.merge.index.4)
    network.pan.E16.5.merge.index.5 = network_extract_construct(src.cellres.pan.E16.5.merge.index.5)
    
    network.pan.E16.5.merge.index.6 = network_extract_construct(src.cellres.pan.E16.5.merge.index.6)
    network.pan.E16.5.merge.index.7 = network_extract_construct(src.cellres.pan.E16.5.merge.index.7)
    network.pan.E16.5.merge.index.8 = network_extract_construct(src.cellres.pan.E16.5.merge.index.8)
    network.pan.E16.5.merge.index.9 = network_extract_construct(src.cellres.pan.E16.5.merge.index.9)
    
    distance.pan.E16.5.merge.index.1 = distance_NW_network_index(src.cellres.pan.E16.5.merge, network.pan.E16.5.merge.index.1)
    distance.pan.E16.5.merge.index.2 = distance_NW_network_index(src.cellres.pan.E16.5.merge, network.pan.E16.5.merge.index.2)
    distance.pan.E16.5.merge.index.3 = distance_NW_network_index(src.cellres.pan.E16.5.merge, network.pan.E16.5.merge.index.3)
    distance.pan.E16.5.merge.index.4 = distance_NW_network_index(src.cellres.pan.E16.5.merge, network.pan.E16.5.merge.index.4)
    distance.pan.E16.5.merge.index.5 = distance_NW_network_index(src.cellres.pan.E16.5.merge, network.pan.E16.5.merge.index.5)
    
    distance.pan.E16.5.merge.index.6 = distance_NW_network_index(src.cellres.pan.E16.5.merge, network.pan.E16.5.merge.index.6)
    distance.pan.E16.5.merge.index.7 = distance_NW_network_index(src.cellres.pan.E16.5.merge, network.pan.E16.5.merge.index.7)
    distance.pan.E16.5.merge.index.8 = distance_NW_network_index(src.cellres.pan.E16.5.merge, network.pan.E16.5.merge.index.8)
    distance.pan.E16.5.merge.index.9 = distance_NW_network_index(src.cellres.pan.E16.5.merge, network.pan.E16.5.merge.index.9)
    
  }
  
  for(i.network in c("src.cellres.pan.E17.5.merge")){
    network.pan.E17.5.merge.index.1 = network_extract_construct(src.cellres.pan.E17.5.merge.index.1)
    network.pan.E17.5.merge.index.2 = network_extract_construct(src.cellres.pan.E17.5.merge.index.2)
    network.pan.E17.5.merge.index.3 = network_extract_construct(src.cellres.pan.E17.5.merge.index.3)
    network.pan.E17.5.merge.index.4 = network_extract_construct(src.cellres.pan.E17.5.merge.index.4)
    network.pan.E17.5.merge.index.5 = network_extract_construct(src.cellres.pan.E17.5.merge.index.5)
    
    distance.pan.E17.5.merge.index.1 = distance_NW_network_index(src.cellres.pan.E17.5.merge, network.pan.E17.5.merge.index.1)
    distance.pan.E17.5.merge.index.2 = distance_NW_network_index(src.cellres.pan.E17.5.merge, network.pan.E17.5.merge.index.2)
    distance.pan.E17.5.merge.index.3 = distance_NW_network_index(src.cellres.pan.E17.5.merge, network.pan.E17.5.merge.index.3)
    distance.pan.E17.5.merge.index.4 = distance_NW_network_index(src.cellres.pan.E17.5.merge, network.pan.E17.5.merge.index.4)
    distance.pan.E17.5.merge.index.5 = distance_NW_network_index(src.cellres.pan.E17.5.merge, network.pan.E17.5.merge.index.5)
    
  }
  
  
  #-- Set lineage-mouse: src.cellres.pan.E15.5.merge
  for(i.src in c("","index.1","index.2","index.3","index.4")){
    if(i.src == ""){
      src.name = "src.cellres.pan.E15.5.merge"
    }else{
      src.name = paste("src.cellres.pan.E15.5.merge.", i.src, sep = "")
    }
    src = get(src.name)
    src$lineage = NA
    for(i.lineage in 1:length(list.mouse.cellres.pan)){
      if(length(names(table(grepl(list.mouse.cellres.pan[i.lineage], src$Source4))))==1){next()}
      src@meta.data[grepl(list.mouse.cellres.pan[i.lineage], src$Source4),]$lineage =
        list.lineage.cellres.pan[i.lineage]
    }
    assign(src.name, src)
  }
  
  #-- Set lineage-mouse: src.cellres.pan.E16.5.merge
  for(i.src in c(# "",
    "index.1","index.2","index.3","index.4","index.5",
    "index.6","index.7","index.8","index.9")){
    if(i.src == ""){
      src.name = "src.cellres.pan.E16.5.merge"
    }else{
      src.name = paste("src.cellres.pan.E16.5.merge.", i.src, sep = "")
    }
    src = get(src.name)
    src$lineage = NA
    for(i.lineage in 1:length(list.mouse.cellres.pan)){
      if(length(names(table(grepl(list.mouse.cellres.pan[i.lineage], src$Source4))))==1){next()}
      src@meta.data[grepl(list.mouse.cellres.pan[i.lineage], src$Source4),]$lineage =
        list.lineage.cellres.pan[i.lineage]
    }
    assign(src.name, src)
  }
  
  #-- Set lineage-mouse: src.cellres.pan.E17.5.merge
  for(i.src in c("","index.1","index.2","index.3","index.4","index.5")){
    if(i.src == ""){
      src.name = "src.cellres.pan.E17.5.merge"
    }else{
      src.name = paste("src.cellres.pan.E17.5.merge.", i.src, sep = "")
    }
    src = get(src.name)
    src$lineage = NA
    for(i.lineage in 1:length(list.mouse.cellres.pan)){
      if(length(names(table(grepl(list.mouse.cellres.pan[i.lineage], src$Source4))))==1){next()}
      src@meta.data[grepl(list.mouse.cellres.pan[i.lineage], src$Source4),]$lineage =
        list.lineage.cellres.pan[i.lineage]
    }
    assign(src.name, src)
  }
  
  
  #-- Network MinPath Metric of KNN
  type = c("Net_MinPath")
  reduction = c("mnn", "pca_integrated")
  
  #--- calculate_tracing_index: E15.5, E16.5 & E17.5
  for(i.time in c("E15.5","E16.5","E17.5")){
    seurat_name = paste("src.cellres.pan.", i.time, ".merge", sep="")
    seurat = get(seurat_name)
    for(i.type in c("Net_MinPath")){
      for(i.reduction in c("mnn", "pca_integrated")){
        seurat = calculate_tracing_index(seurat, 
                                         time = i.time, type = i.type, reduction = i.reduction,
                                         lineage.list = c("Pdx1","Neurog3","Neruog3GFP",
                                                          "Ins1","Gcg","GcgGFP","Sst","Ppy","Ghrl"))
      }
    }
    assign(seurat_name, seurat)
  }
  
  
  #>> add df.lineage.affinity.index
  for(i.time in c("E15.5","E16.5","E17.5")){
    seurat = get(paste("src.cellres.pan.", i.time, ".merge", sep=""))
    df.lineage.affinity.index[[paste(i.time,"_",i.replicate,sep="")]] = 
      seurat@meta.data[
        ,c(list.lineage.affinity.index_mnn,
           list.lineage.affinity.index_pca_integrated)]
  }
}

#> calulate mean and log(x+1)
for(i.time in c("E15.5","E16.5","E17.5")){
  seurat = get(paste("src.cellres.pan.", i.time, ".merge", sep=""))
  
  for(i.index in c(list.lineage.affinity.index_mnn,
                   list.lineage.affinity.index_pca_integrated)){
    seurat@meta.data[,i.index] = log(
      apply(cbind(
        df.lineage.affinity.index[[paste(i.time,"_1",sep="")]][,i.index],
        df.lineage.affinity.index[[paste(i.time,"_2",sep="")]][,i.index],
        df.lineage.affinity.index[[paste(i.time,"_3",sep="")]][,i.index]), 1, mean) + 1
    )
  }
  
  assign(paste("src.cellres.pan.", i.time, ".merge", sep=""), seurat)
}

#--------------------------------------------------------------------------------


src.seu.Yu.all.10x.E15.5 = src.seu.Yu.all.10x[,src.seu.Yu.all.10x$Time%in%"E15.5"]
src.seu.Yu.all.10x.E16.5 = src.seu.Yu.all.10x[,src.seu.Yu.all.10x$Time%in%"E16.5"]
src.seu.Yu.all.10x.E17.5 = src.seu.Yu.all.10x[,src.seu.Yu.all.10x$Time%in%"E17.5"]

src.seu.Yu.all.10x.E15.5 = AddMetaData(
  src.seu.Yu.all.10x.E15.5, 
  metadata = src.cellres.pan.E15.5.merge@meta.data[
    colnames(src.seu.Yu.all.10x.E15.5),
    setdiff(colnames(src.cellres.pan.E15.5.merge@meta.data),
            colnames(src.seu.Yu.all.10x.E15.5@meta.data))])
src.seu.Yu.all.10x.E16.5 = AddMetaData(
  src.seu.Yu.all.10x.E16.5, 
  metadata = src.cellres.pan.E16.5.merge@meta.data[
    colnames(src.seu.Yu.all.10x.E16.5),
    setdiff(colnames(src.cellres.pan.E16.5.merge@meta.data),
            colnames(src.seu.Yu.all.10x.E16.5@meta.data))])
src.seu.Yu.all.10x.E17.5 = AddMetaData(
  src.seu.Yu.all.10x.E17.5, 
  metadata = src.cellres.pan.E17.5.merge@meta.data[
    colnames(src.seu.Yu.all.10x.E17.5),
    setdiff(colnames(src.cellres.pan.E17.5.merge@meta.data),
            colnames(src.seu.Yu.all.10x.E17.5@meta.data))])

#>> MNN was used in Figure S3H
type = c("Net_MinPath")
reduction = c("mnn", "pca_integrated")
lineage_list = list.lineage.cellres.pan
for(i.time in c("E15.5","E16.5","E17.5")){
  for( i in c(1:length(type))){
    for( j in c(1:length(reduction))){
      for( k in c(1:length(lineage_list))){
        for( m in c("tsne","fdl","fdl.rotated","umap")){
          
          library("scales")
          typename = paste("tracing.index", type[i], reduction[j], lineage_list[k], sep = "_")
          
          filename = paste("Simulation/distance_png/", i.time, "_", m, "_", type[i],"_", reduction[j],"_",
                           "Integrated_",typename,"_",".png", sep="")
          srcname = paste("src.seu.Yu.all.10x.", i.time, sep="")
          print(filename)
          seurat = get(srcname)
          
          metadata = cbind(seurat@meta.data, seurat@reductions[[m]]@cell.embeddings[,c(1,2)])
          colnames(metadata) = c(colnames(seurat@meta.data), 'Coord_1', "Coord_2")
          metadata[["tracing.index"]] = metadata[[typename]]
          p_temp = ggplot() +
            geom_point(metadata[metadata[["tracing.index"]]%in%0,],
                       mapping = aes(x=Coord_1, y=Coord_2, color=tracing.index),  size=3) +
            geom_point(metadata[metadata[["tracing.index"]]>0,],
                       mapping = aes(x=Coord_1, y=Coord_2, color=tracing.index),  size=3) +
            scale_colour_gradient(low = "gray", high ="#C70039")+
            theme_void() + 
            theme(legend.position = "none")
          
          png(filename = filename, width = 1000, height = 1000, pointsize = 20)
          print(p_temp)
          dev.off() 
        }
      }
    }
  }
}

#================================================================







