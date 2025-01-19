#===============================================================================
#>> 3.Integration of 10x with Smart-seq2
#===============================================================================

#===========================================================
#>>>> Load data: scRNA-seq of mouse endoderm at 9SS
#===========================================================
#> GSA: CRA003104
#> Li, LC., Wang, X., Xu, ZR. et al. \
#> Single-cell patterning and axis characterization in the murine and human definitive endoderm. \
#> Cell Res 31, 326–344 (2021). \
#> https://doi.org/10.1038/s41422-020-00426-0

src.s9 = UpdateSeuratObject(src.s9)
src.s9$Type = "Sm2"

src.s9$cluster = sapply(src.s9$cluster, as.character)
src.s9$cluster[grep("Foregut.AIP",src.s9$cluster)] = 
  gsub("Foregut.AIP","FG", src.s9$cluster[grep("Foregut.AIP",src.s9$cluster)])
src.s9$cluster[grep("Foregut.Lip",src.s9$cluster)] = 
  gsub("Foregut.Lip","AL", src.s9$cluster[grep("Foregut.Lip",src.s9$cluster)])
src.s9$cluster[grep("Midgut",src.s9$cluster)] = 
  gsub("Midgut","MG", src.s9$cluster[grep("Midgut",src.s9$cluster)])
src.s9$cluster[grep("Hindgut",src.s9$cluster)] = 
  gsub("Hindgut","HG", src.s9$cluster[grep("Hindgut",src.s9$cluster)])

src.s9$organ = sapply(src.s9$organ, as.character)
src.s9$organ[grep("Foregut.AIP",src.s9$organ)] = 
  gsub("Foregut.AIP","FG", src.s9$organ[grep("Foregut.AIP",src.s9$organ)])
src.s9$organ[grep("Foregut.Lip",src.s9$organ)] = 
  gsub("Foregut.Lip","AL", src.s9$organ[grep("Foregut.Lip",src.s9$organ)])
src.s9$organ[grep("Midgut",src.s9$organ)] = 
  gsub("Midgut","MG", src.s9$organ[grep("Midgut",src.s9$organ)])
src.s9$organ[grep("Hindgut",src.s9$organ)] = 
  gsub("Hindgut","HG", src.s9$organ[grep("Hindgut",src.s9$organ)])

src.s9$cluster.extract.v1.1...re = src.s9$cluster
src.s9$cluster.extract.v1.1 = src.s9$organ
#===========================================================


#===========================================================
#>>>> Load src.9ss.integrated: 
#===========================================================
#>>> Load endoderm-origin in supplymental-table.S1-sheet3-10x-v3
src.9ss.integrated$cluster.extract.v1.1...re = src.9ss.integrated$cluster.extract.v1.1
src.9ss.integrated@meta.data[grep("FG",src.9ss.integrated$cluster.extract.v1.1),]$cluster.extract.v1.1...re = 'FG'
src.9ss.integrated@meta.data[grep("AL",src.9ss.integrated$cluster.extract.v1.1),]$cluster.extract.v1.1...re = 'AL'
src.9ss.integrated@meta.data[grep("MG",src.9ss.integrated$cluster.extract.v1.1),]$cluster.extract.v1.1...re = 'MG'
src.9ss.integrated@meta.data[grep("HG",src.9ss.integrated$cluster.extract.v1.1),]$cluster.extract.v1.1...re = 'HG'
src.9ss.integrated$Type = "10X"

pdf("FG1_correct/9ss.type.pdf",12,9)
DimPlot(src.9ss.integrated, group.by = "cluster.extract.v1.1",pt.size = 1.1,
        cols = cluster.endoderm.color.v5) +
  theme_classic2() + p_add +
  guides(color = guide_legend(override.aes = list(size = 5)))
DimPlot(src.9ss.integrated, group.by = "cluster.extract.v1.1...re",pt.size = 1.1,
        cols = cluster.endoderm.color.v3) +
  theme_classic2() + p_add +
  guides(color = guide_legend(override.aes = list(size = 5)))
dev.off()
#===========================================================


#===========================================================
#>>>> CCA-10X&Sm2
#===========================================================
src.1 = src.9ss.integrated; src.2 = src.s9
src.1$index = 1; src.2$index = 2

src.merge = RunCCA(src.1, src.2,
                   features = rownames(src.9ss.integrated@reductions$pca@feature.loadings),
                   num.cc = 30) # reduction: CCA

cca.merge.9ss.integrated.region = 
  by(as.data.frame(src.merge@reductions$cca@cell.embeddings[,1:20]),
     INDICES = paste(src.merge$Type, src.merge$cluster.extract.v1.1...re, sep = "_"),
     FUN = colMeans)
cca.merge.9ss.integrated.region = do.call(cbind, cca.merge.9ss.integrated.region)
cor.cca.merge.9ss.integrated.region = cor(cca.merge.9ss.integrated.region)

pdf("figure.v08.07/10X_Smartseq2_CCA.pdf",7,7)
for(i.plot in c("plot.cca.heatmap")){
  data = cor.cca.merge.9ss.integrated.region
  name_list1 = paste("10X_", c("FG","AL","MG",'HG'), sep = "")
  name_list2 = paste("Sm2_", c("FG","AL","MG",'HG'), sep = "")
  
  # name_list1 = paste("10X_",
  #                    c("FG.1","FG.2","FG.3","FG.4","FG.5","FG.6",
  #                      "AL.1","AL.2",'AL.3',
  #                      "MG.1","MG.2","MG.3",
  #                      "HG.1","HG.2"
  #                    ), sep = "")
  # name_list2 = paste("Sm2_",
  #                    c("FG.1","FG.2","FG.3","FG.4","FG.5",
  #                      "AL.1","AL.2",'AL.3',
  #                      "MG.1","MG.2","MG.3",
  #                      "HG.1","HG.2"
  #                    ), sep = "")
  
  name_list_all = c(name_list1,name_list2)
  
  a = MyHeatmap(as.data.frame(data[name_list_all, name_list_all]),
                type = "raw",
                hc.c.data.type = "raw",
                c.cov.method = "p",
                c.hc.method = "ward.D",
                ColSideColors = cbind(
                  MyName2Col(strsplit2(name_list_all, split = "_")[,2], cluster.endoderm.color.v5)),
                RowSideColors = t(cbind(
                  MyName2Col(strsplit2(name_list_all, split = "_")[,2], cluster.endoderm.color.v5))),
                color.palette = colorRampPalette(c("#5aa7dd","#ebebeb","#ff523f"), 
                                                 space="Lab"),
                #Rowv = "none", Colv = "none", 
                # return.tree = "col",
                ColSideColorsSize = 2, RowSideColorsSize = 2,
                labCol = name_list_all,  labRow = name_list_all,
                margins = c(8,8))
  
  a = MyHeatmap(as.data.frame(data[name_list1, name_list2]),
                type = "raw",
                hc.c.data.type = "raw",
                c.cov.method = "p",
                c.hc.method = "ward.D",
                ColSideColors = cbind(
                  MyName2Col(strsplit2(name_list2, split = "_")[,2], cluster.endoderm.color.v5)),
                RowSideColors = t(cbind(
                  MyName2Col(strsplit2(name_list1, split = "_")[,2], cluster.endoderm.color.v5))),
                color.palette = colorRampPalette(c("#5aa7dd","#ebebeb","#ff523f"), 
                                                 space="Lab"),
                Rowv = "none", Colv = "none", 
                # return.tree = "col",
                ColSideColorsSize = 2, RowSideColorsSize = 2,
                labCol = name_list2,  labRow = name_list1,
                margins = c(8,8))
}
dev.off()

#===========================================================










