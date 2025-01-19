#-------------------------------
#>> 11.3 Ngn3-Cre
#--------------------------------

#>> after pre-process for Smartseq-3
mat_sm3_ngn3cre = src.sm3.update.rawdata_Ngn3Cre@meta.data

pdf("figure.v08.07/ngn3cre_qc.pdf",6,6)
ggplot() +
  geom_vline(xintercept = 6, linetype = "dashed", colour = "black")+
  geom_violin(data = mat_sm3_ngn3cre, 
              mapping = aes(y = "", x = nFeature_RNA/1000), fill = "#50d2f0")+
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=, size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        plot.title = element_blank(),
        aspect.ratio=0.25) +
  xlab("The number of\n detected gene (x 1000)")  + ylab("Tracing\n Lineeag") 

mat_sm3_ngn3cre$cluster = factor(mat_sm3_ngn3cre$cluster,
                                 levels = c("Endoderm","Ectoderm","Yolk sac","Low Quality"))
mat_sm3_ngn3cre[mat_sm3_ngn3cre$cluster%in%NA,]$cluster = "Low Quality"
mat_temp = as.data.frame(table(mat_sm3_ngn3cre$cluster))
mat_temp$Var1 = factor(mat_temp$Var1, levels = c(names(qc.color))) 

pie = ggplot() +
  #geom_vline(xintercept = 2.5, linetype = "dashed", colour = "black")+
  geom_bar(data = mat_temp, 
           mapping = aes(x = "", y = Freq,
                         group = Var1, fill = Var1),
           stat = "identity", position = "stack")+
  coord_polar(theta = "y") +  theme_void() +
  scale_fill_manual(values = qc.color) 
pie1 = pie +
  geom_text(aes(x='', y=mat_temp$order*300,
                label= paste(round(mat_temp$Freq/sum(mat_temp$Freq),3)*100,"%",sep="")))+
  pie;pie1
dev.off()


src.sm3.endoderm_Ngn3Cre = src.sm3.update.rawdata_Ngn3Cre[,src.sm3.update.rawdata_Ngn3Cre$cluster%in%"Endoderm"]
src.sm3.endoderm_Ngn3Cre$lineage  = "Ngn3Cre"

src.endoderm.mg3.merge_Ngn3Cre = merge(
  src.sm3.endoderm_Ngn3Cre,
  src.endoderm.mg3.ext.re)

src.endoderm.mg3.merge_Ngn3Cre$type = 1
src.endoderm.mg3.merge_Ngn3Cre@meta.data[src.endoderm.mg3.merge_Ngn3Cre$lineage%in%NA,]$type = 2
src.endoderm.mg3.merge_Ngn3Cre$batch_type = 
  paste(src.endoderm.mg3.merge_Ngn3Cre$batch,
        src.endoderm.mg3.merge_Ngn3Cre$type, sep = "_")

src.endoderm.mg3.merge_Ngn3Cre = ScaleData(src.endoderm.mg3.merge_Ngn3Cre, 
                                           features = rownames(src.endoderm.mg3.merge_Ngn3Cre),
                                           split.by = "batch_type")
src.endoderm.mg3.merge_Ngn3Cre.selectgene = src.endoderm.mg3.ext.re.selectgene.fin
src.endoderm.mg3.merge_Ngn3Cre = RunPCA(src.endoderm.mg3.merge_Ngn3Cre,
                                        features = src.endoderm.mg3.merge_Ngn3Cre.selectgene)
src.endoderm.mg3.merge_Ngn3Cre = RunUMAP(src.endoderm.mg3.merge_Ngn3Cre,
                                         dims = 1:20,n.neighbors = 50)

src.endoderm.mg3.merge_Ngn3Cre@reductions$umap@cell.embeddings = 
  cbind(src.endoderm.mg3.merge_Ngn3Cre@reductions$umap@cell.embeddings[,2],
        src.endoderm.mg3.merge_Ngn3Cre@reductions$umap@cell.embeddings[,1])
colnames(src.endoderm.mg3.merge_Ngn3Cre@reductions$umap@cell.embeddings) = c("UMAP_1","UMAP_2")

# MNN
#---------------
MNN.res = 
  mnnCorrect(as.matrix(src.endoderm.mg3.merge_Ngn3Cre@assays$RNA@data[src.endoderm.mg3.merge_Ngn3Cre.selectgene, 
                                                                      src.endoderm.mg3.merge_Ngn3Cre$batch_type%in%"1_1"]),
             as.matrix(src.endoderm.mg3.merge_Ngn3Cre@assays$RNA@data[src.endoderm.mg3.merge_Ngn3Cre.selectgene, 
                                                                      src.endoderm.mg3.merge_Ngn3Cre$batch_type%in%"2_1"]),
             as.matrix(src.endoderm.mg3.merge_Ngn3Cre@assays$RNA@data[src.endoderm.mg3.merge_Ngn3Cre.selectgene, 
                                                                      src.endoderm.mg3.merge_Ngn3Cre$batch_type%in%"1_2"]),
             as.matrix(src.endoderm.mg3.merge_Ngn3Cre@assays$RNA@data[src.endoderm.mg3.merge_Ngn3Cre.selectgene, 
                                                                      src.endoderm.mg3.merge_Ngn3Cre$batch_type%in%"2_2"]),
             k = 5,cos.norm.out=F)
src.endoderm.mg3.merge_Ngn3Cre@assays$mnnRNA = CreateAssayObject(data = MNN.res@assays@data$corrected)
src.endoderm.mg3.merge_Ngn3Cre@assays$mnnRNA@key = "mnn_"
src.endoderm.mg3.merge_Ngn3Cre = ScaleData(src.endoderm.mg3.merge_Ngn3Cre, 
                                           rownames(src.endoderm.mg3.merge_Ngn3Cre@assays$mnnRNA),
                                           assay = "mnnRNA")

src.endoderm.mg3.merge_Ngn3Cre@meta.data = src.endoderm.mg3.merge_Ngn3Cre@meta.data[!src.endoderm.mg3.merge_Ngn3Cre$Time%in%NA,]

anchor.mg3.merge_Ngn3Cre = 
  FindTransferAnchors(reference = src.endoderm.mg3.merge_Ngn3Cre[,src.endoderm.mg3.merge_Ngn3Cre$type%in%2],
                      query = src.endoderm.mg3.merge_Ngn3Cre[,src.endoderm.mg3.merge_Ngn3Cre$type%in%1],
                      reference.assay = "mnnRNA",
                      query.assay = "mnnRNA",scale = T,
                      features = src.endoderm.mg3.merge_Ngn3Cre.selectgene)
umap.transfer.mg3.merge_Ngn3Cre = TransferData(anchor.mg3.merge_Ngn3Cre,
                                               t(src.endoderm.mg3.ext.re@reductions$umap@cell.embeddings))

src.endoderm.mg3.merge_Ngn3Cre@reductions$umap_fta = src.endoderm.mg3.merge_Ngn3Cre@reductions$umap
src.endoderm.mg3.merge_Ngn3Cre@reductions$umap_fta@cell.embeddings[colnames(src.endoderm.mg3.ext.re),] = 
  src.endoderm.mg3.ext.re@reductions$umap@cell.embeddings
src.endoderm.mg3.merge_Ngn3Cre@reductions$umap_fta@cell.embeddings[colnames(src.sm3.endoderm_Ngn3Cre),] = 
  as.matrix(t(umap.transfer.mg3.merge_Ngn3Cre@data))
src.endoderm.mg3.merge_Ngn3Cre@reductions$umap_fta@key = "Coord_"
colnames(src.endoderm.mg3.merge_Ngn3Cre@reductions$umap_fta@cell.embeddings) = c("Coord_1","Coord_2")
DimPlot(src.endoderm.mg3.merge_Ngn3Cre, reduction = "umap_fta",
        group.by = "lineage", cols = color.lineage)

a = FNN::knn(
  src.endoderm.mg3.merge_Ngn3Cre@reductions$umap_fta@cell.embeddings[
    rownames(src.endoderm.mg3.merge_Ngn3Cre@meta.data[src.endoderm.mg3.merge_Ngn3Cre$type%in%2,]),],
  src.endoderm.mg3.merge_Ngn3Cre@reductions$umap_fta@cell.embeddings[
    rownames(src.endoderm.mg3.merge_Ngn3Cre@meta.data[src.endoderm.mg3.merge_Ngn3Cre$type%in%1,]),],
  src.endoderm.mg3.merge_Ngn3Cre@meta.data[src.endoderm.mg3.merge_Ngn3Cre$type%in%2,]$cluster.v06.26.re, k = 10)
src.endoderm.mg3.merge_Ngn3Cre$cluster.predict.umap_int.ext.v1.1 = src.endoderm.mg3.merge_Ngn3Cre$cluster.v06.26.re
src.endoderm.mg3.merge_Ngn3Cre@meta.data[src.endoderm.mg3.merge_Ngn3Cre$type%in%1,]$cluster.predict.umap_int.ext.v1.1 = as.character(a)
DimPlot(src.endoderm.mg3.merge_Ngn3Cre, group.by = "cluster.predict.umap_int.ext.v1.1")


p = ggplot()+
  geom_point(data = 
               cbind(src.endoderm.mg3.ext.re@meta.data,
                     src.endoderm.mg3.ext.re@reductions$umap@cell.embeddings),
             mapping = aes(x=UMAP_1, y=UMAP_2), colour="#eeeeee") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_classic2() + p_add
p0 = p + geom_point(data = 
                      cbind(src.endoderm.mg3.merge_Ngn3Cre@meta.data[colnames(src.sm3.endoderm_Ngn3Cre),],
                            src.endoderm.mg3.merge_Ngn3Cre@reductions$umap_fta@cell.embeddings[colnames(src.sm3.endoderm_Ngn3Cre),]),
                    mapping = aes(x=Coord_1, y=Coord_2, color = lineage)) +
  scale_color_manual(values = color.lineage)
p1 = p + geom_point(data = 
                      cbind(src.endoderm.mg3.merge_Ngn3Cre@meta.data[colnames(src.sm3.endoderm_Ngn3Cre),],
                            src.endoderm.mg3.merge_Ngn3Cre@reductions$umap_fta@cell.embeddings[colnames(src.sm3.endoderm_Ngn3Cre),]),
                    mapping = aes(x=Coord_1, y=Coord_2, color = Time)) +
  scale_color_manual(values = colors.time.2)
p2 = p + geom_point(data = 
                      cbind(src.endoderm.mg3.merge_Ngn3Cre@meta.data[colnames(src.sm3.endoderm_Ngn3Cre),],
                            src.endoderm.mg3.merge_Ngn3Cre@reductions$umap_fta@cell.embeddings[colnames(src.sm3.endoderm_Ngn3Cre),]),
                    mapping = aes(x=Coord_1, y=Coord_2, color = cluster.predict.umap_int.ext.v1.1)) +
  scale_color_manual(values = cluster.endoderm.color.v5)

mat = src.endoderm.mg3.merge_Ngn3Cre@meta.data[colnames(src.sm3.endoderm_Ngn3Cre),]
mat$Time = factor(mat$Time,
                  levels = c("27ss","24ss","21ss","18ss"))
time_list = table(src.endoderm.mg3.merge_Ngn3Cre@meta.data[colnames(src.sm3.endoderm_Ngn3Cre),"Time"])
p3 = ggplot() +
  geom_bar(data = mat, 
           mapping = aes(y = Time,
                         group = cluster.predict.umap_int.ext.v1.1, 
                         fill = cluster.predict.umap_int.ext.v1.1),
           stat = "count")+ # position = 'fill')+
  annotate("text", x = rev(time_list)+5, y = c(1:4),
           label = paste("(",rev(time_list),")",sep="")) +
  scale_fill_manual(values = cluster.endoderm.color.v5)+
  theme_classic() +
  scale_x_continuous(position = "top")+
  #scale_y_continuous(position = "right")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        plot.title = element_blank(),
        aspect.ratio=0.3) +
  ylab("Time") + xlab("Cell Type")

pdf("figure.v08.07/ngn3cre_fta.pdf",9,7)
p0;p1;p2;p3;
dev.off()


pdf("figure.v08.07/mg3.marler.ngn3.pdf",5,5)
for(i in c(c("Cdx2","Mnx1","Sox2","Osr1","Pdx1","Pitx2","Osr2","Dpp4",
             "2610528A11Rik","Neurod1","Neurog3","Gcg",'Nkx6-1'))){
  gene = i
  print(gene)
  gene_frame = t(src.endoderm.mg3.ext.re@assays$RNA[gene,cell_permission_mg3])
  colnames(gene_frame) = "gene"
  print(
    ggplot() +
      geom_point(data = cbind(src.endoderm.mg3.ext.re@meta.data[cell_permission_mg3,], gene_frame,
                              src.endoderm.mg3.ext.re@reductions$umap@cell.embeddings[cell_permission_mg3,]),
                 #fdl.src.endoderm.mg3.ext.re[cell_permission_mg3,]),
                 mapping = aes(x = UMAP_1, y = UMAP_2,
                               color = gene))+
      theme_void() + theme(legend.position = "none") +
      ggtitle(label = gene)+
      scale_color_viridis() + p_add)
}
dev.off()


mat_ngn3 = src.endoderm.mg3.ext.re@meta.data[src.endoderm.mg3.ext.re@assays$RNA@data["Neurog3",]>0.5,]
mat_ngn3$Time = factor(mat_ngn3$Time,
                       rev(c("ss12","ss15","ss18","ss21","ss24","ss27")))
mat_ngn3$cluster.v06.26.re = factor(
  mat_ngn3$cluster.v06.26.re,
  rev(c("MG.3","MG.3.A","MG.3.M","MG.3.P","Stomach","DP","Pancreas","EP.1","EP.2","Small.intestine.1","Small.intestine.2")))


pdf("figure.v08.07/ngn3_mg3_mat_sum.pdf",7,7)
ggplot() +
  geom_bar(data = mat_ngn3, 
           mapping = aes(y = Time,
                         group = cluster.v06.26.re, 
                         fill = cluster.v06.26.re),
           stat = "count")+ # position = 'fill')+
  annotate("text", x = table(mat_ngn3$Time)+5, y = c(1:6),
           label = paste("(",table(mat_ngn3$Time),")",sep="")) +
  scale_fill_manual(values = cluster.endoderm.color.v5)+
  theme_classic() +
  scale_x_continuous(position = "top")+
  #scale_y_continuous(position = "right")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        plot.title = element_blank(),
        aspect.ratio=0.3) +
  ylab("Time") + xlab("Cell Type")
dev.off()


