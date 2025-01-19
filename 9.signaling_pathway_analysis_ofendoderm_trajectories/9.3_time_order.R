#===============================================================================
#>>> 9.signaling pathways analysis 
#===============================================================================

#===============================================================================
#>>> 9.3 Time trajectory
#===============================================================================

#==========================================================
#   pathway.endoderm.time.list.add.refine :: 9SS 15SS  (v1)
#==========================================================
celllist.endoderm.time.list = list()
for(i.time in c('ss9','ss12','ss15','ss18','ss21','ss24','ss27')){
  celllist.endoderm.time.list[[i.time]] = list()
}

#-- subset cell:  ------------------------------------------
#-- ss9
celllist.endoderm.time.list$ss9$merge = colnames(src.9ss.integrated)
celllist.endoderm.time.list$ss9$FG = rownames(src.9ss.integrated@meta.data[
  src.9ss.integrated$cluster.v06.26.re..correct..un%in%c("FG.1","FG.2","FG.3","FG.4","FG.5","FG.6"),])
celllist.endoderm.time.list$ss9$AL = rownames(src.9ss.integrated@meta.data[
  src.9ss.integrated$cluster.v06.26.re..correct..un%in%c("AL.1","AL.2","AL.3"),])
celllist.endoderm.time.list$ss9$MG = rownames(src.9ss.integrated@meta.data[
  src.9ss.integrated$cluster.v06.26.re..correct..un%in%c("MG.1","MG.2","MG.3"),])
celllist.endoderm.time.list$ss9$HG = rownames(src.9ss.integrated@meta.data[
  src.9ss.integrated$cluster.v06.26.re..correct..un%in%c("HG.1","HG.2"),])
celllist.endoderm.time.list$ss9$Line1 = rownames(src.9ss.integrated@meta.data[
  src.9ss.integrated$cluster.v06.26.re..correct..un%in%c("FG.1","FG.6",'MG.3',"HG.2"),])
celllist.endoderm.time.list$ss9$Line2 = rownames(src.9ss.integrated@meta.data[
  src.9ss.integrated$cluster.v06.26.re..correct..un%in%c("FG.2","FG.3","MG.1","MG.2","HG.1"),])
celllist.endoderm.time.list$ss9$Line3 = rownames(src.9ss.integrated@meta.data[
  src.9ss.integrated$cluster.v06.26.re..correct..un%in%c("FG.5","FG.4","AL.1","AL.2","AL.3"),])

#-- ss15
pdf("try.pdf", 14, 14)
DimPlot(src.15ss.integrated, reduction = 'umap.rotated', pt.size = 1.5,
        group.by = "cluster.v06.26.re..correct..un", label = T, label.color = "red", label.size = 5,
        cols = cluster.endoderm.color.v5) +
  theme(legend.position = "bottom",
        aspect.ratio = 1)
dev.off()

src.15ss.integrated$cluster.v06.26.re..correct..un =
  src.15ss.integrated$cluster.v06.26.re..correct

celllist.endoderm.time.list$ss15$merge = colnames(src.15ss.integrated)
celllist.endoderm.time.list$ss15$Line1 = rownames(src.15ss.integrated@meta.data[
  src.15ss.integrated$cluster.v06.26.re..correct..un%in%c("Pharynx.organ.2","FG.1","FG.6","Esophagus","MG.3.A","DP","EP.1","MG.3.M","MG.3","MG.3.P"),])
celllist.endoderm.time.list$ss15$Line2 = rownames(src.15ss.integrated@meta.data[
  src.15ss.integrated$cluster.v06.26.re..correct..un%in%c('FG.2',"Pharynx.organ.1","FG.3","Pharynx.organ.4","Stomach","MG.1",
                                                          "Small.intestine.1","Small.intestine.2","MG.2","HG.1","Large.intestine.1",
                                                          "HG.1-Large.intestine.2","Large.intestine.2","HG.2","Large.intestine.3"),])
celllist.endoderm.time.list$ss15$Line3 = rownames(src.15ss.integrated@meta.data[
  src.15ss.integrated$cluster.v06.26.re..correct..un%in%c("FG.5","Thyroid","Pharynx.organ.3","Pharynx.organ.5","Lung","FG.4","FG.4-Lung/Stomach",
                                                          "FG.4-Liver","AL.1/2-Liver","Liver","EHBD","AL.3-Liver","VP","AL.3-EHBD/VP",'AL.3',"AL.3-Small.intestine.1"),])
#-------------------------------

seurat.plsda.pathway.endoderm.time.list = list()
for(i.time in c("ss9","ss15")){
  seurat.plsda.pathway.endoderm.time.list[[i.time]] = list()
  for(j.type in names(celllist.endoderm.time.list[[i.time]])){
    seurat.plsda.pathway.endoderm.time.list[[i.time]][[j.type]] =
      seurat_pathway_pre(pathway.matrix = pathway.endoderm.time.list[[paste(gsub("ss","",i.time), "ss", sep="")]][,celllist.endoderm.time.list[[i.time]][[j.type]]], 
                         seurat.ref = get(paste("src.", gsub("ss","",i.time), "ss.integrated", sep = ""))[,celllist.endoderm.time.list[[i.time]][[j.type]]]) 
  }
}

pdf('try.pdf',8,8)
for(i.time in names(seurat.plsda.pathway.endoderm.time.list)){
  # if(i.time == "ss15"){}else{next()}
  for(i.type in names(seurat.plsda.pathway.endoderm.time.list[[i.time]])){
    seurat.plsda.pathway.endoderm.time.list[[i.time]][[i.type]] = 
      seurat_plsda_pathway(seurat.plsda.pathway.endoderm.time.list[[i.time]][[i.type]],
                           seurat.plsda.pathway.endoderm.time.list[[i.time]][[i.type]]$cluster.v06.26.re..correct..un)[["seurat"]]
  }
}
dev.off()

pathway.plsda.pathway.endoderm.time.list = list()
for(i.type in names(seurat.plsda.pathway.endoderm.time.list)){
  # if(i.type == "ss15"){}else{next()}
  for(i.name in names(seurat.plsda.pathway.endoderm.time.list[[i.type]])){
    pathway.plsda.pathway.endoderm.time.list[[i.type]][[i.name]] =
      plsda_select(seurat.plsda.pathway.endoderm.time.list[[i.type]][[i.name]], comp = 3, threshold = 30)
  }
}

save(pathway.endoderm.time.list,
     pathway.plsda.pathway.endoderm.time.list, 
     file = "Milestones/pathway.endoderm.time.list.summary.Rdata")
save(seurat.plsda.pathway.endoderm.time.list, 
     file = "Milestones/seurat.plsda.pathway.endoderm.time.list.Rdata")

tree.plsda.pathway.endoderm.time.list = list()
for(i.name in names(seurat.plsda.pathway.endoderm.time.list)){
  tree.plsda.pathway.endoderm.time.list[[i.name]] = list()
  pdf(paste("Milestones/try.", i.name, '.pdf', sep = ""), 10, 10)
  for(j.name in names(seurat.plsda.pathway.endoderm.time.list[[i.name]])){
    
    seurat = seurat.plsda.pathway.endoderm.time.list[[i.name]][[j.name]]
    seurat$cluster.temp = seurat$cluster.v06.26.re..correct..un
    gene_list = pathway.plsda.pathway.endoderm.time.list[[i.name]][[j.name]]
    
    cellorder_list = get(paste("list_cluster_order.",gsub("ss","",i.name),"ss",sep=''))
    cell_name = c()
    for(i.cluster in intersect(cellorder_list, unique(seurat$cluster.temp))){
      cell_name = c(cell_name,
                    sample(rownames(seurat@meta.data[seurat$cluster.temp%in%i.cluster,]),
                           size = length(rownames(seurat@meta.data[seurat$cluster.temp%in%i.cluster,])),
                           replace = F))
    }
    
    seurat = ScaleData(seurat[,cell_name], features = rownames(seurat))
    
    data = seurat@assays$Pathway@scale.data[gene_list, cell_name]
    data_re = t(apply(data,1,kernelsmooth))
    rownames(data_re) = rownames(data); colnames(data_re) = colnames(data)
    
    plot(c(1:10), c(1:10))
    text(x=5,y=5,labels=paste(i.type, i.name, j.name, sep = " "), cex = 5)
    
    pathway.heatmap  =
      MyHeatmap(data_re,
                type = "raw",
                ColSideColors = cbind(
                  # MyName2Col(seurat@meta.data[cell_name,]$Time, c(colors.time, colors.time.2)),
                  MyName2Col(seurat@meta.data[cell_name,]$cluster.temp, cluster.endoderm.color.v5)
                ),
                # RowSideColors = t(cbind(
                #   MyName2Col(names(gene_list), colors.geneset)
                # )),
                color.palette = colorRampPalette(c("#5aa7dd","#EBEBEB","red"), space="Lab"),
                ColSideColorsSize = 3,
                RowSideColorsSize = 2,
                #Rowv = "none", 
                Colv = "none",
                return.tree = "row",
                labCol="none", graph = T)
    
    tree.plsda.pathway.endoderm.time.list[[i.name]][[j.name]] = as.dendrogram(pathway.heatmap) 
  }
  dev.off()
}
save(tree.plsda.pathway.endoderm.time.list,
     file = "Milestones/tree.plsda.pathway.endoderm.time.list.Rdata")


princurve.withorder = function(i.coord, j.class, j.classorder){
  j.coord.list = list()
  for (i in 1:length(j.classorder)) {
    j.coord.list[[i]] = i.coord[which(j.class == j.classorder[i]),]
  }
  j.gravity_center.list = lapply(j.coord.list, colMeans)
  j.gravity_center = do.call("rbind", j.gravity_center.list)
  
  princurve.j = princurve::principal_curve(
    x = i.coord,
    start = j.gravity_center,
    # maxit=1, stretch = 1, spar = 1.25,
    smoother = "smooth.spline")
  return(princurve.j)
}

princurve.endoderm.time.instances = list()
for(i.time in c('ss9','ss12','ss15','ss18','ss21','ss24','ss27')){
  princurve.endoderm.time.instances[[i.time]] = list()
}

#-- ss9 -----
princurve.endoderm.time.instances$ss9$Line1 = princurve::principal_curve(
  x = cbind(src.9ss.integrated@reductions$umap@cell.embeddings[,1],
            src.9ss.integrated@reductions$umap@cell.embeddings[,2])[
              celllist.endoderm.time.list$ss9$Line1,],
  smoother = "smooth.spline")
princurve.endoderm.time.instances$ss9$Line2 = princurve::principal_curve(
  x = cbind(src.9ss.integrated@reductions$umap@cell.embeddings[,1],
            src.9ss.integrated@reductions$umap@cell.embeddings[,2])[
              intersect(celllist.endoderm.time.list$ss9$Line2,
                        rownames(src.9ss.integrated@meta.data[src.9ss.integrated$cluster.v06.26.re..correct..un%in%c(
                          "FG.2","FG.3","MG.1","MG.2","HG.1"),])),],
  smoother = "smooth.spline")
princurve.endoderm.time.instances$ss9$Line3 = princurve::principal_curve(
  x = cbind(src.9ss.integrated@reductions$umap@cell.embeddings[,1],
            src.9ss.integrated@reductions$umap@cell.embeddings[,2])[
              celllist.endoderm.time.list$ss9$Line3,],
  smoother = "smooth.spline")

for(i.type in c("Line1", "Line2", "Line3")){
  celllist.endoderm.time.list$ss9[[paste(i.type,".order",sep = "")]] = rownames(
    princurve.endoderm.time.instances$ss9[[paste(i.type,sep = "")]]$s[
      order(princurve.endoderm.time.instances$ss9[[paste(i.type,sep = "")]]$lambda),])
}

celllist.endoderm.time.list$ss9$Line3.order = rev(celllist.endoderm.time.list$ss9$Line3.order)


plot(x = rep(1:length(celllist.endoderm.time.list$ss9[[paste("Line3",".order",sep = "")]])),
     y = rep(1:length(celllist.endoderm.time.list$ss9[[paste("Line3",".order",sep = "")]])),
     col = cluster.endoderm.color.v5[src.9ss.integrated@meta.data[
       celllist.endoderm.time.list$ss9[[paste("Line3",".order",sep = "")]],]$cluster.v06.26.re..correct..un])
#--------

#-- ss15 -----
princurve.endoderm.time.instances$ss15$Line1 = princurve::principal_curve(
  x = cbind(src.15ss.integrated@reductions$umap.rotated@cell.embeddings[,1],
            src.15ss.integrated@reductions$umap.rotated@cell.embeddings[,2])[
              intersect(celllist.endoderm.time.list$ss15$Line1,
                        rownames(src.15ss.integrated@meta.data[
                          src.15ss.integrated@reductions$umap.rotated@cell.embeddings[,2]>(-5),])),],
  smoother = "smooth.spline")
princurve.endoderm.time.instances$ss15$Line2 = princurve::principal_curve(
  x = cbind(src.15ss.integrated@reductions$umap.rotated@cell.embeddings[,1],
            src.15ss.integrated@reductions$umap.rotated@cell.embeddings[,2])[
              celllist.endoderm.time.list$ss15$Line2,],
  smoother = "smooth.spline")
princurve.endoderm.time.instances$ss15$Line3 = princurve::principal_curve(
  x = cbind(src.15ss.integrated@reductions$umap.rotated@cell.embeddings[,1],
            src.15ss.integrated@reductions$umap.rotated@cell.embeddings[,2])[
              celllist.endoderm.time.list$ss15$Line3,],
  smoother = "smooth.spline")

for(i.type in c("Line1", "Line2", "Line3")){
  celllist.endoderm.time.list$ss15[[paste(i.type,".order",sep = "")]] = rownames(
    princurve.endoderm.time.instances$ss15[[paste(i.type,sep = "")]]$s[
      order(princurve.endoderm.time.instances$ss15[[paste(i.type,sep = "")]]$lambda),])
}
celllist.endoderm.time.list$ss15$Line3.order = rev(celllist.endoderm.time.list$ss15$Line3.order)

plot(x = rep(1:length(celllist.endoderm.time.list$ss15[[paste("Line3",".order",sep = "")]])),
     y = rep(1:length(celllist.endoderm.time.list$ss15[[paste("Line3",".order",sep = "")]])),
     col = cluster.endoderm.color.v5[src.15ss.integrated@meta.data[
       celllist.endoderm.time.list$ss15[[paste("Line3",".order",sep = "")]],]$cluster.v06.26.re..correct..un])
#--------

save(princurve.endoderm.time.instances, file = "Milestones/princurve.endoderm.time.instances.Rdata")
save(celllist.endoderm.time.list, file = "Milestones/celllist.endoderm.time.list.Rdata")

tree.plsda.pathway.endoderm.time.list = list()
for(i.name in names(seurat.plsda.pathway.endoderm.time.list)){
  tree.plsda.pathway.endoderm.time.list[[i.name]] = list()
  pdf(paste("Milestones/try.", i.name, '.pdf', sep = ""), 10, 10)
  for(j.name in names(seurat.plsda.pathway.endoderm.time.list[[i.name]])){
    
    seurat = seurat.plsda.pathway.endoderm.time.list[[i.name]][[j.name]]
    seurat$cluster.temp = seurat$cluster.v06.26.re..correct..un
    gene_list = pathway.plsda.pathway.endoderm.time.list[[i.name]][[j.name]]
    
    cellorder_list = get(paste("list_cluster_order.",gsub("ss","",i.name),"ss",sep=''))
    cell_name = c()
    for(i.cluster in intersect(cellorder_list, unique(seurat$cluster.temp))){
      cell_name = c(cell_name,
                    sample(rownames(seurat@meta.data[seurat$cluster.temp%in%i.cluster,]),
                           size = length(rownames(seurat@meta.data[seurat$cluster.temp%in%i.cluster,])),
                           replace = F))
    }
    
    seurat = ScaleData(seurat[,cell_name], features = rownames(seurat))  
    data = seurat@assays$Pathway@scale.data[gene_list, cell_name]
    data_re = t(apply(data,1,kernelsmooth))
    rownames(data_re) = rownames(data); colnames(data_re) = colnames(data)
    
    plot(c(1:10), c(1:10))
    text(x=5,y=5,labels=paste(i.type, i.name, j.name, sep = " "), cex = 5)
    
    pathway.heatmap  =
      MyHeatmap(data_re,
                type = "raw",
                ColSideColors = cbind(
                  # MyName2Col(seurat@meta.data[cell_name,]$Time, c(colors.time, colors.time.2)),
                  MyName2Col(seurat@meta.data[cell_name,]$cluster.temp, cluster.endoderm.color.v5)
                ),
                # RowSideColors = t(cbind(
                #   MyName2Col(names(gene_list), colors.geneset)
                # )),
                color.palette = colorRampPalette(c("#5aa7dd","#EBEBEB","red"), space="Lab"),
                ColSideColorsSize = 3,
                RowSideColorsSize = 2,
                # Rowv = "none", 
                Colv = "none",
                return.tree = "row",
                labCol="none", graph = T)
    
    tree.plsda.pathway.endoderm.time.list[[i.name]][[j.name]] = as.dendrogram(pathway.heatmap) 
    
    if((i.name == "ss27" & j.name %in% c("MG", "FG.cluster1", "FG.cluster2", "FG.cluster3",'AL.cluster1'))|
       (i.name == "ss9" & j.name %in% c("Line1", "Line2", "Line3"))|
       (i.name == "ss15" & j.name %in% c("Line1", "Line2", "Line3"))){
      
      cell_name = celllist.endoderm.time.list[[i.name]][[paste(j.name,".order",sep = "")]]
      pathway.heatmap  =
        MyHeatmap(data_re[c(labels(as.dendrogram(pathway.heatmap)[[1]]),
                            labels(as.dendrogram(pathway.heatmap)[[2]])), cell_name],
                  type = "raw",
                  ColSideColors = cbind(
                    # MyName2Col(seurat@meta.data[cell_name,]$Time, c(colors.time, colors.time.2)),
                    MyName2Col(seurat@meta.data[cell_name,]$cluster.temp, cluster.endoderm.color.v5)
                  ),
                  # RowSideColors = t(cbind(
                  #   MyName2Col(names(gene_list), colors.geneset)
                  # )),
                  color.palette = colorRampPalette(c("#5aa7dd","#EBEBEB","red"), space="Lab"),
                  ColSideColorsSize = 2,
                  RowSideColorsSize = 2,
                  Rowv = "none", 
                  Colv = "none",
                  # return.tree = "row",
                  labCol="none", graph = T)
    }
  }
  dev.off()
}
save(tree.plsda.pathway.endoderm.time.list,
     file = "Milestones/tree.plsda.pathway.endoderm.time.list.Rdata")


gene.tree.plsda.pathway.endoderm.time.list = list()
for(i.type in names(tree.plsda.pathway.endoderm.time.list)){
  gene.tree.plsda.pathway.endoderm.time.list[[i.type]] = list()
  for(i.name in names(tree.plsda.pathway.endoderm.time.list[[i.type]])){
    # gene.tree.plsda.pathway.endoderm.time.list[[i.type]][[i.name]] = c(NA)
  }
}

#-- Set endoderm pathways
#-------------------------------------------------------------------------------
#-- ss9 Line1
gene.tree.plsda.pathway.endoderm.time.list$ss9$Line1 = c(
  labels(tree.plsda.pathway.endoderm.time.list$ss9$Line1[[2]][[1]][[1]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss9$Line1[[2]][[1]][[2]][[1]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss9$Line1[[1]][[2]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss9$Line1[[1]][[1]][[2]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss9$Line1[[1]][[1]][[1]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss9$Line1[[2]][[2]][[1]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss9$Line1[[2]][[2]][[2]]))
names(gene.tree.plsda.pathway.endoderm.time.list$ss9$Line1) = c(
  rep(4, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss9$Line1[[2]][[1]][[1]])))),
  rep(3, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss9$Line1[[2]][[1]][[2]][[1]])))),
  rep(2, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss9$Line1[[1]][[2]]),
                  labels(tree.plsda.pathway.endoderm.time.list$ss9$Line1[[1]][[1]][[2]]),
                  labels(tree.plsda.pathway.endoderm.time.list$ss9$Line1[[1]][[1]][[1]])))),
  rep(5, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss9$Line1[[2]][[2]][[1]])))),
  rep(7, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss9$Line1[[2]][[2]][[2]])))))
#-- ss9 Line2
gene.tree.plsda.pathway.endoderm.time.list$ss9$Line2 = c(
  labels(tree.plsda.pathway.endoderm.time.list$ss9$Line2[[2]][[2]][[2]][[2]][[1]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss9$Line2[[2]][[1]][[1]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss9$Line2[[1]][[2]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss9$Line2[[1]][[1]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss9$Line2[[2]][[2]][[1]]),
  rev(labels(tree.plsda.pathway.endoderm.time.list$ss9$Line2[[2]][[2]][[2]][[2]][[2]])),
  rev(labels(tree.plsda.pathway.endoderm.time.list$ss9$Line2[[2]][[2]][[2]][[1]])),
  rev(labels(tree.plsda.pathway.endoderm.time.list$ss9$Line2[[2]][[1]][[2]])))
names(gene.tree.plsda.pathway.endoderm.time.list$ss9$Line2) = c(
  rep(4, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss9$Line2[[2]][[2]][[2]][[2]][[1]]),
                  labels(tree.plsda.pathway.endoderm.time.list$ss9$Line2[[2]][[1]][[1]])))),
  rep(3, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss9$Line2[[1]][[2]])))),
  rep(2, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss9$Line2[[1]][[1]])))),
  rep(5, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss9$Line2[[2]][[2]][[1]]),
                  labels(tree.plsda.pathway.endoderm.time.list$ss9$Line2[[2]][[2]][[2]][[2]][[2]])))),
  rep(11, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss9$Line2[[2]][[2]][[2]][[1]])))),
  rep(7, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss9$Line2[[2]][[1]][[2]])))))
#-- ss9 Line3
gene.tree.plsda.pathway.endoderm.time.list$ss9$Line3 = c(
  labels(tree.plsda.pathway.endoderm.time.list$ss9$Line3[[1]][[1]][[1]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss9$Line3[[1]][[1]][[2]]),
  rev(labels(tree.plsda.pathway.endoderm.time.list$ss9$Line3[[2]][[1]][[2]][[1]])),
  rev(labels(tree.plsda.pathway.endoderm.time.list$ss9$Line3[[2]][[1]][[1]])),
  # labels(tree.plsda.pathway.endoderm.time.list$ss9$Line3[[2]][[2]][[1]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss9$Line3[[1]][[2]][[2]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss9$Line3[[2]][[2]][[2]][[2]]))
names(gene.tree.plsda.pathway.endoderm.time.list$ss9$Line3) = c(
  rep(4, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss9$Line3[[1]][[1]][[1]])))),
  rep(3, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss9$Line3[[1]][[1]][[2]])))),
  rep(2, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss9$Line3[[2]][[1]][[2]][[1]])))),
  rep(5, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss9$Line3[[2]][[1]][[1]])))),
  rep(7, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss9$Line3[[1]][[2]][[2]]),
                  labels(tree.plsda.pathway.endoderm.time.list$ss9$Line3[[2]][[2]][[2]][[2]])))))

#-- ss15 Line1
gene.tree.plsda.pathway.endoderm.time.list$ss15$Line1 = c(
  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line1[[1]][[2]][[1]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line1[[1]][[2]][[2]][[1]][[1]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line1[[1]][[2]][[2]][[1]][[2]][[2]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line1[[1]][[2]][[2]][[1]][[2]][[1]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line1[[1]][[2]][[2]][[2]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line1[[1]][[1]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line1[[2]][[1]]),
  rev(labels(tree.plsda.pathway.endoderm.time.list$ss15$Line1[[2]][[2]][[2]])),
  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line1[[2]][[2]][[1]]))
names(gene.tree.plsda.pathway.endoderm.time.list$ss15$Line1) = c(
  rep(4, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss15$Line1[[1]][[2]][[1]]),
                  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line1[[1]][[2]][[2]][[1]][[1]]),
                  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line1[[1]][[2]][[2]][[1]][[2]][[2]])))),
  rep(3, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss15$Line1[[1]][[2]][[2]][[1]][[2]][[1]]),
                  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line1[[1]][[2]][[2]][[2]])))),
  rep(2, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss15$Line1[[1]][[1]])))),
  rep(5, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss15$Line1[[2]][[1]])))),
  rep(7, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss15$Line1[[2]][[2]])))))
#-- ss15 Line2


gene.tree.plsda.pathway.endoderm.time.list$ss15$Line2 = c(
  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line2[[1]][[1]][[1]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line2[[1]][[2]][[1]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line2[[1]][[2]][[2]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line2[[2]][[1]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line2[[2]][[2]][[2]][[2]][[1]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line2[[2]][[2]][[2]][[1]]),
  rev(c(labels(tree.plsda.pathway.endoderm.time.list$ss15$Line2[[2]][[2]][[1]]),
        labels(tree.plsda.pathway.endoderm.time.list$ss15$Line2[[2]][[2]][[2]][[2]][[2]]))))
names(gene.tree.plsda.pathway.endoderm.time.list$ss15$Line2) = c(
  rep(4, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss15$Line2[[1]][[1]][[1]]),
                  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line2[[1]][[2]][[1]])))),
  rep(3, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss15$Line2[[1]][[2]][[2]])))),
  rep(2, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss15$Line2[[2]][[1]])))),
  rep(5, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss15$Line2[[2]][[2]][[1]]),
                  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line2[[2]][[2]][[2]][[2]][[2]])))),
  rep(7, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss15$Line2[[2]][[2]][[2]][[2]][[1]]),
                  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line2[[2]][[2]][[2]][[1]])))))


#-- ss15 Line3
gene.tree.plsda.pathway.endoderm.time.list$ss15$Line3 = c(
  rev(labels(tree.plsda.pathway.endoderm.time.list$ss15$Line3[[1]][[1]])),
  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line3[[2]][[2]][[2]][[2]][[2]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line3[[2]][[2]][[2]][[1]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line3[[2]][[2]][[2]][[2]][[1]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line3[[2]][[1]][[2]]),
  # labels(tree.plsda.pathway.endoderm.time.list$ss15$Line3[[2]][[1]][[1]][[1]]),
  # labels(tree.plsda.pathway.endoderm.time.list$ss15$Line3[[1]][[2]][[1]][[2]][[2]]),
  # labels(tree.plsda.pathway.endoderm.time.list$ss15$Line3[[1]][[2]][[1]]),
  # labels(tree.plsda.pathway.endoderm.time.list$ss15$Line3[[1]][[2]][[1]][[2]][[1]]),
  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line3[[2]][[1]][[1]][[2]]))
names(gene.tree.plsda.pathway.endoderm.time.list$ss15$Line3) = c(
  rep(4, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss15$Line3[[1]][[1]])))),
  rep(3, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss15$Line3[[1]][[1]][[2]])))),
  rep(2, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss15$Line3[[2]][[2]][[2]][[1]]),
                  labels(tree.plsda.pathway.endoderm.time.list$ss15$Line3[[2]][[2]][[2]][[2]][[1]])))),
  rep(5, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss15$Line3[[2]][[1]][[2]])))),
  # rep(11, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss15$Line3[[2]][[1]][[1]][[1]]),
  #                 labels(tree.plsda.pathway.endoderm.time.list$ss15$Line3[[1]][[2]][[1]][[2]][[2]]),
  #                 labels(tree.plsda.pathway.endoderm.time.list$ss15$Line3[[1]][[2]][[1]]),
  #                 labels(tree.plsda.pathway.endoderm.time.list$ss15$Line3[[1]][[2]][[1]][[2]][[1]])))),
  rep(7, length(c(labels(tree.plsda.pathway.endoderm.time.list$ss15$Line3[[2]][[1]][[1]][[2]])))))

#-------------------------------------------------------------

save(gene.tree.plsda.pathway.endoderm.time.list, 
     file = "Milestones/gene.tree.plsda.pathway.endoderm.time.list.Rdata")

for(i.name in names(seurat.plsda.pathway.endoderm.time.list)){
  # if(i.name %in% c("ss15")){}else{next()}
  # pdf(paste("Milestones/try.order.", i.name, '.pdf', sep = ""), 10, 10)
  pdf(paste("trajectory_signal/try.", i.name, '.pdf', sep = ""), 8, 8)
  for(j.name in names(seurat.plsda.pathway.endoderm.time.list[[i.name]])){
    
    seurat = seurat.plsda.pathway.endoderm.time.list[[i.name]][[j.name]]
    seurat$cluster.temp = seurat$cluster.v06.26.re..correct..un
    
    if(j.name %in% names(gene.tree.plsda.pathway.endoderm.time.list[[i.name]])){
      gene_list = gene.tree.plsda.pathway.endoderm.time.list[[i.name]][[j.name]]
    }else{next()}
    
    # cellorder_list = get(paste("list_cluster_order.",gsub("ss","",i.name),"ss",sep=''))
    cell_name = celllist.endoderm.time.list[[i.name]][[paste(j.name, ".order", sep = "")]]
    
    
    seurat = ScaleData(seurat[,cell_name], features = rownames(seurat))
    data = seurat@assays$Pathway@scale.data[gene_list, cell_name]
    data_re = t(apply(data,1,kernelsmooth))
    rownames(data_re) = rownames(data); colnames(data_re) = colnames(data)
    
    plot(c(1:10), c(1:10))
    text(x=5,y=5,labels=paste(i.type, i.name, j.name, sep = " "), cex = 5)
    
    pathway.heatmap  =
      MyHeatmap(data_re,
                type = "raw",
                ColSideColors = cbind(
                  MyName2Col(seurat@meta.data[cell_name,]$cluster.temp, cluster.endoderm.color.v5)
                ),
                RowSideColors = t(cbind(
                  MyName2Col(names(gene_list), colors.geneset)
                )),
                color.palette = colorRampPalette(c("#5aa7dd","#EBEBEB","red"), space="Lab"),
                ColSideColorsSize = 2,
                RowSideColorsSize = 2,
                Rowv = "none", 
                Colv = "none",
                # return.tree = "row",
                margins = c(10,10),
                labCol="none", graph = T)
  }
  dev.off()
}



#-- gene.tree.plsda.pathway.endoderm.time.list.filter
gene.tree.plsda.pathway.endoderm.time.list.filter = list()
for(i.type in names(gene.tree.plsda.pathway.endoderm.time.list)){
  gene.tree.plsda.pathway.endoderm.time.list.filter[[i.type]] = list()
  
  for(j.type in names(gene.tree.plsda.pathway.endoderm.time.list[[i.type]])){
    print(paste(i.type,j.type))
    
    type.src = "cluster.v06.26.re..correct..un"
    src.plsda = seurat.plsda.pathway.endoderm.time.list[[i.type]][[j.type]]
    gene.plsda = gene.tree.plsda.pathway.endoderm.time.list[[i.type]][[j.type]]
    
    gene.anova = anova.test(tpm = src.plsda@assays$Pathway@data[gene.plsda, colnames(src.plsda)],
                            variable = as.character(src.plsda@meta.data[,type.src]))
    gene.anova.scale = anova.test(tpm = src.plsda@assays$Pathway@scale.data[gene.plsda, colnames(src.plsda)],
                                  variable = as.character(src.plsda@meta.data[,type.src]))
    
    gene.tree.plsda.pathway.endoderm.time.list.filter[[i.type]][[j.type]][["gene"]] = gene.plsda
    gene.tree.plsda.pathway.endoderm.time.list.filter[[i.type]][[j.type]][["p.val.data"]] = gene.anova
    gene.tree.plsda.pathway.endoderm.time.list.filter[[i.type]][[j.type]][["p.val.scaledata"]] = gene.anova.scale
    
  }
}
save(gene.tree.plsda.pathway.endoderm.time.list.filter,
     file = "Milestones/gene.tree.plsda.pathway.endoderm.time.list.filter.Rdata")

anova.curve.gene.tree.plsda.pathway.endoderm.time.list.filter = list()
anova.df.gene.tree.plsda.pathway.endoderm.time.list.filter = cbind(c(1:50))
for(i.type in names(gene.tree.plsda.pathway.endoderm.time.list.filter)){
  for(j.type in names(gene.tree.plsda.pathway.endoderm.time.list.filter[[i.type]])){
    
    gene = gene.tree.plsda.pathway.endoderm.time.list.filter[[i.type]][[j.type]][["gene"]]
    pvalue.data = gene.tree.plsda.pathway.endoderm.time.list.filter[[i.type]][[j.type]][["p.val.data"]]
    num = c()
    for(i in c(1:50)){
      num = c(num, length(gene[pvalue.data<10^(-i)]))
      if(i==1 & length(gene[pvalue.data<10^(-i)])==0){print(paste(i.type,j.type))}
    }
    names(num) = c(1:50)
    anova.curve.gene.tree.plsda.pathway.endoderm.time.list.filter[[i.type]][[j.type]] = num 
    anova.df.gene.tree.plsda.pathway.endoderm.time.list.filter = cbind(
      anova.df.gene.tree.plsda.pathway.endoderm.time.list.filter, num)
    
  }
}
save(anova.curve.gene.tree.plsda.pathway.endoderm.time.list.filter,
     anova.df.gene.tree.plsda.pathway.endoderm.time.list.filter,
     file = "Milestones/anova.df.gene.tree.plsda.pathway.endoderm.time.list.filter.summary.Rdata")

for(i.plot in c("gene.tree.plsda.pathway.endoderm.time.list.filter")){
  df.plot = melt(anova.df.gene.tree.plsda.pathway.endoderm.time.list.filter[,2:ncol(anova.df.gene.tree.plsda.pathway.endoderm.time.list.filter)])
  list.index = c()
  for(i in c(1:(ncol(anova.df.gene.tree.plsda.pathway.endoderm.time.list.filter)-1))){list.index = c(list.index, rep(i,50))}
  df.plot$Var2 = list.index
  
  pdf("Milestones/anova.curve.pathways.time.pdf",6,6)
  ggplot(data = df.plot) + 
    geom_line(mapping = aes(x=Var1, y=value, color=factor(Var2)))+
    scale_color_manual(values = unique(c(as.character(color.cluster2.new),
                                         as.character(color.time.new)))) +
    theme_classic() +
    geom_vline(xintercept=c(10,20,30),lty=4,col="black",lwd=0.8) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # axis.ticks =  element_blank(),
          # axis.text = element_blank(),
          legend.position = "none",
          # plot.title = element_blank(),
          axis.title.x = element_text(color="black", size=20, face="bold"),
          axis.title.y = element_text(color="black", size=20, face="bold"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.line.x = element_line(linetype=1, color="black", size=1.5),
          axis.line.y = element_line(linetype=1, color="black", size=1.5),
          axis.ticks.x = element_line(linetype=1, color="black", size=1.5),
          axis.ticks.y = element_line(linetype=1, color="black", size=1.5),
          axis.ticks.length  = unit(0.2, "cm"),
          aspect.ratio=1)+
    xlab("-log10(p.value)") + ylab("Detected pathways")
  dev.off()
  
  df.plot.50 = df.plot[(df.plot$value>45)&(df.plot$value<55),]
  hist(df.plot.50$Var1, breaks = 25)
  
  pdf("Milestones/anova.curve.P50.pathways.time.pdf",6,6)
  ggplot(df.plot.50) + 
    geom_histogram(mapping = aes(x=Var1),bins = 25,colour="black",fill="#eeeeee")+
    theme_classic() +
    geom_vline(xintercept=c(10,20,30),lty=4,col="black",lwd=0.8) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # axis.ticks =  element_blank(),
          # axis.text = element_blank(),
          legend.position = "none",
          # plot.title = element_blank(),
          axis.title.x = element_text(color="black", size=20, face="bold"),
          axis.title.y = element_text(color="black", size=20, face="bold"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.line.x = element_line(linetype=1, color="black", size=1.5),
          axis.line.y = element_line(linetype=1, color="black", size=1.5),
          axis.ticks.x = element_line(linetype=1, color="black", size=1.5),
          axis.ticks.y = element_line(linetype=1, color="black", size=1.5),
          axis.ticks.length  = unit(0.2, "cm"),
          aspect.ratio=1)+
    xlab("-log10(p.value)") + ylab("Frequency")
  dev.off()
}

for(i.type in names(gene.tree.plsda.pathway.endoderm.time.list.filter)){
  for(j.type in names(gene.tree.plsda.pathway.endoderm.time.list.filter[[i.type]])){
    
    gene = gene.tree.plsda.pathway.endoderm.time.list.filter[[i.type]][[j.type]][["gene"]]
    pvalue.data = gene.tree.plsda.pathway.endoderm.time.list.filter[[i.type]][[j.type]][["p.val.data"]]
    gene.rev = names(gene)
    names(gene.rev) = gene
    
    #-- p10 & p20 (Value range)
    #---------------------------- 
    gene.tree.plsda.pathway.endoderm.time.list.filter[[i.type]][[j.type]][["gene.p10"]] = gene[pvalue.data<10^(-10)]
    gene.tree.plsda.pathway.endoderm.time.list.filter[[i.type]][[j.type]][["gene.p20"]] = gene[pvalue.data<10^(-20)]
    
    if(length(gene)==0){print(paste(i.type, j.type, k.type));next()}
    if(length(gene[pvalue.data<10^(-10)])/length(gene) < 0.5){
      gene.tree.plsda.pathway.endoderm.time.list.filter[[i.type]][[j.type]][["gene.set"]] = "gene.p10"
    }else{
      gene.tree.plsda.pathway.endoderm.time.list.filter[[i.type]][[j.type]][["gene.set"]] = "gene.p20"
    }
    #---------------------------- 
    
    #-- p20 with top30 & top50
    #---------------------------- 
    gene.p10.top30 = intersect(gene, names(sort(pvalue.data)[1:min(30, length(gene[pvalue.data<10^(-20)]))]))
    names(gene.p10.top30) = gene.rev[gene.p10.top30]
    
    gene.p10.top50 = intersect(gene, names(sort(pvalue.data)[1:min(50, length(gene[pvalue.data<10^(-20)]))]))
    names(gene.p10.top50) = gene.rev[gene.p10.top50]
    
    gene.tree.plsda.pathway.endoderm.time.list.filter[[i.type]][[j.type]][["gene.p20.top30"]] = gene.p10.top30
    gene.tree.plsda.pathway.endoderm.time.list.filter[[i.type]][[j.type]][["gene.p20.top50"]] = gene.p10.top50
    #---------------------------- 
    
    #-- p20 with top30 & top50 stratified
    #---------------------------- 
    gene.p10 = gene[pvalue.data<10^(-20)]
    type.gene.p10.top30 = round(table(names(gene.p10)) / sum(table(gene.p10)) * 30, 0)
    type.gene.p10.top50 = round(table(names(gene.p10)) / sum(table(gene.p10)) * 50, 0)
    type.order.list = c(4,3,2,5,11,7)
    
    gene.p10.top30.stratified = c()
    for(i.num in intersect(type.order.list, names(type.gene.p10.top30))){
      gene.p10.i.num = gene.p10[names(gene.p10)%in%i.num]
      gene.p10.top30.stratified = c(gene.p10.top30.stratified,
                                    names(sort(pvalue.data[gene.p10.i.num])[1:max(1, type.gene.p10.top30[i.num])]))
    }
    names(gene.p10.top30.stratified) = gene.rev[gene.p10.top30.stratified] 
    
    gene.p10.top50.stratified = c()
    for(i.num in intersect(type.order.list, names(type.gene.p10.top50))){
      gene.p10.i.num = gene.p10[names(gene.p10)%in%i.num]
      gene.p10.top50.stratified = c(gene.p10.top50.stratified,
                                    names(sort(pvalue.data[gene.p10.i.num])[1:max(1, type.gene.p10.top50[i.num])]))
    }
    names(gene.p10.top50.stratified) = gene.rev[gene.p10.top50.stratified]
    
    
    gene.tree.plsda.pathway.endoderm.time.list.filter[[i.type]][[j.type]][["gene.p20.top30.stratified"]] = gene.p10.top30.stratified
    gene.tree.plsda.pathway.endoderm.time.list.filter[[i.type]][[j.type]][["gene.p20.top50.stratified"]] = gene.p10.top50.stratified
    #---------------------------- 
    
  }
}

#-- tree.plot :: gene.tree.plsda.pathway.endoderm.time.list.filter
for(i.name in names(seurat.plsda.pathway.endoderm.time.list)){
  pdf(paste("trajectory_signal/try.order.", i.name, '.pdf', sep = ""), 8, 8)
  for(j.name in names(seurat.plsda.pathway.endoderm.time.list[[i.name]])){
    
    seurat = seurat.plsda.pathway.endoderm.time.list[[i.name]][[j.name]]
    seurat$cluster.temp = seurat$cluster.v06.26.re..correct..un
    
    if(j.name %in% names(gene.tree.plsda.pathway.endoderm.time.list[[i.name]])){
      # gene_list = gene.tree.plsda.pathway.endoderm.time.list[[i.name]][[j.name]]
      gene_list = gene.tree.plsda.pathway.endoderm.time.list.filter[[i.name]][[j.name]][["gene.p20.top30.stratified"]]
      gene_list = gene_list[is.na(gene_list)==F]
    }else{next()}
    
    print(paste(i.name, j.name))
    
    cell_name = celllist.endoderm.time.list[[i.name]][[paste(j.name, ".order", sep = "")]]
    if(i.name=="ss27" & j.name == "FG.cluster2"){cell_name = rev(cell_name)}
    
    seurat = ScaleData(seurat[,cell_name], features = rownames(seurat))
    data = seurat@assays$Pathway@scale.data[gene_list, cell_name]
    data_re = t(apply(data,1,kernelsmooth))
    rownames(data_re) = rownames(data); colnames(data_re) = colnames(data)
    
    plot(c(1:10), c(1:10))
    text(x=5,y=5,labels=paste(i.type, i.name, j.name, sep = " "), cex = 5)
    
    pathway.heatmap  =
      MyHeatmap(data_re,
                type = "raw",
                ColSideColors = cbind(MyName2Col(seurat@meta.data[cell_name,]$cluster.temp, cluster.endoderm.color.v5)),
                RowSideColors = t(cbind(MyName2Col(names(gene_list), colors.geneset))),
                color.palette = colorRampPalette(c("#5aa7dd","#EBEBEB","red"), space="Lab"),
                ColSideColorsSize = 2,
                RowSideColorsSize = 2,
                Rowv = "none", 
                Colv = "none",
                # return.tree = "row",
                margins = c(10,10),
                labCol="none", graph = T)
  }
  dev.off()
}


df.summary.plsda.pathway.endoderm.time.list.filter = rbind(c("Types","Subtypes","Clusters","Pathways"))
colnames(df.summary.plsda.pathway.endoderm.time.list.filter) = c("Types","Subtypes","Clusters","Pathways")
for(i.type in names(gene.tree.plsda.pathway.endoderm.time.list.filter)){
  for(j.type in names(gene.tree.plsda.pathway.endoderm.time.list.filter[[i.type]])){
    gene.temp = gene.tree.plsda.pathway.endoderm.time.list.filter[[i.type]][[j.type]]$gene.p20.top30.stratified
    gene.temp = gene.temp[is.na(gene.temp)==F]
    df.temp = cbind.data.frame(
      rep(i.type, length(gene.temp)),
      rep(j.type, length(gene.temp)),
      names(gene.temp), gene.temp)
    colnames(df.temp) = c("Types","Subtypes","Clusters","Pathways")
    df.summary.plsda.pathway.endoderm.time.list.filter = rbind(df.summary.plsda.pathway.endoderm.time.list.filter, df.temp)
  }
}
df.summary.plsda.pathway.endoderm.time.list.filter = as.data.frame(df.summary.plsda.pathway.endoderm.time.list.filter)

#==========================================================



#==========================================================
#   pathway.endoderm.time.list.add.refine :: 12SS
#==========================================================
pathway.endoderm.time_add.list = list()
for(i.time in c("ss12")){
  pathway.endoderm.time_add.list[[i.time]] = pathway_calculation(
    get(paste("src.",gsub("ss","",i.time),"ss.integrated",sep = "")))
}

celllist.endoderm.time_add.list = list()
for(i.celllist in c("celllist.endoderm.time_add.list")){
  for(i.time in c("ss12")){
    celllist.endoderm.time_add.list[[i.time]] = list()
  }
  
  DimPlot(src.12ss.integrated, group.by = "cluster.v06.26.re..correct..un",
          cols = cluster.endoderm.color.v5)
  celllist.endoderm.time_add.list$ss12$merge = colnames(src.12ss.integrated)
  celllist.endoderm.time_add.list$ss12$Line1 = rownames(src.12ss.integrated@meta.data[
    src.12ss.integrated$cluster.v06.26.re..correct..un%in%c("Pharynx.organ.2","FG.1","FG.6","Esophagus","MG.3.A","DP","EP.1",
                                                            "MG.3.M","MG.3","MG.3.P","HG.2","Large.intestine.3"),])
  celllist.endoderm.time_add.list$ss12$Line2 = rownames(src.12ss.integrated@meta.data[
    src.12ss.integrated$cluster.v06.26.re..correct..un%in%c('FG.2',"Pharynx.organ.1","FG.3","Pharynx.organ.4","Stomach","MG.1",
                                                            "Small.intestine.1","Small.intestine.2","MG.2","HG.1","Large.intestine.1",
                                                            "HG.1-Large.intestine.2","Large.intestine.2"),])
  celllist.endoderm.time_add.list$ss12$Line3 = rownames(src.12ss.integrated@meta.data[
    src.12ss.integrated$cluster.v06.26.re..correct..un%in%c("FG.5","Thyroid","Pharynx.organ.3","Pharynx.organ.5","Lung","FG.4","FG.4-Lung/Stomach",
                                                            "FG.4-Liver","AL.1/2-Liver","Liver","EHBD","AL.3-Liver","VP","AL.3-EHBD/VP",'AL.3',"AL.3-Small.intestine.1"),])
}

seurat.plsda.pathway.endoderm.time_add.list = list()
for(i.time in c(names(pathway.endoderm.time_add.list))){
  seurat.plsda.pathway.endoderm.time_add.list[[i.time]] = list()
  for(j.type in names(celllist.endoderm.time_add.list[[i.time]])){
    seurat.plsda.pathway.endoderm.time_add.list[[i.time]][[j.type]] =
      seurat_pathway_pre(pathway.matrix = pathway.endoderm.time_add.list[[i.time]][,celllist.endoderm.time_add.list[[i.time]][[j.type]]], 
                         seurat.ref = get(paste("src.", gsub("ss","",i.time), "ss.integrated", sep = ""))[,celllist.endoderm.time_add.list[[i.time]][[j.type]]]) 
  }
}

#-- add pathway-plsda
for(i.time in names(seurat.plsda.pathway.endoderm.time_add.list)){
  for(i.type in names(seurat.plsda.pathway.endoderm.time_add.list[[i.time]])){
    seurat.plsda.pathway.endoderm.time_add.list[[i.time]][[i.type]] = 
      seurat_plsda_pathway(seurat.plsda.pathway.endoderm.time_add.list[[i.time]][[i.type]],
                           seurat.plsda.pathway.endoderm.time_add.list[[i.time]][[i.type]]$cluster.v06.26.re..correct..un)[["seurat"]]
  }
}

princurve.endoderm.time_add.instances = list()
for(i.curve_plot in c("princurve.endoderm.time_add.instances")){
  
  for(i.time in c("ss12","ss18","ss21","ss24")){
    princurve.endoderm.time_add.instances[[i.time]] = list()
  }
  
  #-- 12ss
  #--------------------------
  princurve.endoderm.time_add.instances$ss12$Line1 = princurve::principal_curve(
    x = cbind(src.12ss.integrated@reductions$umap@cell.embeddings[,1],
              src.12ss.integrated@reductions$umap@cell.embeddings[,2])[
                celllist.endoderm.time_add.list$ss12$Line1,],
    smoother = "smooth.spline")
  princurve.endoderm.time_add.instances$ss12$Line2 = princurve::principal_curve(
    x = cbind(src.12ss.integrated@reductions$umap@cell.embeddings[,1],
              src.12ss.integrated@reductions$umap@cell.embeddings[,2])[
                celllist.endoderm.time_add.list$ss12$Line2,],
    smoother = "smooth.spline")
  princurve.endoderm.time_add.instances$ss12$Line3 = princurve::principal_curve(
    x = cbind(src.12ss.integrated@reductions$umap@cell.embeddings[,1],
              src.12ss.integrated@reductions$umap@cell.embeddings[,2])[
                celllist.endoderm.time_add.list$ss12$Line3,],
    smoother = "smooth.spline")
  
  for(i.type in c("Line1","Line2","Line3")){
    celllist.endoderm.time_add.list$ss12[[paste(i.type,".order",sep = "")]] = rownames(
      princurve.endoderm.time_add.instances$ss12[[paste(i.type,sep = "")]]$s[
        order(princurve.endoderm.time_add.instances$ss12[[paste(i.type,sep = "")]]$lambda),])
  }
  #--------------------------
  

}

#-- test
for(i.test in c("princurve.endoderm.time_add.instances")){
  
  #--  12ss
  #--------------------------
  p12 = DimPlot(src.12ss.integrated, reduction = "umap", pt.size = 1.5,
                group.by = "Time") + 
    scale_color_manual(values = c("gray")) + # cluster.endoderm.color.v5) +
    theme_void()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          plot.title = element_blank(),
          aspect.ratio=1) +
    geom_path(data = as.data.frame(princurve.endoderm.time_add.instances$ss12$Line1$s[
      celllist.endoderm.time_add.list$ss12$Line1.order,]), 
      mapping = aes(x = V1, y = V2), size=3.5, col=I("#87afcc")) +
    geom_path(data = as.data.frame(princurve.endoderm.time_add.instances$ss12$Line3$s[
      celllist.endoderm.time_add.list$ss12$Line3.order,]), 
      mapping = aes(x = V1, y = V2), size=3.5, col=I("#E7AE27")) +
    geom_path(data = as.data.frame(princurve.endoderm.time_add.instances$ss12$Line2$s[
      celllist.endoderm.time_add.list$ss12$Line2.order,]), 
      mapping = aes(x = V1, y = V2), size=3.5, col=I("#FB867A"))
  #--------------------------
  
  for(i.p in c("p12")){
    png(filename = paste("trajectory_signal/Pattern_code.", i.p, ".png", sep = ""),
        width = 1000,height = 1000,pointsize = 20)
    print(get(i.p))
    dev.off()
  }
}

#-- correct: cell order and reverse
for(i.cellorderlist in c("celllist.endoderm.time_add.list")){
  celllist.endoderm.time_add.list$ss12$Line1.order = rev(celllist.endoderm.time_add.list$ss12$Line1.order)
  celllist.endoderm.time_add.list$ss12$Line2.order = rev(celllist.endoderm.time_add.list$ss12$Line2.order)
  celllist.endoderm.time_add.list$ss12$Line3.order = rev(celllist.endoderm.time_add.list$ss12$Line3.order)
}


#-- Add-new seurat with plsda (step1)
for(i.time in c(names(pathway.endoderm.time_add.list))){
  for(j.type in names(celllist.endoderm.time_add.list[[i.time]])){
    # if(i.time == "ss21" & j.type == "FG.cluster3.order"){}else{next()}
    if(grepl("order",j.type)==F){next()}
    seurat.plsda.pathway.endoderm.time_add.list[[i.time]][[j.type]] =
      seurat_pathway_pre(pathway.matrix = pathway.endoderm.time_add.list[[i.time]][,celllist.endoderm.time_add.list[[i.time]][[j.type]]], 
                         seurat.ref = get(paste("src.", gsub("ss","",i.time), "ss.integrated", sep = ""))[,celllist.endoderm.time_add.list[[i.time]][[j.type]]]) 
  }
}

#-- Add-new seurat with plsda (step2)
for(i.time in names(seurat.plsda.pathway.endoderm.time_add.list)){
  for(i.type in names(seurat.plsda.pathway.endoderm.time_add.list[[i.time]])){
    if(grepl("order",j.type)==F){next()}
    seurat.plsda.pathway.endoderm.time_add.list[[i.time]][[i.type]] = 
      seurat_plsda_pathway(seurat.plsda.pathway.endoderm.time_add.list[[i.time]][[i.type]],
                           seurat.plsda.pathway.endoderm.time_add.list[[i.time]][[i.type]]$cluster.v06.26.re..correct..un)[["seurat"]]
  }
}

pathway.plsda.pathway.endoderm.time_add.list = list()
for(i.type in names(seurat.plsda.pathway.endoderm.time_add.list)){
  for(i.name in names(seurat.plsda.pathway.endoderm.time_add.list[[i.type]])){
    pathway.plsda.pathway.endoderm.time_add.list[[i.type]][[i.name]] =
      plsda_select(seurat.plsda.pathway.endoderm.time_add.list[[i.type]][[i.name]], comp = 3, threshold = 30)
  }
}

tree.plsda.pathway.endoderm.time_add.list = list()
for(i.name in names(seurat.plsda.pathway.endoderm.time_add.list)){
  tree.plsda.pathway.endoderm.time_add.list[[i.name]] = list()
  pdf(paste("Milestones/try.", i.name, '.pdf', sep = ""), 10, 10)
  
  for(j.name in names(seurat.plsda.pathway.endoderm.time_add.list[[i.name]])){
    if(grepl("order",j.type)==F){next()}
    
    seurat = seurat.plsda.pathway.endoderm.time_add.list[[i.name]][[j.name]]
    seurat$cluster.temp = seurat$cluster.v06.26.re..correct..un
    gene_list = pathway.plsda.pathway.endoderm.time_add.list[[i.name]][[j.name]]
    
    cellorder_list = celllist.endoderm.time_add.list[[i.name]][[j.name]]
    cell_name = intersect(cellorder_list, colnames(seurat))
    
    seurat = ScaleData(seurat[,cell_name], features = rownames(seurat))
    
    data = seurat@assays$Pathway@scale.data[gene_list, cell_name]
    data_re = t(apply(data,1,kernelsmooth))
    rownames(data_re) = rownames(data); colnames(data_re) = colnames(data)
    
    plot(c(1:10), c(1:10))
    text(x=5,y=5,labels=paste(i.type, i.name, j.name, sep = " "), cex = 5)
    
    pathway.heatmap  =
      MyHeatmap(data_re,
                type = "raw",
                ColSideColors = cbind(
                  MyName2Col(seurat@meta.data[cell_name,]$cluster.temp, cluster.endoderm.color.v5)),
                color.palette = colorRampPalette(c("#5aa7dd","#EBEBEB","red"), space="Lab"),
                ColSideColorsSize = 3,
                RowSideColorsSize = 2,
                # Rowv = "none", 
                # Colv = "none",
                return.tree = "row",
                labCol="none", graph = T)
    
    tree.plsda.pathway.endoderm.time_add.list[[i.name]][[j.name]] = as.dendrogram(pathway.heatmap) 
  }
  dev.off()
}

gene.tree.plsda.pathway.endoderm.time_add.list = list()
for(i.tree.plsda in c("gene.tree.plsda.pathway.endoderm.list.refine")){
  #-- Add-tree
  for(i.type in names(tree.plsda.pathway.endoderm.time_add.list)){
    gene.tree.plsda.pathway.endoderm.time_add.list[[i.type]] = list()
    for(j.type in names(tree.plsda.pathway.endoderm.time_add.list[[i.type]])){
      gene.tree.plsda.pathway.endoderm.time_add.list[[i.type]][[j.type]] = NA
    }
  }
  
  #-- SS12: Line1, Line2, Line3
  #----------------------------------
  gene.tree.plsda.pathway.endoderm.time_add.list$ss12$Line1.order = c(
    labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line1.order[[2]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line1.order[[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line1.order[[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line1.order[[1]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line1.order[[1]][[2]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.time_add.list$ss12$Line1.order) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line1.order[[2]][[2]][[2]])))),
    rep(3, length(c(labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line1.order[[2]][[2]][[1]])))),
    rep(5, length(c(labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line1.order[[2]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line1.order[[1]][[2]][[2]]),
                    labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line1.order[[1]][[2]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.time_add.list$ss12$Line2.order = c(
    labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line2.order[[2]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line2.order[[2]][[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line2.order[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line2.order[[1]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line2.order[[1]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line2.order[[2]][[2]][[1]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.time_add.list$ss12$Line2.order) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line2.order[[2]][[2]][[2]])))),
    rep(3, length(c(labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line2.order[[2]][[1]][[2]])))),
    rep(2, length(c(labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line2.order[[1]][[1]])))),
    rep(5, length(c(labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line2.order[[1]][[2]][[2]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line2.order[[1]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line2.order[[2]][[2]][[1]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.time_add.list$ss12$Line3.order = c(
    labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line3.order[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line3.order[[2]][[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line3.order[[2]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line3.order[[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line3.order[[1]][[2]]))
  names(gene.tree.plsda.pathway.endoderm.time_add.list$ss12$Line3.order) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line3.order[[1]][[1]])))),
    rep(3, length(c(labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line3.order[[2]][[1]][[2]])))),
    rep(5, length(c(labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line3.order[[2]][[2]][[2]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line3.order[[2]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.time_add.list$ss12$Line3.order[[1]][[2]])))))
  #----------------------------------
}

#-- Filter by anova
#--------------------------------------
gene.tree.plsda.pathway.endoderm.time_add.list.filter = list()
for(i.type in names(gene.tree.plsda.pathway.endoderm.time_add.list)){
  gene.tree.plsda.pathway.endoderm.time_add.list.filter[[i.type]] = list()
  for(j.type in names(gene.tree.plsda.pathway.endoderm.time_add.list[[i.type]])){
    if(grepl("order", j.type)==F){next()}
    print(paste(i.type,j.type))
    
    src.plsda = seurat.plsda.pathway.endoderm.time_add.list[[i.type]][[j.type]]
    gene.plsda = gene.tree.plsda.pathway.endoderm.time_add.list[[i.type]][[j.type]]
    gene.anova = anova.test(tpm = src.plsda@assays$Pathway@data[gene.plsda, colnames(src.plsda)],
                            variable = src.plsda$cluster.v06.26.re..correct..un)
    gene.anova.scale = anova.test(tpm = src.plsda@assays$Pathway@scale.data[gene.plsda, colnames(src.plsda)],
                                  variable = src.plsda$cluster.v06.26.re..correct..un)
    
    gene.tree.plsda.pathway.endoderm.time_add.list.filter[[i.type]][[j.type]][["gene"]] = gene.plsda
    gene.tree.plsda.pathway.endoderm.time_add.list.filter[[i.type]][[j.type]][["p.val.data"]] = gene.anova
    gene.tree.plsda.pathway.endoderm.time_add.list.filter[[i.type]][[j.type]][["p.val.scaledata"]] = gene.anova.scale
  }
}

anova.curve.gene.tree.plsda.pathway.endoderm.time_add.list.filter = list()
anova.df.gene.tree.plsda.pathway.endoderm.time_add.list.filter = cbind(c(1:50))
for(i.type in names(gene.tree.plsda.pathway.endoderm.time_add.list.filter)){
  for(j.type in names(gene.tree.plsda.pathway.endoderm.time_add.list.filter[[i.type]])){
    gene = gene.tree.plsda.pathway.endoderm.time_add.list.filter[[i.type]][[j.type]][["gene"]]
    pvalue.data = gene.tree.plsda.pathway.endoderm.time_add.list.filter[[i.type]][[j.type]][["p.val.data"]]
    num = c()
    for(i in c(1:50)){
      num = c(num, length(gene[pvalue.data<10^(-i)]))
      if(i==1 & length(gene[pvalue.data<10^(-i)])==0){print(paste(i.type,j.type))}
    }
    names(num) = c(1:50)
    anova.curve.gene.tree.plsda.pathway.endoderm.time_add.list.filter[[i.type]][[j.type]] = num 
    anova.df.gene.tree.plsda.pathway.endoderm.time_add.list.filter = cbind(
      anova.df.gene.tree.plsda.pathway.endoderm.time_add.list.filter, num)
  }
}

for(i.plot in c("gene.tree.plsda.pathway.endoderm.time_add.list.filter")){
  df.plot = melt(anova.df.gene.tree.plsda.pathway.endoderm.time_add.list.filter[,2:ncol(anova.df.gene.tree.plsda.pathway.endoderm.time_add.list.filter)])
  list.index = c()
  for(i in c(1:(ncol(anova.df.gene.tree.plsda.pathway.endoderm.time_add.list.filter)-1))){list.index = c(list.index, rep(i,50))}
  df.plot$Var2 = list.index
  
  pdf("Milestones/anova.curve.pathways.pdf",6,6)
  ggplot(data = df.plot) + 
    geom_line(mapping = aes(x=Var1, y=value, color=factor(Var2)))+
    scale_color_manual(values = unique(c(as.character(color.cluster2.new),
                                         as.character(color.time.new)))) +
    theme_classic() +
    geom_vline(xintercept=c(10,20,30),lty=4,col="black",lwd=0.8) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # axis.ticks =  element_blank(),
          # axis.text = element_blank(),
          legend.position = "none",
          # plot.title = element_blank(),
          axis.title.x = element_text(color="black", size=20, face="bold"),
          axis.title.y = element_text(color="black", size=20, face="bold"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.line.x = element_line(linetype=1, color="black", size=1.5),
          axis.line.y = element_line(linetype=1, color="black", size=1.5),
          axis.ticks.x = element_line(linetype=1, color="black", size=1.5),
          axis.ticks.y = element_line(linetype=1, color="black", size=1.5),
          axis.ticks.length  = unit(0.2, "cm"),
          aspect.ratio=1)+
    xlab("-log10(p.value)") + ylab("Detected pathways")
  dev.off()
  
  df.plot.50 = df.plot[(df.plot$value>25)&(df.plot$value<35),]
  hist(df.plot.50$Var1, breaks = 40)
  
  pdf("Milestones/anova.curve.P50.pathways.pdf",6,6)
  ggplot(df.plot.50) + 
    geom_histogram(mapping = aes(x=Var1),bins = 40,colour="black",fill="#eeeeee")+
    theme_classic() +
    geom_vline(xintercept=c(10,20,30),lty=4,col="black",lwd=0.8) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # axis.ticks =  element_blank(),
          # axis.text = element_blank(),
          legend.position = "none",
          # plot.title = element_blank(),
          axis.title.x = element_text(color="black", size=20, face="bold"),
          axis.title.y = element_text(color="black", size=20, face="bold"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.line.x = element_line(linetype=1, color="black", size=1.5),
          axis.line.y = element_line(linetype=1, color="black", size=1.5),
          axis.ticks.x = element_line(linetype=1, color="black", size=1.5),
          axis.ticks.y = element_line(linetype=1, color="black", size=1.5),
          axis.ticks.length  = unit(0.2, "cm"),
          aspect.ratio=1)+
    xlab("-log10(p.value)") + ylab("Frequency")
  dev.off()
}

for(i.type in names(gene.tree.plsda.pathway.endoderm.time_add.list.filter)){
  for(j.type in names(gene.tree.plsda.pathway.endoderm.time_add.list.filter[[i.type]])){
    
    gene = gene.tree.plsda.pathway.endoderm.time_add.list.filter[[i.type]][[j.type]][["gene"]]
    pvalue.data = gene.tree.plsda.pathway.endoderm.time_add.list.filter[[i.type]][[j.type]][["p.val.data"]]
    gene.rev = names(gene)
    names(gene.rev) = gene
    
    #-- p10 & p20 (Value range)
    #---------------------------- 
    gene.tree.plsda.pathway.endoderm.time_add.list.filter[[i.type]][[j.type]][["gene.p10"]] = gene[pvalue.data<10^(-10)]
    gene.tree.plsda.pathway.endoderm.time_add.list.filter[[i.type]][[j.type]][["gene.p20"]] = gene[pvalue.data<10^(-20)]
    
    if(length(gene)==0){print(paste(i.type, j.type, k.type));next()}
    if(length(gene[pvalue.data<10^(-10)])/length(gene) < 0.5){
      gene.tree.plsda.pathway.endoderm.time_add.list.filter[[i.type]][[j.type]][["gene.set"]] = "gene.p10"
    }else{
      gene.tree.plsda.pathway.endoderm.time_add.list.filter[[i.type]][[j.type]][["gene.set"]] = "gene.p20"
    }
    #---------------------------- 
    
    #-- p10 with top30 & top50
    #---------------------------- 
    gene.p10.top30 = intersect(gene, names(sort(pvalue.data)[1:min(30, length(gene[pvalue.data<10^(-10)]))]))
    names(gene.p10.top30) = gene.rev[gene.p10.top30]
    
    gene.p10.top50 = intersect(gene, names(sort(pvalue.data)[1:min(50, length(gene[pvalue.data<10^(-10)]))]))
    names(gene.p10.top50) = gene.rev[gene.p10.top50]
    
    gene.tree.plsda.pathway.endoderm.time_add.list.filter[[i.type]][[j.type]][["gene.p10.top30"]] = gene.p10.top30
    gene.tree.plsda.pathway.endoderm.time_add.list.filter[[i.type]][[j.type]][["gene.p10.top50"]] = gene.p10.top50
    #---------------------------- 
    
    #-- p10 with top30 & top50 stratified
    #---------------------------- 
    gene.p10 = gene[pvalue.data<10^(-10)]
    type.gene.p10.top30 = round(table(names(gene.p10)) / sum(table(gene.p10)) * 30, 0)
    type.gene.p10.top50 = round(table(names(gene.p10)) / sum(table(gene.p10)) * 50, 0)
    type.order.list = c(4,3,2,5,7)
    
    gene.p10.top30.stratified = c()
    for(i.num in intersect(type.order.list, names(type.gene.p10.top30))){
      gene.p10.i.num = gene.p10[names(gene.p10)%in%i.num]
      gene.p10.top30.stratified = c(gene.p10.top30.stratified,
                                    names(sort(pvalue.data[gene.p10.i.num])[1:max(1, type.gene.p10.top30[i.num])]))
    }
    names(gene.p10.top30.stratified) = gene.rev[gene.p10.top30.stratified] 
    
    gene.p10.top50.stratified = c()
    for(i.num in intersect(type.order.list, names(type.gene.p10.top50))){
      gene.p10.i.num = gene.p10[names(gene.p10)%in%i.num]
      gene.p10.top50.stratified = c(gene.p10.top50.stratified,
                                    names(sort(pvalue.data[gene.p10.i.num])[1:max(1, type.gene.p10.top50[i.num])]))
    }
    names(gene.p10.top50.stratified) = gene.rev[gene.p10.top50.stratified]
    
    
    gene.tree.plsda.pathway.endoderm.time_add.list.filter[[i.type]][[j.type]][["gene.p10.top30.stratified"]] = gene.p10.top30.stratified
    gene.tree.plsda.pathway.endoderm.time_add.list.filter[[i.type]][[j.type]][["gene.p10.top50.stratified"]] = gene.p10.top50.stratified
    #---------------------------- 
    
  }
}
#--------------------------------------

save(pathway.endoderm.time_add.list,
     seurat.plsda.pathway.endoderm.time_add.list,
     celllist.endoderm.time_add.list,
     princurve.endoderm.time_add.instances,
     pathway.plsda.pathway.endoderm.time_add.list,
     tree.plsda.pathway.endoderm.time_add.list,
     gene.tree.plsda.pathway.endoderm.time_add.list,
     gene.tree.plsda.pathway.endoderm.time_add.list.filter,
     file = "trajectory_signal/seurat.plsda.pathway.endoderm.time_add.list.parameter.Rdata")

#-- fin-plot
for(i.name in names(seurat.plsda.pathway.endoderm.time_add.list)){
  pdf(paste("trajectory_signal/try.", i.name, '.pdf', sep = ""), 8, 8)
  for(j.name in names(seurat.plsda.pathway.endoderm.time_add.list[[i.name]])){
    if(grepl("order",j.name)==F){next()}
    seurat = seurat.plsda.pathway.endoderm.time_add.list[[i.name]][[j.name]]
    seurat$cluster.temp = seurat$cluster.v06.26.re..correct..un
    
    cell_name = intersect(celllist.endoderm.time_add.list[[i.name]][[j.name]], colnames(seurat))
    gene_list = gene.tree.plsda.pathway.endoderm.time_add.list.filter[[i.name]][[j.name]][["gene.p10.top30.stratified"]]
    gene_list = gene_list[is.na(gene_list)==F]
    print(paste(i.name, j.name))
    
    data = seurat@assays$Pathway@scale.data[gene_list, cell_name]
    data_re = t(apply(data,1,kernelsmooth))
    rownames(data_re) = rownames(data); colnames(data_re) = colnames(data)
    
    pathway.heatmap  =
      MyHeatmap(data_re,
                type = "raw",
                ColSideColors = cbind(
                  MyName2Col(seurat@meta.data[cell_name,]$cluster.temp, cluster.endoderm.color.v5)
                ),
                RowSideColors = t(cbind(
                  MyName2Col(names(gene_list), colors.geneset)
                )),
                color.palette = colorRampPalette(c("#5aa7dd","#EBEBEB","red"), space="Lab"),
                ColSideColorsSize = 2,
                RowSideColorsSize = 2,
                Rowv = "none", 
                Colv = "none",
                # return.tree = "row",
                margins = c(10,10),
                labCol="none", graph = T)
  }
  dev.off()
}

#==========================================================


#==========================================================
#   pathway.endoderm.time.list.add.refine :: 15SS(re) 18SS 21SS 24SS 27SS
#==========================================================

pathway.endoderm.time_add_new.list = list()
for(i.time in c('ss15',"ss18","ss21","ss24","ss27")){
  pathway.endoderm.time_add_new.list[[i.time]] = pathway_calculation(
    get(paste("src.",gsub("ss","",i.time),"ss.integrated",sep = "")))
}

celllist.endoderm.time_add_new.list = list()
for(i.celllist in c("celllist.endoderm.time_add_new.list")){
  for(i.time in c('ss15',"ss18","ss21","ss24","ss27")){
    celllist.endoderm.time_add_new.list[[i.time]] = list()
  }
  
  DimPlot(src.15ss.integrated, group.by = "cluster.v06.26.re..correct..un",
          cols = cluster.endoderm.color.v5)
  celllist.endoderm.time_add_new.list$ss15$Line3 = rownames(src.15ss.integrated@meta.data[
    src.15ss.integrated$cluster.v06.26.re..correct..un%in%c("FG.5","Thyroid","Pharynx.organ.3","Pharynx.organ.5","Lung","FG.4",
                                                            "FG.4-Lung/Stomach","FG.4-Liver"),])
  
  DimPlot(src.18ss.integrated, group.by = "cluster.v06.26.re..correct..un",
          cols = cluster.endoderm.color.v5)
  celllist.endoderm.time_add_new.list$ss18$MG = rownames(src.18ss.integrated@meta.data[
    src.18ss.integrated$cluster.v06.26.re..correct..un%in%c("Pharynx.organ.1","FG.2",
                                                            "Pharynx.organ.4","Lung","FG.3","FG.4-Lung/Stomach",
                                                            "Pharynx.organ.2","Esophagus","FG.1","FG.6","Stomach",
                                                            "Thyroid","Pharynx.organ.3","Pharynx.organ.5","Lung",
                                                            "Small.intestine.1","Small.intestine.2",
                                                            "Large.intestine.1","Large.intestine.2","Large.intestine.3",
                                                            "MG.3.A","MG.3.M","DP","MG.3.P","MG.1","MG.2",
                                                            "HG.1","HG.2","HG.1-Large.intestine.2"),])
  
  DimPlot(src.21ss.integrated, group.by = "cluster.v06.26.re..correct..un",
          cols = cluster.endoderm.color.v5)
  celllist.endoderm.time_add_new.list$ss21$MG = rownames(src.21ss.integrated@meta.data[
    src.21ss.integrated$cluster.v06.26.re..correct..un%in%c("Pharynx.organ.1","FG.2",
                                                            "Pharynx.organ.4","Lung","FG.3","FG.4-Lung/Stomach",
                                                            "Pharynx.organ.2","Esophagus","FG.1","FG.6","Stomach",
                                                            "Thyroid","Pharynx.organ.3","Pharynx.organ.5","Lung",
                                                            "Small.intestine.1","Small.intestine.2",
                                                            "Large.intestine.1","Large.intestine.2","Large.intestine.3",
                                                            "MG.3.A","MG.3.M","DP","MG.3.P","MG.1","MG.2",
                                                            "HG.1","HG.2","HG.1-Large.intestine.2"),])
  
  DimPlot(src.24ss.integrated, group.by = "cluster.v06.26.re..correct..un",
          cols = cluster.endoderm.color.v5)
  celllist.endoderm.time_add_new.list$ss24$MG = rownames(src.24ss.integrated@meta.data[
    src.24ss.integrated$cluster.v06.26.re..correct..un%in%c("Pharynx.organ.1","FG.2",
                                                            "Pharynx.organ.4","Lung","FG.3","FG.4-Lung/Stomach",
                                                            "Pharynx.organ.2","Esophagus","FG.1","FG.6","Stomach",
                                                            "Thyroid","Pharynx.organ.3","Pharynx.organ.5","Lung",
                                                            "Small.intestine.1","Small.intestine.2",
                                                            "Large.intestine.1","Large.intestine.2","Large.intestine.3",
                                                            "MG.3.A","MG.3.M","MG.3.P","MG.1","MG.2",
                                                            "HG.1","HG.2","HG.1-Large.intestine.2"),])
  
  DimPlot(src.27ss.integrated, group.by = "cluster.v06.26.re..correct..un",
          cols = cluster.endoderm.color.v5)
  celllist.endoderm.time_add_new.list$ss27$MG = rownames(src.27ss.integrated@meta.data[
    src.27ss.integrated$cluster.v06.26.re..correct..un%in%c("Pharynx.organ.1","FG.2",
                                                            "Pharynx.organ.4","Lung","FG.3","FG.4-Lung/Stomach",
                                                            "Pharynx.organ.2","Esophagus","FG.1","FG.6","Stomach",
                                                            "Thyroid","Pharynx.organ.3","Pharynx.organ.5","Lung",
                                                            "Small.intestine.1","Small.intestine.2",
                                                            "Large.intestine.1","Large.intestine.2","Large.intestine.3",
                                                            "MG.3.A","MG.3.M","MG.3.P","MG.1","MG.2",
                                                            "HG.1","HG.2","HG.1-Large.intestine.2"),])
  
}

seurat.plsda.pathway.endoderm.time_add_new.list = list()
for(i.time in c(names(pathway.endoderm.time_add_new.list))){
  seurat.plsda.pathway.endoderm.time_add_new.list[[i.time]] = list()
  for(j.type in names(celllist.endoderm.time_add_new.list[[i.time]])){
    seurat.plsda.pathway.endoderm.time_add_new.list[[i.time]][[j.type]] =
      seurat_pathway_pre(pathway.matrix = pathway.endoderm.time_add_new.list[[i.time]][,celllist.endoderm.time_add_new.list[[i.time]][[j.type]]], 
                         seurat.ref = get(paste("src.", gsub("ss","",i.time), "ss.integrated", sep = ""))[,celllist.endoderm.time_add_new.list[[i.time]][[j.type]]]) 
  }
}

#-- add pathway-plsda
for(i.time in names(seurat.plsda.pathway.endoderm.time_add_new.list)){
  for(i.type in names(seurat.plsda.pathway.endoderm.time_add_new.list[[i.time]])){
    seurat.plsda.pathway.endoderm.time_add_new.list[[i.time]][[i.type]] = 
      seurat_plsda_pathway(seurat.plsda.pathway.endoderm.time_add_new.list[[i.time]][[i.type]],
                           seurat.plsda.pathway.endoderm.time_add_new.list[[i.time]][[i.type]]$cluster.v06.26.re..correct..un)[["seurat"]]
  }
}

princurve.endoderm.time_add_new.instances = list()
for(i.curve_plot in c("princurve.endoderm.time_add_new.instances")){
  
  for(i.time in c("ss15","ss18","ss21","ss24","ss27")){
    princurve.endoderm.time_add_new.instances[[i.time]] = list()
  }
  
  #-- 15ss
  #--------------------------
  princurve.endoderm.time_add_new.instances$ss15$Line3 = princurve::principal_curve(
    x = cbind(src.15ss.integrated@reductions$umap@cell.embeddings[,1],
              src.15ss.integrated@reductions$umap@cell.embeddings[,2])[
                celllist.endoderm.time_add_new.list$ss15$Line3,],
    smoother = "smooth.spline")
  
  for(i.type in c("Line3")){
    celllist.endoderm.time_add_new.list$ss15[[paste(i.type,".order",sep = "")]] = rownames(
      princurve.endoderm.time_add_new.instances$ss15[[paste(i.type,sep = "")]]$s[
        order(princurve.endoderm.time_add_new.instances$ss15[[paste(i.type,sep = "")]]$lambda),])
  }
  #--------------------------
  
  #-- 18ss
  #--------------------------
  princurve.endoderm.time_add_new.instances$ss18$MG = princurve::principal_curve(
    x = cbind(src.18ss.integrated@reductions$umap.rotated@cell.embeddings[,1],
              src.18ss.integrated@reductions$umap.rotated@cell.embeddings[,2])[
                celllist.endoderm.time_add_new.list$ss18$MG,],
    smoother = "smooth.spline")
  
  for(i.type in c("MG")){
    celllist.endoderm.time_add_new.list$ss18[[paste(i.type,".order",sep = "")]] = rownames(
      princurve.endoderm.time_add_new.instances$ss18[[paste(i.type,sep = "")]]$s[
        order(princurve.endoderm.time_add_new.instances$ss18[[paste(i.type,sep = "")]]$lambda),])
  }
  #--------------------------
  
  #-- 21ss
  #--------------------------
  princurve.endoderm.time_add_new.instances$ss21$MG = princurve::principal_curve(
    x = cbind(src.21ss.integrated@reductions$umap@cell.embeddings[,1],
              src.21ss.integrated@reductions$umap@cell.embeddings[,2])[
                celllist.endoderm.time_add_new.list$ss21$MG,],
    smoother = "smooth.spline")
  
  for(i.type in c("MG")){
    celllist.endoderm.time_add_new.list$ss21[[paste(i.type,".order",sep = "")]] = rownames(
      princurve.endoderm.time_add_new.instances$ss21[[paste(i.type,sep = "")]]$s[
        order(princurve.endoderm.time_add_new.instances$ss21[[paste(i.type,sep = "")]]$lambda),])
  }
  #--------------------------
  
  #-- 24ss
  #--------------------------
  princurve.endoderm.time_add_new.instances$ss24$MG = princurve::principal_curve(
    x = cbind(src.24ss.integrated@reductions$umap@cell.embeddings[,1],
              src.24ss.integrated@reductions$umap@cell.embeddings[,2])[
                celllist.endoderm.time_add_new.list$ss24$MG,],
    smoother = "smooth.spline")
  
  for(i.type in c("MG")){
    celllist.endoderm.time_add_new.list$ss24[[paste(i.type,".order",sep = "")]] = rownames(
      princurve.endoderm.time_add_new.instances$ss24[[paste(i.type,sep = "")]]$s[
        order(princurve.endoderm.time_add_new.instances$ss24[[paste(i.type,sep = "")]]$lambda),])
  }
  #--------------------------
  
  #-- 27ss
  #--------------------------
  princurve.endoderm.time_add_new.instances$ss27$MG = princurve::principal_curve(
    x = cbind(src.27ss.integrated@reductions$umap@cell.embeddings[,1],
              src.27ss.integrated@reductions$umap@cell.embeddings[,2])[
                celllist.endoderm.time_add_new.list$ss27$MG,],
    smoother = "smooth.spline")
  
  for(i.type in c("MG")){
    celllist.endoderm.time_add_new.list$ss27[[paste(i.type,".order",sep = "")]] = rownames(
      princurve.endoderm.time_add_new.instances$ss27[[paste(i.type,sep = "")]]$s[
        order(princurve.endoderm.time_add_new.instances$ss27[[paste(i.type,sep = "")]]$lambda),])
  }
  #--------------------------
  
}

#-- test
for(i.test in c("princurve.endoderm.time_add_new.instances")){
  
  #--  12ss :: time_add
  #--------------------------
  p12 = DimPlot(src.12ss.integrated, reduction = "umap", pt.size = 1.5,
                group.by = "Time") + 
    scale_color_manual(values = c("gray")) + # cluster.endoderm.color.v5) +
    theme_void()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          plot.title = element_blank(),
          aspect.ratio=1) +
    geom_path(data = as.data.frame(princurve.endoderm.time_add.instances$ss12$Line1$s[
      celllist.endoderm.time_add.list$ss12$Line1.order,]), 
      mapping = aes(x = V1, y = V2), size=5, col=I("#87afcc")) +
    geom_path(data = as.data.frame(princurve.endoderm.time_add.instances$ss12$Line3$s[
      celllist.endoderm.time_add.list$ss12$Line3.order,]), 
      mapping = aes(x = V1, y = V2), size=5, col=I("#E7AE27")) +
    geom_path(data = as.data.frame(princurve.endoderm.time_add.instances$ss12$Line2$s[
      celllist.endoderm.time_add.list$ss12$Line2.order,]), 
      mapping = aes(x = V1, y = V2), size=5, col=I("#FB867A"))
  #--------------------------
  
  #--  18ss :: time_add_new
  #--------------------------
  p18 = (DimPlot(src.18ss.integrated, 
                 reduction = "umap.rotated", pt.size = 1.5,
                 group.by = "Time") + 
           scale_color_manual(values = c("gray")) + # cluster.endoderm.color.v5) +
           theme_void()+
           theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 legend.position = "none",
                 plot.title = element_blank(),
                 aspect.ratio=1) +
           geom_path(data = as.data.frame(princurve.endoderm.time_add_new.instances$ss18$MG$s[
             celllist.endoderm.time_add_new.list$ss18$MG.order,]), 
             mapping = aes(x = V1, y = V2), size=5, col=I("#9C93D7")))
  #--------------------------
  
  #--  21ss :: time_add_new
  #--------------------------
  p21 = (DimPlot(src.21ss.integrated, 
                 reduction = "umap", pt.size = 1.5,
                 group.by = "Time") + 
           scale_color_manual(values = c("gray")) + # cluster.endoderm.color.v5) +
           theme_void()+
           theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 legend.position = "none",
                 plot.title = element_blank(),
                 aspect.ratio=1) +
           geom_path(data = as.data.frame(princurve.endoderm.time_add_new.instances$ss21$MG$s[
             celllist.endoderm.time_add_new.list$ss21$MG.order,]), 
             mapping = aes(x = V1, y = V2), size=5, col=I("#9C93D7")))
  #--------------------------
  
  #--  24ss :: time_add_new
  #--------------------------
  p24 = (DimPlot(src.24ss.integrated, 
                 reduction = "umap", pt.size = 1.5,
                 group.by = "Time") + 
           scale_color_manual(values = c("gray")) + # cluster.endoderm.color.v5) +
           theme_void()+
           theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 legend.position = "none",
                 plot.title = element_blank(),
                 aspect.ratio=1) +
           geom_path(data = as.data.frame(princurve.endoderm.time_add_new.instances$ss24$MG$s[
             celllist.endoderm.time_add_new.list$ss24$MG.order,]), 
             mapping = aes(x = V1, y = V2), size=5, col=I("#9C93D7")))
  #--------------------------
  
  #---plot-old ------
  #-- 9ss :: time
  #--------------------------
  p9 = (DimPlot(src.9ss.integrated, reduction = "umap", pt.size = 1.5,
                group.by = "Time") + 
          scale_color_manual(values = c("gray")) + # cluster.endoderm.color.v5) +
          theme_void()+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.position = "none",
                plot.title = element_blank(),
                aspect.ratio=1) +
          geom_path(data = as.data.frame(princurve.endoderm.time.instances$ss9$Line1$s[
            celllist.endoderm.time.list$ss9$Line1.order,]), 
            mapping = aes(x = V1, y = V2), size=5, col=I("#87afcc")) +
          geom_path(data = as.data.frame(princurve.endoderm.time.instances$ss9$Line3$s[
            celllist.endoderm.time.list$ss9$Line3.order,]), 
            mapping = aes(x = V1, y = V2), size=5, col=I("#E7AE27")) +
          geom_path(data = as.data.frame(princurve.endoderm.time.instances$ss9$Line2$s[
            celllist.endoderm.time.list$ss9$Line2.order,]), 
            mapping = aes(x = V1, y = V2), size=5, col=I("#FB867A")))
  #--------------------------
  
  #-- 15ss :: time + time_add_new
  #--------------------------
  p15 = (DimPlot(src.15ss.integrated, reduction = "umap.rotated", pt.size = 1.5,
                 group.by = "Time") + 
           scale_color_manual(values = c("gray")) + # cluster.endoderm.color.v5) +
           theme_void()+
           theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 legend.position = "none",
                 plot.title = element_blank(),
                 aspect.ratio=1) +
           geom_path(data = as.data.frame(princurve.endoderm.time.instances$ss15$Line1$s[
             celllist.endoderm.time.list$ss15$Line1.order,]), 
             mapping = aes(x = V1, y = V2), size=5, col=I("#87afcc")) +
           geom_path(data = as.data.frame(princurve.endoderm.time.instances$ss15$Line2$s[
             celllist.endoderm.time.list$ss15$Line2.order,]), 
             mapping = aes(x = V1, y = V2), size=5, col=I("#FB867A")) +
           geom_path(data = as.data.frame(princurve.endoderm.time_add_new.instances$ss15$Line3$s[
             celllist.endoderm.time_add_new.list$ss15$Line3.order,]), 
             mapping = aes(x = V1, y = V2), size=5, col=I("#E7AE27")))
  #--------------------------
  
  #-- 27ss :: time_add_new
  #--------------------------
  p27 = (DimPlot(src.27ss.integrated, reduction = "umap", pt.size = 1.5,
                 group.by = "Time") + 
           scale_color_manual(values = c("gray")) + # cluster.endoderm.color.v5) +
           theme_void()+
           theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 legend.position = "none",
                 plot.title = element_blank(),
                 aspect.ratio=1) +
           geom_path(data = as.data.frame(cbind(
             princurve.endoderm.time_add_new.instances$ss27$MG$s[,1],
             princurve.endoderm.time_add_new.instances$ss27$MG$s[,2])[
               celllist.endoderm.time_add_new.list$ss27$MG.order,]), 
             mapping = aes(x = V1, y = V2), size=5, col=I("#9C93D7")))
  #--------------------------
  
  for(i.p in c("p9","p12","p15","p18","p21","p24","p27")){
    png(filename = paste("trajectory_signal/Pattern_code.new.", i.p, ".png", sep = ""),
        width = 1000,height = 1000,pointsize = 20)
    print(get(i.p))
    dev.off()
  }
}

#-- correct: cell order and reverse
for(i.cellorderlist in c("celllist.endoderm.time_add_new.list")){
  # celllist.endoderm.time_add_new.list$ss15$Line3.order = rev(celllist.endoderm.time_add_new.list$ss12$Line3.order)
  celllist.endoderm.time_add_new.list$ss18$MG.order = rev(celllist.endoderm.time_add_new.list$ss18$MG.order)
  celllist.endoderm.time_add_new.list$ss21$MG.order = rev(celllist.endoderm.time_add_new.list$ss21$MG.order)
  celllist.endoderm.time_add_new.list$ss24$MG.order = rev(celllist.endoderm.time_add_new.list$ss24$MG.order)
  celllist.endoderm.time_add_new.list$ss27$MG.order = rev(celllist.endoderm.time_add_new.list$ss27$MG.order)
  
  src.check = src.27ss.integrated
  cell.check = celllist.endoderm.time_add_new.list$ss27$MG.order
  plot(x = 1:length(cell.check), y = 1:length(cell.check),
       col = cluster.endoderm.color.v5[src.check@meta.data[cell.check,]$cluster.v06.26.re..correct..un])
  rm(src.check); gc()
}

#-- Add-new seurat with plsda (step1)
for(i.time in c(names(pathway.endoderm.time_add_new.list))){
  for(j.type in names(celllist.endoderm.time_add_new.list[[i.time]])){
    if(grepl("order",j.type)==F){next()}
    seurat.plsda.pathway.endoderm.time_add_new.list[[i.time]][[j.type]] =
      seurat_pathway_pre(pathway.matrix = pathway.endoderm.time_add_new.list[[i.time]][,celllist.endoderm.time_add_new.list[[i.time]][[j.type]]], 
                         seurat.ref = get(paste("src.", gsub("ss","",i.time), "ss.integrated", sep = ""))[,celllist.endoderm.time_add_new.list[[i.time]][[j.type]]]) 
  }
}

#-- Add-new seurat with plsda (step2)
for(i.time in names(seurat.plsda.pathway.endoderm.time_add_new.list)){
  for(i.type in names(seurat.plsda.pathway.endoderm.time_add_new.list[[i.time]])){
    if(grepl("order",j.type)==F){next()}
    seurat.plsda.pathway.endoderm.time_add_new.list[[i.time]][[i.type]] = 
      seurat_plsda_pathway(seurat.plsda.pathway.endoderm.time_add_new.list[[i.time]][[i.type]],
                           seurat.plsda.pathway.endoderm.time_add_new.list[[i.time]][[i.type]]$cluster.v06.26.re..correct..un)[["seurat"]]
  }
}

pathway.plsda.pathway.endoderm.time_add_new.list = list()
for(i.type in names(seurat.plsda.pathway.endoderm.time_add_new.list)){
  for(i.name in names(seurat.plsda.pathway.endoderm.time_add_new.list[[i.type]])){
    pathway.plsda.pathway.endoderm.time_add_new.list[[i.type]][[i.name]] =
      plsda_select(seurat.plsda.pathway.endoderm.time_add_new.list[[i.type]][[i.name]], comp = 3, threshold = 30)
  }
}

tree.plsda.pathway.endoderm.time_add_new.list = list()
for(i.name in names(seurat.plsda.pathway.endoderm.time_add_new.list)){
  tree.plsda.pathway.endoderm.time_add_new.list[[i.name]] = list()
  pdf(paste("Milestones/try.time_add_new.", i.name, '.pdf', sep = ""), 10, 10)
  
  for(j.name in names(seurat.plsda.pathway.endoderm.time_add_new.list[[i.name]])){
    if(grepl("order",j.type)==F){next()}
    
    seurat = seurat.plsda.pathway.endoderm.time_add_new.list[[i.name]][[j.name]]
    seurat$cluster.temp = seurat$cluster.v06.26.re..correct..un
    gene_list = pathway.plsda.pathway.endoderm.time_add_new.list[[i.name]][[j.name]]
    
    cellorder_list = celllist.endoderm.time_add_new.list[[i.name]][[j.name]]
    cell_name = intersect(cellorder_list, colnames(seurat))
    
    seurat = ScaleData(seurat[,cell_name], features = rownames(seurat))
    
    data = seurat@assays$Pathway@scale.data[gene_list, cell_name]
    data_re = t(apply(data,1,kernelsmooth))
    rownames(data_re) = rownames(data); colnames(data_re) = colnames(data)
    
    plot(c(1:10), c(1:10))
    text(x=5,y=5,labels=paste(i.type, i.name, j.name, sep = " "), cex = 5)
    
    pathway.heatmap  =
      MyHeatmap(data_re,
                type = "raw",
                ColSideColors = cbind(
                  MyName2Col(seurat@meta.data[cell_name,]$cluster.temp, cluster.endoderm.color.v5)),
                color.palette = colorRampPalette(c("#5aa7dd","#EBEBEB","red"), space="Lab"),
                ColSideColorsSize = 3,
                RowSideColorsSize = 2,
                # Rowv = "none", 
                # Colv = "none",
                return.tree = "row",
                labCol="none", graph = T)
    
    tree.plsda.pathway.endoderm.time_add_new.list[[i.name]][[j.name]] = as.dendrogram(pathway.heatmap) 
  }
  dev.off()
}

gene.tree.plsda.pathway.endoderm.time_add_new.list = list()
for(i.tree.plsda in c("gene.tree.plsda.pathway.endoderm.list.refine")){
  #-- Add-tree
  for(i.type in names(tree.plsda.pathway.endoderm.time_add_new.list)){
    gene.tree.plsda.pathway.endoderm.time_add_new.list[[i.type]] = list()
    for(j.type in names(tree.plsda.pathway.endoderm.time_add_new.list[[i.type]])){
      gene.tree.plsda.pathway.endoderm.time_add_new.list[[i.type]][[j.type]] = NA
    }
  }
  
  #-- SS15: Line3
  #----------------------------------
  gene.tree.plsda.pathway.endoderm.time_add_new.list$ss15$Line3.order = c(
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss15$Line3.order[[2]][[1]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss15$Line3.order[[2]][[1]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss15$Line3.order[[2]][[1]][[1]]),
    rev(labels(tree.plsda.pathway.endoderm.time_add_new.list$ss15$Line3.order[[1]][[2]])),
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss15$Line3.order[[2]][[2]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.time_add_new.list$ss15$Line3.order) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.time_add_new.list$ss15$Line3.order[[2]][[1]][[2]][[2]]),
                    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss15$Line3.order[[2]][[1]][[2]][[1]])))),
    rep(3, length(c(labels(tree.plsda.pathway.endoderm.time_add_new.list$ss15$Line3.order[[2]][[1]][[1]])))),
    rep(5, length(c(labels(tree.plsda.pathway.endoderm.time_add_new.list$ss15$Line3.order[[1]][[2]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.time_add_new.list$ss15$Line3.order[[2]][[2]][[1]])))))
  #----------------------------------
  
  #-- SS18: MG.order
  #----------------------------------
  gene.tree.plsda.pathway.endoderm.time_add_new.list$ss18$MG.order = c(
    setdiff(labels(tree.plsda.pathway.endoderm.time_add_new.list$ss18$MG.order[[2]][[1]]),
            c()),
    
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss18$MG.order[[2]][[2]][[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss18$MG.order[[2]][[2]][[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss18$MG.order[[2]][[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss18$MG.order[[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss18$MG.order[[1]][[1]][[2]]))
  names(gene.tree.plsda.pathway.endoderm.time_add_new.list$ss18$MG.order) = c(
    rep(4, length(c(setdiff(labels(tree.plsda.pathway.endoderm.time_add_new.list$ss18$MG.order[[2]][[1]]),
                            c())))),
    rep(3, length(c(labels(tree.plsda.pathway.endoderm.time_add_new.list$ss18$MG.order[[2]][[2]][[1]][[2]]),
                    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss18$MG.order[[2]][[2]][[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss18$MG.order[[2]][[2]][[2]][[1]])))),
    rep(5, length(c(labels(tree.plsda.pathway.endoderm.time_add_new.list$ss18$MG.order[[1]][[2]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.time_add_new.list$ss18$MG.order[[1]][[1]][[2]])))))
  #----------------------------------
  
  #-- SS21: MG.order
  #----------------------------------
  gene.tree.plsda.pathway.endoderm.time_add_new.list$ss21$MG.order = c(
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss21$MG.order[[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss21$MG.order[[2]][[2]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss21$MG.order[[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss21$MG.order[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss21$MG.order[[1]][[2]][[2]]))
  names(gene.tree.plsda.pathway.endoderm.time_add_new.list$ss21$MG.order) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.time_add_new.list$ss21$MG.order[[2]][[1]])))),
    rep(3, length(c(labels(tree.plsda.pathway.endoderm.time_add_new.list$ss21$MG.order[[2]][[2]][[2]][[2]])))),
    rep(5, length(c(labels(tree.plsda.pathway.endoderm.time_add_new.list$ss21$MG.order[[2]][[2]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.time_add_new.list$ss21$MG.order[[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss21$MG.order[[1]][[2]][[2]])))))
  #----------------------------------
  
  #-- SS24: MG.order
  #----------------------------------
  gene.tree.plsda.pathway.endoderm.time_add_new.list$ss24$MG.order = c(
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss24$MG.order[[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss24$MG.order[[2]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss24$MG.order[[2]][[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss24$MG.order[[1]]))
  
  names(gene.tree.plsda.pathway.endoderm.time_add_new.list$ss24$MG.order) = c(
    rep(4, length(c(setdiff(labels(tree.plsda.pathway.endoderm.time_add_new.list$ss24$MG.order[[2]][[2]][[1]]),c())))),
    rep(3, length(c(labels(tree.plsda.pathway.endoderm.time_add_new.list$ss24$MG.order[[2]][[2]][[2]])))),
    rep(5, length(c(labels(tree.plsda.pathway.endoderm.time_add_new.list$ss24$MG.order[[2]][[1]][[2]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.time_add_new.list$ss24$MG.order[[1]])))))
  #----------------------------------
  
  #-- SS27: MG.order
  #----------------------------------
  gene.tree.plsda.pathway.endoderm.time_add_new.list$ss27$MG.order = c(
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss27$MG.order[[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss27$MG.order[[2]][[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss27$MG.order[[2]][[1]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss27$MG.order[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss27$MG.order[[1]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss27$MG.order[[1]][[2]][[2]]))
  
  names(gene.tree.plsda.pathway.endoderm.time_add_new.list$ss27$MG.order) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.time_add_new.list$ss27$MG.order[[2]][[2]])))),
    rep(3, length(c(labels(tree.plsda.pathway.endoderm.time_add_new.list$ss27$MG.order[[2]][[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss27$MG.order[[2]][[1]][[2]][[2]])))),
    rep(5, length(c(labels(tree.plsda.pathway.endoderm.time_add_new.list$ss27$MG.order[[1]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.time_add_new.list$ss27$MG.order[[1]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.time_add_new.list$ss27$MG.order[[1]][[2]][[2]])))))
  #----------------------------------
  
}

#-- Filter by anova
#--------------------------------------
gene.tree.plsda.pathway.endoderm.time_add_new.list.filter = list()
for(i.type in names(gene.tree.plsda.pathway.endoderm.time_add_new.list)){
  gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]] = list()
  for(j.type in names(gene.tree.plsda.pathway.endoderm.time_add_new.list[[i.type]])){
    if(grepl("order", j.type)==F){next()}
    print(paste(i.type,j.type))
    
    src.plsda = seurat.plsda.pathway.endoderm.time_add_new.list[[i.type]][[j.type]]
    gene.plsda = gene.tree.plsda.pathway.endoderm.time_add_new.list[[i.type]][[j.type]]
    gene.anova = anova.test(tpm = src.plsda@assays$Pathway@data[gene.plsda, colnames(src.plsda)],
                            variable = src.plsda$cluster.v06.26.re..correct..un)
    gene.anova.scale = anova.test(tpm = src.plsda@assays$Pathway@scale.data[gene.plsda, colnames(src.plsda)],
                                  variable = src.plsda$cluster.v06.26.re..correct..un)
    
    gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]][[j.type]][["gene"]] = gene.plsda
    gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]][[j.type]][["p.val.data"]] = gene.anova
    gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]][[j.type]][["p.val.scaledata"]] = gene.anova.scale
  }
}

anova.curve.gene.tree.plsda.pathway.endoderm.time_add_new.list.filter = list()
anova.df.gene.tree.plsda.pathway.endoderm.time_add_new.list.filter = cbind(c(1:50))
for(i.type in names(gene.tree.plsda.pathway.endoderm.time_add_new.list.filter)){
  for(j.type in names(gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]])){
    gene = gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]][[j.type]][["gene"]]
    pvalue.data = gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]][[j.type]][["p.val.data"]]
    num = c()
    for(i in c(1:50)){
      num = c(num, length(gene[pvalue.data<10^(-i)]))
      if(i==1 & length(gene[pvalue.data<10^(-i)])==0){print(paste(i.type,j.type))}
    }
    names(num) = c(1:50)
    anova.curve.gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]][[j.type]] = num 
    anova.df.gene.tree.plsda.pathway.endoderm.time_add_new.list.filter = cbind(
      anova.df.gene.tree.plsda.pathway.endoderm.time_add_new.list.filter, num)
  }
}
for(i.plot in c("gene.tree.plsda.pathway.endoderm.time_add_new.list.filter")){
  df.plot = melt(anova.df.gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[,2:ncol(anova.df.gene.tree.plsda.pathway.endoderm.time_add_new.list.filter)])
  list.index = c()
  for(i in c(1:(ncol(anova.df.gene.tree.plsda.pathway.endoderm.time_add_new.list.filter)-1))){list.index = c(list.index, rep(i,50))}
  df.plot$Var2 = list.index
  
  pdf("Milestones/anova.curve.pathways.pdf",6,6)
  ggplot(data = df.plot) + 
    geom_line(mapping = aes(x=Var1, y=value, color=factor(Var2)))+
    scale_color_manual(values = unique(c(as.character(color.cluster2.new),
                                         as.character(color.time.new)))) +
    theme_classic() +
    geom_vline(xintercept=c(10,20,30),lty=4,col="black",lwd=0.8) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # axis.ticks =  element_blank(),
          # axis.text = element_blank(),
          legend.position = "none",
          # plot.title = element_blank(),
          axis.title.x = element_text(color="black", size=20, face="bold"),
          axis.title.y = element_text(color="black", size=20, face="bold"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.line.x = element_line(linetype=1, color="black", size=1.5),
          axis.line.y = element_line(linetype=1, color="black", size=1.5),
          axis.ticks.x = element_line(linetype=1, color="black", size=1.5),
          axis.ticks.y = element_line(linetype=1, color="black", size=1.5),
          axis.ticks.length  = unit(0.2, "cm"),
          aspect.ratio=1)+
    xlab("-log10(p.value)") + ylab("Detected pathways")
  dev.off()
  
  df.plot.50 = df.plot[(df.plot$value>25)&(df.plot$value<35),]
  hist(df.plot.50$Var1, breaks = 40)
  
  pdf("Milestones/anova.curve.P50.pathways.pdf",6,6)
  ggplot(df.plot.50) + 
    geom_histogram(mapping = aes(x=Var1),bins = 40,colour="black",fill="#eeeeee")+
    theme_classic() +
    geom_vline(xintercept=c(10,20,30),lty=4,col="black",lwd=0.8) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # axis.ticks =  element_blank(),
          # axis.text = element_blank(),
          legend.position = "none",
          # plot.title = element_blank(),
          axis.title.x = element_text(color="black", size=20, face="bold"),
          axis.title.y = element_text(color="black", size=20, face="bold"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.line.x = element_line(linetype=1, color="black", size=1.5),
          axis.line.y = element_line(linetype=1, color="black", size=1.5),
          axis.ticks.x = element_line(linetype=1, color="black", size=1.5),
          axis.ticks.y = element_line(linetype=1, color="black", size=1.5),
          axis.ticks.length  = unit(0.2, "cm"),
          aspect.ratio=1)+
    xlab("-log10(p.value)") + ylab("Frequency")
  dev.off()
}

for(i.type in names(gene.tree.plsda.pathway.endoderm.time_add_new.list.filter)){
  for(j.type in names(gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]])){
    
    gene = gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]][[j.type]][["gene"]]
    pvalue.data = gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]][[j.type]][["p.val.data"]]
    gene.rev = names(gene)
    names(gene.rev) = gene
    
    #-- p10 & p20 (Value range)
    #---------------------------- 
    gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]][[j.type]][["gene.p10"]] = gene[pvalue.data<10^(-10)]
    gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]][[j.type]][["gene.p20"]] = gene[pvalue.data<10^(-20)]
    
    if(length(gene)==0){print(paste(i.type, j.type, k.type));next()}
    if(length(gene[pvalue.data<10^(-10)])/length(gene) < 0.5){
      gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]][[j.type]][["gene.set"]] = "gene.p10"
    }else{
      gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]][[j.type]][["gene.set"]] = "gene.p20"
    }
    #---------------------------- 
    
    #-- p10 with top30 & top50
    #---------------------------- 
    gene.p10.top30 = intersect(gene, names(sort(pvalue.data)[1:min(30, length(gene[pvalue.data<10^(-10)]))]))
    names(gene.p10.top30) = gene.rev[gene.p10.top30]
    
    gene.p10.top50 = intersect(gene, names(sort(pvalue.data)[1:min(50, length(gene[pvalue.data<10^(-10)]))]))
    names(gene.p10.top50) = gene.rev[gene.p10.top50]
    
    gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]][[j.type]][["gene.p10.top30"]] = gene.p10.top30
    gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]][[j.type]][["gene.p10.top50"]] = gene.p10.top50
    #---------------------------- 
    
    #-- p10 with top30 & top50 stratified
    #---------------------------- 
    gene.p10 = gene[pvalue.data<10^(-10)]
    type.gene.p10.top30 = round(table(names(gene.p10)) / sum(table(gene.p10)) * 30, 0)
    type.gene.p10.top50 = round(table(names(gene.p10)) / sum(table(gene.p10)) * 50, 0)
    type.order.list = c(4,3,2,5,7)
    
    gene.p10.top30.stratified = c()
    for(i.num in intersect(type.order.list, names(type.gene.p10.top30))){
      gene.p10.i.num = gene.p10[names(gene.p10)%in%i.num]
      gene.p10.top30.stratified = c(gene.p10.top30.stratified,
                                    names(sort(pvalue.data[gene.p10.i.num])[1:max(1, type.gene.p10.top30[i.num])]))
    }
    names(gene.p10.top30.stratified) = gene.rev[gene.p10.top30.stratified] 
    
    gene.p10.top50.stratified = c()
    for(i.num in intersect(type.order.list, names(type.gene.p10.top50))){
      gene.p10.i.num = gene.p10[names(gene.p10)%in%i.num]
      gene.p10.top50.stratified = c(gene.p10.top50.stratified,
                                    names(sort(pvalue.data[gene.p10.i.num])[1:max(1, type.gene.p10.top50[i.num])]))
    }
    names(gene.p10.top50.stratified) = gene.rev[gene.p10.top50.stratified]
    
    
    gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]][[j.type]][["gene.p10.top30.stratified"]] = gene.p10.top30.stratified
    gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]][[j.type]][["gene.p10.top50.stratified"]] = gene.p10.top50.stratified
    #---------------------------- 
    
  }
}
#--------------------------------------

save(pathway.endoderm.time_add_new.list,
     seurat.plsda.pathway.endoderm.time_add_new.list,
     celllist.endoderm.time_add_new.list,
     princurve.endoderm.time_add_new.instances,
     pathway.plsda.pathway.endoderm.time_add_new.list,
     tree.plsda.pathway.endoderm.time_add_new.list,
     gene.tree.plsda.pathway.endoderm.time_add_new.list,
     gene.tree.plsda.pathway.endoderm.time_add_new.list.filter,
     file = "trajectory_signal/seurat.plsda.pathway.endoderm.time_add_new.list.parameter.Rdata")

#-- fin-plot
for(i.name in names(seurat.plsda.pathway.endoderm.time_add_new.list)){
  pdf(paste("trajectory_signal/try.time_add_new.top30.", i.name, '.pdf', sep = ""), 8, 8)
  for(j.name in names(seurat.plsda.pathway.endoderm.time_add_new.list[[i.name]])){
    if(grepl("order",j.name)==F){next()}
    seurat = seurat.plsda.pathway.endoderm.time_add_new.list[[i.name]][[j.name]]
    seurat$cluster.temp = seurat$cluster.v06.26.re..correct..un
    
    cell_name = intersect(celllist.endoderm.time_add_new.list[[i.name]][[j.name]], colnames(seurat))
    if(i.name %in% c("ss18","ss21","ss24", "ss27")){
      cell_name = intersect(celllist.endoderm.time_add_new.list[[i.name]][[j.name]], 
                            rownames(seurat@meta.data[!seurat$cluster.temp%in%c("Thyroid","FG.5","Pharynx.organ.3"),]))
    }
    
    gene_list = gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.name]][[j.name]][["gene.p10.top30.stratified"]]
    gene_list = gene_list[is.na(gene_list)==F]
    print(paste(i.name, j.name))
    
    data = seurat@assays$Pathway@scale.data[gene_list, cell_name]
    data_re = t(apply(data,1,kernelsmooth))
    rownames(data_re) = rownames(data); colnames(data_re) = colnames(data)
    
    pathway.heatmap  =
      MyHeatmap(data_re,
                type = "raw",
                ColSideColors = cbind(
                  MyName2Col(seurat@meta.data[cell_name,]$cluster.temp, cluster.endoderm.color.v5)
                ),
                RowSideColors = t(cbind(
                  MyName2Col(names(gene_list), colors.geneset)
                )),
                color.palette = colorRampPalette(c("#5aa7dd","#EBEBEB","red"), space="Lab"),
                ColSideColorsSize = 2,
                RowSideColorsSize = 2,
                Rowv = "none", 
                Colv = "none",
                # return.tree = "row",
                margins = c(10,10),
                labCol="none", graph = T)
  }
  dev.off()
}

gene.tree.plsda.pathway.endoderm.time_add_new.list.filter$ss21$MG.order$gene.p10.top30.stratified.compress = c(
  "FGFR_compress",
  gene.tree.plsda.pathway.endoderm.time_add_new.list.filter$ss21$MG.order$gene.p10.top30.stratified[
    names(gene.tree.plsda.pathway.endoderm.time_add_new.list.filter$ss21$MG.order$gene.p10.top30.stratified)==3],
  gene.tree.plsda.pathway.endoderm.time_add_new.list.filter$ss21$MG.order$gene.p10.top30.stratified[
    names(gene.tree.plsda.pathway.endoderm.time_add_new.list.filter$ss21$MG.order$gene.p10.top30.stratified)==5],
  gene.tree.plsda.pathway.endoderm.time_add_new.list.filter$ss21$MG.order$gene.p10.top30.stratified[
    names(gene.tree.plsda.pathway.endoderm.time_add_new.list.filter$ss21$MG.order$gene.p10.top30.stratified)==7]
)
names(gene.tree.plsda.pathway.endoderm.time_add_new.list.filter$ss21$MG.order$gene.p10.top30.stratified.compress) = c(
  c(4),
  rep(3, length(c(gene.tree.plsda.pathway.endoderm.time_add_new.list.filter$ss21$MG.order$gene.p10.top30.stratified[
    names(gene.tree.plsda.pathway.endoderm.time_add_new.list.filter$ss21$MG.order$gene.p10.top30.stratified)==3]))),
  rep(5, length(c(gene.tree.plsda.pathway.endoderm.time_add_new.list.filter$ss21$MG.order$gene.p10.top30.stratified[
    names(gene.tree.plsda.pathway.endoderm.time_add_new.list.filter$ss21$MG.order$gene.p10.top30.stratified)==5]))),
  rep(7, length(c(gene.tree.plsda.pathway.endoderm.time_add_new.list.filter$ss21$MG.order$gene.p10.top30.stratified[
    names(gene.tree.plsda.pathway.endoderm.time_add_new.list.filter$ss21$MG.order$gene.p10.top30.stratified)==7])))
)

for(i.name in c("ss21")){
  pdf(paste("trajectory_signal/try.time_add_new.top30.compress.", i.name, '.pdf', sep = ""), 8, 8)
  for(j.name in c("MG.order")){
    if(grepl("order",j.name)==F){next()}
    
    seurat = seurat.plsda.pathway.endoderm.time_add_new.list[[i.name]][[j.name]]
    seurat$cluster.temp = seurat$cluster.v06.26.re..correct..un
    
    cell_name = intersect(celllist.endoderm.time_add_new.list[[i.name]][[j.name]], colnames(seurat))
    if(i.name %in% c("ss18","ss21","ss24","ss27")){
      cell_name = intersect(celllist.endoderm.time_add_new.list[[i.name]][[j.name]], 
                            rownames(seurat@meta.data[!seurat$cluster.temp%in%c("Thyroid","FG.5","Pharynx.organ.3"),]))
    }
    
    gene_list = gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.name]][[j.name]][["gene.p10.top30.stratified.compress"]]
    gene_list = gene_list[is.na(gene_list)==F]
    print(paste(i.name, j.name))
    
    data = seurat@assays$Pathway@scale.data[
      gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.name]][[j.name]][["gene.p10.top30.stratified"]],
      cell_name]
    data = rbind(data[intersect(gene_list, rownames(data)),],
                 apply(data[
                   gene.tree.plsda.pathway.endoderm.time_add_new.list.filter$ss21$MG.order$gene.p10.top30.stratified[
                     names(gene.tree.plsda.pathway.endoderm.time_add_new.list.filter$ss21$MG.order$gene.p10.top30.stratified)==4],
                 ], 2, mean))
    data_re = t(apply(data,1,kernelsmooth))
    rownames(data_re) = c(intersect(gene_list, rownames(data)), "FGFR_compress") 
    colnames(data_re) = colnames(data)
    
    pathway.heatmap  =
      MyHeatmap(data_re[gene_list,],
                type = "raw",
                ColSideColors = cbind(
                  MyName2Col(seurat@meta.data[cell_name,]$cluster.temp, cluster.endoderm.color.v5)
                ),
                RowSideColors = t(cbind(
                  MyName2Col(names(gene_list), colors.geneset)
                )),
                color.palette = colorRampPalette(c("#5aa7dd","#EBEBEB","red"), space="Lab"),
                ColSideColorsSize = 2,
                RowSideColorsSize = 2,
                Rowv = "none", 
                Colv = "none",
                # return.tree = "row",
                margins = c(10,10),
                labCol="none", graph = T)
    
  }
  dev.off()
}


#-- top 30
df.summary.plsda.pathway.endoderm.time_add_new.list.filter = rbind(c("Types","Subtypes","Clusters","Pathways"))
colnames(df.summary.plsda.pathway.endoderm.time_add_new.list.filter) = c("Types","Subtypes","Clusters","Pathways")
for(i.type in names(gene.tree.plsda.pathway.endoderm.time_add_new.list.filter)){
  for(j.type in names(gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]])){
    gene.temp = gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]][[j.type]]$gene.p20.top30.stratified
    gene.temp = gene.temp[is.na(gene.temp)==F]
    df.temp = cbind.data.frame(
      rep(i.type, length(gene.temp)),
      rep(j.type, length(gene.temp)),
      names(gene.temp), gene.temp)
    colnames(df.temp) = c("Types","Subtypes","Clusters","Pathways")
    df.summary.plsda.pathway.endoderm.time_add_new.list.filter = rbind(df.summary.plsda.pathway.endoderm.time_add_new.list.filter, df.temp)
  }
}
df.summary.plsda.pathway.endoderm.time_add_new.list.filter = as.data.frame(df.summary.plsda.pathway.endoderm.time_add_new.list.filter)

write.csv(df.summary.plsda.pathway.endoderm.time_add_new.list.filter,
          file = "trajectory_signal/df.summary.plsda.pathway.endoderm.time_add_new.list.filter.csv")


#-- top 30
df.summary.plsda.pathway.endoderm.time_add_new.list.filter = rbind(c("Types","Subtypes","Clusters","Pathways"))
colnames(df.summary.plsda.pathway.endoderm.time_add_new.list.filter) = c("Types","Subtypes","Clusters","Pathways")
for(i.type in names(gene.tree.plsda.pathway.endoderm.time_add_new.list.filter)){
  for(j.type in names(gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]])){
    gene.temp = gene.tree.plsda.pathway.endoderm.time_add_new.list.filter[[i.type]][[j.type]]$gene.p10.top30.stratified
    gene.temp = gene.temp[is.na(gene.temp)==F]
    df.temp = cbind.data.frame(
      rep(i.type, length(gene.temp)),
      rep(j.type, length(gene.temp)),
      names(gene.temp), gene.temp)
    colnames(df.temp) = c("Types","Subtypes","Clusters","Pathways")
    df.summary.plsda.pathway.endoderm.time_add_new.list.filter = rbind(df.summary.plsda.pathway.endoderm.time_add_new.list.filter, df.temp)
  }
}
df.summary.plsda.pathway.endoderm.time_add_new.list.filter = as.data.frame(df.summary.plsda.pathway.endoderm.time_add_new.list.filter)

write.csv(df.summary.plsda.pathway.endoderm.time_add_new.list.filter,
          file = "trajectory_signal/df.summary.plsda.pathway.endoderm.time_add_new.list.filter.top50.csv")
#==========================================================


#==========================================================
#   use 3d-umap to correct cell-order::  9SS 12SS 15SS 
#==========================================================
src.12ss.integrated.3d = RunUMAP(src.12ss.integrated, reduction = "mnnpca", 
                                 dims = 1:30, min.dist = 0.3,
                                 n.neighbors = 100, n.components = 3)

src.15ss.integrated.3d = RunUMAP(src.15ss.integrated, reduction = "mnnpca", 
                                 dims = 1:30, min.dist = 0.3,
                                 n.neighbors = 100, n.components = 3)

save(src.12ss.integrated.3d, src.15ss.integrated.3d,
     file = "trajectory_signal/umap.integrated.3d.Rdata")

rm(src.12ss.integrated.3d,
   src.15ss.integrated.3d)

#---  ss12
for(i.rotated in c("src.12ss.integrated")){
  src.12ss.integrated.3d = RunUMAP(src.12ss.integrated.3d, reduction = "mnnpca", 
                                   dims = 1:30, n.neighbors = 100, n.components = 3)
  
  data = cbind(src.12ss.integrated.3d@meta.data,
               src.12ss.integrated.3d@reductions$umap@cell.embeddings)
  colnames(data) = gsub("umap","Coord",colnames(data))
  
  open3d() 
  view = par3d(family="arial", cex=20, font=1)
  par3d(userMatrix = view)
  
  plot3d(data[,c("Coord_1","Coord_2","Coord_3")],
         col=cluster.endoderm.color.v5[data$cluster.v06.26.re..correct..un],
         lwd = 1, type="p", size=5, axes= F,
         xlab = "", ylab = "", zlab = "",
         ticktype = "detailed")
  
  axes3d(
    col = "black",
    lwd = 1,
    family = "serif", font = 1,
    marklen = 0,
    marklen.rel = F)
  #axes3d(c("x","y","z")) 
  grid3d(c("x","y","z"),n = 8)
  #ticktype = "detailed")
  close3d()
  
  view_data = as.matrix(dput(par3d("userMatrix")))
  view.umap.src.12ss.integrated.3d = view_data
  src.12ss.integrated.3d@reductions$umap.rotated = src.12ss.integrated.3d@reductions$umap
  data.temp =  view.umap.src.12ss.integrated.3d[1:3,1:3] %*% 
    t(as.matrix(src.12ss.integrated.3d@reductions$umap@cell.embeddings)); data.temp = t(data.temp)
  colnames(data.temp) = paste("UMAP_", c(1:3), sep = "")
  src.12ss.integrated.3d@reductions$umap.rotated@cell.embeddings = data.temp
  
  umap.rotated.src.12ss.integrated.3d = src.12ss.integrated.3d@reductions$umap.rotated
}
src.12ss.integrated@reductions$umap.3d = umap.rotated.src.12ss.integrated.3d
src.12ss.integrated@reductions$umap.3d@key = "UMAP_"
src.12ss.integrated@reductions$umap.3d@cell.embeddings = cbind(
  -src.12ss.integrated@reductions$umap.3d@cell.embeddings[,1],
  src.12ss.integrated@reductions$umap.3d@cell.embeddings[,c(2,3)])
colnames(src.12ss.integrated@reductions$umap.3d@cell.embeddings) = paste("UMAP_", 1:3, sep = "")

DimPlot(src.12ss.integrated, reduction = "umap.3d", group.by = "Time", dims = c(1,2))

princurve.endoderm.time_add_new.instances$ss12$Line1 = princurve::principal_curve(
  x = cbind(src.12ss.integrated@reductions$umap.3d@cell.embeddings[,1],
            src.12ss.integrated@reductions$umap.3d@cell.embeddings[,2])[
              celllist.endoderm.time_add.list$ss12$Line1,],
  smoother = "smooth.spline")
princurve.endoderm.time_add_new.instances$ss12$Line2 = princurve::principal_curve(
  x = cbind(src.12ss.integrated@reductions$umap.3d@cell.embeddings[,1],
            src.12ss.integrated@reductions$umap.3d@cell.embeddings[,2])[
              celllist.endoderm.time_add.list$ss12$Line2,],
  smoother = "smooth.spline")
princurve.endoderm.time_add_new.instances$ss12$Line3 = princurve::principal_curve(
  x = cbind(src.12ss.integrated@reductions$umap.3d@cell.embeddings[,1],
            src.12ss.integrated@reductions$umap.3d@cell.embeddings[,2])[
              celllist.endoderm.time_add.list$ss12$Line3,],
  smoother = "smooth.spline")


p12 = DimPlot(src.12ss.integrated, reduction = "umap.3d", pt.size = 1.5,
              dims = c(1,2), group.by = "Time") + 
  scale_color_manual(values = c("gray")) + # cluster.endoderm.color.v5) +
  theme_void()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_blank(),
        aspect.ratio=1) +
  geom_path(data = as.data.frame(princurve.endoderm.time_add_new.instances$ss12$Line1$s[
    order(princurve.endoderm.time_add_new.instances$ss12$Line1$lambda),]), 
    mapping = aes(x = V1, y = V2), size=5, col=I("#87afcc")) +
  geom_path(data = as.data.frame(princurve.endoderm.time_add_new.instances$ss12$Line3$s[
    order(princurve.endoderm.time_add_new.instances$ss12$Line3$lambda),]), 
    mapping = aes(x = V1, y = V2), size=5, col=I("#E7AE27")) +
  geom_path(data = as.data.frame(princurve.endoderm.time_add_new.instances$ss12$Line2$s[
    order(princurve.endoderm.time_add_new.instances$ss12$Line2$lambda),]), 
    mapping = aes(x = V1, y = V2), size=5, col=I("#FB867A"))

png(filename = paste("trajectory_signal/Pattern_code.new.", "ss12", ".png", sep = ""),
    width = 1000,height = 1000,pointsize = 20)
print(p12)
dev.off()



#---  ss15
for(i.rotated in c("src.12ss.integrated")){
  src.15ss.integrated.3d = RunUMAP(src.15ss.integrated.3d, reduction = "mnnpca", 
                                   dims = 1:30, n.neighbors = 100, n.components = 3)
  
  data = cbind(src.15ss.integrated.3d@meta.data,
               src.15ss.integrated.3d@reductions$umap@cell.embeddings)
  colnames(data) = gsub("umap","Coord",colnames(data))
  
  open3d() 
  view = par3d(family="arial", cex=20, font=1)
  par3d(userMatrix = view)
  
  plot3d(data[,c("Coord_1","Coord_2","Coord_3")],
         col=cluster.endoderm.color.v5[data$cluster.v06.26.re..correct..un],
         lwd = 1, type="p", size=5, axes= F,
         xlab = "", ylab = "", zlab = "",
         ticktype = "detailed")
  
  axes3d(
    col = "black",#
    family = "serif", font = 1,
    marklen = 0,
    marklen.rel = F)
  #axes3d(c("x","y","z")) 
  grid3d(c("x","y","z"),n = 8)
  # ticktype = "detailed")
  close3d()
  
  view_data = as.matrix(dput(par3d("userMatrix")))
  view.umap.src.15ss.integrated.3d = view_data
  src.15ss.integrated.3d@reductions$umap.rotated = src.15ss.integrated.3d@reductions$umap
  data.temp =  view.umap.src.15ss.integrated.3d[1:3,1:3] %*% 
    t(as.matrix(src.15ss.integrated.3d@reductions$umap@cell.embeddings)); data.temp = t(data.temp)
  colnames(data.temp) = paste("UMAP_", c(1:3), sep = "")
  src.15ss.integrated.3d@reductions$umap.rotated@cell.embeddings = data.temp
  
  umap.rotated.src.15ss.integrated.3d = src.15ss.integrated.3d@reductions$umap.rotated
  save(view.umap.src.15ss.integrated.3d,
       umap.rotated.src.15ss.integrated.3d,
       file = "Simulation/umap.rotated.src.15ss.integrated.3d.summary.Rdata")
}


src.15ss.integrated@reductions$umap.3d = umap.rotated.src.15ss.integrated.3d
src.15ss.integrated@reductions$umap.3d@key = "UMAP_"
colnames(src.15ss.integrated@reductions$umap.3d@cell.embeddings) = paste("UMAP_", 1:3, sep = "")

DimPlot(src.15ss.integrated, reduction = "umap.3d", group.by = "Time", dims = c(1,2))

princurve.endoderm.time_add_new.instances$ss15$Line1 = princurve::principal_curve(
  x = cbind(src.15ss.integrated@reductions$umap.3d@cell.embeddings[,1],
            src.15ss.integrated@reductions$umap.3d@cell.embeddings[,2])[
              setdiff(celllist.endoderm.time.list$ss15$Line1.order,
                      rownames(src.15ss.integrated@reductions$umap.3d@cell.embeddings[
                        src.15ss.integrated@reductions$umap.3d@cell.embeddings[,1]>2,])),],
  smoother = "smooth.spline")
princurve.endoderm.time_add_new.instances$ss15$Line2 = princurve::principal_curve(
  x = cbind(src.15ss.integrated@reductions$umap.3d@cell.embeddings[,1],
            src.15ss.integrated@reductions$umap.3d@cell.embeddings[,2])[
              setdiff(celllist.endoderm.time.list$ss15$Line2.order,
                      rownames(src.15ss.integrated@reductions$umap.3d@cell.embeddings[
                        src.15ss.integrated@reductions$umap.3d@cell.embeddings[,2]>3.5,])),],
  smoother = "smooth.spline")
princurve.endoderm.time_add_new.instances$ss15$Line3 = princurve::principal_curve(
  x = cbind(src.15ss.integrated@reductions$umap.3d@cell.embeddings[,1],
            src.15ss.integrated@reductions$umap.3d@cell.embeddings[,2])[
              setdiff(celllist.endoderm.time_add_new.list$ss15$Line3.order,
                      rownames(src.15ss.integrated@reductions$umap.3d@cell.embeddings[
                        src.15ss.integrated@reductions$umap.3d@cell.embeddings[,1]>2|
                          src.15ss.integrated@reductions$umap.3d@cell.embeddings[,2]<1,])),],
  smoother = "smooth.spline")


p15 = DimPlot(src.15ss.integrated, reduction = "umap.3d", pt.size = 1.5,
              dims = c(1,2), group.by = "Time") + 
  scale_color_manual(values = c("gray")) + # cluster.endoderm.color.v5) +
  theme_void()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_blank(),
        aspect.ratio=1) +
  geom_path(data = as.data.frame(princurve.endoderm.time_add_new.instances$ss15$Line1$s[
    order(princurve.endoderm.time_add_new.instances$ss15$Line1$lambda),]), 
    mapping = aes(x = V1, y = V2), size=5, col=I("#87afcc")) +
  geom_path(data = as.data.frame(princurve.endoderm.time_add_new.instances$ss15$Line3$s[
    order(princurve.endoderm.time_add_new.instances$ss15$Line3$lambda),]), 
    mapping = aes(x = V1, y = V2), size=5, col=I("#E7AE27")) +
  geom_path(data = as.data.frame(princurve.endoderm.time_add_new.instances$ss15$Line2$s[
    order(princurve.endoderm.time_add_new.instances$ss15$Line2$lambda),]), 
    mapping = aes(x = V1, y = V2), size=5, col=I("#FB867A"))
p15
png(filename = paste("trajectory_signal/Pattern_code.new.", "ss15", ".png", sep = ""),
    width = 1000,height = 1000,pointsize = 20)
print(p15)
dev.off()

save(src.12ss.integrated,  file = "trajectory_signal/src.12ss.integrated.Rdata")
save(src.15ss.integrated,  file = "trajectory_signal/src.15ss.integrated.Rdata")


#---  ss9
src.9ss.integrated.3d = RunUMAP(src.9ss.integrated, reduction = "pca", 
                                dims = 1:30, min.dist = 0.3,
                                n.neighbors = 100, n.components = 3)
save(src.9ss.integrated.3d, file = "trajectory_signal/umap.integrated.3d.Rdata")
rm(src.9ss.integrated.3d)

for(i.rotated in c("src.9ss.integrated")){
  src.9ss.integrated.3d = RunUMAP(src.9ss.integrated.3d, reduction = "mnnpca", 
                                  dims = 1:30, n.neighbors = 100, n.components = 3)
  
  data = cbind(src.9ss.integrated.3d@meta.data,
               src.9ss.integrated.3d@reductions$umap@cell.embeddings)
  colnames(data) = gsub("UMAP","Coord",colnames(data))
  
  open3d() 
  view = par3d(family="arial", cex=20, font=1)
  par3d(userMatrix = view)
  
  plot3d(data[,c("Coord_1","Coord_2","Coord_3")],
         col=cluster.endoderm.color.v5[data$cluster.v06.26.re..correct..un],
         lwd = 1, type="p", size=5, axes= F,
         xlab = "", ylab = "", zlab = "",
         ticktype = "detailed")
  
  axes3d(
    col = "black",
    lwd = 1,# 
    family = "serif", font = 1,
    marklen = 0,
    marklen.rel = F)
  #axes3d(c("x","y","z")) # 
  grid3d(c("x","y","z"),n = 8) #
  # ticktype = "detailed")
  close3d()
  view_data = as.matrix(dput(par3d("userMatrix")))
  view.umap.src.9ss.integrated.3d = view_data
  src.9ss.integrated.3d@reductions$umap.rotated = src.9ss.integrated.3d@reductions$umap
  data.temp =  view.umap.src.9ss.integrated.3d[1:3,1:3] %*% 
    t(as.matrix(src.9ss.integrated.3d@reductions$umap@cell.embeddings)); data.temp = t(data.temp)
  colnames(data.temp) = paste("UMAP_", c(1:3), sep = "")
  src.9ss.integrated.3d@reductions$umap.rotated@cell.embeddings = data.temp
  umap.rotated.src.9ss.integrated.3d = src.9ss.integrated.3d@reductions$umap.rotated
  DimPlot(src.9ss.integrated.3d, reduction = "umap.rotated", pt.size = 2.5)+
    theme(aspect.ratio = 1)
  
  save(view.umap.src.9ss.integrated.3d,
       umap.rotated.src.9ss.integrated.3d,
       file = "Simulation/umap.rotated.src.9ss.integrated.3d.summary.Rdata")
}


src.9ss.integrated@reductions$umap.3d = umap.rotated.src.9ss.integrated.3d
src.9ss.integrated@reductions$umap.3d@key = "UMAP_"
colnames(src.9ss.integrated@reductions$umap.3d@cell.embeddings) = paste("UMAP_", 1:3, sep = "")
DimPlot(src.9ss.integrated, reduction = "umap.3d", group.by = "Time", dims = c(1,2))

princurve.endoderm.time_add_new.instances$ss9$Line1 = princurve::principal_curve(
  x = cbind(src.9ss.integrated@reductions$umap.3d@cell.embeddings[,1],
            src.9ss.integrated@reductions$umap.3d@cell.embeddings[,2])[
              celllist.endoderm.time.list$ss9$Line1,],
  smoother = "smooth.spline")
princurve.endoderm.time_add_new.instances$ss9$Line2 = princurve::principal_curve(
  x = cbind(src.9ss.integrated@reductions$umap.3d@cell.embeddings[,1],
            src.9ss.integrated@reductions$umap.3d@cell.embeddings[,2])[
              celllist.endoderm.time.list$ss9$Line2,],
  smoother = "smooth.spline")
princurve.endoderm.time_add_new.instances$ss9$Line3 = princurve::principal_curve(
  x = cbind(src.9ss.integrated@reductions$umap.3d@cell.embeddings[,1],
            src.9ss.integrated@reductions$umap.3d@cell.embeddings[,2])[
              celllist.endoderm.time.list$ss9$Line3,],
  smoother = "smooth.spline")


p9 = DimPlot(src.9ss.integrated, reduction = "umap.3d", pt.size = 1.5,
             dims = c(1,2), group.by = "Time") + 
  scale_color_manual(values = c("gray")) + # cluster.endoderm.color.v5) +
  theme_void()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_blank(),
        aspect.ratio=1) +
  geom_path(data = as.data.frame(princurve.endoderm.time_add_new.instances$ss9$Line1$s[
    order(princurve.endoderm.time_add_new.instances$ss9$Line1$lambda),]), 
    mapping = aes(x = V1, y = V2), size=5, col=I("#87afcc")) +
  geom_path(data = as.data.frame(princurve.endoderm.time_add_new.instances$ss9$Line3$s[
    order(princurve.endoderm.time_add_new.instances$ss9$Line3$lambda),]), 
    mapping = aes(x = V1, y = V2), size=5, col=I("#E7AE27")) +
  geom_path(data = as.data.frame(princurve.endoderm.time_add_new.instances$ss9$Line2$s[
    order(princurve.endoderm.time_add_new.instances$ss9$Line2$lambda),]), 
    mapping = aes(x = V1, y = V2), size=5, col=I("#FB867A"))

png(filename = paste("trajectory_signal/Pattern_code.new.", "ss9", ".png", sep = ""),
    width = 1000,height = 1000,pointsize = 20)
print(p9)
dev.off()


celllist.endoderm.time.list$ss9$Line1.order.refine = rownames(
  princurve.endoderm.time_add_new.instances$ss9$Line1$s[
    order(princurve.endoderm.time_add_new.instances$ss9$Line1$lambda),])
celllist.endoderm.time.list$ss9$Line2.order.refine = rownames(
  princurve.endoderm.time_add_new.instances$ss9$Line2$s[
    order(princurve.endoderm.time_add_new.instances$ss9$Line2$lambda),])
celllist.endoderm.time.list$ss9$Line3.order.refine = rownames(
  princurve.endoderm.time_add_new.instances$ss9$Line3$s[
    order(princurve.endoderm.time_add_new.instances$ss9$Line3$lambda),])

gene.tree.plsda.pathway.endoderm.time.list.filter$ss9$Line1.order.refine = 
  gene.tree.plsda.pathway.endoderm.time.list.filter$ss9$Line1
gene.tree.plsda.pathway.endoderm.time.list.filter$ss9$Line2.order.refine = 
  gene.tree.plsda.pathway.endoderm.time.list.filter$ss9$Line2
gene.tree.plsda.pathway.endoderm.time.list.filter$ss9$Line3.order.refine = 
  gene.tree.plsda.pathway.endoderm.time.list.filter$ss9$Line3


for(i.name in names(seurat.plsda.pathway.endoderm.time.list)){
  if(i.name == "ss9"){}else{next()}
  
  pdf(paste("trajectory_signal/try.order.refine.", i.name, '.pdf', sep = ""), 8, 8)
  for(j.name in names(seurat.plsda.pathway.endoderm.time.list[[i.name]])){
    if(j.name %in% c("Line1","Line2","Line3")){}else{next()}
    seurat = seurat.plsda.pathway.endoderm.time.list[[i.name]][[j.name]]
    seurat$cluster.temp = seurat$cluster.v06.26.re..correct..un
    
    if(paste(j.name,".order.refine",sep = "") %in% names(gene.tree.plsda.pathway.endoderm.time.list.filter[[i.name]])){
      # gene_list = gene.tree.plsda.pathway.endoderm.time.list[[i.name]][[j.name]]
      gene_list = gene.tree.plsda.pathway.endoderm.time.list.filter[[i.name]][[paste(j.name,".order.refine",sep = "")]][["gene.p20.top30.stratified"]]
      gene_list = gene_list[is.na(gene_list)==F]
    }else{next()}
    
    print(paste(i.name, j.name))
    cell_name = celllist.endoderm.time.list[[i.name]][[paste(j.name,".order.refine",sep = "")]]
    
    seurat = ScaleData(seurat[,cell_name], features = rownames(seurat))
    data = seurat@assays$Pathway@scale.data[gene_list, cell_name]
    data_re = t(apply(data,1,kernelsmooth))
    rownames(data_re) = rownames(data); colnames(data_re) = colnames(data)
    
    plot(c(1:10), c(1:10))
    text(x=5,y=5,labels=paste(i.type, i.name, j.name, sep = " "), cex = 5)
    
    pathway.heatmap  =
      MyHeatmap(data_re,
                type = "raw",
                ColSideColors = cbind(MyName2Col(seurat@meta.data[cell_name,]$cluster.temp, cluster.endoderm.color.v5)),
                RowSideColors = t(cbind(MyName2Col(names(gene_list), colors.geneset))),
                color.palette = colorRampPalette(c("#5aa7dd","#EBEBEB","red"), space="Lab"),
                ColSideColorsSize = 2,
                RowSideColorsSize = 2,
                Rowv = "none", 
                Colv = "none",
                # return.tree = "row",
                margins = c(10,10),
                labCol="none", graph = T)
  }
  dev.off()
}
#==========================================================


#===============================================================================
