#===============================================================================
#>> 8.optimal transport
#>>>> R-part
#===============================================================================

#==========================================================
#>>>> (8.1) to loom
#==========================================================
library(loomR)
gi=read.csv("~/genome/mm10/mm10.gene.inf.merge.10X.v1.csv",stringsAsFactors = F)
rownames(gi)=gi$Symbol2
gi$Accession=gi$Symbol2
gi$Gene=gi$Symbol2
View(gi)

loomR::create(filename = "endoderm.loom",data = src.endoderm@assays$mnnRNA@data)
lfile=loomR::connect("endoderm.loom",mode = "r+")
lfile$add.row.attribute(gi[rownames(src.endoderm),
                           c("EnsemblGeneID","Accession","Gene","Symbol2","GeneFullName","Synonyms",
                             "GeneLength","ChromosomeName","Strand","GeneType","TranscriptCount",
                             "GCContent","PlasmaMembrane","Nucleus","CellCycle",
                             "TF","EpiFactor")],overwrite = T)

metadata=src.endoderm@meta.data
metadata[is.na(metadata)]="/"
lfile$add.col.attribute(metadata,overwrite = T)
lfile$close_all()
DefaultAssay(src.endoderm)="RNA"
src.endoderm[["mnnRNA"]]=NULL
plot3d(src.endoderm@reductions$umap@cell.embeddings[,1],
       src.endoderm@reductions$umap@cell.embeddings[,2],
       src.endoderm@reductions$umap@cell.embeddings[,3],
       col = MyName2Col(src.endoderm$Time,
                        time.color),
       size=3,
       xlab = "", ylab = "", zlab = "",
       box=F,
       axes = FALSE)

save(mnnlogtpm.src.12ss.integrated,
     mnnlogtpm.src.15ss.integrated,
     mnnlogtpm.src.18ss.integrated,
     mnnlogtpm.src.21ss.integrated,
     mnnlogtpm.src.24ss.integrated,
     mnnlogtpm.src.27ss.integrated,
     file = "mnncorrected.logtpm.RData")

rm(mnnlogtpm.src.12ss.integrated,
   mnnlogtpm.src.15ss.integrated,
   mnnlogtpm.src.18ss.integrated,
   mnnlogtpm.src.21ss.integrated,
   mnnlogtpm.src.24ss.integrated,
   mnnlogtpm.src.27ss.integrated)
save(src.endoderm,src.endoderm.selectgene,
     src.endoderm.rowtree,mnnlogtpm,
     file = "endoderm.merge.RData")

MyWriteTable(src.endoderm.selectgene,
             "endoderm.merge.selectgene.tab",
             col.names = F)

MyWriteTable(c(labels(src.endoderm.rowtree[[2]][[1]]),
               labels(src.endoderm.rowtree[[2]][[2]][[2]][[2]])),
             "endoderm.merge.cellcycle.gene.tab",
             col.names = F)

endoderm.merge.cellcycle.score= colMeans(
  ScaleData(src.endoderm,
            features = c(labels(src.endoderm.rowtree[[2]][[1]]),
                         labels(src.endoderm.rowtree[[2]][[2]][[2]][[2]])))@assays$RNA@scale.data)
MyWriteTable(data.frame(Cell.cycle=endoderm.merge.cellcycle.score),
             "endoderm.merge.cellcycle.score.tab",
             row.names = T,col.names = T)
#==========================================================

#==========================================================
#>>>> (8.2) calculate OT transition matirx using wot in [Python] 
#==========================================================
#==========================================================

#==========================================================
#>>>> (8.3) OT prediction
#========================================================== 

#>> Read CSV
#---------------------------------------
library(readr)
trans_matrix_9_12 <- read_csv("~/Bioinformatic/project_20221224_endoderm.refine/figure_v6.05/figure.v08.07/ot_refine/OT new/trans_matrix_9_12.csv", col_names = F)
trans_matrix_12_15 <- read_csv("~/Bioinformatic/project_20221224_endoderm.refine/figure_v6.05/figure.v08.07/ot_refine/OT new/trans_matrix_12_15.csv", col_names = F)
trans_matrix_15_18 <- read_csv("~/Bioinformatic/project_20221224_endoderm.refine/figure_v6.05/figure.v08.07/ot_refine/OT new/trans_matrix_15_18.csv", col_names = F)
trans_matrix_18_21 <- read_csv("~/Bioinformatic/project_20221224_endoderm.refine/figure_v6.05/figure.v08.07/ot_refine/OT new/trans_matrix_18_21.csv", col_names = F)
trans_matrix_21_24 <- read_csv("~/Bioinformatic/project_20221224_endoderm.refine/figure_v6.05/figure.v08.07/ot_refine/OT new/trans_matrix_21_24.csv", col_names = F)
trans_matrix_24_27 <- read_csv("~/Bioinformatic/project_20221224_endoderm.refine/figure_v6.05/figure.v08.07/ot_refine/OT new/trans_matrix_24_27.csv", col_names = F)

trans_list = c("trans_matrix_9_12","trans_matrix_12_15","trans_matrix_15_18",
               "trans_matrix_18_21","trans_matrix_21_24","trans_matrix_24_27")
cell_count = 0
for( mat_re in trans_list){
  cell_count = cell_count +nrow(get(mat_re))
  if(mat_re == "trans_matrix_24_27"){
    cell_count = cell_count + ncol(get(mat_re))
  }
}
print(cell_count)
#---------------------------------------

#>> Set OT function
#---------------------------------------
Softmax = function(exp, filter = F, filter_value = 0.75, 
                   filter_mean = F, plot = F,
                   out_filter_value =F){
  ncols  = ncol(exp)
  nrows  = nrow(exp)
  outmatrix = matrix(0, ncol = ncols, nrow = nrows)
  threshold_matrix = matrix(0, ncol = nrows, nrow = 1)
  
  for(i in c(1:ncols)){
    outmatrix[,i] = sapply(exp[,i], FUN = exp) / sum(sapply(exp[,i], FUN = exp))
  }
  
  if(filter == T){
    for(j in c(1:nrows)){
      x = outmatrix[j,]
      y = ecdf(x)
      z = summary(y)
      
      #Set threshold value
      ##########################
      if(filter_value==0.50){
        threshold_value = z[3]}
      else if(filter_value>0.50&filter_value<0.75){
        test_list = seq(z[3],z[5],by=0.0001)
        filter_list = y(test_list)
        filter_list = abs(filter_list-filter_value)
        threshold_value = 
          test_list[which(filter_list==min(filter_list), arr.in=TRUE)]
        if(length(threshold_value)!=1){
          threshold_value = threshold_value[1]
        }
      }
      else if(filter_value==0.75){threshold_value = z[5]}
      else if(filter_value>0.75&filter_value<1){
        test_list = seq(z[5],z[6],by=0.0001)
        filter_list = y(test_list)
        filter_list = abs(filter_list-filter_value)
        threshold_value = 
          test_list[which(filter_list==min(filter_list), arr.in=TRUE)]
        if(length(threshold_value)!=1){
          threshold_value = threshold_value[1]
        }
      }
      else{
        return("ValureError: Filter value must be between 0.50 and 1")}
      
      if(filter_mean == T){
        threshold_value = z[4] 
      }
      
      if(plot==T){
        hist(x, breaks = 1000)
        abline(v = threshold_value, lty=2,col = "red")
        labs(title = rownames(exp[,i]),
             subtitle = paste("Threshold: ",filter_value,sep=""))
        
        plot(ecdf(x),verticals = TRUE, do.p = FALSE)
        abline(v = threshold_value, lty=2,col = "red")
        abline(h = filter_value, lty=2,col = "red")
        labs(title = rownames(exp[,i]),
             subtitle = paste("Threshold: ",filter_value,sep=""))
      }
      
      ##########################
      
      outmatrix[j,] = ifelse(x>=threshold_value,x,0) # 最终录入原score还是1
      threshold_matrix[1,j] = threshold_value
    }
  }
  
  if(out_filter_value ==T){
    return(threshold_matrix)
  }
  
  return(outmatrix)
}

OT_process = function(rev = F, softmax=T, softvalue=c(), softfilter=T){
  
  # rev: T. 27ss to 9ss; F. 9ss to 27ss
  # trans matirx: row (T-1) to col (T)
  trans_list = c("trans_matrix_9_12","trans_matrix_12_15","trans_matrix_15_18",
                 "trans_matrix_18_21","trans_matrix_21_24","trans_matrix_24_27")
  trans_list_rev = rev(trans_list)
  if(rev==T){
    trans_list_use = trans_list_rev
  }else{
    trans_list_use = trans_list
  }
  names(trans_list_use) = c(1:length(trans_list_use))
  
  # Set metadata time start
  metadata_9ss_start = src.9ss.integrated@meta.data
  metadata_27ss_start = src.27ss.integrated@meta.data
  if(rev==T){
    metadata_start = metadata_27ss_start
  }else{
    metadata_start = metadata_9ss_start
  }
  
  # Set ident
  metadata_start$ident = metadata_start$cluster.v06.26.re..correct..un
  
  # Set softmax filter threshold
  if(length(softvalue) == 0 ){
    softvalue = rep(0.92,6)
  }
  
  # Set Index-0: rev T col, F row 
  if(rev==T){
    cluster_m.index0 = 
      matrix(0, 
             nrow = length(unique(metadata_start$ident)),
             ncol = ncol(get(trans_list_use[1])),
             dimnames = list(unique(metadata_start$ident),
                             colnames(get(trans_list_use[1]))))
    
    for (cell in colnames( get(trans_list_use[1]) )) {
      cluster_m.index0[metadata_start[cell,]$ident, cell]=1
    }
    
  }else{
    cluster_m.index0 = 
      matrix(0, 
             nrow = length(unique(metadata_start$ident)),
             ncol = nrow(get(trans_list_use[1])),
             dimnames = list(unique(metadata_start$ident),
                             rownames(get(trans_list_use[1]))))
    
    for (cell in rownames( get(trans_list_use[1]) )) {
      cluster_m.index0[metadata_start[cell,]$ident, cell]=1
    }
  }
  
  # Set Index-1-6: rev T col, F row 
  for(i in c(1:6)){
    cluster_m.time0 = get(paste("cluster_m.index",as.character(i-1),sep=""))
    # cluster_m.time1 = get(paste("cluster_m.index",str(i),sep=""))
    
    if(rev==T){
      cluster_m.time1 = cluster_m.time0 %*% t(as.matrix(get(trans_list_use[i])))
    }else{
      cluster_m.time1 = cluster_m.time0 %*% as.matrix(get(trans_list_use[i]))
    }
    
    if(softmax==T){
      cluster_m.time1 = Softmax(cluster_m.time1, 
                                filter = softfilter, 
                                filter_value = softvalue[i])
    }
    
    if(rev==T){
      colnames(cluster_m.time1) = rownames(get(trans_list_use[i]))
    }else{
      colnames(cluster_m.time1) = colnames(get(trans_list_use[i]))
    }
    
    rownames(cluster_m.time1) = rownames(cluster_m.time0)
    assign(paste("cluster_m.index",as.character(i),sep=""), cluster_m.time1)
  }
  
  ckuster_m.index = list(index0 = cluster_m.index0,
                         index1 = cluster_m.index1,
                         index2 = cluster_m.index2,
                         index3 = cluster_m.index3,
                         index4 = cluster_m.index4,
                         index5 = cluster_m.index5,
                         index6 = cluster_m.index6)
  
  return(ckuster_m.index)
}

OT_soft_output = function(ot_list, filter = T, softvalue = c()){
  # ot_list = ot_list 
  
  ot_list = ot_list_forward
  ot_list_index = names(ot_list)
  names(ot_list_index) = c(1:length(ot_list_index))
  
  # Set softmax filter threshold
  if(length(softvalue) == 0 ){
    softvalue = rep(0.92,6)
  }
  
  for(i in names(ot_list_index)){
    
    cluste_s.index = ot_list[ot_list_index[i]][[1]]
    cluste_s.index = Softmax(cluste_s.index, filter = filter, 
                             filter_value = softvalue[i])
    
    assign(paste("cluster_s.index",as.character(i-1),sep=""),
           cluste_s.index)
  }
  
  ot_list_res = list(index0 = cluster_s.index0,
                     index1 = cluster_s.index1,
                     index2 = cluster_s.index2,
                     index3 = cluster_s.index3,
                     index4 = cluster_s.index4,
                     index5 = cluster_s.index5,
                     index6 = cluster_s.index6)
  
  return(ot_list_res)
}


softmax_small = function(exp){
  ncols  = ncol(exp)
  nrows  = nrow(exp)
  outmatrix = matrix(0, ncol = ncols, nrow = nrows)
  threshold_matrix = matrix(0, ncol = nrows, nrow = 1)
  for(i in c(1:ncols)){
    scale_list = exp[,i] / (1e-10)
    outmatrix[,i] = sapply(scale_list, FUN = exp) / sum(sapply(scale_list, FUN = exp))
  }
  return(outmatrix)
}

softmax_small_add = function(exp, scale.value = 1e-3){
  ncols  = ncol(exp)
  nrows  = nrow(exp)
  outmatrix = matrix(0, ncol = ncols, nrow = nrows)
  threshold_matrix = matrix(0, ncol = nrows, nrow = 1)
  for(i in c(1:ncols)){
    scale_list = exp[,i] / (scale.value)
    outmatrix[,i] = sapply(scale_list, FUN = exp) / sum(sapply(scale_list, FUN = exp))
  }
  return(outmatrix)
}

ot_series_calculation = function(start = 9, end = 27,
                                 index_start = "cluster.v06.26.re..correct..un",
                                 index_end = "cluster.v06.26.re..correct..un"){
  
  # creat trans_matrix
  i = min(start,end)
  j = max(start,end)
  if(i == start){
    index_i = index_start
    index_j = index_end
  }else{
    index_j = index_start
    index_i = index_end
  }
  
  if(i<j){ # i <= j
    m = i+3
    trans_matrix = as.matrix(get(paste("trans_matrix_",i,"_",m, sep = "")))
    
    while (m<j) {
      trans_matrix = trans_matrix %*% as.matrix(get(paste("trans_matrix_",m,"_",m+3, sep = "")))
      m = m+3
    }
    
  }else{ # i = j
    trans_matrix = as.matrix(get(paste("trans_matrix_",i,"_",m, sep = "")))
  }
  
  metadata_i = get(paste("src.",i,"ss.integrated",sep=""))@meta.data
  metadata_j = get(paste("src.",j,"ss.integrated",sep=""))@meta.data
  colname = rownames(metadata_i[metadata_i$nFeature_RNA>2500, ])
  rowname = rownames(metadata_j[metadata_j$nFeature_RNA>2500, ])
  trans_matrix_re = trans_matrix[colname, rowname]
  
  trans_summary = by(data = trans_matrix_re,
                     INDICES = metadata_i[colname, index_i],
                     FUN = colSums)
  trans_summary = do.call(rbind,trans_summary)
  
  trans_summary = t(by(data = t(trans_summary),
                       INDICES = metadata_j[rowname, index_j],
                       FUN = colSums))
  trans_summary = do.call(cbind,trans_summary)
  colnames(trans_summary) = names(table(metadata_j[rowname, index_j]))
  
  trans_summary.mean = trans_summary / 
    sqrt(table(metadata_i[colname, index_i]) %*%
           t(table(metadata_j[rowname, index_j])))
  
  trans_summary.softmax.S2E = softmax_small(t(trans_summary.mean))
  trans_summary.softmax.E2S = softmax_small(trans_summary.mean)
  # trans_summary.mean
  # rownames(trans_summary.softmax) = colnames(trans_summary.mean)
  # colnames(trans_summary.softmax) = rownames(trans_summary.mean)
  trans_index = list(raw_matrix = trans_matrix,
                     use_matrix = trans_matrix_re,
                     raw_summary = trans_summary,
                     norm = trans_summary.mean,
                     norm_softmax.S2E = trans_summary.softmax.S2E,
                     norm_softmax.E2S = trans_summary.softmax.E2S)
  
  return(trans_index)
}
#---------------------------------------

#>> Set cellorder
#> Need to add cell type information 
#> Check [Supplymental Table S1. sheet3. 10x-v3 celltypes]
#> Use "cluster.v06.26.re..correct..un" to store
#---------------------------------------
celllist.order.endoderm = c(
  "FG.1","Pharynx.organ.2","FG.2","Pharynx.organ.1","FG.3","Pharynx.organ.4",
  "FG.4","FG.4-Lung/Stomach","FG.4-Liver","Pharynx.organ.5","Lung","FG.5","Pharynx.organ.3","Thyroid","FG.6","Esophagus",
  "AL.1","AL.2", "Liver","AL.3",'AL.3-Small.intestine.1',"AL.3-EHBD/VP","AL.3-Liver","EHBD","VP","DP","EP","EP.1","EP.2",
  "MG.1","Stomach","MG.2","MG.3","MG.3.M","MG.3.P","Small.intestine.1","Small.intestine.2",
  "HG.1","Large.intestine.1","HG.1-Large.intestine.2","Large.intestine.2","HG.2","Large.intestine.3")

celllist.order.9ss = c("FG.1","FG.2","FG.3","FG.4","FG.5","FG.6",
                       "AL.1","AL.2","AL.3","MG.1","MG.2","MG.3","HG.1","HG.2")

celllist.order.12ss = c(
  "FG.1","FG.2","FG.3","FG.4","FG.4-Lung/Stomach","FG.4-Liver","FG.5","FG.6",
  "AL.1","AL.2", # "AL.3",
  'AL.3-Small.intestine.1',"AL.3-EHBD/VP","AL.3-Liver",
  "MG.1","MG.2", # "MG.3",
  "MG.3.A","MG.3.P",
  "HG.1","HG.1-Large.intestine.2","HG.2")

celllist.order.15ss = c(
  "FG.1","Pharynx.organ.2","FG.2","FG.3","Pharynx.organ.4",
  "FG.4","FG.4-Lung/Stomach","FG.4-Liver","Pharynx.organ.5",
  "FG.5","Pharynx.organ.3","Thyroid","FG.6","Esophagus",
  "AL.1/2-Liver","Liver",'AL.3-Small.intestine.1',"AL.3-EHBD/VP","AL.3-Liver",
  "MG.1","MG.2","MG.3.A","MG.3.P","Small.intestine.1","Small.intestine.2",
  "HG.1","HG.1-Large.intestine.2","HG.2")

celllist.order.18ss = c(
  "FG.1","Pharynx.organ.2","FG.2","Pharynx.organ.1","Pharynx.organ.4",
  "FG.4-Lung/Stomach","FG.4-Liver","Pharynx.organ.5","Lung",
  "FG.5","Pharynx.organ.3","Thyroid","FG.6","Esophagus",
  "AL.1/2-Liver","Liver","AL.3-EHBD/VP","AL.3-Liver",
  "EHBD","VP","DP","EP.1","EP.2","Stomach",
  "MG.3.A","MG.3.P","Small.intestine.1","Small.intestine.2",
  "HG.1","Large.intestine.1","HG.1-Large.intestine.2","Large.intestine.2","HG.2","Large.intestine.3")

celllist.order.21ss = c(
  "FG.1","Pharynx.organ.2","FG.2","Pharynx.organ.1","Pharynx.organ.4",
  "FG.4-Lung/Stomach","Pharynx.organ.5","Lung",
  "Pharynx.organ.3","Thyroid","Esophagus",
  "AL.1/2-Liver","Liver","EHBD","VP","DP","EP.1","EP.2","Stomach",
  "MG.3.P","Small.intestine.1","Small.intestine.2",
  "Large.intestine.1","Large.intestine.2","HG.2","Large.intestine.3")

celllist.order.24ss = c(
  "Pharynx.organ.2","Pharynx.organ.1",
  "Pharynx.organ.4","Pharynx.organ.5","Pharynx.organ.3",
  "Thyroid","Lung","Esophagus","Stomach",
  "Liver","EHBD","VP","DP","EP.1","EP.2",
  "Small.intestine.1","Small.intestine.2",
  "Large.intestine.1","Large.intestine.2","Large.intestine.3")

celllist.order.27ss = c(
  "Pharynx.organ.2","Pharynx.organ.1",
  "Pharynx.organ.4","Pharynx.organ.5","Pharynx.organ.3",
  "Thyroid","Lung","Esophagus","Stomach",
  "Liver","EHBD","VP","DP","EP.1","EP.2",
  "Small.intestine.1","Small.intestine.2",
  "Large.intestine.1","Large.intestine.2","Large.intestine.3")

celllist.order.27ss.mergeEP = c(
  "Pharynx.organ.2","Pharynx.organ.1",
  "Pharynx.organ.4","Pharynx.organ.5","Pharynx.organ.3",
  "Thyroid","Lung","Esophagus","Stomach",
  "Liver","EHBD","VP","DP","EP",
  "Small.intestine.1","Small.intestine.2",
  "Large.intestine.1","Large.intestine.2","Large.intestine.3")
#---------------------------------------


#> Merge EP.1 (PEP) / EP.2 (alpha) to EP
src.27ss.integrated$cluster.v06.26.re..correct..un.mergeEP = 
  src.27ss.integrated$cluster.v06.26.re..correct..un
src.27ss.integrated@meta.data[
  src.27ss.integrated$cluster.v06.26.re..correct..un.mergeEP%in%c("EP.1","EP.2"),]$cluster.v06.26.re..correct..un.mergeEP = "EP"

ot_series_9_27 = ot_series_calculation(start = 9, end = 27,
                                       index_start = "cluster.v06.26.re..correct..un",
                                       index_end = "cluster.v06.26.re..correct..un.mergeEP")

ot_series_9_27 = ot_series_calculation(start = 9, end = 27,
                                       index_start = "cluster.v06.26.re..correct..un",
                                       index_end = "cluster.v06.26.re..correct..un.mergeEP")

save(ot_series_9_27, file = "_revise_240801/OT_remake/ot_series_9_27.Rdata")

#------------------------------------------------
#>>> In Figure.S5O, we perfomrm: 
#> (1) from 9SS to 27SS
#> (2) from 27SS to 9SS
#------------------------------------------------
pdf("figure.v08.07/ot_refine/ot.series.pdf",7,7)
MyHeatmap(ot_series_9_27$norm_softmax.E2S[index_9, index_27],
          type = "row.relat",
          labRow = index_9,
          labCol = index_27,
          Colv = "none", Rowv = "none",
          color.palette = colorRampPalette(c("#eeeeee","#F64E60","#C70039"), # "#900C3F","#581845"), 
                                           space="Lab"),
          margins = c(10,10),
          key.title = "norm_softmax row relat 9ss to 27ss")

MyHeatmap(ot_series_9_27$norm_softmax.S2E[index_27, index_9],
          type = "row.relat",
          labRow = index_27,
          labCol = index_9,
          Colv = "none", Rowv = "none",
          color.palette = colorRampPalette(c("#eeeeee","#F64E60","#C70039"), # "#900C3F","#581845"), 
                                           space="Lab"),
          margins = c(10,10),
          key.title = "norm_softmax row.relat 27ss to 9ss")
dev.off()
#------------------------------------------------

#-------------------------------------------------------------------------
#> You can also perform OT computation in: 
#> (1) adjacent times
#> (2) from 9SS to other time points
#> (3) from 27SS to other time points 
#-------------------------------------------------------------------------
#->>> adjacent
ot_series_9_12 = ot_series_calculation(start = 9, end = 12)
ot_series_12_15 = ot_series_calculation(start = 12, end = 15)
ot_series_15_18 = ot_series_calculation(start = 15, end = 18)
ot_series_18_21 = ot_series_calculation(start = 18, end = 21)
ot_series_21_24 = ot_series_calculation(start = 21, end = 24)
ot_series_24_27 = ot_series_calculation(start = 24, end = 27)

ot_series_9_12$norm_softmax.E2S = softmax_small_add(ot_series_9_12$norm[celllist.order.9ss, celllist.order.12ss], scale.value = 1e-3)
ot_series_12_15$norm_softmax.E2S = softmax_small_add(ot_series_12_15$norm[celllist.order.12ss, celllist.order.15ss], scale.value = 1e-3)
ot_series_15_18$norm_softmax.E2S = softmax_small_add(ot_series_15_18$norm[celllist.order.15ss, celllist.order.18ss], scale.value = 1e-3)
ot_series_18_21$norm_softmax.E2S = softmax_small_add(ot_series_18_21$norm[celllist.order.18ss, celllist.order.21ss], scale.value = (10^(-2.75)))
ot_series_21_24$norm_softmax.E2S = softmax_small_add(ot_series_21_24$norm[celllist.order.21ss, celllist.order.24ss], scale.value = (10^(-2.75)))
ot_series_24_27$norm_softmax.E2S = softmax_small_add(ot_series_24_27$norm[celllist.order.24ss, celllist.order.27ss], scale.value = 1e-3)

ot_series_9_12$norm_softmax.S2E = softmax_small_add(t(ot_series_9_12$norm[celllist.order.9ss, celllist.order.12ss]), scale.value = 1e-3)
ot_series_12_15$norm_softmax.S2E = softmax_small_add(t(ot_series_12_15$norm[celllist.order.12ss, celllist.order.15ss]), scale.value = 1e-3)
ot_series_15_18$norm_softmax.S2E = softmax_small_add(t(ot_series_15_18$norm[celllist.order.15ss, celllist.order.18ss]), scale.value = 1e-3)
ot_series_18_21$norm_softmax.S2E = softmax_small_add(t(ot_series_18_21$norm[celllist.order.18ss, celllist.order.21ss]), scale.value = (10^(-2.75)))
ot_series_21_24$norm_softmax.S2E = softmax_small_add(t(ot_series_21_24$norm[celllist.order.21ss, celllist.order.24ss]), scale.value = (10^(-2.75)))
ot_series_24_27$norm_softmax.S2E = softmax_small_add(t(ot_series_24_27$norm[celllist.order.24ss, celllist.order.27ss]), scale.value = 1e-3)


#->>> Start-9SS
ot_series_9_12 = ot_series_calculation(start = 9, end = 12)
ot_series_9_15 = ot_series_calculation(start = 9, end = 15)
ot_series_9_18 = ot_series_calculation(start = 9, end = 18)
ot_series_9_21 = ot_series_calculation(start = 9, end = 21)
ot_series_9_24 = ot_series_calculation(start = 9, end = 24)

ot_series_9_12$norm_softmax.E2S = softmax_small_add(ot_series_9_12$norm[celllist.order.9ss, celllist.order.12ss], scale.value = 10^(-4.5))
ot_series_9_15$norm_softmax.E2S = softmax_small_add(ot_series_9_15$norm[celllist.order.9ss, celllist.order.15ss], scale.value = 10^(-4.5))
ot_series_9_18$norm_softmax.E2S = softmax_small_add(ot_series_9_18$norm[celllist.order.9ss, celllist.order.18ss], scale.value = 10^(-6))
ot_series_9_21$norm_softmax.E2S = softmax_small_add(ot_series_9_21$norm[celllist.order.9ss, celllist.order.21ss], scale.value = 10^(-7.5))
ot_series_9_24$norm_softmax.E2S = softmax_small_add(ot_series_9_24$norm[celllist.order.9ss, celllist.order.24ss], scale.value = 10^(-9))
ot_series_9_27$norm_softmax.E2S = softmax_small_add(ot_series_9_27$norm[celllist.order.9ss, celllist.order.27ss], scale.value = 10^(-10))

ot_series_9_12$norm_softmax.S2E = softmax_small_add(t(ot_series_9_12$norm[celllist.order.9ss, celllist.order.12ss]), scale.value = 10^(-4.5))
ot_series_9_15$norm_softmax.S2E = softmax_small_add(t(ot_series_9_15$norm[celllist.order.9ss, celllist.order.15ss]), scale.value = 10^(-4.5))
ot_series_9_18$norm_softmax.S2E = softmax_small_add(t(ot_series_9_18$norm[celllist.order.9ss, celllist.order.18ss]), scale.value = 10^(-6))
ot_series_9_21$norm_softmax.S2E = softmax_small_add(t(ot_series_9_21$norm[celllist.order.9ss, celllist.order.21ss]), scale.value = 10^(-7.5))
ot_series_9_24$norm_softmax.S2E = softmax_small_add(t(ot_series_9_24$norm[celllist.order.9ss, celllist.order.24ss]), scale.value = 10^(-9))

#->>> End-27SS
ot_series_12_27 = ot_series_calculation(start = 12, end = 27)
ot_series_15_27 = ot_series_calculation(start = 15, end = 27)
ot_series_18_27 = ot_series_calculation(start = 18, end = 27)
ot_series_21_27 = ot_series_calculation(start = 21, end = 27)
ot_series_24_27 = ot_series_calculation(start = 24, end = 27)

ot_series_12_27$norm_softmax.E2S = softmax_small_add(ot_series_12_27$norm[celllist.order.12ss, celllist.order.27ss], scale.value = 10^(-8.5))
ot_series_15_27$norm_softmax.E2S = softmax_small_add(ot_series_15_27$norm[celllist.order.15ss, celllist.order.27ss], scale.value = 10^(-7.5))
ot_series_18_27$norm_softmax.E2S = softmax_small_add(ot_series_18_27$norm[celllist.order.18ss, celllist.order.27ss], scale.value = 10^(-6))
ot_series_21_27$norm_softmax.E2S = softmax_small_add(ot_series_21_27$norm[celllist.order.21ss, celllist.order.27ss], scale.value = 10^(-4))
ot_series_24_27$norm_softmax.E2S = softmax_small_add(ot_series_24_27$norm[celllist.order.24ss, celllist.order.27ss], scale.value = 10^(-4))

ot_series_12_27$norm_softmax.S2E = softmax_small_add(t(ot_series_12_27$norm[celllist.order.12ss, celllist.order.27ss]), scale.value = 10^(-8.5))
ot_series_15_27$norm_softmax.S2E = softmax_small_add(t(ot_series_15_27$norm[celllist.order.15ss, celllist.order.27ss]), scale.value = 10^(-7.5))
ot_series_18_27$norm_softmax.S2E = softmax_small_add(t(ot_series_18_27$norm[celllist.order.18ss, celllist.order.27ss]), scale.value = 10^(-6))
ot_series_21_27$norm_softmax.S2E = softmax_small_add(t(ot_series_21_27$norm[celllist.order.21ss, celllist.order.27ss]), scale.value = 10^(-4))
ot_series_24_27$norm_softmax.S2E = softmax_small_add(t(ot_series_24_27$norm[celllist.order.24ss, celllist.order.27ss]), scale.value = 10^(-4))


#-- adjacent
for(i.ot_series in c(c("ot_series_9_12", "ot_series_12_15", "ot_series_15_18",
                       "ot_series_18_21", "ot_series_21_24", "ot_series_24_27"))){
  print(i.ot_series)
  
  ot_series = get(i.ot_series)
  index.i = strsplit(i.ot_series, "_")[[1]][3]
  index.j = strsplit(i.ot_series, "_")[[1]][4]
  
  celllist.order.i = get(paste("celllist.order.", index.i, "ss",sep=""))
  celllist.order.j = get(paste("celllist.order.", index.j, "ss",sep=""))
  
  if(i.ot_series == "ot_series_9_27"){
    celllist.order.i = get(paste("celllist.order.", index.i, "ss",sep=""))
    celllist.order.j = get(paste("celllist.order.", index.j, "ss.mergeEP",sep=""))
  }
  
  rownames(ot_series$norm_softmax.E2S) = celllist.order.i
  colnames(ot_series$norm_softmax.E2S) = celllist.order.j
  
  rownames(ot_series$norm_softmax.S2E) = celllist.order.j
  colnames(ot_series$norm_softmax.S2E) = celllist.order.i
  
  pdf(paste("_revise_240801/OT_remake/", i.ot_series, ".pdf", sep = ""), 7,7)
  MyHeatmap(ot_series$norm_softmax.E2S[celllist.order.i, celllist.order.j],
            type = "row.relat",
            labRow = celllist.order.i,
            labCol = celllist.order.j,
            Colv = "none", Rowv = "none",
            ColSideColorsSize = 1.5,
            RowSideColorsSize = 1.5,
            ColSideColors = cbind(MyName2Col(celllist.order.j, cluster.endoderm.color.v5)),
            RowSideColors = t(cbind(MyName2Col(celllist.order.i, cluster.endoderm.color.v5))),
            color.palette = colorRampPalette(c("#eeeeee","#F64E60","#C70039"), # "#900C3F","#581845"), 
                                             space="Lab"),
            margins = c(10,10),
            key.title = paste("norm_softmax row relat ", index.i, "ss to ", index.j, "ss", sep = ""))
  dev.off()
}
#-- forkward
for(i.ot_series in c(c("ot_series_9_12", "ot_series_9_15", "ot_series_9_18",
                       "ot_series_9_21", "ot_series_9_24"))){
  print(i.ot_series)
  
  ot_series = get(i.ot_series)
  index.i = strsplit(i.ot_series, "_")[[1]][3]
  index.j = strsplit(i.ot_series, "_")[[1]][4]
  
  celllist.order.i = get(paste("celllist.order.", index.i, "ss",sep=""))
  celllist.order.j = get(paste("celllist.order.", index.j, "ss",sep=""))
  
  if(i.ot_series == "ot_series_9_27"){
    celllist.order.i = get(paste("celllist.order.", index.i, "ss",sep=""))
    celllist.order.j = get(paste("celllist.order.", index.j, "ss.mergeEP",sep=""))
  }
  
  rownames(ot_series$norm_softmax.E2S) = celllist.order.i
  colnames(ot_series$norm_softmax.E2S) = celllist.order.j
  
  rownames(ot_series$norm_softmax.S2E) = celllist.order.j
  colnames(ot_series$norm_softmax.S2E) = celllist.order.i
  
  pdf(paste("_revise_240801/OT_remake/", i.ot_series, ".pdf", sep = ""), 7,7)
  MyHeatmap(ot_series$norm_softmax.E2S[celllist.order.i, celllist.order.j],
            type = "row.relat",
            labRow = celllist.order.i,
            labCol = celllist.order.j,
            Colv = "none", Rowv = "none",
            ColSideColorsSize = 1.5,
            RowSideColorsSize = 1.5,
            ColSideColors = cbind(MyName2Col(celllist.order.j, cluster.endoderm.color.v5)),
            RowSideColors = t(cbind(MyName2Col(celllist.order.i, cluster.endoderm.color.v5))),
            color.palette = colorRampPalette(c("#eeeeee","#F64E60","#C70039"), # "#900C3F","#581845"), 
                                             space="Lab"),
            margins = c(10,10),
            key.title = paste("norm_softmax row relat ", index.i, "ss to ", index.j, "ss", sep = ""))
  dev.off()
}
#-- backward
for(i.ot_series in c(c("ot_series_12_27", "ot_series_15_27",
                       "ot_series_18_27", "ot_series_21_27", "ot_series_24_27"))){
  print(i.ot_series)
  
  ot_series = get(i.ot_series)
  index.i = strsplit(i.ot_series, "_")[[1]][3]
  index.j = strsplit(i.ot_series, "_")[[1]][4]
  
  celllist.order.i = intersect(get(paste("celllist.order.", index.i, "ss",sep="")), rownames(ot_series$norm))
  celllist.order.j = intersect(get(paste("celllist.order.", index.j, "ss",sep="")), colnames(ot_series$norm))
  if(i.ot_series == "ot_series_9_27"){
    celllist.order.i = intersect(get(paste("celllist.order.", index.i, "ss",sep="")), rownames(ot_series$norm))
    celllist.order.j = intersect(get(paste("celllist.order.", index.j, "ss.mergeEP",sep="")), colnames(ot_series$norm))
  }
  
  rownames(ot_series$norm_softmax.E2S) = celllist.order.i
  colnames(ot_series$norm_softmax.E2S) = celllist.order.j
  rownames(ot_series$norm_softmax.S2E) = celllist.order.j
  colnames(ot_series$norm_softmax.S2E) = celllist.order.i
  
  pdf(paste("_revise_240801/OT_remake/", i.ot_series, ".pdf", sep = ""), 7,7)
  MyHeatmap(ot_series$norm_softmax.S2E[celllist.order.j, celllist.order.i],
            type = "row.relat",
            labRow = celllist.order.i,
            labCol = celllist.order.j,
            Colv = "none", Rowv = "none",
            ColSideColorsSize = 1.5,
            RowSideColorsSize = 1.5,
            RowSideColors = t(cbind(MyName2Col(celllist.order.j, cluster.endoderm.color.v5))),
            ColSideColors = cbind(MyName2Col(celllist.order.i, cluster.endoderm.color.v5)),
            color.palette = colorRampPalette(c("#eeeeee","#F64E60","#C70039"), # "#900C3F","#581845"), 
                                             space="Lab"),
            margins = c(10,10),
            key.title = paste("norm_softmax row relat ", index.i, "ss to ", index.j, "ss", sep = ""))
  dev.off()
}
#-------------------------------------------------------------------------


#==========================================================









