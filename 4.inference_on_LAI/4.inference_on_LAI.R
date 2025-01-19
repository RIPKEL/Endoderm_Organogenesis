#===============================================================================
#>> 4.Inference on LAI
#===============================================================================
library(batchelor)
library(progress)
library(progressr)
library(scran)
library(Seurat)
library(SeuratObject)
library(future)
library(furrr)
library(doParallel)
library(igraph)
library(scales)

#===============================================================================
#>>>>> 4.1 batch_process_integration
#===============================================================================
#>> Need: 
#> (1) 10x-v3: src.xxss.integrated [9-27]
#> (2) Smart-seq3: src.sm3.merge
batch_process_integration = function(time){
  timess = paste(time,"ss",sep="")
  sstime = paste("ss",time,sep="")
  print(timess);print(sstime)
  
  src.sm3.merge.timess = src.sm3.merge[,src.sm3.merge$Time%in%timess]
  
  src.timess.integrated = get(paste("src.",timess,".integrated",sep=""))
  src.timess.integrated.merge.selectgene =
    get(paste("src.",timess,".integrated.merge.selectgene",sep=""))
  
  
  src.timess.integrated.merge = merge(src.timess.integrated, src.sm3.merge.timess)
  src.timess.integrated.merge$batch_Time = paste(src.timess.integrated.merge$batch,
                                                 src.timess.integrated.merge$Time,sep="_")
  
  src.timess.integrated.merge = ScaleData(src.timess.integrated.merge, 
                                          features = rownames(src.timess.integrated.merge),
                                          split.by = "batch_Time")
  
  src.timess.integrated.selectgene = rownames(src.timess.integrated@reductions$pca@feature.loadings)
  src.timess.integrated.merge.selectgene = src.timess.integrated.selectgene
  
  src.timess.integrated.merge = RunPCA(src.timess.integrated.merge,
                                       features = src.timess.integrated.merge.selectgene)
  src.timess.integrated.merge = RunUMAP(src.timess.integrated.merge,
                                        dims = 1:20,n.neighbors = 50)
  
  src.timess.integrated.merge@reductions$umap@cell.embeddings = 
    cbind(src.timess.integrated.merge@reductions$umap@cell.embeddings[,2],
          src.timess.integrated.merge@reductions$umap@cell.embeddings[,1])
  colnames(src.timess.integrated.merge@reductions$umap@cell.embeddings) = c("UMAP_1","UMAP_2")
  
  # MNN
  #---------------
  MNN.res = 
    mnnCorrect(as.matrix(src.timess.integrated.merge@assays$RNA@data[src.timess.integrated.merge.selectgene, src.timess.integrated.merge$Time%in%sstime]),
               as.matrix(src.timess.integrated.merge@assays$RNA@data[src.timess.integrated.merge.selectgene, src.timess.integrated.merge$Time%in%timess]),
               k = 5,cos.norm.out=F)
  src.timess.integrated.merge@assays$mnnRNA = CreateAssayObject(data = MNN.res@assays@data$corrected)
  src.timess.integrated.merge@assays$mnnRNA@key = "mnn_"
  src.timess.integrated.merge = ScaleData(src.timess.integrated.merge, 
                                          rownames(src.timess.integrated.merge@assays$mnnRNA),
                                          assay = "mnnRNA")
  
  src.timess.integrated.merge@meta.data = src.timess.integrated.merge@meta.data[!src.timess.integrated.merge$Time%in%NA,]
  
  # FTA
  #---------------
  anchor.timess.integrated = 
    FindTransferAnchors(reference = src.timess.integrated.merge[,src.timess.integrated.merge$Time%in%sstime],
                        query = src.timess.integrated.merge[,src.timess.integrated.merge$Time%in%timess],
                        reference.assay = "mnnRNA",
                        query.assay = "mnnRNA",scale = T,
                        features = src.timess.integrated.selectgene)
  umap.transfer.timess.integrated = TransferData(anchor.timess.integrated,
                                                 t(src.timess.integrated@reductions$umap@cell.embeddings))
  
  src.timess.integrated.merge@reductions$umap_fta = src.timess.integrated.merge@reductions$umap
  src.timess.integrated.merge@reductions$umap_fta@cell.embeddings[colnames(src.timess.integrated),] = 
    src.timess.integrated@reductions$umap@cell.embeddings
  src.timess.integrated.merge@reductions$umap_fta@cell.embeddings[colnames(
    src.timess.integrated.merge[,src.timess.integrated.merge$Time%in%timess]),] = 
    as.matrix(t(umap.transfer.timess.integrated@data))
  
  src.timess.integrated.merge@reductions$umap_fta@key = "Coord_"
  colnames(src.timess.integrated.merge@reductions$umap_fta@cell.embeddings) = c("Coord_1","Coord_2")
  
  # FIA
  #---------------
  anchor.timess.integrated.merge = FindIntegrationAnchors(c(src.timess.integrated.merge[,src.timess.integrated.merge$Time%in%timess],
                                                            src.timess.integrated.merge[,src.timess.integrated.merge$Time%in%sstime]),
                                                          reduction = "cca",
                                                          anchor.features = src.timess.integrated.merge.selectgene ,
                                                          dims = 1:30)
  
  src.timess.integrated.merge.re = IntegrateData(anchorset = anchor.timess.integrated.merge , dims = 1:30)
  DefaultAssay(src.timess.integrated.merge.re) = "integrated"
  src.timess.integrated.merge.re = ScaleData(src.timess.integrated.merge.re, features = rownames(src.timess.integrated.merge.re))
  src.timess.integrated.merge.re = RunPCA(src.timess.integrated.merge.re, features = rownames(src.timess.integrated@reductions$pca@feature.loadings))
  src.timess.integrated.merge.re = FindNeighbors(src.timess.integrated.merge.re, dims = 1:30)
  src.timess.integrated.merge.re = FindClusters(src.timess.integrated.merge.re, resolution = 2)
  src.timess.integrated.merge.re = FindClusters(src.timess.integrated.merge.re, resolution = 3)
  src.timess.integrated.merge.re = FindClusters(src.timess.integrated.merge.re, resolution = 4)
  
  DimPlot(src.timess.integrated.merge.re, reduction = "pca",group.by = "Time")
  DimPlot(src.timess.integrated.merge.re, reduction = "pca",group.by = "cluster.revise.re.v1.30.re", cols = cluster.endoderm.color.v3, pt.size = 2)
  
  src.timess.integrated.merge.re = RunUMAP(src.timess.integrated.merge.re, dims = 1:20, n.neighbors = 100)
  
  src.timess.integrated.merge@reductions$umap_integrated = src.timess.integrated.merge.re@reductions$umap
  src.timess.integrated.merge@reductions$pca_integrated = src.timess.integrated.merge.re@reductions$pca
  
  src.timess.integrated.merge@reductions$umap_integrated@assay.used = "RNA"
  src.timess.integrated.merge@reductions$pca_integrated@assay.used = "RNA"
  DimPlot(src.timess.integrated.merge, group.by = "lineage", reduction = "umap_integrated",
          cols = color.lineage, na.value = "#eeeeee")
  rm(src.timess.integrated.merge.re)
  
  src.timess.integrated.merge = 
    src.timess.integrated.merge[,src.timess.integrated.merge$nFeature_RNA>2500]
  
  return(src.timess.integrated.merge) 
}

#> Set integrated.merge.selectgene as integrated.selectgene
for(i.time in c("9ss","12ss","15ss","18ss","21ss","24ss","27ss")){
  assign(paste("src.",timess,".integrated.merge.selectgene",sep=""),
         get(paste("src.",timess,".integrated.selectgene",sep="")))
}

src.9ss.integrated.merge = batch_process_integration("9")
src.12ss.integrated.merge = batch_process_integration("12")
src.15ss.integrated.merge = batch_process_integration("15")
src.18ss.integrated.merge = batch_process_integration("18")
src.21ss.integrated.merge = batch_process_integration("21")
src.24ss.integrated.merge = batch_process_integration("24")
src.27ss.integrated.merge = batch_process_integration("27")

#>> Add UMAP based on MNN
for(i.time in c("9ss","12ss","15ss","18ss","21ss","24ss","27ss")){
  seurat = get(paste("src.",timess,".integrated.merge",sep=""))
  seurat.selectgene = get(paste("src.",timess,".integrated.merge.selectgene",sep="")) 
  
  if(i.time == "9ss"){
    fMNN.res =  fastMNN(as.matrix(seurat@assays$RNA@data[seurat.selectgene, seurat$Time%in%paste(gsub("ss","",i.time),"ss",sep="")]),
                        as.matrix(seurat@assays$RNA@data[seurat.selectgene, seurat$Time%in%i.time]),
                        d=50,
                        subset.row = seurat.selectgene)
  }else{
    fMNN.res =  fastMNN(as.matrix(seurat@assays$RNA@data[seurat.selectgene, seurat$Time%in%paste(gsub("ss","",i.time),"ss",sep="")]),
                        as.matrix(seurat@assays$RNA@data[seurat.selectgene, seurat$Time%in%i.time&
                                                           seurat$batch%in%1]),
                        as.matrix(seurat@assays$RNA@data[seurat.selectgene, seurat$Time%in%i.time&
                                                           seurat$batch%in%2]),
                        d=50,
                        subset.row = seurat.selectgene)
  }
  
  seurat@reductions$mnn = seurat@reductions$umap
  seurat@reductions$mnn@cell.embeddings = fMNN.res@assays@data$reconstructed@seed@components
  seurat@reductions$mnn@feature.loadings = fMNN.res@assays@data$reconstructed@seed@rotation
  colnames(seurat@reductions$mnn@cell.embeddings) = paste("MNN",1:50,sep="_")
  colnames(seurat@reductions$mnn@feature.loadings) = paste("MNN",1:50,sep="_")
  seurat@reductions$mnn@key = "MNN_"
  
  seurat.re = get(paste("src.",timess,".integrated.merge",sep=""))
  seurat.re@reductions$mnn = seurat@reductions$mnn
  seurat.re@reductions$umap_mnn = seurat@reductions$umap
  
  assign(paste("src.",timess,".integrated.merge",sep=""), seurat.re)
  rm(seurat, seurat.re); gc()
}

#>> Add UMAP based on Scale
for(i.time in c("12ss","15ss","18ss","21ss","24ss","27ss")){
  seurat = get(paste("src.",timess,".integrated.merge",sep=""))
  seurat.selectgene = get(paste("src.",timess,".integrated.merge.selectgene",sep="")) 
  
  seurat = ScaleData(seurat, features = rownames(seurat), split.by = "batch_Time")
  seurat = RunPCA(seurat, features = seurat.selectgene)
  seurat = RunUMAP(seurat, dims = 1:20,n.neighbors = 50)
  # DimPlot(seurat, group.by = "lineage", reduction = "umap")
  
  seurat@reductions$umap@cell.embeddings = cbind(
    seurat@reductions$umap@cell.embeddings[,2],
    -seurat@reductions$umap@cell.embeddings[,1])
  colnames(seurat@reductions$umap@cell.embeddings) = c(
    paste(seurat@reductions$umap@key,1,sep=""),
    paste(seurat@reductions$umap@key,2,sep=""))
  assign(paste("src.",timess,".integrated.merge",sep=""), seurat)
}


pdf("integrated/integrated.summary.pdf",9,18)
for(seurat in c("src.9ss.integrated.merge",
                "src.12ss.integrated.merge",
                "src.15ss.integrated.merge",
                "src.18ss.integrated.merge",
                "src.21ss.integrated.merge",
                "src.24ss.integrated.merge",
                "src.27ss.integrated.merge")){
  seurat_name = gsub("src.","",gsub("ss.integrated.merge","",seurat))
  
  seurat = get(seurat)
  for(reduction in c("umap","umap_mnn",
                     "umap_fta","umap_integrated")){
    assign(paste(reduction,seurat_name,sep = "_"),
           DimPlot(seurat, group.by = "lineage", reduction = reduction,
                   cols = color.lineage, na.value = "#eeeeee", pt.size = 1.2) +
             theme_void() + p_add + 
             theme(legend.position = "none",
                   title = element_blank()))
  }
}
dev.off()


#>>> set rotated umap for ss18 (selected)
embedding_rotated = function(seurat, reduction = "umap",theta= pi/2){
  rotation_matrix = matrix(c(cos(theta), sin(theta), 
                             -sin(theta), cos(theta)), ncol = 2, byrow = TRUE)
  reduction_rotated = paste(reduction, ".rotated", sep = "")
  
  seurat[[reduction_rotated]] = seurat[[reduction]]
  seurat[[reduction_rotated]]@cell.embeddings =
    seurat[[reduction]]@cell.embeddings %*% rotation_matrix
  colnames(seurat[[reduction_rotated]]@cell.embeddings) = c("UMAP_1","UMAP_2")
  seurat[[reduction_rotated]]@key = "UMAP_"
  return(seurat)
}

src.18ss.integrated.merge = embedding_rotated(src.18ss.integrated.merge, reduction = "umap_integrated", theta = pi/7)

#>>> set colnames for reduction matrix
set_umap = function(seurat, reduction){
  col = ncol(seurat[["umap"]]@cell.embeddings)
  seurat[["umap"]]@cell.embeddings = seurat[[reduction]]@cell.embeddings[,1:col]
  colnames(seurat[["umap"]]@cell.embeddings) = paste("UMAP_", c(1:col), sep="")
  seurat[["umap"]]@key = "UMAP_"
  return(seurat)
}

integrated_reduce = hash::hash(keys=Time)
integrated_reduce["ss9"] = c("umap_integrated")
integrated_reduce["ss12"] = c("umap_integrated")
integrated_reduce["ss15"] = c("umap_integrated")
integrated_reduce["ss18"] = c("umap_integrated.rotated")
integrated_reduce["ss21"] = c("umap_integrated")
integrated_reduce["ss24"] = c("umap_integrated")
integrated_reduce["ss27"] = c("umap_integrated")

for(time in Time){
  seurat_name = paste("src.",gsub("ss","",time),"ss.integrated.merge", sep="")
  reduction = hash::values(integrated_reduce, keys=time)[1]
  seurat = get(seurat_name)
  seurat = set_umap(seurat, reduction = reduction)
  assign(seurat_name, seurat)
  if(time == "ss18"){
    seurat_name = paste("src.",gsub("ss","",time),"ss.integrated.merge", sep="")
    reduction = "umap_integrated.rotated"
    seurat = get(seurat_name)
    seurat = set_umap(seurat, reduction = reduction)
    assign(seurat_name, seurat)
  }
}


#===============================================================================


#===============================================================================
#>>>>> 4.2 set network distance and calculate Lineage affinity index
#===============================================================================
#> 
#> We attempted to construct distance matrices \
#> using Euclidean and network shortest paths to compute LAI, \  
#> finding the latter more effective, which was not discussed in the original text.
#> 
#> When executing the following code, you can disregard the part where type is set to 'Euclidean'.
#> 

distance_euclidean_index = function(seurat, dims=1:20){
  
  metric_euclidean = list(
    pca = dist(seurat@reductions$pca@cell.embeddings[,dims], method = "euclidean"),
    mnn = dist(seurat@reductions$mnn@cell.embeddings[,dims], method = "euclidean"),
    pca_integrated = dist(seurat@reductions$pca_integrated@cell.embeddings[,dims], method = "euclidean"),
    mnn_pca_integrated  = dist(seurat@reductions$mnn_pca_integrated@cell.embeddings[,dims], method = "euclidean"))
  
  metric_euclidean$pca = as.matrix(metric_euclidean$pca)
  metric_euclidean$mnn = as.matrix(metric_euclidean$mnn)
  metric_euclidean$pca_integrated = as.matrix(metric_euclidean$pca_integrated)
  metric_euclidean$mnn_pca_integrated = as.matrix(metric_euclidean$mnn_pca_integrated)
  
  return(metric_euclidean )
}

network_extract_construct = function(seurat, reduction.list = c('mnn',"pca_integrated")){
  
  for (reduction in reduction.list) {
    seurat = FindNeighbors(seurat, reduction = reduction, dims = 1:30)
    
    graph_name = ifelse( reduction %in% c("pca_integrated"), "RNA", "mnnRNA")
    
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

distance_network_index = function(seurat){
  
  library(future)
  library(furrr)
  library(doParallel)
  library(Seurat)
  library(igraph)
  library(progress)
  
  num_cores <- 18
  
  # Set up future plan to use multiple cores
  plan(multisession, workers = num_cores)
  
  DefaultAssay(seurat) <- "RNA"
  
  sm3_start <- table(seurat$source_batch)["refer"] + 1
  
  # Define function to calculate distance
  calculate_distance <- function(graph, cor_list, vi) {
    distances <- foreach(vj = 1:length(cor_list), .combine = 'c') %dopar% {
      if (vi <= vj) {
        distance <- 0
      } else if (vi < sm3_start) {
        distance <- 0
      } else if (vj >= sm3_start) {
        distance <- 0
      } else {
        shortest_paths_details <- get.shortest.paths(
          graph, from = cor_list[vi], to = cor_list[vj],
          output = "both", weights = E(graph)$weight)
        
        node_path <- shortest_paths_details$vpath[[1]]
        edge_path <- shortest_paths_details$epath[[1]]
        
        distance <- (length(node_path) - 1) - sum(E(graph)[edge_path]$weight)
      }
      distance
    }
    return(distances)
  }
  
  pb_outer <- progress_bar$new(total = 2, format = "[:bar] :percent :eta")
  
  for (reduction in c("mnn", "pca_integrated")) {
    seurat <- FindNeighbors(seurat, reduction = reduction, dims = 1:30)
    
    cor.data <- as.matrix(seurat@graphs[[names(seurat@graphs)[grep("snn", names(seurat@graphs))]]])
    cor.data <- as(cor.data, "dgCMatrix")
    diag(cor.data) <- 0
    cor_list <- rownames(cor.data)
    graph <- igraph::graph.adjacency(cor.data, mode = "undirected", weighted = TRUE)
    
    pb_inner <- progress_bar$new(total = length(cor_list), format = "  - [:bar] :percent :eta")
    
    distances <- future_map_dfr(1:length(cor_list), ~{
      pb_inner$tick()  # Update inner progress bar
      
      calculate_distance(graph, cor_list, .x)
    })
    
    pb_inner$terminate()  # Close inner progress bar
    
    assign(paste('distances_', reduction, sep=""), distances)
    
    pb_outer$tick()  # Update outer progress bar
  }
  
  metric_minpath <- list(
    pca_integrated = distances_pca_integrated,
    mnn = distances_mnn)
  
  metric_minpath$mnn <- as.matrix(metric_minpath$mnn)
  metric_minpath$pca_integrated <- as.matrix(metric_minpath$pca_integrated)
  
  return(metric_minpath)
}

distance_NW_network_index <- function(seurat, network_index, 
                                      reduction.list = c('mnn',"pca_integrated","mnn_pca_integrated")) {
  
  network_index_list <- rownames(network_index[[names(network_index)[1]]])
  
  refer_list <- intersect(network_index_list, rownames(seurat[seurat$source_batch %in% "refer", ]))
  query_list <- intersect(network_index_list, rownames(seurat[seurat$source_batch %in% "query", ]))
  
  plan(multisession, workers = 20)  
  
  future_values <- list()
  
  for (reduction in reduction.list) {
    
    future_values[[reduction]] <- future({
      cor.data <- as.matrix(network_index[[reduction]])
      # cor.data
      cor.data[cor.data != 0] = 1 / cor.data[cor.data != 0] # distance pendlty
      
      
      graph <- igraph::graph.adjacency(cor.data, mode = "undirected", weighted = TRUE)
      
      path_distance <- matrix(0, nrow = length(refer_list), ncol = length(query_list))
      rownames(path_distance) <- refer_list
      colnames(path_distance) <- query_list
      
      cat(reduction, "Start\n")
      
      for (i in seq_along(query_list)) {
        shortest_paths_details <- igraph::get.shortest.paths(
          graph, from = query_list[i], to = refer_list,
          output = "both", weights = E(graph)$weight
        )
        
        path_distance[, i] <- sapply(shortest_paths_details$epath, function(epath) { length(epath) })
        
        cat(" ", i)
      }
      
      cat("\n")
      
      return(path_distance)
    })
  }
  
  metric_minpath <- future::value(future_values)
  
  return(metric_minpath)
}

batch_process_integration_index = function(seurat, time){
  
  timess = paste(gsub("ss","",time),"ss",sep="")
  
  sstime = paste("ss",gsub("ss","",time),sep="")
  
  seurat.selectgene = get(paste("src.",timess,".integrated.merge.selectgene",sep=""))
  
  seurat = ScaleData(seurat, features = rownames(seurat), split.by = "source_batch")
  seurat = RunPCA(seurat, features = seurat.selectgene)
  seurat = RunUMAP(seurat, dims = 1:30,n.neighbors = 100)
  
  seurat@reductions$umap@cell.embeddings = 
    cbind(seurat@reductions$umap@cell.embeddings[,2],
          seurat@reductions$umap@cell.embeddings[,1])
  colnames(seurat@reductions$umap@cell.embeddings) = c("UMAP_1","UMAP_2")
  
  
  # MNN on Counts
  #---------------
  MNN.res = 
    mnnCorrect(as.matrix(seurat@assays$RNA@data[seurat.selectgene, seurat$Time%in%sstime]),
               as.matrix(seurat@assays$RNA@data[seurat.selectgene, seurat$Time%in%timess]),
               k = 5,cos.norm.out=F)
  seurat@assays$mnnRNA = CreateAssayObject(data = MNN.res@assays@data$corrected)
  seurat@assays$mnnRNA@key = "mnn_"
  seurat = ScaleData(seurat, 
                     rownames(seurat@assays$mnnRNA),
                     assay = "mnnRNA")
  
  seurat@meta.data = seurat@meta.data[!seurat$Time%in%NA,]
  
  # MNN on PCA
  #---------------
  fMNN.res =  fastMNN(as.matrix(seurat@assays$RNA@data[seurat.selectgene, seurat$Time%in%timess]),
                      as.matrix(seurat@assays$RNA@data[seurat.selectgene, seurat$Time%in%sstime&
                                                         seurat$batch%in%1]),
                      as.matrix(seurat@assays$RNA@data[seurat.selectgene, seurat$Time%in%sstime&
                                                         seurat$batch%in%2]),
                      d=50,
                      subset.row = seurat.selectgene)
  seurat@reductions$mnn = seurat@reductions$umap
  seurat@reductions$mnn@cell.embeddings = fMNN.res@assays@data$reconstructed@seed@components
  seurat@reductions$mnn@feature.loadings = fMNN.res@assays@data$reconstructed@seed@rotation
  colnames(seurat@reductions$mnn@cell.embeddings) = paste("MNN",1:50,sep="_")
  colnames(seurat@reductions$mnn@feature.loadings) = paste("MNN",1:50,sep="_")
  seurat@reductions$mnn@key = "MNN_"
  
  # Finde integrated on RNA
  #--------------------------
  anchor.timess.integrated.merge = FindIntegrationAnchors(c(seurat[,seurat$Time%in%timess],
                                                            seurat[,seurat$Time%in%sstime]),
                                                          reduction = "cca",
                                                          assay = c("RNA","RNA"),
                                                          anchor.features = seurat.selectgene,
                                                          dims = 1:30)
  
  seurat.re = IntegrateData(anchorset = anchor.timess.integrated.merge , dims = 1:30)
  DefaultAssay(seurat.re) = "integrated"
  seurat.re = ScaleData(seurat.re, features = rownames(seurat.re))
  seurat.re = RunPCA(seurat.re, features = seurat.selectgene)
  seurat.re = RunUMAP(seurat.re, dims = 1:30, n.neighbors = 100)
  
  seurat@reductions$umap_integrated = seurat.re@reductions$umap
  seurat@reductions$pca_integrated = seurat.re@reductions$pca
  
  seurat@reductions$umap_integrated@assay.used = "RNA"
  seurat@reductions$pca_integrated@assay.used = "RNA"
  
  rm(seurat.re,  anchor.timess.integrated.merge)
  
  # Finde integrated on mnnRNA
  #--------------------------
  anchor.timess.integrated.merge = FindIntegrationAnchors(c(seurat[,seurat$Time%in%timess],
                                                            seurat[,seurat$Time%in%sstime]),
                                                          reduction = "cca",
                                                          assay = c("mnnRNA","mnnRNA"),
                                                          anchor.features = seurat.selectgene,
                                                          dims = 1:30)
  
  seurat.re = IntegrateData(anchorset = anchor.timess.integrated.merge , dims = 1:30)
  DefaultAssay(seurat.re) = "integrated"
  seurat.re = ScaleData(seurat.re, features = rownames(seurat.re))
  seurat.re = RunPCA(seurat.re, features = seurat.selectgene)
  seurat.re = RunUMAP(seurat.re, dims = 1:30, n.neighbors = 100)
  
  seurat@reductions$mnn_umap_integrated = seurat.re@reductions$umap
  seurat@reductions$mnn_pca_integrated = seurat.re@reductions$pca
  seurat@reductions$mnn_umap_integrated@assay.used = "mnnRNA"
  seurat@reductions$mnn_pca_integrated@assay.used = "mnnRNA"
  
  rm(seurat.re,  anchor.timess.integrated.merge)
  
  return(seurat)
}

calculate_tracing_index = function(seurat, time, type = c("Euclidean", "Net_MinPath"), reduction = "mnn"){
  
  seurat.merge = seurat # integrated.merge
  
  if(type == "Euclidean"){
    
    if(time %in% c("12ss","15ss","18ss","21ss","24ss")){
      distance_index1 = get(paste("distance_euclidean_",time,"_index1", sep=""))
      distance_index2 = get(paste("distance_euclidean_",time,"_index2", sep=""))
      distance_index3 = get(paste("distance_euclidean_",time,"_index3", sep=""))
    }else{
      distance_index1 = get(paste("distance_euclidean_",time,"_index1.re", sep=""))
    }
    
  }else if(type == "Net_MinPath"){
    
    if(time %in% c("12ss","15ss","18ss","21ss","24ss")){
      distance_index1 = get(paste("distance_NW_network_",time,"_index1", sep=""))
      distance_index2 = get(paste("distance_NW_network_",time,"_index2", sep=""))
      distance_index3 = get(paste("distance_NW_network_",time,"_index3", sep=""))
    }else{
      distance_index1 = get(paste("distance_NW_network_",time,"_index1", sep=""))
    }
    
  }else{
    return("Change type !")
  }
  
  if(time %in% c("12ss","15ss","18ss","21ss","24ss")){
    seurat_index1 = get(paste("src.",time,".integrated.merge_index1", sep=""))
    seurat_index2 = get(paste("src.",time,".integrated.merge_index2", sep=""))
    seurat_index3 = get(paste("src.",time,".integrated.merge_index3", sep=""))
  }else{
    seurat_index1 = get(paste("src.",time,".integrated.merge_index1", sep=""))
  }
  
  
  index_type = paste("tracing.index_",type,"_",reduction, sep="")
  cat(index_type, '\n')
  
  lineage_list = c("Pax9","Sox2","Nepn",'Wnt5b',"Nkx2_3","Hhex",'Pdx1',"Mnx1","Mnx1GFP")
  
  if(time %in% c("12ss","15ss","18ss","21ss","24ss")){
    seurat_index_list = c("seurat_index1","seurat_index2","seurat_index3")
  }else{
    seurat_index_list = c("seurat_index1")
  }
  
  
  if(type == "Euclidean"){
    
    for(seurat_temp_index in seurat_index_list){
      
      seurat_temp = get(seurat_temp_index)
      refer_list = rownames(seurat_temp@meta.data[seurat_temp$source_batch%in%"refer",])
      nwn_data = get(paste("distance", gsub("seurat","", seurat_temp_index), sep=""))
      nwn_data = nwn_data[[reduction]]
      
      for(lineage in lineage_list){
        
        name_list = rownames(seurat_temp@meta.data[seurat_temp$lineage%in%lineage,])
        seurat_temp[[paste(index_type, lineage, sep="_")]] = 0
        
        for(cell in name_list){
          if(cell %in% c(black_list_summary)){
            next()
          }
          
          cell_list = names(sort(nwn_data[cell, refer_list], decreasing = F)[1:25])
          
          seurat_temp@meta.data[cell_list, paste(index_type, lineage, sep="_")] = 
            seurat_temp@meta.data[cell_list, paste(index_type, lineage, sep="_")] + 1
          
        }
        
      }
      assign(seurat_temp_index, seurat_temp)
    }
    
    
  }else if(type == "Net_MinPath"){
    
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
    
  }else{
    return("Change type !")
  }
  
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

lineage_list = c("Pax9","Sox2","Nepn",'Wnt5b',"Nkx2_3","Hhex",'Pdx1',"Mnx1","Mnx1GFP")

list.netwrok.index.index = c(
  paste("tracing.index_Euclidean_mnn_", lineage_list, sep = ""),
  paste("tracing.index_Euclidean_pca_integrated_", lineage_list, sep = ""),
  paste("tracing.index_Net_MinPath_mnn_", lineage_list, sep = ""),
  paste("tracing.index_Net_MinPath_integrated_", lineage_list, sep = ""))

df.netwrok.index = list()
for(i.replicate in c(1,2,3)){
  #>>> shuffle and size-balance
  for(i.define.batch in c("12ss")){
    table(src.12ss.integrated.merge$source_batch)["refer"] 
    table(src.12ss.integrated.merge$source_batch)["query"]
    query_12ss_num = table(src.12ss.integrated.merge$source_batch)["query"]
    refer_12ss_num = table(src.12ss.integrated.merge$source_batch)["refer"]
    
    # shuffle
    src.12ss.integrated.merge$index_refer_raw = NA
    src.12ss.integrated.merge@meta.data[
      src.12ss.integrated.merge$source_batch%in%"refer",]$index_refer_raw = 
      sample(x = c(1: table(src.12ss.integrated.merge$source_batch)["refer"]),
             size = table(src.12ss.integrated.merge$source_batch)["refer"])
    
    
    # sample by batch size
    size_12ss = floor(refer_12ss_num / (query_12ss_num * 1.5))
    print(size_12ss)
    refer_12ss_num_1 = floor(refer_12ss_num/3)
    refer_12ss_num_2 = floor(refer_12ss_num*2/3)
    
    src.12ss.integrated.merge$index_refer = NA
    src.12ss.integrated.merge@meta.data[
      src.12ss.integrated.merge$index_refer_raw%in%c(1:refer_12ss_num_1),]$index_refer = 1
    src.12ss.integrated.merge@meta.data[
      src.12ss.integrated.merge$index_refer_raw%in%c((refer_12ss_num_1+1):(refer_12ss_num_2)),]$index_refer = 2
    src.12ss.integrated.merge@meta.data[
      src.12ss.integrated.merge$index_refer_raw%in%c((refer_12ss_num_2+1):(refer_12ss_num)),]$index_refer = 3
    table(src.12ss.integrated.merge$index_refer)
    
    
    # set sample for integrated
    src.12ss.integrated.merge_index1 = merge(
      src.12ss.integrated.merge[,colnames(src.12ss.integrated.merge)%in%rownames(
        src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$index_refer%in%1,])],
      src.12ss.integrated.merge[,colnames(src.12ss.integrated.merge)%in%rownames(
        src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$source_batch%in%"query",])])
    src.12ss.integrated.merge_index2 = merge(
      src.12ss.integrated.merge[,colnames(src.12ss.integrated.merge)%in%rownames(
        src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$index_refer%in%2,])],
      src.12ss.integrated.merge[,colnames(src.12ss.integrated.merge)%in%rownames(
        src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$source_batch%in%"query",])])
    src.12ss.integrated.merge_index3 = merge(
      src.12ss.integrated.merge[,colnames(src.12ss.integrated.merge)%in%rownames(
        src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$index_refer%in%3,])],
      src.12ss.integrated.merge[,colnames(src.12ss.integrated.merge)%in%rownames(
        src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$source_batch%in%"query",])])
    
    #---- integrated for batch
    src.12ss.integrated.merge_index1 = batch_process_integration_index(src.12ss.integrated.merge_index1, "12ss")
    src.12ss.integrated.merge_index2 = batch_process_integration_index(src.12ss.integrated.merge_index2, "12ss")
    src.12ss.integrated.merge_index3 = batch_process_integration_index(src.12ss.integrated.merge_index3, "12ss")
    
  }
  
  for(i.define.batch in c("15ss")){
    table(src.15ss.integrated.merge$source_batch)["refer"] 
    table(src.15ss.integrated.merge$source_batch)["query"]
    query_15ss_num = table(src.15ss.integrated.merge$source_batch)["query"]
    refer_15ss_num = table(src.15ss.integrated.merge$source_batch)["refer"]
    
    # shuffle
    src.15ss.integrated.merge$index_refer_raw = NA
    src.15ss.integrated.merge@meta.data[
      src.15ss.integrated.merge$source_batch%in%"refer",]$index_refer_raw = 
      sample(x = c(1: table(src.15ss.integrated.merge$source_batch)["refer"]),
             size = table(src.15ss.integrated.merge$source_batch)["refer"])
    
    
    # sample by batch size
    size_15ss = floor(refer_15ss_num / (query_15ss_num * 1.5))
    print(size_15ss)
    refer_15ss_num_1 = floor(refer_15ss_num/3)
    refer_15ss_num_2 = floor(refer_15ss_num*2/3)
    
    src.15ss.integrated.merge$index_refer = NA
    src.15ss.integrated.merge@meta.data[
      src.15ss.integrated.merge$index_refer_raw%in%c(1:refer_15ss_num_1),]$index_refer = 1
    src.15ss.integrated.merge@meta.data[
      src.15ss.integrated.merge$index_refer_raw%in%c((refer_15ss_num_1+1):(refer_15ss_num_2)),]$index_refer = 2
    src.15ss.integrated.merge@meta.data[
      src.15ss.integrated.merge$index_refer_raw%in%c((refer_15ss_num_2+1):(refer_15ss_num)),]$index_refer = 3
    table(src.15ss.integrated.merge$index_refer)
    
    
    # set sample for integrated
    src.15ss.integrated.merge_index1 = merge(
      src.15ss.integrated.merge[,colnames(src.15ss.integrated.merge)%in%rownames(
        src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$index_refer%in%1,])],
      src.15ss.integrated.merge[,colnames(src.15ss.integrated.merge)%in%rownames(
        src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$source_batch%in%"query",])])
    src.15ss.integrated.merge_index2 = merge(
      src.15ss.integrated.merge[,colnames(src.15ss.integrated.merge)%in%rownames(
        src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$index_refer%in%2,])],
      src.15ss.integrated.merge[,colnames(src.15ss.integrated.merge)%in%rownames(
        src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$source_batch%in%"query",])])
    src.15ss.integrated.merge_index3 = merge(
      src.15ss.integrated.merge[,colnames(src.15ss.integrated.merge)%in%rownames(
        src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$index_refer%in%3,])],
      src.15ss.integrated.merge[,colnames(src.15ss.integrated.merge)%in%rownames(
        src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$source_batch%in%"query",])])
    
    #---- integrated for batch
    src.15ss.integrated.merge_index1 = batch_process_integration_index(src.15ss.integrated.merge_index1, "15ss")
    src.15ss.integrated.merge_index2 = batch_process_integration_index(src.15ss.integrated.merge_index2, "15ss")
    src.15ss.integrated.merge_index3 = batch_process_integration_index(src.15ss.integrated.merge_index3, "15ss")
  }
  
  for(i.define.batch in c("18ss")){
    able(src.18ss.integrated.merge$source_batch)["refer"] 
    table(src.18ss.integrated.merge$source_batch)["query"]
    query_18ss_num = table(src.18ss.integrated.merge$source_batch)["query"]
    refer_18ss_num = table(src.18ss.integrated.merge$source_batch)["refer"]
    
    # shuffle
    src.18ss.integrated.merge$index_refer_raw = NA
    src.18ss.integrated.merge@meta.data[
      src.18ss.integrated.merge$source_batch%in%"refer",]$index_refer_raw = 
      sample(x = c(1: table(src.18ss.integrated.merge$source_batch)["refer"]),
             size = table(src.18ss.integrated.merge$source_batch)["refer"])
    
    
    # sample by batch size
    size_18ss = floor(refer_18ss_num / (query_18ss_num * 1.5))
    print(size_18ss)
    refer_18ss_num_1 = floor(refer_18ss_num/3)
    refer_18ss_num_2 = floor(refer_18ss_num*2/3)
    
    src.18ss.integrated.merge$index_refer = NA
    src.18ss.integrated.merge@meta.data[
      src.18ss.integrated.merge$index_refer_raw%in%c(1:refer_18ss_num_1),]$index_refer = 1
    src.18ss.integrated.merge@meta.data[
      src.18ss.integrated.merge$index_refer_raw%in%c((refer_18ss_num_1+1):(refer_18ss_num_2)),]$index_refer = 2
    src.18ss.integrated.merge@meta.data[
      src.18ss.integrated.merge$index_refer_raw%in%c((refer_18ss_num_2+1):(refer_18ss_num)),]$index_refer = 3
    table(src.18ss.integrated.merge$index_refer)
    
    
    # set sample for integrated
    src.18ss.integrated.merge_index1 = merge(
      src.18ss.integrated.merge[,colnames(src.18ss.integrated.merge)%in%rownames(
        src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$index_refer%in%1,])],
      src.18ss.integrated.merge[,colnames(src.18ss.integrated.merge)%in%rownames(
        src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$source_batch%in%"query",])])
    src.18ss.integrated.merge_index2 = merge(
      src.18ss.integrated.merge[,colnames(src.18ss.integrated.merge)%in%rownames(
        src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$index_refer%in%2,])],
      src.18ss.integrated.merge[,colnames(src.18ss.integrated.merge)%in%rownames(
        src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$source_batch%in%"query",])])
    src.18ss.integrated.merge_index3 = merge(
      src.18ss.integrated.merge[,colnames(src.18ss.integrated.merge)%in%rownames(
        src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$index_refer%in%3,])],
      src.18ss.integrated.merge[,colnames(src.18ss.integrated.merge)%in%rownames(
        src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$source_batch%in%"query",])])
    
    
    #---- integrated for batch
    src.18ss.integrated.merge_index1 = batch_process_integration_index(src.18ss.integrated.merge_index1, "18ss")
    src.18ss.integrated.merge_index2 = batch_process_integration_index(src.18ss.integrated.merge_index2, "18ss")
    src.18ss.integrated.merge_index3 = batch_process_integration_index(src.18ss.integrated.merge_index3, "18ss")
    
  }
  
  for(i.define.batch in c("21ss")){
    table(src.21ss.integrated.merge$source_batch)["refer"] 
    table(src.21ss.integrated.merge$source_batch)["query"]
    query_21ss_num = table(src.21ss.integrated.merge$source_batch)["query"]
    refer_21ss_num = table(src.21ss.integrated.merge$source_batch)["refer"]
    
    # shuffle
    src.21ss.integrated.merge$index_refer_raw = NA
    src.21ss.integrated.merge@meta.data[
      src.21ss.integrated.merge$source_batch%in%"refer",]$index_refer_raw = 
      sample(x = c(1: table(src.21ss.integrated.merge$source_batch)["refer"]),
             size = table(src.21ss.integrated.merge$source_batch)["refer"])
    
    
    # sample by batch size
    size_21ss = floor(refer_21ss_num / (query_21ss_num * 1.5))
    print(size_21ss)
    refer_21ss_num_1 = floor(refer_21ss_num/3)
    refer_21ss_num_2 = floor(refer_21ss_num*2/3)
    
    src.21ss.integrated.merge$index_refer = NA
    src.21ss.integrated.merge@meta.data[
      src.21ss.integrated.merge$index_refer_raw%in%c(1:refer_21ss_num_1),]$index_refer = 1
    src.21ss.integrated.merge@meta.data[
      src.21ss.integrated.merge$index_refer_raw%in%c((refer_21ss_num_1+1):(refer_21ss_num_2)),]$index_refer = 2
    src.21ss.integrated.merge@meta.data[
      src.21ss.integrated.merge$index_refer_raw%in%c((refer_21ss_num_2+1):(refer_21ss_num)),]$index_refer = 3
    table(src.21ss.integrated.merge$index_refer)
    
    
    # set sample for integrated
    src.21ss.integrated.merge_index1 = merge(
      src.21ss.integrated.merge[,colnames(src.21ss.integrated.merge)%in%rownames(
        src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$index_refer%in%1,])],
      src.21ss.integrated.merge[,colnames(src.21ss.integrated.merge)%in%rownames(
        src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$source_batch%in%"query",])])
    src.21ss.integrated.merge_index2 = merge(
      src.21ss.integrated.merge[,colnames(src.21ss.integrated.merge)%in%rownames(
        src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$index_refer%in%2,])],
      src.21ss.integrated.merge[,colnames(src.21ss.integrated.merge)%in%rownames(
        src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$source_batch%in%"query",])])
    src.21ss.integrated.merge_index3 = merge(
      src.21ss.integrated.merge[,colnames(src.21ss.integrated.merge)%in%rownames(
        src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$index_refer%in%3,])],
      src.21ss.integrated.merge[,colnames(src.21ss.integrated.merge)%in%rownames(
        src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$source_batch%in%"query",])])
    
    
    #---- integrated for batch
    src.21ss.integrated.merge_index1 = batch_process_integration_index(src.21ss.integrated.merge_index1, "21ss")
    src.21ss.integrated.merge_index2 = batch_process_integration_index(src.21ss.integrated.merge_index2, "21ss")
    src.21ss.integrated.merge_index3 = batch_process_integration_index(src.21ss.integrated.merge_index3, "21ss")
  }
  
  for(i.define.batch in c("24ss")){
    table(src.24ss.integrated.merge$source_batch)["refer"] 
    table(src.24ss.integrated.merge$source_batch)["query"]
    query_24ss_num = table(src.24ss.integrated.merge$source_batch)["query"]
    refer_24ss_num = table(src.24ss.integrated.merge$source_batch)["refer"]
    
    # shuffle
    src.24ss.integrated.merge$index_refer_raw = NA
    src.24ss.integrated.merge@meta.data[
      src.24ss.integrated.merge$source_batch%in%"refer",]$index_refer_raw = 
      sample(x = c(1: table(src.24ss.integrated.merge$source_batch)["refer"]),
             size = table(src.24ss.integrated.merge$source_batch)["refer"])
    
    
    # sample by batch size
    size_24ss = floor(refer_24ss_num / (query_24ss_num * 1.5))
    print(size_24ss)
    refer_24ss_num_1 = floor(refer_24ss_num/3)
    refer_24ss_num_2 = floor(refer_24ss_num*2/3)
    
    src.24ss.integrated.merge$index_refer = NA
    src.24ss.integrated.merge@meta.data[
      src.24ss.integrated.merge$index_refer_raw%in%c(1:refer_24ss_num_1),]$index_refer = 1
    src.24ss.integrated.merge@meta.data[
      src.24ss.integrated.merge$index_refer_raw%in%c((refer_24ss_num_1+1):(refer_24ss_num_2)),]$index_refer = 2
    src.24ss.integrated.merge@meta.data[
      src.24ss.integrated.merge$index_refer_raw%in%c((refer_24ss_num_2+1):(refer_24ss_num)),]$index_refer = 3
    table(src.24ss.integrated.merge$index_refer)
    
    
    # set sample for integrated
    src.24ss.integrated.merge_index1 = merge(
      src.24ss.integrated.merge[,colnames(src.24ss.integrated.merge)%in%rownames(
        src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$index_refer%in%1,])],
      src.24ss.integrated.merge[,colnames(src.24ss.integrated.merge)%in%rownames(
        src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$source_batch%in%"query",])])
    src.24ss.integrated.merge_index2 = merge(
      src.24ss.integrated.merge[,colnames(src.24ss.integrated.merge)%in%rownames(
        src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$index_refer%in%2,])],
      src.24ss.integrated.merge[,colnames(src.24ss.integrated.merge)%in%rownames(
        src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$source_batch%in%"query",])])
    src.24ss.integrated.merge_index3 = merge(
      src.24ss.integrated.merge[,colnames(src.24ss.integrated.merge)%in%rownames(
        src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$index_refer%in%3,])],
      src.24ss.integrated.merge[,colnames(src.24ss.integrated.merge)%in%rownames(
        src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$source_batch%in%"query",])])
    
    
    #---- integrated for batch
    src.24ss.integrated.merge_index1 = batch_process_integration_index(src.24ss.integrated.merge_index1, "24ss")
    src.24ss.integrated.merge_index2 = batch_process_integration_index(src.24ss.integrated.merge_index2, "24ss")
    src.24ss.integrated.merge_index3 = batch_process_integration_index(src.24ss.integrated.merge_index3, "24ss")
  }
  
  for(i.define.batch in c("27ss")){
    table(src.27ss.integrated.merge$source_batch)["refer"] 
    table(src.27ss.integrated.merge$source_batch)["query"]
    query_27ss_num = table(src.27ss.integrated.merge$source_batch)["query"]
    refer_27ss_num = table(src.27ss.integrated.merge$source_batch)["refer"]
    
    # shuffle
    src.27ss.integrated.merge$index_refer_raw = NA
    src.27ss.integrated.merge@meta.data[
      src.27ss.integrated.merge$source_batch%in%"refer",]$index_refer_raw = 
      sample(x = c(1: table(src.27ss.integrated.merge$source_batch)["refer"]),
             size = table(src.27ss.integrated.merge$source_batch)["refer"])
    
    src.27ss.integrated.merge$index_refer = NA
    src.27ss.integrated.merge@meta.data[colnames(src.27ss.integrated),]$index_refer = 1
    
    # set sample for integrated
    src.27ss.integrated.merge_index1 = src.27ss.integrated.merge
    src.27ss.integrated.merge.selectgene = rownames(src.27ss.integrated@reductions$pca@feature.loadings)
    src.27ss.integrated.merge_index1@reductions$mnn_pca_integrated = 
      src.27ss.integrated.merge_index1@reductions$mnnpca_integrated
    
    #---- integrated for batch
    src.27ss.integrated.merge_index1 = batch_process_integration_index(src.27ss.integrated.merge_index1, "27ss")
  }
  
  
  #>>> batch_process_integration for each index
  for(seurat in c(# "src.12ss.integrated.merge_index1",
    "src.12ss.integrated.merge_index2","src.12ss.integrated.merge_index3",
    "src.15ss.integrated.merge_index1","src.15ss.integrated.merge_index2","src.15ss.integrated.merge_index3",
    "src.18ss.integrated.merge_index1","src.18ss.integrated.merge_index2","src.18ss.integrated.merge_index3",
    "src.21ss.integrated.merge_index1","src.21ss.integrated.merge_index2","src.21ss.integrated.merge_index3",
    "src.24ss.integrated.merge_index1","src.24ss.integrated.merge_index2","src.24ss.integrated.merge_index3",
    "src.27ss.integrated.merge_index1")){
    
    seurat.re = get(seurat)
    time = gsub(".*src\\.([0-9]+)ss.*", "\\1ss", seurat)
    
    seurat.re = batch_process_integration_index(seurat.re, time)
    assign(seurat, seurat.re)
  }
  
  #>>> extract network (try)
  for(i.time in c("12ss","15ss","18ss","21ss","24ss","27ss")){
    if(i.time == "27ss"){
      assign(paste("network_", i.time, "_index1", sep = ""), 
             network_extract_construct(get(paste("src.", i.time, ".integrated.merge_index1",sep = ""))))
    }else{
      assign(paste("network_", i.time, "_index1", sep = ""), 
             network_extract_construct(get(paste("src.", i.time, ".integrated.merge_index1",sep = ""))))
      assign(paste("network_", i.time, "_index2", sep = ""), 
             network_extract_construct(get(paste("src.", i.time, ".integrated.merge_index2",sep = ""))))
      assign(paste("network_", i.time, "_index3", sep = ""), 
             network_extract_construct(get(paste("src.", i.time, ".integrated.merge_index3",sep = ""))))
    }
  }
  
  #>>> set metadata:
  for(i.time in c("9ss","12ss","15ss","18ss","21ss","24ss","27ss")){
    assign(paste("metadata_src.", i.time, ".integrated.merge",sep = ""),
           get(paste("src.", i.time, ".integrated.merge",sep = ""))@meta.data)
  }
  
  #>> distance_NW_network for index
  #-----------------------------
  distance_NW_network_12ss_index1 = distance_NW_network_index(metadata_src.12ss.integrated.merge, network_12ss_index1)
  distance_NW_network_12ss_index2 = distance_NW_network_index(metadata_src.12ss.integrated.merge, network_12ss_index2)
  distance_NW_network_12ss_index3 = distance_NW_network_index(metadata_src.12ss.integrated.merge, network_12ss_index3)
  
  distance_NW_network_15ss_index1 = distance_NW_network_index(src.15ss.integrated.merge, network_15ss_index1)
  distance_NW_network_15ss_index2 = distance_NW_network_index(src.15ss.integrated.merge, network_15ss_index2)
  distance_NW_network_15ss_index3 = distance_NW_network_index(src.15ss.integrated.merge, network_15ss_index3)
  
  distance_NW_network_18ss_index1 = distance_NW_network_index(metadata_src.18ss.integrated.merge, network_18ss_index1)
  distance_NW_network_18ss_index2 = distance_NW_network_index(metadata_src.18ss.integrated.merge, network_18ss_index2)
  distance_NW_network_18ss_index3 = distance_NW_network_index(metadata_src.18ss.integrated.merge, network_18ss_index3)
  
  distance_NW_network_21ss_index1 = distance_NW_network_index(metadata_src.21ss.integrated.merge, network_21ss_index1)
  distance_NW_network_21ss_index2 = distance_NW_network_index(metadata_src.21ss.integrated.merge, network_21ss_index2)
  distance_NW_network_21ss_index3 = distance_NW_network_index(metadata_src.21ss.integrated.merge, network_21ss_index3)
  
  distance_NW_network_24ss_index1 = distance_NW_network_index(metadata_src.24ss.integrated.merge, network_24ss_index1)
  distance_NW_network_24ss_index2 = distance_NW_network_index(metadata_src.24ss.integrated.merge, network_24ss_index2)
  distance_NW_network_24ss_index3 = distance_NW_network_index(metadata_src.24ss.integrated.merge, network_24ss_index3)
  
  distance_NW_network_27ss_index1 = distance_NW_network_index(metadata_src.27ss.integrated.merge, network_27ss_index1)
  #-----------------------------
  
  type = c("Euclidean", "Net_MinPath")
  reduction = c("mnn", "pca_integrated")
  
  #--- calculate_tracing_index
  for(time in c("ss12","ss15","ss18","ss21","ss24","ss27")){
    # if(time %in% c("ss9","ss12","ss15","ss18","ss21","ss24","ss27")){next()}
    seurat_name = paste("src.",gsub("ss","",time),"ss.integrated.merge", sep="")
    seurat = get(seurat_name)
    for(type in c("Euclidean", "Net_MinPath")){
      for(reduction in c("mnn", "pca_integrated")){
        seurat = calculate_tracing_index(
          seurat, 
          time = paste(gsub("ss","",time),"ss",sep=""), 
          type = type, 
          reduction = reduction)
      }
    }
    assign(seurat_name, seurat)
  }
  
  for(time in c("ss12","ss15","ss18","ss21","ss24","ss27")){
    seurat = get(paste("src.",gsub("ss","",time),"ss.integrated.merge", sep=""))
    df.netwrok.index[[paste(i.time,"_",i.replicate,sep="")]] = 
      seurat@meta.data[, list.netwrok.index.index]
  }
  
}

#> calulate mean and log(x+1) for df.netwrok.index
for(time in c('12ss','15ss','18ss','21ss','24ss','27ss')){
  seurat = get(paste("src.",gsub("ss","",time),"ss.integrated.merge", sep=""))
  for(i.index in c(list.netwrok.index.index)){
    seurat@meta.data[,i.index] = log(
      apply(cbind(
        df.netwrok.index[[paste(i.time,"_1",sep="")]][,i.index],
        df.netwrok.index[[paste(i.time,"_2",sep="")]][,i.index],
        df.netwrok.index[[paste(i.time,"_3",sep="")]][,i.index]), 1, mean) + 1
    )
  }
  assign(paste("src.",gsub("ss","",time),"ss.integrated.merge", sep=""), seurat)
}

type = c("Euclidean", "Net_MinPath")
reduction = c("mnn", "pca_integrated", "mnn_pca_integrated")
lineage_list = c("Pax9","Sox2","Nepn",'Wnt5b',"Nkx2_3","Hhex",'Pdx1',"Mnx1","Mnx1GFP")

#>> plot-pattern 
for(time in c('12ss','15ss','18ss','21ss','24ss','27ss')){
  for( i in c(1:length(type))){
    for( j in c(1:length(reduction))){
      
      for( k in c(1:length(lineage_list))){
        
        library("scales")
        
        typename = paste("tracing.index", type[i], reduction[j], lineage_list[k], sep = "_")
        filename = paste("integrated/distance_png/",type[i],"/", reduction[j],"/",
                         "Integrated_",typename,"_",time,".png", sep="")
        srcname = paste("src.",time,".integrated.merge", sep = "")
        print(filename)
        
        seurat = get(srcname)
        metadata = cbind(seurat@meta.data, seurat@reductions$umap@cell.embeddings)
        metadata[["tracing.index"]] = metadata[[typename]]
        
        p_temp = ggplot() +
          geom_point(metadata[metadata[["tracing.index"]]%in%0,],
                     mapping = aes(x=UMAP_1, y =UMAP_2, color=tracing.index),  size=3) +
          geom_point(metadata[metadata[["tracing.index"]]>0,],
                     mapping = aes(x=UMAP_1, y =UMAP_2, color=tracing.index),  size=3) +
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
#===============================================================================


#===============================================================================
#>>>>> 4.3 define tracing code
#===============================================================================
for(i.define.time in c("9ss")){
  #>> set cell type as tracing code for 9SS
  src.9ss.integrated$cluster.extract.v1.1 = src.9ss.integrated$cluster.endoderm
}

for(i.define.time in c("12ss")){
  src.12ss.integrated.merge@reductions$umap_integrated@cell.embeddings = cbind(
    src.12ss.integrated.merge@reductions$umap_integrated@cell.embeddings[,1],
    -src.12ss.integrated.merge@reductions$umap_integrated@cell.embeddings[,2])
  colnames(src.12ss.integrated.merge@reductions$umap_integrated@cell.embeddings) = c("UMAP_1",'UMAP_2')
  
  DimPlot(src.12ss.integrated.merge, 
          group.by = "lineage", reduction = "umap_integrated",
          cols = color.lineage, na.value = "#eeeeee")
  
  src.12ss.integrated = FindNeighbors(src.12ss.integrated, reduction = "pca", dims = 1:30)
  src.12ss.integrated = FindClusters(src.12ss.integrated, resolution = 5)
  src.12ss.integrated = FindClusters(src.12ss.integrated, resolution = 10)
  
  src.12ss.integrated.merge$cluster_snn_res.5 = NA
  src.12ss.integrated.merge@meta.data[rownames(src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$Time%in%"ss12",]),]$cluster_snn_res.5 = 
    src.12ss.integrated@meta.data[rownames(src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$Time%in%"ss12",]),]$RNA_snn_res.5 
  src.12ss.integrated.merge$cluster_snn_res.10 = NA
  src.12ss.integrated.merge@meta.data[rownames(src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$Time%in%"ss12",]),]$cluster_snn_res.10 = 
    src.12ss.integrated@meta.data[rownames(src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$Time%in%"ss12",]),]$RNA_snn_res.10 
  
  src.12ss.integrated.merge$cluster.extract.v1.1 = NA
  src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$cluster_snn_res.5%in%c(32),]$cluster.extract.v1.1 = "FG.2"
  src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$cluster_snn_res.5%in%c(7,24),]$cluster.extract.v1.1 = "FG.2"
  src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$cluster_snn_res.5%in%c(21,28),]$cluster.extract.v1.1 = "FG.3"
  src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$cluster_snn_res.5%in%c(11,40),]$cluster.extract.v1.1 = "FG.4"
  src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$cluster_snn_res.5%in%c(2,26),]$cluster.extract.v1.1 = "FG.5"
  src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$cluster_snn_res.5%in%c(27),]$cluster.extract.v1.1 = "FG.6"
  src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$cluster_snn_res.5%in%c(20,13,0)&
                                        src.12ss.integrated$RNA_snn_res.10%in%c(3,16,80,113),]$cluster.extract.v1.1 = "AL.1" 
  src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$cluster_snn_res.5%in%c(20,13,0)&
                                        !src.12ss.integrated$RNA_snn_res.10%in%c(3,16,80,113),]$cluster.extract.v1.1 = "AL.2" 
  src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$cluster_snn_res.5%in%c(22,31,33,38),]$cluster.extract.v1.1 = "AL.3" 
  src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$cluster_snn_res.5%in%c(8,17,30,35),]$cluster.extract.v1.1 = "MG.1"
  src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$cluster_snn_res.5%in%c(3,12,29,25,19),]$cluster.extract.v1.1 = "MG.2"
  src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$cluster_snn_res.5%in%c(1,4,36,14,34),]$cluster.extract.v1.1 = "MG.3" 
  src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$cluster_snn_res.5%in%c(23,9,17,35),]$cluster.extract.v1.1 = "HG.1" 
  src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$cluster_snn_res.5%in%c(37,15,39),]$cluster.extract.v1.1 = "HG.2" 

  #>> need more resolution to distinguish or cover by different LAI patterns
  a = FNN::knn(src.12ss.integrated.merge@reductions$umap_integrated@cell.embeddings[
    rownames(src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$Time%in%"ss12"&
                                                   !src.12ss.integrated.merge$cluster_snn_res.5%in%c(5,6,10,16,18),]),],
    src.12ss.integrated.merge@reductions$umap_integrated@cell.embeddings[
      rownames(src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$Time%in%"ss12"&
                                                     src.12ss.integrated.merge$cluster_snn_res.5%in%c(5,6,10,16,18),]),],
    src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$Time%in%"ss12"&
                                          !src.12ss.integrated.merge$cluster_snn_res.5%in%c(5,6,10,16,18),]$cluster.extract.v1.1, k = 10)
  src.12ss.integrated.merge@meta.data[src.12ss.integrated.merge$cluster_snn_res.5%in%c(5,6,10,16,18),]$cluster.extract.v1.1 = as.character(a)
  
}

for(i.define.time in c("15ss")){
  src.15ss.integrated.merge@reductions$umap_integrated@cell.embeddings = cbind(
    src.15ss.integrated.merge@reductions$umap_integrated@cell.embeddings[,1],
    -src.15ss.integrated.merge@reductions$umap_integrated@cell.embeddings[,2])
  colnames(src.15ss.integrated.merge@reductions$umap_integrated@cell.embeddings) = c("UMAP_1",'UMAP_2')
  
  DimPlot(src.15ss.integrated.merge, 
          group.by = "lineage", reduction = "umap_integrated",
          cols = color.lineage, na.value = "#eeeeee")
  
  src.15ss.integrated = FindNeighbors(src.15ss.integrated, reduction = "pca", dims = 1:30)
  src.15ss.integrated = FindClusters(src.15ss.integrated, resolution = 5)
  src.15ss.integrated.merge$cluster_snn_res.5 = NA
  src.15ss.integrated.merge@meta.data[rownames(src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$Time%in%"ss15",]),]$cluster_snn_res.5 = 
    src.15ss.integrated@meta.data[rownames(src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$Time%in%"ss15",]),]$RNA_snn_res.5 
  
  
  src.15ss.integrated.merge$cluster.extract.v1.1 = NA
  src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$cluster_snn_res.5%in%c(23),]$cluster.extract.v1.1 = "FG.1"
  src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$cluster_snn_res.5%in%c(6,14),]$cluster.extract.v1.1 = "FG.2"
  src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$cluster_snn_res.5%in%c(21,37),]$cluster.extract.v1.1 = "FG.3"
  src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$cluster_snn_res.5%in%c(13,31),]$cluster.extract.v1.1 = "FG.4"
  src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$cluster_snn_res.5%in%c(15,34),]$cluster.extract.v1.1 = "FG.5"
  src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$cluster_snn_res.5%in%c(30,40),]$cluster.extract.v1.1 = "FG.6"
  src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$cluster_snn_res.5%in%c(2,11,16,22,26),]$cluster.extract.v1.1 = "FG.4-AL.1/2/3" 
  src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$cluster_snn_res.5%in%c(25,28,29,35),]$cluster.extract.v1.1 = "AL.3" 
  src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$cluster_snn_res.5%in%c(0,33),]$cluster.extract.v1.1 = "MG.1"
  src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$cluster_snn_res.5%in%c(1,3,7,20,41,17,19,42),]$cluster.extract.v1.1 = "MG.2"
  src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$cluster_snn_res.5%in%c(8,10,39,12,43),]$cluster.extract.v1.1 = "MG.3"
  src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$cluster_snn_res.5%in%c(9,24,36),]$cluster.extract.v1.1 = "HG.1"
  src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$cluster_snn_res.5%in%c(4,38,44),]$cluster.extract.v1.1 = "HG.2"
  

  #>> need more resolution to distinguish or cover by different LAI patterns
  a = FNN::knn(src.15ss.integrated.merge@reductions$umap_integrated@cell.embeddings[
    rownames(src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$Time%in%"ss15"&
                                                   !src.15ss.integrated.merge$cluster_snn_res.5%in%c(5,18,27,32),]),],
    src.15ss.integrated.merge@reductions$umap_integrated@cell.embeddings[
      rownames(src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$Time%in%"ss15"&
                                                     src.15ss.integrated.merge$cluster_snn_res.5%in%c(5,18,27,32),]),],
    src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$Time%in%"ss15"&
                                          !src.15ss.integrated.merge$cluster_snn_res.5%in%c(5,18,27,32),]$cluster.extract.v1.1, k = 10)
  src.15ss.integrated.merge@meta.data[src.15ss.integrated.merge$cluster_snn_res.5%in%c(5,18,27,32),]$cluster.extract.v1.1 = as.character(a)
}

for(i.define.time in c("18ss")){
  DimPlot(src.18ss.integrated.merge, 
          group.by = "lineage", reduction = "umap_integrated",
          cols = color.lineage, na.value = "#eeeeee")
  DimPlot(src.18ss.integrated.merge, 
          group.by = "cluster.extract.v1.0", reduction = "umap_integrated",
          cols = cluster.endoderm.color.v5, na.value = "#eeeeee")
  
  src.18ss.integrated = FindNeighbors(src.18ss.integrated, reduction = "pca", dims = 1:30)
  src.18ss.integrated = FindClusters(src.18ss.integrated, resolution = 5)
  src.18ss.integrated.merge$cluster_snn_res.5 = NA
  src.18ss.integrated.merge@meta.data[rownames(src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$Time%in%"ss18",]),]$cluster_snn_res.5 = 
    src.18ss.integrated@meta.data[rownames(src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$Time%in%"ss18",]),]$RNA_snn_res.5 
  

  src.18ss.integrated.merge$cluster.extract.v1.1 = NA
  src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$cluster_snn_res.5%in%c(29),]$cluster.extract.v1.1 = "FG.1"
  src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$cluster_snn_res.5%in%c(12,0,20,34,35,42),]$cluster.extract.v1.1 = "FG.2"
  src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$cluster_snn_res.5%in%c(16,26,44),]$cluster.extract.v1.1 = "FG.5"
  src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$cluster_snn_res.5%in%c(40),]$cluster.extract.v1.1 = "FG.6"
  src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$cluster_snn_res.5%in%c(7,8,9,23),]$cluster.extract.v1.1 = "FG.4-AL.1/2/3" # 16?
  src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$cluster_snn_res.5%in%c(6,13,39,43),]$cluster.extract.v1.1 = "HG.1"
  src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$cluster_snn_res.5%in%c(5,18,25),]$cluster.extract.v1.1 = "HG.2"
  
  src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$cluster_snn_res.5%in%c(35,42),]$cluster.extract.v1.1 = "FG.3"
  src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$cluster_snn_res.5%in%c(22,24,37),]$cluster.extract.v1.1 = "FG.3/4"
  src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$cluster_snn_res.5%in%c(14,32),]$cluster.extract.v1.1 = "FG.4"
  src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$cluster_snn_res.5%in%c(10),]$cluster.extract.v1.1 = "AL.3"
  src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$cluster_snn_res.5%in%c(2,28),]$cluster.extract.v1.1 = "AL.3-MG.2" # 12?
  
  src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$cluster_snn_res.5%in%c(19,27,30),]$cluster.extract.v1.1 = "MG.1"
  src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$cluster_snn_res.5%in%c(11,17,31,36),]$cluster.extract.v1.1 = "MG.2"
  src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$cluster_snn_res.5%in%c(3,21,33,41),]$cluster.extract.v1.1 = "MG.3"

  
  #>> need more resolution to distinguish or cover by different LAI patterns
  a = FNN::knn(src.18ss.integrated.merge@reductions$umap_integrated@cell.embeddings[
    rownames(src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$Time%in%"ss18"&
                                                   !src.18ss.integrated.merge$cluster_snn_res.5%in%c(1,4,15,38),]),],
    src.18ss.integrated.merge@reductions$umap_integrated@cell.embeddings[
      rownames(src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$Time%in%"ss18"&
                                                     src.18ss.integrated.merge$cluster_snn_res.5%in%c(1,4,15,38),]),],
    src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$Time%in%"ss18"&
                                          !src.18ss.integrated.merge$cluster_snn_res.5%in%c(1,4,15,38),]$cluster.extract.v1.1, k = 10)
  src.18ss.integrated.merge@meta.data[src.18ss.integrated.merge$cluster_snn_res.5%in%c(1,4,15,38),]$cluster.extract.v1.1 = as.character(a)
}

for(i.define.time in c("21ss")){
  src.21ss.integrated.merge@reductions$umap_integrated@cell.embeddings = cbind(
    src.21ss.integrated.merge@reductions$umap_integrated@cell.embeddings[,1],
    -src.21ss.integrated.merge@reductions$umap_integrated@cell.embeddings[,2])
  colnames(src.21ss.integrated.merge@reductions$umap_integrated@cell.embeddings) = c("UMAP_1",'UMAP_2')
  
  DimPlot(src.21ss.integrated.merge, 
          group.by = "lineage", reduction = "umap_integrated",
          cols = color.lineage, na.value = "#eeeeee")
  
  src.21ss.integrated.merge$cluster_snn_res.5 = NA
  src.21ss.integrated.merge@meta.data[rownames(src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$Time%in%"ss21",]),]$cluster_snn_res.5 = 
    src.21ss.integrated@meta.data[rownames(src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$Time%in%"ss21",]),]$RNA_snn_res.5 
  
  src.21ss.integrated.merge$cluster.extract.v1.1 = NA
  src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$cluster_snn_res.5%in%c(27,35),]$cluster.extract.v1.1 = "FG.5"
  src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$cluster_snn_res.5%in%c(33,41,14),]$cluster.extract.v1.1 = "FG.2"
  src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$cluster_snn_res.5%in%c(34,26),]$cluster.extract.v1.1 = "FG.3"
  src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$cluster_snn_res.5%in%c(23),]$cluster.extract.v1.1 = "AL.3"
  src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$cluster_snn_res.5%in%c(0,10,18,22,43),]$cluster.extract.v1.1 = "FG.4-AL.1/2/3" 
  src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$cluster_snn_res.5%in%c(11,29,44),]$cluster.extract.v1.1 = "AL.3-MG.2" 
  
  src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$cluster_snn_res.5%in%c(15,19,25),]$cluster.extract.v1.1 = "FG.3/4"
  src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$cluster_snn_res.5%in%c(3,4,39,40),]$cluster.extract.v1.1 = "FG.4"
  src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$cluster_snn_res.5%in%c(24),]$cluster.extract.v1.1 = "FG.1"
  src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$cluster_snn_res.5%in%c(20),]$cluster.extract.v1.1 = "FG.6"
  
  src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$cluster_snn_res.5%in%c(13,36),]$cluster.extract.v1.1 = "MG.1"
  src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$cluster_snn_res.5%in%c(2,5,6,17,37),]$cluster.extract.v1.1 = "MG.2"
  src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$cluster_snn_res.5%in%c(8,12,16,38,42),]$cluster.extract.v1.1 = "MG.3" 
  src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$cluster_snn_res.5%in%c(28,30),]$cluster.extract.v1.1 = "HG.1"
  src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$cluster_snn_res.5%in%c(21,31),]$cluster.extract.v1.1 = "HG.2" # 2?
  
  
  #>> need more resolution to distinguish or cover by different LAI patterns
  a = FNN::knn(src.21ss.integrated.merge@reductions$umap_integrated@cell.embeddings[
    rownames(src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$Time%in%"ss21"&
                                                   !src.21ss.integrated.merge$cluster_snn_res.5%in%c(1,7,9,32),]),],
    src.21ss.integrated.merge@reductions$umap_integrated@cell.embeddings[
      rownames(src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$Time%in%"ss21"&
                                                     src.21ss.integrated.merge$cluster_snn_res.5%in%c(1,7,9,32),]),],
    src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$Time%in%"ss21"&
                                          !src.21ss.integrated.merge$cluster_snn_res.5%in%c(1,7,9,32),]$cluster.extract.v1.1, k = 10)
  src.21ss.integrated.merge@meta.data[src.21ss.integrated.merge$cluster_snn_res.5%in%c(1,7,9,32),]$cluster.extract.v1.1 = as.character(a)
  
}

for(i.define.time in c("24ss")){
  src.24ss.integrated.merge@reductions$umap_integrated@cell.embeddings = cbind(
    -src.24ss.integrated.merge@reductions$umap_integrated@cell.embeddings[,1],
    -src.24ss.integrated.merge@reductions$umap_integrated@cell.embeddings[,2])
  colnames(src.24ss.integrated.merge@reductions$umap_integrated@cell.embeddings) = c("UMAP_1",'UMAP_2')
  
  DimPlot(src.24ss.integrated.merge, 
          group.by = "lineage", reduction = "umap_integrated",
          cols = color.lineage, na.value = "#eeeeee")
  DimPlot(src.24ss.integrated.merge, 
          group.by = "cluster.extract.v1.0", reduction = "umap_integrated",
          cols = cluster.endoderm.color.v5, na.value = "#eeeeee")
  
  src.24ss.integrated.merge$cluster_snn_res.5 = NA
  src.24ss.integrated.merge@meta.data[rownames(src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$Time%in%"ss24",]),]$cluster_snn_res.5 = 
    src.24ss.integrated@meta.data[rownames(src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$Time%in%"ss24",]),]$RNA_snn_res.5
  
  
  src.24ss.integrated.merge$cluster.extract.v1.1 = NA
  src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$cluster_snn_res.5%in%c(34,35),]$cluster.extract.v1.1 = "FG.5"
  src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$cluster_snn_res.5%in%c(14),]$cluster.extract.v1.1 = "FG.6"
  src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$cluster_snn_res.5%in%c(27),]$cluster.extract.v1.1 = "FG.1"
  src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$cluster_snn_res.5%in%c(3,38),]$cluster.extract.v1.1 = "FG.3"
  src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$cluster_snn_res.5%in%c(7,13,23,24,30),]$cluster.extract.v1.1 = "FG.2"
  src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$cluster_snn_res.5%in%c(22,1),]$cluster.extract.v1.1 = "FG.3/4"
  src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$cluster_snn_res.5%in%c(16),]$cluster.extract.v1.1 = "FG.4" # 7
  
  src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$cluster_snn_res.5%in%c(0,10,31),]$cluster.extract.v1.1 = "FG.4-MG.1/3"
  src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$cluster_snn_res.5%in%c(9,42),]$cluster.extract.v1.1 = "MG.3" 
  src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$cluster_snn_res.5%in%c(21),]$cluster.extract.v1.1 = "MG.1/2/3" 
  src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$cluster_snn_res.5%in%c(2,4,39),]$cluster.extract.v1.1 = "MG.2/3-HG.1" 
  src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$cluster_snn_res.5%in%c(5,25,43),]$cluster.extract.v1.1 = "MG.2/3"
  
  src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$cluster_snn_res.5%in%c(12,15,18,32,40),]$cluster.extract.v1.1 = "FG.4-AL.1/2/3"
  src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$cluster_snn_res.5%in%c(17,36),]$cluster.extract.v1.1 = "HG.1/2"
  src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$cluster_snn_res.5%in%c(20),]$cluster.extract.v1.1 = "AL.3-MG.2"
  src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$cluster_snn_res.5%in%c(29),]$cluster.extract.v1.1 = "HG.1"
  src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$cluster_snn_res.5%in%c(11,33,37),]$cluster.extract.v1.1 = "HG.2" # 20
  src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$cluster_snn_res.5%in%c(6),]$cluster.extract.v1.1 = "AL.3"
  
 
  #>> need more resolution to distinguish or cover by different LAI patterns
  a = FNN::knn(src.24ss.integrated.merge@reductions$umap_integrated@cell.embeddings[
    rownames(src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$Time%in%"ss24"&
                                                   !src.24ss.integrated.merge$cluster_snn_res.5%in%c(8,19,26,28),]),],
    src.24ss.integrated.merge@reductions$umap_integrated@cell.embeddings[
      rownames(src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$Time%in%"ss24"&
                                                     src.24ss.integrated.merge$cluster_snn_res.5%in%c(8,19,26,28),]),],
    src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$Time%in%"ss24"&
                                          !src.24ss.integrated.merge$cluster_snn_res.5%in%c(8,19,26,28),]$cluster.extract.v1.1, k = 10)
  src.24ss.integrated.merge@meta.data[src.24ss.integrated.merge$cluster_snn_res.5%in%c(8,19,26,28),]$cluster.extract.v1.1 = as.character(a)
}

for(i.define.time in c("27ss")){
  src.27ss.integrated.merge@reductions$umap_integrated@cell.embeddings = cbind(
    src.27ss.integrated.merge@reductions$umap_integrated@cell.embeddings[,1],
    -src.27ss.integrated.merge@reductions$umap_integrated@cell.embeddings[,2])
  colnames(src.27ss.integrated.merge@reductions$umap_integrated@cell.embeddings) = c("UMAP_1",'UMAP_2')
  
  DimPlot(src.27ss.integrated.merge, 
          group.by = "lineage", reduction = "umap_integrated",
          cols = color.lineage, na.value = "#eeeeee")
  
  src.27ss.integrated = FindNeighbors(src.27ss.integrated, reduction = "pca", dims = 1:30)
  src.27ss.integrated = FindClusters(src.27ss.integrated, resolution = 5)
  src.27ss.integrated.merge$cluster_snn_res.5 = NA
  src.27ss.integrated.merge@meta.data[
    rownames(src.27ss.integrated.merge@meta.data[src.27ss.integrated.merge$Time%in%"ss27",]),]$cluster_snn_res.5 = 
    src.27ss.integrated@meta.data[
      rownames(src.27ss.integrated.merge@meta.data[src.27ss.integrated.merge$Time%in%"ss27",]),]$RNA_snn_res.5 
  
  src.27ss.integrated.merge$cluster.extract.v1.1 = NA
  src.27ss.integrated.merge@meta.data[src.27ss.integrated.merge$cluster_snn_res.5%in%c(5),]$cluster.extract.v1.1 = "FG.5"
  src.27ss.integrated.merge@meta.data[src.27ss.integrated.merge$cluster_snn_res.5%in%c(15),]$cluster.extract.v1.1 = "FG.3"
  src.27ss.integrated.merge@meta.data[src.27ss.integrated.merge$cluster_snn_res.5%in%c(10,12,14,30,34),]$cluster.extract.v1.1 = "FG.2"
  src.27ss.integrated.merge@meta.data[src.27ss.integrated.merge$cluster_snn_res.5%in%c(4,16,17,26),]$cluster.extract.v1.1 = "FG.3/4"
  src.27ss.integrated.merge@meta.data[src.27ss.integrated.merge$cluster_snn_res.5%in%c(9),]$cluster.extract.v1.1 = "FG.6"
  src.27ss.integrated.merge@meta.data[src.27ss.integrated.merge$cluster_snn_res.5%in%c(29),]$cluster.extract.v1.1 = "FG.1"
  src.27ss.integrated.merge@meta.data[src.27ss.integrated.merge$cluster_snn_res.5%in%c(2,20,27,28),]$cluster.extract.v1.1 = "FG.4-AL.1/2/3"
  src.27ss.integrated.merge@meta.data[src.27ss.integrated.merge$cluster_snn_res.5%in%c(18,33),]$cluster.extract.v1.1 = "MG.3"
  src.27ss.integrated.merge@meta.data[src.27ss.integrated.merge$cluster_snn_res.5%in%c(24),]$cluster.extract.v1.1 = "AL.3"
  src.27ss.integrated.merge@meta.data[src.27ss.integrated.merge$cluster_snn_res.5%in%c(0,3,31),]$cluster.extract.v1.1 = "MG.2/3-HG.1"
  src.27ss.integrated.merge@meta.data[src.27ss.integrated.merge$cluster_snn_res.5%in%c(6,7,11),]$cluster.extract.v1.1 = "FG.4-MG.1/3"
  src.27ss.integrated.merge@meta.data[src.27ss.integrated.merge$cluster_snn_res.5%in%c(1,21,22,23),]$cluster.extract.v1.1 = "AL.3-MG.1/2/3"
  src.27ss.integrated.merge@meta.data[src.27ss.integrated.merge$cluster_snn_res.5%in%c(19,25),]$cluster.extract.v1.1 = "HG.1/2"
  src.27ss.integrated.merge@meta.data[src.27ss.integrated.merge$cluster_snn_res.5%in%c(8,32),]$cluster.extract.v1.1 = "HG.1"
  src.27ss.integrated.merge@meta.data[src.27ss.integrated.merge$cluster_snn_res.5%in%c(13),]$cluster.extract.v1.1 = "HG.2"
}
#===============================================================================


#===============================================================================
#>>>>> 4.4 extract endoderm trajectories
#===============================================================================
src.endoderm$cluster.extract.v1.1 = NA
src.endoderm@meta.data[paste("ss9_",colnames(src.9ss.integrated),sep=""),]$cluster.extract.v1.1 = src.9ss.integrated$cluster.extract.v1.1
src.endoderm@meta.data[paste("ss12_",colnames(src.12ss.integrated),sep=""),]$cluster.extract.v1.1 = src.12ss.integrated$cluster.extract.v1.1
src.endoderm@meta.data[paste("ss15_",colnames(src.15ss.integrated),sep=""),]$cluster.extract.v1.1 = src.15ss.integrated$cluster.extract.v1.1
src.endoderm@meta.data[paste("ss18_",colnames(src.18ss.integrated),sep=""),]$cluster.extract.v1.1 = src.18ss.integrated$cluster.extract.v1.1
src.endoderm@meta.data[paste("ss21_",colnames(src.21ss.integrated),sep=""),]$cluster.extract.v1.1 = src.21ss.integrated$cluster.extract.v1.1
src.endoderm@meta.data[paste("ss24_",colnames(src.24ss.integrated),sep=""),]$cluster.extract.v1.1 = src.24ss.integrated$cluster.extract.v1.1
src.endoderm@meta.data[paste("ss27_",colnames(src.27ss.integrated),sep=""),]$cluster.extract.v1.1 = src.27ss.integrated$cluster.extract.v1.1

celllist_ext_v1.1_FG.1 = colnames(src.endoderm[,src.endoderm$cluster.extract.v1.1%in%c("FG.1")])
celllist_ext_v1.1_FG.2 = colnames(src.endoderm[,src.endoderm$cluster.extract.v1.1%in%c("FG.2")])
celllist_ext_v1.1_FG.3 = colnames(src.endoderm[,src.endoderm$cluster.extract.v1.1%in%c("FG.3","FG.3/4")])
celllist_ext_v1.1_FG.4 = colnames(src.endoderm[,src.endoderm$cluster.extract.v1.1%in%c("FG.4","FG.3/4","FG.4-AL.1/2/3","FG.4-MG.1/3")])
celllist_ext_v1.1_FG.5 = colnames(src.endoderm[,src.endoderm$cluster.extract.v1.1%in%c("FG.5")])
celllist_ext_v1.1_FG.6 = colnames(src.endoderm[,src.endoderm$cluster.extract.v1.1%in%c("FG.6")])

celllist_ext_v1.1_AL.1 = colnames(src.endoderm[,src.endoderm$cluster.extract.v1.1%in%c("AL.1","FG.4-AL.1/2/3")])
celllist_ext_v1.1_AL.2 = colnames(src.endoderm[,src.endoderm$cluster.extract.v1.1%in%c("AL.2","FG.4-AL.1/2/3")])
celllist_ext_v1.1_AL.3 = colnames(src.endoderm[,src.endoderm$cluster.extract.v1.1%in%c("AL.3","AL.3-MG.2", "AL.3-MG.1/2/3","FG.4-AL.1/2/3")])

celllist_ext_v1.1_MG.1 = colnames(src.endoderm[,src.endoderm$cluster.extract.v1.1%in%c("MG.1","AL.3-MG.1/2/3","FG.4-MG.1/3","MG.1/2/3")])
celllist_ext_v1.1_MG.2 = colnames(src.endoderm[,src.endoderm$cluster.extract.v1.1%in%c("MG.2","AL.3-MG.2","MG.1/2/3","MG.2/3","AL.3-MG.1/2/3","MG.2/3-HG.1")])
celllist_ext_v1.1_MG.3 = colnames(src.endoderm[,src.endoderm$cluster.extract.v1.1%in%c("MG.3","MG.2/3","FG.4-MG.1/3","MG.1/2/3","AL.3-MG.1/2/3","MG.2/3-HG.1")])

celllist_ext_v1.1_HG.1 = colnames(src.endoderm[,src.endoderm$cluster.extract.v1.1%in%c("HG.1","HG.1/2","MG.2/3-HG.1")])
celllist_ext_v1.1_HG.2 = colnames(src.endoderm[,src.endoderm$cluster.extract.v1.1%in%c("HG.2","HG.1/2")])



src.endoderm.fg1.ext_v1.1 = src.endoderm[,colnames(src.endoderm)%in%celllist_ext_v1.1_FG.1]
src.endoderm.fg2.ext_v1.1 = src.endoderm[,colnames(src.endoderm)%in%celllist_ext_v1.1_FG.2]
src.endoderm.fg3.ext_v1.1 = src.endoderm[,colnames(src.endoderm)%in%celllist_ext_v1.1_FG.3]
src.endoderm.fg4.ext_v1.1 = src.endoderm[,colnames(src.endoderm)%in%celllist_ext_v1.1_FG.4]
src.endoderm.fg5.ext_v1.1 = src.endoderm[,colnames(src.endoderm)%in%celllist_ext_v1.1_FG.5]
src.endoderm.fg6.ext_v1.1 = src.endoderm[,colnames(src.endoderm)%in%celllist_ext_v1.1_FG.6]

src.endoderm.al1.ext_v1.1 = src.endoderm[,colnames(src.endoderm)%in%celllist_ext_v1.1_AL.1]
src.endoderm.al2.ext_v1.1 = src.endoderm[,colnames(src.endoderm)%in%celllist_ext_v1.1_AL.2]
src.endoderm.al3.ext_v1.1 = src.endoderm[,colnames(src.endoderm)%in%celllist_ext_v1.1_AL.3]

src.endoderm.mg1.ext_v1.1 = src.endoderm[,colnames(src.endoderm)%in%celllist_ext_v1.1_MG.1]
src.endoderm.mg2.ext_v1.1 = src.endoderm[,colnames(src.endoderm)%in%celllist_ext_v1.1_MG.2]
src.endoderm.mg3.ext_v1.1 = src.endoderm[,colnames(src.endoderm)%in%celllist_ext_v1.1_MG.3]

src.endoderm.hg1.ext_v1.1 = src.endoderm[,colnames(src.endoderm)%in%celllist_ext_v1.1_HG.1]
src.endoderm.hg2.ext_v1.1 = src.endoderm[,colnames(src.endoderm)%in%celllist_ext_v1.1_HG.2]

#===============================================================================













