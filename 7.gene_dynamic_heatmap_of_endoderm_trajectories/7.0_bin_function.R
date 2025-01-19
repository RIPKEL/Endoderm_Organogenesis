# correction
#-----------------------
seurat_mnn = function(seurat_name, seurat.selectgene, 
                      dims = 1:30, n.neighbors = 100, n.components=2){
  seurat = get(seurat_name)
  fMNN.res =  fastMNN(
    as.matrix(seurat@assays$RNA@data[seurat.selectgene, seurat$Phase%in%"G1"&seurat$batch%in%1]),
    as.matrix(seurat@assays$RNA@data[seurat.selectgene, seurat$Phase%in%"G1"&seurat$batch%in%2]),
    as.matrix(seurat@assays$RNA@data[seurat.selectgene, (!seurat$Phase%in%"G1")&seurat$batch%in%1]),
    as.matrix(seurat@assays$RNA@data[seurat.selectgene, (!seurat$Phase%in%"G1")&seurat$batch%in%2]),
    d=50,
    subset.row = seurat.selectgene)
  seurat@reductions$mnn = seurat@reductions$umap
  seurat@reductions$mnn@cell.embeddings = fMNN.res@assays@data$reconstructed@seed@components
  seurat@reductions$mnn@feature.loadings = fMNN.res@assays@data$reconstructed@seed@rotation
  colnames(seurat@reductions$mnn@cell.embeddings) = paste("MNN",1:50,sep="_")
  colnames(seurat@reductions$mnn@feature.loadings) = paste("MNN",1:50,sep="_")
  seurat@reductions$mnn@key = "MNN_"
  
  seurat_re = RunUMAP(seurat, reduction = "mnn", 
                      dims = dims, n.neighbors = n.neighbors, n.components = n.components)
  seurat@reductions$umap_mnn = seurat_re@reductions$umap
  assign(seurat_name, seurat)
}

seurat_int = function(seurat, seurat.selectgene, 
                      dims = 1:30, n.neighbors = 100, n.components=2){
  # Considering Batch Cellcycle
  anchor.seurat = FindIntegrationAnchors(c(seurat[,(seurat$Phase%in%"G1")&seurat$batch%in%c(2)],
                                           seurat[,(!seurat$Phase%in%"G1")&seurat$batch%in%c(2)],
                                           seurat[,(seurat$Phase%in%"G1")&seurat$batch%in%c(1)],
                                           seurat[,(!seurat$Phase%in%"G1")&seurat$batch%in%c(1)]),
                                         reduction = "rpca",  
                                         normalization.method = "LogNormalize",
                                         anchor.features = seurat.selectgene,
                                         dims = 1:30)
  
  seurat.re = IntegrateData(anchorset = anchor.seurat , dims = 1:30)
  
  DefaultAssay(seurat.re) = "integrated"
  seurat.re = ScaleData(seurat.re, features = rownames(seurat.re))
  seurat.re = RunPCA(seurat.re, features = seurat.selectgene)
  seurat.re = RunUMAP(seurat.re, dims = dims, 
                      n.neighbors = n.neighbors, n.components = n.components)
  print("Integrated-PCA+UMAP Done!")
  
  seurat.re = seurat_fdl(seurat.re, red = "pca", assay="integrated")
  print("FDL-Integrated Done!")
  
  seurat@assays$integrated = seurat.re@assays$integrated
  seurat@reductions$umap_integrated = seurat.re@reductions$umap
  seurat@reductions$umap_integrated@cell.embeddings = 
    seurat@reductions$umap_integrated@cell.embeddings[colnames(seurat),]
  
  seurat@reductions$pca_integrated = seurat.re@reductions$pca
  seurat@reductions$pca_integrated@cell.embeddings = 
    seurat@reductions$pca_integrated@cell.embeddings[colnames(seurat),]
  
  seurat@reductions$fdl_integrated = seurat.re@reductions$fdl
  seurat@reductions$fdl_integrated@cell.embeddings = 
    seurat@reductions$fdl_integrated@cell.embeddings[colnames(seurat),]
  DefaultAssay(seurat) = "RNA"
  
  rm(seurat.re)
  
  return(seurat)
}

seurat_fdl = function(seurat, red, assay){
  
  seurat = FindNeighbors(seurat, reduction = red, dims = 1:30, assay = assay)
  fdl.seurat = igraph::layout_with_fr(graph.adjacency(seurat@graphs[[paste(assay,"_nn",sep="")]]), dim=3)
  rownames(fdl.seurat) = colnames(seurat)
  colnames(fdl.seurat) = c("FDL_1","FDL_2","FDL_3")
  
  seurat[[paste("fdl_",red,sep="")]] = seurat@reductions$pca
  seurat[[paste("fdl_",red,sep="")]]@cell.embeddings = fdl.seurat
  seurat[[paste("fdl_",red,sep="")]]@key = "FDL_" 
  
  return(seurat)
}


batch_process_integration_Re = 
  function(cluster_names,
           seurat=NA, set_src=F, set_src.refine=F,
           set_red_refer = F,
           red_refer = "umap_integrated"){
  
  cluster_names = cluster_names
  cluster =  gsub("\\.","",tolower(cluster_names))
  cluster_list = t( hash::values(endoderm_transform, keys=cluster_names))
  print(cluster); print(cluster_list)
  
  if(set_red_refer){
    red_refer = red_refer
  }else{
    red_refer = hash::values(endoderm_embedding, keys=cluster_names)
  }

  if(set_src.refine==T){
    seurat = seurat
    # Tech
    seurat$Source_tech = NA
    seurat@meta.data[seurat$Time%in%paste("ss",3*c(3:9),sep=""),]$Source_tech = "refer"
    seurat@meta.data[seurat$Time%in%paste(3*c(3:9),"ss",sep=""),]$Source_tech = "query"
  }else{
    if(set_src==F){
      src.10x.refer = get(paste("src.",cluster,".integrated", sep=""))
      src.10x.refer@meta.data[,"Source_tech"] = "refer"
      
      src.sm3.query = src.sm3.merge[,src.sm3.merge$cluster.extract.v1.1_define%in%cluster_list]
      src.sm3.query@meta.data[,"Source_tech"] = "query" 
      
      red_refer = ifelse("umap_integrated" %in% names(src.10x.refer@reductions),
                         "umap_integrated", "umap_mnn")
      seurat = merge(src.10x.refer, src.sm3.query)
      
    }else{
      src.10x.refer = get(paste("src.",cluster,".integrated", sep=""))
      
      seurat = seurat
      # Tech
      seurat$Source_tech = NA
      seurat@meta.data[seurat$Time%in%paste("ss",3*c(3:9),sep=""),]$Source_tech = "refer"
      seurat@meta.data[seurat$Time%in%paste(3*c(3:9),"ss",sep=""),]$Source_tech = "query"
    }
  }
  
  # Phase
  seurat = CellCycleScoring(seurat, s.features = s_genes, g2m.features = g2m_genes)
  # Batch
  seurat$batch_tech = paste(seurat$batch, seurat$Phase, 
                            seurat$Source_tech, sep="_")
  print("Set batch_tech done!")
  
  seurat = ScaleData(seurat, 
                     features = rownames(seurat),
                     split.by = "batch_tech")
  seurat.selectgene = get(paste("src.",cluster,".integrated.selectgene", sep=""))
  
  seurat = RunPCA(seurat, features = seurat.selectgene)
  seurat = RunUMAP(seurat, dims = 1:20, n.neighbors = 100)
  print("Based processing done!")
  
  cell_sort_raw = colnames(seurat)
  cell_refer = colnames(seurat[,seurat$Source_tech%in%"refer"])
  cell_query = colnames(seurat[,seurat$Source_tech%in%"query"])
  
  # MNN on TPM
  #---------------
  if(cluster_names == "FG.2"){
    
    MNN.res = 
      mnnCorrect(
        as.matrix(seurat@assays$RNA@data[seurat.selectgene, 
                                         seurat$Source_tech%in%"refer"]),
        as.matrix(seurat@assays$RNA@data[seurat.selectgene,
                                         seurat$Source_tech%in%"query"]),
        k = 20, cos.norm.out=F)
    
  }else{
    
    MNN.res = 
      mnnCorrect(
        as.matrix(seurat@assays$RNA@data[seurat.selectgene, 
                                         seurat$Source_tech%in%"refer"&
                                           seurat$batch%in%1]),
        as.matrix(seurat@assays$RNA@data[seurat.selectgene, 
                                         seurat$Source_tech%in%"refer"&
                                           seurat$batch%in%2]),
        as.matrix(seurat@assays$RNA@data[seurat.selectgene,
                                         seurat$Source_tech%in%"query"]),
        k = 20, cos.norm.out=F)
  }
  
  
  seurat@assays$mnnRNA = CreateAssayObject(data = MNN.res@assays@data$corrected[,cell_sort_raw])
  seurat@assays$mnnRNA@key = "mnn_"
  seurat = ScaleData(seurat,  rownames(seurat@assays$mnnRNA), assay = "mnnRNA")
  
  seurat_re = RunPCA(seurat, features = rownames(seurat@assays$mnnRNA@data), assay = "mnnRNA")
  seurat[["mnn_pca"]] = seurat_re[["pca"]]; rm(seurat_re)
  seurat@reductions$mnn_pca@key = "Coord_"
  colnames(seurat[["mnn_pca"]]@cell.embeddings) = paste("Coord_",c(1:50),sep="")
  
  seurat_re = RunUMAP(seurat, dims = 1:20, reduction = 'mnn_pca', assay = "mnnRNA")
  seurat[["mnn_umap"]] = seurat_re[["umap"]]; rm(seurat_re)
  seurat@reductions$mnn_umap@key = "Coord_"
  colnames(seurat[["mnn_umap"]]@cell.embeddings) = paste("Coord_",c(1:2),sep="")
  #---------------
  print("MNN processing done!")
  
  # Fast-MNN: umap on mnn
  #---------------------
  if(cluster_names == "FG.2"){
    fMNN.res =  fastMNN(as.matrix(seurat@assays$RNA@data[seurat.selectgene, 
                                                         seurat$Source_tech%in%"refer"]),
                        as.matrix(seurat@assays$RNA@data[seurat.selectgene,
                                                         seurat$Source_tech%in%"query"]),
                        d=50,
                        subset.row = seurat.selectgene)
  }else{
    fMNN.res =  fastMNN(as.matrix(seurat@assays$RNA@data[seurat.selectgene, 
                                                         seurat$Source_tech%in%"refer"&
                                                           seurat$batch%in%1]),
                        as.matrix(seurat@assays$RNA@data[seurat.selectgene, 
                                                         seurat$Source_tech%in%"refer"&
                                                           seurat$batch%in%2]),
                        as.matrix(seurat@assays$RNA@data[seurat.selectgene,
                                                         seurat$Source_tech%in%"query"]),
                        d = 50,
                        subset.row = seurat.selectgene)
  }
  
  seurat@reductions$mnn = seurat@reductions$umap
  seurat@reductions$mnn@cell.embeddings = 
    fMNN.res@assays@data$reconstructed@seed@components[cell_sort_raw,]
  seurat@reductions$mnn@feature.loadings = 
    fMNN.res@assays@data$reconstructed@seed@rotation
  colnames(seurat@reductions$mnn@cell.embeddings) = paste("MNN",1:50,sep="_")
  colnames(seurat@reductions$mnn@feature.loadings) = paste("MNN",1:50,sep="_")
  seurat@reductions$mnn@key = "MNN_"
  
  seurat_re = RunUMAP(seurat, dims = 1:20, n.neighbors = 100, reduction = "mnn")
  seurat[["umap_mnn"]] = seurat_re[["umap"]]; rm(seurat_re)
  seurat@reductions$umap_mnn@key = "Coord_"
  colnames(seurat[["umap_mnn"]]@cell.embeddings) = paste("Coord_",c(1:2),sep="")
  #---------------------
  print("Fast-MNN processing done!")
  
  # FTA on Dif-Assay
  #---------------
  for(assay in c(#"RNA",
    "mnnRNA")){
    anchor.integrated = 
      FindTransferAnchors(reference = seurat[,cell_refer],
                          query = seurat[,cell_query],
                          reference.assay =  assay,
                          query.assay =  assay, scale = T,
                          features = seurat.selectgene)
    
    umap_embedding = src.10x.refer[[red_refer]]@cell.embeddings[cell_refer,c(1:2)]
    umap.transfer = TransferData(anchor.integrated, t(umap_embedding))
    assay_name = gsub('RNA',"", gsub("mnn","mnn_",assay))
    seurat[[paste(assay_name, "umap_fta", sep="")]] =  seurat[["umap"]]
    seurat[[paste(assay_name, "umap_fta", sep="")]]@cell.embeddings[cell_refer,] = umap_embedding
    seurat[[paste(assay_name, "umap_fta", sep="")]]@cell.embeddings[cell_query,] = 
      as.matrix(t(umap.transfer@data))
  }
  #---------------
  print("Find-Transfer-anchor processing done!")
  
  
  # FIA on Dif-Assay
  #---------------
  for(assay in c("RNA","mnnRNA")){
    anchor.integrated.merge = 
      FindIntegrationAnchors(c(seurat[,cell_refer],
                               seurat[,cell_query]),
                             assay = c(assay,assay), 
                             reduction = "cca", 
                             anchor.features = seurat.selectgene ,
                             dims = 1:30)
    
    seurat.re = IntegrateData(anchorset = anchor.integrated.merge, dims = 1:30)
    DefaultAssay(seurat.re) = "integrated"
    seurat.re = ScaleData(seurat.re, features = rownames(seurat.re))
    seurat.re = RunPCA(seurat.re, features = seurat.selectgene)
    seurat.re = RunUMAP(seurat.re, dims = 1:20, n.neighbors = 100)
    
    assay_name_re = gsub('RNA',"", gsub("mnn_","mnn",assay_name))
    
    seurat.re[[paste(assay_name_re,"pca_integrated",sep="")]] = seurat.re[["pca"]]
    seurat.re[[paste(assay_name_re,"umap_integrated",sep="")]] = seurat.re[["umap"]]
    
    seurat[[paste(assay_name_re,"pca_integrated",sep="")]] = seurat[["pca"]]
    seurat[[paste(assay_name_re,"pca_integrated",sep="")]]@cell.embeddings[colnames(seurat),] = 
      seurat.re[[paste(assay_name_re,"pca_integrated",sep="")]]@cell.embeddings[colnames(seurat),] 
    
    seurat[[paste(assay_name_re,"umap_integrated",sep="")]] = seurat[["umap"]]
    seurat[[paste(assay_name_re,"umap_integrated",sep="")]]@cell.embeddings[colnames(seurat),] = 
      seurat.re[[paste(assay_name_re,"umap_integrated",sep="")]]@cell.embeddings[colnames(seurat),] 
    rm(seurat.re)
  }
  #---------------
  print("Find-integrated-anchor processing done!")
  
  seurat = seurat[,seurat$nFeature_RNA>2500]
  
  return(seurat) 
}

batch_process_tracing =
  function(cluster_names,
           Set_gene = F, setgene = NA,
           batch_effect=NULL){
  
  cluster =  gsub("\\.","",tolower(cluster_names))
  print(cluster)
  
  seurat = get(paste("src.",cluster,".tracing", sep=""))
  seurat = FindVariableFeatures(seurat, nfeatures = 2000)
  
  # Set Cell cycle Score & Scale
  seurat = CellCycleScoring(seurat, 
                            s.features = s_genes, g2m.features = g2m_genes)
  seurat = ScaleData(seurat, split.by = "Phase",
                     features = rownames(seurat))
  
  # Set Batch Effect
  seurat$batch_effect = 
    paste(seurat$SeqDate,seurat$Time,seurat$lineage,sep="_")
  
  cell_filter = colnames(seurat)
  
  if( length(batch_effect) != 0){
    cell_filter = 
      rownames(seurat@meta.data[!seurat$batch_effect%in%batch_effect,])
  }
  
  # Set Gene
  seurat.filtergene = 
    Myfilter(as.matrix(seurat@assays$RNA@data[,cell_filter]),
             gene = seurat@assays$RNA@var.features,
             pearson.threshold = 0.15,partner.threshold = 5,
             bottom.dispersion.interval = 0.1)
  
  seurat.selectgene = 
    get(paste("src.",cluster,".integrated.selectgene", sep=""))
  
  seurat.uniongene = 
    unique(c(seurat.filtergene,seurat.selectgene))
  
  seurat.setgene = seurat.uniongene
  
  
  # PCA-UMAP
  seurat = RunPCA(seurat, features = seurat.setgene)
  print("PCA Done!")
  
  seurat = RunUMAP(seurat, dims = 1:30, reduction = "pca")
  print("UMAP Done!")
  
  # PCA-FDL
  seurat_fdl = function(seurat, red, assay){
    
    seurat = FindNeighbors(seurat, reduction = red, 
                           dims = 1:30, assay = assay)
    
    fdl.seurat = layout_with_fr(graph.adjacency(seurat@graphs[[paste(assay,"_nn",sep="")]]), dim=3)
    rownames(fdl.seurat) = colnames(seurat)
    colnames(fdl.seurat) = c("FDL_1","FDL_2","FDL_3")
    
    seurat[[paste("fdl_",red,sep="")]] = seurat@reductions$pca
    seurat[[paste("fdl_",red,sep="")]]@cell.embeddings = fdl.seurat
    seurat[[paste("fdl_",red,sep="")]]@key = "FDL_" 
    
    return(seurat)
  }
  
  seurat = seurat_fdl(seurat, red = "pca", assay="RNA")
  print("FDL-PCA Done!")
  return(seurat)
  
  
  #-------------------------------------------------------
  # Considering Batch
  seurat$batch_define = seurat$batch
  seurat@meta.data[!seurat$batch%in%c(1,2),]$batch_define = 1
  
  seurat_int = 
    function(seurat, seurat.selectgene, cluster_names,
             dims = 1:30, n.neighbors = 100, n.components=2){
      # Considering Batch Cellcycle
      
      cell_phase = table(seurat$Phase)
      
      if( 2* cell_phase["G1"] < (cell_phase["G2M"]+cell_phase['S'])){
        Phase = c("G1","G2M")
      }else{
        Phase = c("G1")
      }
      
      
      if(cluster_names%in%c('NG.1', "HG.2")){
        Phase = c("G1","G2M")
      }
      
      
      anchor.seurat = FindIntegrationAnchors(c(seurat[,(seurat$Phase%in%Phase)&
                                                        seurat$batch_define%in%c(1,2)],
                                               seurat[,(!seurat$Phase%in%Phase)&
                                                        seurat$batch_define%in%c(1,2)]),
                                             reduction = "cca",
                                             anchor.features = seurat.selectgene,
                                             dims = 1:30)
      
      seurat.re = IntegrateData(anchorset = anchor.seurat, dims = 1:30)
      
      DefaultAssay(seurat.re) = "integrated"
      seurat.re = ScaleData(seurat.re, features = rownames(seurat.re))
      seurat.re = RunPCA(seurat.re, features = seurat.selectgene)
      seurat.re = RunUMAP(seurat.re, dims = dims, 
                          n.neighbors = n.neighbors, n.components = n.components)
      print("Integrated-PCA+UMAP Done!")
      
      seurat.re = seurat_fdl(seurat.re, red = "pca", assay="integrated")
      print("FDL-Integrated Done!")
      
      seurat@assays$integrated = seurat.re@assays$integrated
      seurat@reductions$umap_integrated = seurat.re@reductions$umap
      seurat@reductions$umap_integrated@cell.embeddings = 
        seurat@reductions$umap_integrated@cell.embeddings[colnames(seurat),]
      
      seurat@reductions$pca_integrated = seurat.re@reductions$pca
      seurat@reductions$pca_integrated@cell.embeddings = 
        seurat@reductions$pca_integrated@cell.embeddings[colnames(seurat),]
      
      seurat@reductions$fdl_integrated = seurat.re@reductions$fdl
      seurat@reductions$fdl_integrated@cell.embeddings = 
        seurat@reductions$fdl_integrated@cell.embeddings[colnames(seurat),]
      DefaultAssay(seurat) = "RNA"
      
      rm(seurat.re)
      
      return(seurat)
    }
  seura = seurat_int(seurat, seurat.setgene)
  
  return(seurat)
  }


label_KNN_learn_Re =
  function(seurat, 
           reduction, label, 
           group, 
           ref_group="refer", 
           que_group="query"){
  
  cell_refer = colnames(seurat[,seurat[[group]]==ref_group])
  cell_query = colnames(seurat[,seurat[[group]]==que_group])
  
  for(red in reduction){
    embedding_refer = seurat[[red]]@cell.embeddings[cell_refer, ]
    embedding_query = seurat[[red]]@cell.embeddings[cell_query, ]
    label_refer = seurat@meta.data[cell_refer, label]
    
    label_learning = 
      FNN::knn(embedding_refer, embedding_query,
               label_refer, k = 10)
    
    seurat@meta.data[,paste(label,red,sep="_")] = NA
    seurat@meta.data[cell_refer, paste(label,red,sep="_")] =
      seurat@meta.data[cell_refer, label] 
    seurat@meta.data[cell_query, paste(label,red,sep="_")] =
      as.character(label_learning)
  }
  
  return(seurat)
  }


integration_fta = function(seurat, cell_refer, cell_query, 
                           seurat.selectgene, embeddings){
  
  for(assay in c("RNA", "mnnRNA")){
    anchor.integrated = 
      FindTransferAnchors(reference = seurat[,cell_refer],
                          query = seurat[,cell_query],
                          reference.assay =  assay,
                          query.assay =  assay, scale = T,
                          features = seurat.selectgene)
    
    umap_embedding = embeddings[cell_refer,c(1:2)]
    umap.transfer = TransferData(anchor.integrated, t(umap_embedding))
    assay_name = gsub('RNA',"", gsub("mnn","mnn_",assay))
    seurat[[paste(assay_name, "umap_fta", sep="")]] =  seurat[["umap"]]
    seurat[[paste(assay_name, "umap_fta", sep="")]]@cell.embeddings[cell_refer,] = umap_embedding
    seurat[[paste(assay_name, "umap_fta", sep="")]]@cell.embeddings[cell_query,] = 
      as.matrix(t(umap.transfer@data))
  }
  
  #---------------
  print("Find-Transfer-anchor processing done!")
  
  return(seurat)
}
  


