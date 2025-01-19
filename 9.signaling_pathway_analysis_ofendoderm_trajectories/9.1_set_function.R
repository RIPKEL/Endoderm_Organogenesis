#===============================================================================
#>>> 9.signaling pathways analysis 
#===============================================================================

#>>>> These function for running signal pathway analysis is essential.
#>>>> Additionally, matrices with Reactome/KEGG pathway information are required,
#> which can be downloaded from the official websites of the respective databases.

#===============================================================================
#>>> 9.1 Set Functions
#===============================================================================
kernelsmooth = function(x, kern=gaussian.kernel(), norm=TRUE) {
  # how many rows/cols of zeroes are used to pad.
  width <- length(kern)
  pad <- floor(width / 2)
  
  # record the width and height the input data matrix
  x_w <- length(x)
  
  # Are we normalizing the kernel?
  if (norm == TRUE) {
    k <- kern / sum(abs(kern))
  } else {
    k <- kern
  }
  
  # pad all around the matrix an equal width of zeros
  x_pad <- t(ptw::padzeros(data=x, nzeros=pad, side="both"))
  
  # Pre-allocate the final (smoothed) data matrix
  s <- vector(length = x_w)
  
  # Pre-allocate a temporary matrix for the iterative calculations
  temp <- vector(length = width)
  
  # Loop through the data to apply the kernel.
  for (col in 1:x_w ) {
    temp <- x_pad[col:(col + width - 1)]
    s[col] <-  sum(k * temp)
  }
  
  # return the smoothed data
  return(s)
}

gaussian.kernel = function(x=50,sigma=4){
  A=vector(length = x)
  for (i in 1:x) {
    A[i]=exp(-((i-x/2-0.5)^2)/sigma^2)/(2*pi*sigma)
  }
  return(A)
}

Mywilcox = function(tpm, variable){
  p.val=c()
  for (gene in rownames(tpm)) {
    p.val=c(p.val,wilcox.test(as.numeric(tpm[gene,variable==names(table(variable))[1]]),
                              as.numeric(tpm[gene,variable==names(table(variable))[2]]))$p.value)
  }
  names(p.val)=rownames(tpm)
  p.val=sort(p.val)
  return(p.val)
}

MyKEGGpathway2gene = function(pathway,species="mmu"){
  url=paste0("https://rest.kegg.jp/get/",pathway)
  form=strsplit(RCurl::getURL(url),split = "\n")[[1]]
  gene_start=grep("^GENE",form)
  
  if (length(gene_start)!=0) {
    
    gene_list=c()
    gene_list=c(gene_list,gsub(" ","",stringr::str_extract(form[gene_start],"\\d+ ")))
    i=1
    while(substr(form[gene_start+i],1,1)==" "){
      gene_list=c(gene_list,gsub(" ","",stringr::str_extract(form[gene_start+i],"\\d+ ")))
      i=i+1
    }
    
    if(species=="mmu"){
      library(org.Mm.eg.db) 
      symbollist=AnnotationDbi::select(org.Mm.eg.db, keys = gene_list, 
                                       columns = c("SYMBOL","ENSEMBL"), keytype = "ENTREZID")
    }
    if(species=="hsa"){
      library(org.Hs.eg.db)
      symbollist=AnnotationDbi::select(org.Mm.eg.db, keys = gene_list, 
                                       columns = c("SYMBOL","ENSEMBL"), keytype = "ENTREZID")
    }
    return(symbollist)
  }else{return(NULL)}
}

seurat_pathway_pre = function(pathway.matrix, seurat.ref){
  seurat.ref = CellCycleScoring(seurat.ref, s.features = s_genes, g2m.features = g2m_genes, assay = "RNA")
  
  seurat.pathway = CreateSeuratObject(pathway.matrix, assay = "Pathway")
  seurat.pathway@meta.data = 
    cbind(seurat.pathway@meta.data,
          seurat.ref@meta.data[colnames(seurat.pathway),])
  
  seurat.pathway = ScaleData(seurat.pathway, rownames(seurat.pathway))
  seurat.pathway = FindVariableFeatures(seurat.pathway, nfeatures = 150)
  VariableFeaturePlot(seurat.pathway)
  
  #-- PCA & UMAP
  seurat.pathway = RunPCA(seurat.pathway)
  seurat.pathway = RunUMAP(seurat.pathway,dims = 1:5)
  
  #-- FDL
  seurat.pathway = FindNeighbors(seurat.pathway,reduction = "pca",dims = 1:30)
  fdl.seurat.pathway = igraph::layout_with_fr(igraph::graph.adjacency(seurat.pathway@graphs[["Pathway_nn"]]),dim=3)
  rownames(fdl.seurat.pathway) = colnames(seurat.pathway)
  colnames(fdl.seurat.pathway) = c("FDL_1","FDL_2","FDL_3")
  seurat.pathway@reductions$fdl = seurat.pathway@reductions$pca
  seurat.pathway@reductions$fdl@cell.embeddings = fdl.seurat.pathway
  seurat.pathway@reductions$fdl@key = "FDL_"
  return(seurat.pathway)
}

seurat_plsda_pathway = function(seurat.pathway, type, assay = "Pathway"){
  # library(mixOmics)
  data1 = seurat.pathway@assays[[assay]]@data
  # group = seurat.pathway@meta.data
  X = t(data1)
  Y = type
  plsda.datatm = plsda(X, Y, ncomp = 3)
  plotIndiv(plsda.datatm, ind.names = FALSE, 
            legend=TRUE, ellipse = TRUE, title="sPLS-DA - final result")
  background = background.predict(plsda.datatm, comp.predicted=2 ,dist = "max.dist") 
  #plotVar(plsda.datatm) 
  plotIndiv(plsda.datatm, comp = 1:2, 
            ind.names = FALSE, title = "Maximum distance",
            legend = TRUE,  background = background,ellipse = TRUE)
  auc.plsda = auroc(plsda.datatm)
  seurat.pathway@reductions$plsda =seurat.pathway@reductions$pca
  seurat.pathway@reductions$plsda@cell.embeddings = plsda.datatm$variates$X
  seurat.pathway@reductions$plsda@feature.loadings = plsda.datatm$loadings$X
  colnames(seurat.pathway@reductions$plsda@cell.embeddings) = c("Comp_1","Comp_2","Comp_3")
  seurat.pathway@reductions$plsda@key = "Comp_"
  result = list()
  result[["seurat"]] = seurat.pathway
  result[["auc"]] = auc.plsda
  return(result)
}

plsda_select =  function(seurat, comp = 2, threshold = 25){
  plsda_X = seurat@reductions$plsda@feature.loadings
  
  ncol_X = comp
  plsda1_select_fin = c()
  
  for(i in c(1:ncol_X)){
    plsda1_select = c(
      plsda_X[,i][order(plsda_X[,i])][1:threshold],
      plsda_X[,i][rev(order(plsda_X[,i]))][1:threshold])
    plsda1_select = names(
      plsda1_select[order(plsda1_select)])
    plsda1_select_fin = 
      c(plsda1_select_fin, plsda1_select)
  }
  
  plsda1_select_fin = unique(plsda1_select_fin)
  
  return(plsda1_select_fin)
}

anova.test = function(tpm, variable){
  p.val=c()
  for(i.tpm.gene in rownames(tpm)){
    data.anova = as.data.frame(cbind(tpm[i.tpm.gene,], variable))
    data.anova = data.anova[(is.na(data.anova[,1])==F)&(is.na(data.anova[,2])==F),] 
    colnames(data.anova) = c("feature", "ident")
    summary.anova = summary(aov(feature~ident, data.anova))
    char.anova = str_split(as.character(summary.anova), "`")[[1]][9]
    for(i.char in c("*\\("," ","=","c","NA",",",")")){
      char.anova = gsub(i.char, "", char.anova)
    }
    p.val=c(p.val, char.anova)
  }
  p.val = as.numeric(p.val)
  names(p.val)=rownames(tpm)
  
  names(p.val)=rownames(tpm)
  return(p.val)
}
#===============================================================================





