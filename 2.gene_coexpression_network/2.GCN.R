#===============================================================================
#>> 2.gene co-expression network
#===============================================================================


#===============================================================================
#>> 2.1 9SS cell type
#===============================================================================

for(i.cell_type in c("sum")){
  src.9ss.integrated.sum = src.27ss.integrated

  src.9ss.integrated.sum = FindVariableFeatures(src.9ss.integrated.sum,nfeatures = 2000)
  src.9ss.integrated.sum.gene =   
    Myfilter(as.matrix(src.9ss.integrated.sum@assays$RNA@data),
             gene = unique(c(src.9ss.integrated.sum @assays$RNA@var.features,
                             rownames(src.9ss.integrated.sum@reductions$pca@feature.loadings))),
             pearson.threshold = 0.05,  partner.threshold = 5,
             bottom.dispersion.interval = 0.1)
  gene_list_p1 = src.9ss.integrated.sum.gene # rownames(src.9ss.integrated@reductions$pca@feature.loadings)
  names(gene_list_p1) = rep(1, length(gene_list_p1))
  
  pdf("FG1_correct/try_sum.pdf",7,7)
  src.9ss.integrated.sum.rowtree =  
    MyHeatmap(as.matrix(src.9ss.integrated.sum@assays$RNA@data[gene_list_p1, ]),
              type = "row.relat",
              ColSideColorsSize = 3,
              RowSideColorsSize = 2,
              #Rowv = "none", 
              #Colv = "none",
              return.tree = "row",
              labCol="none", graph = T)
  dev.off()
  
  tree_src.9ss.integrated.sum.rowtree = as.dendrogram(src.9ss.integrated.sum.rowtree)
  #>> rm-cell-cycle
  gene_src.9ss.integrated.sum.rowtree = setdiff(
    gene_list_p1,
    c(labels(tree_src.9ss.integrated.sum.rowtree[[2]][[1]]),
      labels(tree_src.9ss.integrated.sum.rowtree[[2]][[2]][[2]][[1]])))
  
  pdf("figure.v08.07/FG1_correct/try_sum.pdf",7,7)
  src.9ss.integrated.sum.rowtree.1 =  
    MyHeatmap(as.matrix(src.9ss.integrated.sum@assays$RNA@data[gene_src.9ss.integrated.sum.rowtree, ]),
              type = "row.relat",
              ColSideColorsSize = 3,
              RowSideColorsSize = 2,
              #Rowv = "none", 
              #Colv = "none",
              return.tree = "row",
              labCol="none", graph = T)
  dev.off()
  tree_src.9ss.integrated.sum.rowtree.1 = as.dendrogram(src.9ss.integrated.sum.rowtree.1)
  #>> rm-cell-cycle
  gene_src.9ss.integrated.sum.rowtree.1 = setdiff(
    c(labels(tree_src.9ss.integrated.sum.rowtree.1[[1]]),
      labels(tree_src.9ss.integrated.sum.rowtree.1[[2]])),
    c(labels(tree_src.9ss.integrated.sum.rowtree.1[[1]][[1]]),
      labels(tree_src.9ss.integrated.sum.rowtree.1[[1]][[2]][[1]][[1]][[1]])))

  pdf("figure.v08.07/FG1_correct/try_sum.pdf",7,7)
  src.9ss.integrated.sum.rowtree.2 =  
    MyHeatmap(as.matrix(src.9ss.integrated.sum@assays$RNA@data[gene_src.9ss.integrated.sum.rowtree.1, ]),
              type = "row.relat",
              ColSideColorsSize = 3,
              RowSideColorsSize = 2,
              #Rowv = "none", 
              #Colv = "none",
              return.tree = "row",
              labCol="none", graph = T)
  dev.off()
  
  tree_src.9ss.integrated.sum.rowtree.2 = as.dendrogram(src.9ss.integrated.sum.rowtree.2)
  #>> save-unique-pattern
  gene_src.9ss.integrated.sum.rowtree.2 = c(
    labels(tree_src.9ss.integrated.sum.rowtree.2[[1]][[2]][[1]]), 
    labels(tree_src.9ss.integrated.sum.rowtree.2[[2]][[2]][[1]][[1]][[2]]), 
    labels(tree_src.9ss.integrated.sum.rowtree.2[[2]][[1]][[2]][[2]]), 
    labels(tree_src.9ss.integrated.sum.rowtree.2[[2]][[2]][[2]]))
  
  logistic = function(x, threshold = 0.1){
    x.log = 1 / (1 + exp(-x))
    x[x < threshold] <- 0
    return(x)} 
  
  seurat = src.9ss.integrated.sum
  gene_list_sum = gene_src.9ss.integrated.sum.rowtree.2; p = 0.2
  #----------------
  pearson.pathway = WGCNA::cor(t(seurat@assays$RNA@data[gene_list_sum, ]), method = "p")
  pearson.pathway = logistic(pearson.pathway, threshold = 0.1)
  pearson.pathway[is.na(pearson.pathway)]=0
  pearson.pathway[pearson.pathway<p]=0
  # pearson.pathway[pearson.pathway<0.58]=0  # Used in Background
  
  graph.pathway = igraph::graph.adjacency(pearson.pathway,mode = "undirected",weighted = T)
  graph.pathway = igraph::simplify(graph.pathway,remove.multiple = T,remove.loops = T)
  graph.pathway = igraph::delete.vertices(
    graph.pathway,
    v = igraph::V(graph.pathway)[igraph::clusters(graph.pathway)$membership%in%which(
      igraph::clusters(graph.pathway)$csize<=5)])
  degree.Foregut.pathway=igraph::degree(graph.pathway)
  layout.pathway = igraph::layout_nicely(graph.pathway)
  
  graph.cluster <- igraph::cluster_walktrap(graph.pathway,steps = 5)
  gene.cluster <- graph.cluster$membership
  names(gene.cluster) <- graph.cluster$names
  set.seed(3)
  
  pdf("gcn_9ss_sum.pdf",7,7)
  plot.igraph(graph.pathway,
              layout = layout.pathway,
              vertex.size=6,
              label.cex =6,
              vertex.label = NA, 
              vertex.label.color = "black",
              vertex.color = colors.num[gene.cluster])
  dev.off()
  #----------------
  
  gene_list_sum.re = c(
    gene.cluster[gene.cluster %in% c(1)],
    gene.cluster[gene.cluster %in% c(2)],
    gene.cluster[gene.cluster %in% c(3)],
    gene.cluster[gene.cluster %in% c(4)])
  gene_list_sum.re[gene_list_sum.re%in%c(1)]="FG"
  gene_list_sum.re[gene_list_sum.re%in%c(2)]="AL"
  gene_list_sum.re[gene_list_sum.re%in%c(3)]="MG"
  gene_list_sum.re[gene_list_sum.re%in%c(4)]="HG"
  
  gene_list_sum = names(gene_list_sum)
  names(gene_list_sum) = gene_list_sum
  
  pdf("heatmap_sum_9ss.pdf",7,7)
  #-------------
  data = src.9ss.integrated.sum@assays$RNA@data
  a = MyHeatmap(
    as.matrix(data[gene_list_sum,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    c.cov.method = "p",
    c.hc.method = "ward.D",
    RowSideColors = t(cbind(
      MyName2Col(
        names(gene_list_sum),
        cluster.endoderm.color.v5
      ))),
    # Rowv = "none", Colv = "none", 
    return.tree = "col",
    ColSideColorsSize = 2, RowSideColorsSize = 2,
    #labCol = name_list_all,  labRow = name_list_all,
    margins = c(8,8))
  #-------------
  dev.off()
  
  sum.coltree_1 = as.dendrogram(a) # cell
  cell_list_sum = c(
    as.dendrogram(labels(sum.coltree_1[[2]])),
    as.dendrogram(labels(sum.coltree_1[[1]][[2]][[2]])),
    as.dendrogram(labels(sum.coltree_1[[1]][[2]][[1]])),
    as.dendrogram(labels(sum.coltree_1[[1]][[1]])))
  names(cell_list_sum) = c(
    rep("FG", length(c(as.dendrogram(labels(sum.coltree_1[[2]]))))),
    rep("AL", length(c(as.dendrogram(labels(sum.coltree_1[[1]][[2]][[2]]))))),
    rep("MG", length(c(as.dendrogram(labels(sum.coltree_1[[1]][[2]][[1]]))))),
    rep("HG", length(c(as.dendrogram(labels(sum.coltree_1[[1]][[1]]))))))
  
  src.9ss.integrated_sum$cluster.v06.26.re = NA
  src.9ss.integrated_sum@meta.data[cell_list_sum,]$cluster.v06.26.re = names(cell_list_sum)
  src.9ss.integrated$cluster.v06.26.re = src.9ss.integrated_sum$cluster.v06.26.re
}
src.9ss.integrated$cluster.v06.26.re = src.9ss.integrated_sum$cluster.v06.26.re

#>>> For sub-clusters
for(i.cell_type in c("FG")){
  src.9ss.integrated.fg = src.9ss.integrated[,src.9ss.integrated$cluster.v06.26.re%in%c("FG")]
  src.9ss.integrated.fg = FindVariableFeatures(src.9ss.integrated.fg,nfeatures = 2000)
  src.9ss.integrated.fg.filtergene =   
    Myfilter(as.matrix(src.9ss.integrated.fg@assays$RNA@data),
             gene = unique(c(src.9ss.integrated.fg @assays$RNA@var.features,
                             rownames(src.9ss.integrated.fg@reductions$pca@feature.loadings))),
             pearson.threshold = 0.15,partner.threshold = 5,
             bottom.dispersion.interval = 0.1)
  
  
  seurat = src.9ss.integrated.fg
  gene_list_sum = gsrc.9ss.integrated.fg.filtergene; p = 0.55
  #----------------
  pearson.pathway = WGCNA::cor(t(seurat@assays$RNA@data[gene_list_sum,]), method = "p")
  pearson.pathway = logistic(pearson.pathway, threshold = 0.25)
  pearson.pathway[is.na(pearson.pathway)]=0
  pearson.pathway[pearson.pathway<p]=0
  # pearson.pathway[pearson.pathway<0.58]=0  # Used in Background
  
  graph.pathway = graph.adjacency(pearson.pathway,mode = "undirected",weighted = T)
  graph.pathway = igraph::simplify(graph.pathway,remove.multiple = T,remove.loops = T)
  graph.pathway = delete.vertices(
    graph.pathway,
    v = V(graph.pathway)[clusters(graph.pathway)$membership%in%which(clusters(graph.pathway)$csize<=10)])
  degree.Foregut.pathway=igraph::degree(graph.pathway)
  layout.pathway = layout_nicely(graph.pathway)
  
  graph.cluster <- cluster_walktrap(graph.pathway,steps = 5)
  gene.cluster <- graph.cluster$membership
  names(gene.cluster) <- graph.cluster$names
  set.seed(3)
  
  pdf("gcn_9ss_fg.pdf",7,7)
  plot.igraph(graph.pathway,
              layout = layout.pathway,
              vertex.size=6,
              label.cex =6,
              vertex.label = NA, 
              vertex.label.color = "black",
              vertex.color = colors.num[gene.cluster]
              # vertex.color = color.p3_gcn[gene_list_sum_p3_gcn])
  )
  dev.off()
  #----------------
  
  gene_list_fg.re = c(
    gene.cluster[gene.cluster %in% c(1)],
    gene.cluster[gene.cluster %in% c(2)],
    gene.cluster[gene.cluster %in% c(3)],
    gene.cluster[gene.cluster %in% c(4)],
    gene.cluster[gene.cluster %in% c(5)])
  gene_list_fg.re[gene_list_fg.re%in%c(1)]="FG.1"
  gene_list_fg.re[gene_list_fg.re%in%c(2)]="FG.2"
  gene_list_fg.re[gene_list_fg.re%in%c(3)]="FG.3"
  gene_list_fg.re[gene_list_fg.re%in%c(4)]="FG.4"
  gene_list_fg.re[gene_list_fg.re%in%c(5)]="FG.5"
  
  color.p3_gcn = c("#9970D7","#D190D8","#D1A955","#A8C2D5","#B6DABE")
  names(color.p3_gcn) = c("FG.1",'FG.2',"FG.3",'FG.4','FG.5')
  
  gene_list_fg = names(gene_list_fg)
  names(gene_list_fg) = gene_list_fg
  
  pdf("heatmap_fg_9ss.pdf",7,7)
  #-------------
  data = src.9ss.integrated.fg@assays$RNA@data
  a = MyHeatmap(
    as.matrix(data[gene_list_fg,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    c.cov.method = "p",
    c.hc.method = "ward.D",
    RowSideColors = t(cbind(
      MyName2Col(
        names(gene_list_fg),
        cluster.endoderm.color.v5
      ))),
    # Rowv = "none", Colv = "none", 
    return.tree = "col",
    ColSideColorsSize = 2, RowSideColorsSize = 2,
    #labCol = name_list_all,  labRow = name_list_all,
    margins = c(8,8))
  #-------------
  dev.off()
  
  fg.coltree_1 = as.dendrogram(a) # cell
  cell_list_fg = c(
    as.dendrogram(labels(fg.coltree_1[[2]][[1]][[1]][[1]])),
    as.dendrogram(labels(fg.coltree_1[[2]][[1]][[2]])),
    as.dendrogram(labels(fg.coltree_1[[2]][[2]])),
    as.dendrogram(labels(fg.coltree_1[[1]][[1]])),
    as.dendrogram(labels(fg.coltree_1[[1]][[2]])),
    as.dendrogram(labels(fg.coltree_1[[2]][[1]][[1]][[1]])))
  names(cell_list_fg) = c(
    rep("FG.1", length(c(as.dendrogram(labels(fg.coltree_1[[2]][[1]][[1]][[1]]))))),
    rep("FG.2", length(c(as.dendrogram(labels(fg.coltree_1[[2]][[1]][[2]]))))),
    rep("FG.3", length(c(as.dendrogram(labels(fg.coltree_1[[2]][[2]]))))),
    rep("FG.4", length(c(as.dendrogram(labels(fg.coltree_1[[1]][[1]]))))),
    rep("FG.5", length(c(as.dendrogram(labels(fg.coltree_1[[1]][[2]]))))),
    rep("FG.6", length(c(as.dendrogram(labels(fg.coltree_1[[2]][[1]][[1]][[1]]))))))
  
  src.9ss.integrated_fg$cluster.v06.26.re..correct = NA
  src.9ss.integrated_fg@meta.data[cell_list_fg,]$cluster.v06.26.re..correct = names(cell_list_fg)
}

for(i.cell_type in c("AL")){
  src.9ss.integrated.al = src.9ss.integrated[,src.9ss.integrated$cluster.v06.26.re%in%c("AL")]
  src.9ss.integrated.al = FindVariableFeatures(src.9ss.integrated.al,nfeatures = 2000)
  src.9ss.integrated.al.filtergene =   
    Myfilter(as.matrix(src.9ss.integrated.al@assays$RNA@data),
             gene = unique(c(src.9ss.integrated.al @assays$RNA@var.features,
                             rownames(src.9ss.integrated.al@reductions$pca@feature.loadings))),
             pearson.threshold = 0.15,partner.threshold = 5,
             bottom.dispersion.interval = 0.1)
  
  
  seurat = src.9ss.integrated.al
  gene_list_sum = gsrc.9ss.integrated.al.filtergene; p = 0.52
  #----------------
  pearson.pathway = WGCNA::cor(t(seurat@assays$RNA@data[gene_list_sum,]), method = "p")
  pearson.pathway = logistic(pearson.pathway, threshold = 0.25)
  pearson.pathway[is.na(pearson.pathway)]=0
  pearson.pathway[pearson.pathway<p]=0
  # pearson.pathway[pearson.pathway<0.58]=0  # Used in Background
  
  graph.pathway = graph.adjacency(pearson.pathway,mode = "undirected",weighted = T)
  graph.pathway = igraph::simplify(graph.pathway,remove.multiple = T,remove.loops = T)
  graph.pathway = delete.vertices(
    graph.pathway,
    v = V(graph.pathway)[clusters(graph.pathway)$membership%in%which(clusters(graph.pathway)$csize<=10)])
  degree.Foregut.pathway=igraph::degree(graph.pathway)
  layout.pathway = layout_nicely(graph.pathway)
  
  graph.cluster <- cluster_walktrap(graph.pathway,steps = 5)
  gene.cluster <- graph.cluster$membership
  names(gene.cluster) <- graph.cluster$names
  set.seed(3)
  
  pdf("gcn_9ss_al.pdf",7,7)
  plot.igraph(graph.pathway,
              layout = layout.pathway,
              vertex.size=6,
              label.cex =6,
              vertex.label = NA, 
              vertex.label.color = "black",
              vertex.color = colors.num[gene.cluster])
  dev.off()
  #----------------
  
  gene_list_al.re = c(
    gene.cluster[gene.cluster %in% c(1)],
    gene.cluster[gene.cluster %in% c(2,3)])
  gene_list_al.re[gene_list_al.re%in%c(1)]="AL.1" # overlap?
  gene_list_al.re[gene_list_al.re%in%c(2,3)]="AL.3"
  
  pdf("figure.v08.07/FG1_correct/try_al.re.pdf",7,7)
  a =  MyHeatmap(as.matrix(src.9ss.integrated.al@assays$RNA@data[names(gene_list_al.re), ]),
                 type = "log.row.relat",
                 RowSideColors = t(cbind(
                   MyName2Col(gene_list_al.re, cluster.endoderm.color.v5))),
                 ColSideColorsSize = 3,
                 RowSideColorsSize = 2,
                 # Rowv = "none", 
                 # Colv = "none",
                 # return.tree = "none",
                 # return.tree = "row", # gene_tree
                 return.tree = "col", # col_tree
                 labCol="none", graph = T)
  dev.off()
  
  al.coltree_1 = as.dendrogram(a) # cell
  cell_list_al = c(
    as.dendrogram(labels(al.coltree_1[[2]][[2]])),
    as.dendrogram(labels(al.coltree_1[[2]][[1]])),
    as.dendrogram(labels(al.coltree_1[[1]])))
  names(cell_list_al) = c(
    rep("AL.1", length(c(as.dendrogram(labels(al.coltree_1[[2]][[2]]))))),
    rep("AL.2", length(c(as.dendrogram(labels(al.coltree_1[[2]][[1]]))))),
    rep("AL.3", length(c(as.dendrogram(labels(al.coltree_1[[1]]))))))
  
  src.9ss.integrated_al$cluster.v06.26.re..correct = NA
  src.9ss.integrated_al@meta.data[cell_list_al,]$cluster.v06.26.re..correct = names(cell_list_al)
  
  
  pdf("figure.v08.07/FG1_correct/try_al.re.pdf",7,7)
  a =  MyHeatmap(as.matrix(src.9ss.integrated.al@assays$RNA@data[names(gene_list_al.re), ]),
                 type = "row.relat",
                 ColSideColors = cbind(
                   MyName2Col(src.9ss.integrated.al$cluster.v06.26.re..correct, 
                              cluster.endoderm.color.v5)),
                 RowSideColors = t(cbind(
                   MyName2Col(gene_list_al.re, cluster.endoderm.color.v5))),
                 ColSideColorsSize = 3,
                 RowSideColorsSize = 2,
                 # Rowv = "none", 
                 # Colv = "none",
                 # return.tree = "none",
                 return.tree = "row", # gene_tree
                 # return.tree = "col", # col_tree
                 labCol="none", graph = T)
  dev.off()
  al.rowtree_1 = as.dendrogram(a)
  
  gene_list_al3.temp = c(
    labels(gene_tree[[2]][[1]][[2]]),
    labels(gene_tree[[2]][[2]][[2]][[1]]),
    labels(gene_tree[[2]][[2]][[1]]),
    labels(gene_tree[[2]][[1]][[1]])) # gene_list_al.re.1: AL.3
  
  #>> analysis as PMID: 33106598 (Fig.3C)
  #>> order by expression
  #----------- 
  gene_test = gene_list_p3_re[names(gene_list_p3_re)=="AL.3"]
  data_test = src.9ss.integrated.al@assays$RNA@data[gene_test,]
  colnames(data_test) = colnames(src.9ss.integrated.al)
  data_test = by(data = t(data_test),
                 INDICES = src.9ss.integrated.al@meta.data[
                   colnames(data_test),"cluster.v06.26.re..correct"],
                 FUN = colSums)
  data_test_summary = do.call(cbind,data_test)
  data_test_summary = data_test_summary / 
    sqrt(table(gene_test) %*%
           t(table(src.9ss.integrated.al$cluster.v06.26.re..correct)))
  data_test_summary = as.data.frame(data_test_summary)
  data_test_summary$gene = rownames(data_test_summary)
  #----------
  
  gene_list_al = c(
    labels(al.rowtree_1[[1]][[2]][[2]]), # gene_list_al.re.1: AL.1
    labels(al.rowtree_1[[1]][[2]][[1]]), # gene_list_al.re.1: AL.1
    labels(al.rowtree_1[[1]][[1]]), # gene_list_al.re.1: AL.1 -> AL.2
    data_test_summary$gene[order(data_test_summary$AL.2,decreasing = T)] # gene_list_al.re.1: AL.3
    )
  names(gene_list_al) = c(
    rep("AL.1", length(c(labels(al.rowtree_1[[1]][[2]][[2]]), # gene_list_al.re.1: AL.1
                         labels(al.rowtree_1[[1]][[2]][[1]])))),
    rep("AL.2", length(c(labels(al.rowtree_1[[1]][[1]])))),
    rep("AL.3.1", length(intersect(
      data_test_summary$gene[order(data_test_summary$AL.2,decreasing = T)],
      data_test_summary$gene[data_test_summary$AL.2>=2.5]))),
    rep("AL.3.2", length(intersect(
      data_test_summary$gene[order(data_test_summary$AL.2,decreasing = T)],
      data_test_summary$gene[data_test_summary$AL.2<2.5]))))
  
  color.p3_gcn =c("#9970D7","#D190D8","#D1A955","#A8C2D5","#B6DABE")
  names(color.p3_gcn) = c("AL.1","AL.1.2","AL.2","AL.3.1","AL.3.2")
  
  #>> re-GCN
  seurat = src.9ss.integrated.al
  gene_list_sum = gene_list_al; p = 0
  #gene_list_sum = gene_list_p3_gcn; p = 0.55
  #----------------
  pearson.pathway = WGCNA::cor(t(seurat@assays$RNA@data[gene_list_sum,]), method = "p")
  pearson.pathway = logistic(pearson.pathway, threshold = 0)
  
  pearson.pathway[is.na(pearson.pathway)]=0
  pearson.pathway[pearson.pathway<p]=0
  # pearson.pathway[pearson.pathway<0.58]=0  # Used in Background
  graph.pathway = igraph::graph.adjacency(pearson.pathway,mode = "undirected",weighted = T)
  graph.pathway = igraph::simplify(graph.pathway,remove.multiple = T,remove.loops = T)
  
  graph.pathway = igraph::delete.vertices(
    graph.pathway,
    v = igraph::V(graph.pathway)[igraph::clusters(graph.pathway)$membership%in%which(igraph::clusters(graph.pathway)$csize<=10)])
  degree.Foregut.pathway=igraph::degree(graph.pathway)
  
  layout.pathway = igraph::layout_nicely(graph.pathway)
  graph.cluster <- igraph::cluster_walktrap(graph.pathway,steps = 10)
  gene.cluster <- graph.cluster$membership
  names(gene.cluster) <- graph.cluster$names
  set.seed(2)
  #----------------
  
  merge_gene_name = function(query,refer){
      
      genr = rep(NA,length(query))
      names(genr) = query
      
      for(i in query){
        genr[query] = ifelse(
          query%in%refer,
          names(refer[refer%in%query]),
          NA)
      }
      
      return(genr)
      
    }
  
  pdf("figure.v08.07/FG1_correct/10X_AL.GCN.pdf",9,7)
  igraph::plot.igraph(graph.pathway,
                      layout = layout.pathway,
                      vertex.size=7,
                      label.cex =8,
                      vertex.label = NA,
                      vertex.label.color = "black",
                      vertex.color = color.p3_gcn[merge_gene_name(
                        names(gene.cluster), gene_list_sum)])
  dev.off()
  
}

for(i.cell_type in c("MG")){
  src.9ss.integrated.mg = src.9ss.integrated[,src.9ss.integrated$cluster.v06.26.re%in%c("MG")]
  src.9ss.integrated.mg = FindVariableFeatures(src.9ss.integrated.mg,nfeatures = 2000)
  src.9ss.integrated.mg.filtergene =   
    Myfilter(as.matrix(src.9ss.integrated.mg@assays$RNA@data),
             gene = unique(c(src.9ss.integrated.mg @assays$RNA@var.features,
                             rownames(src.9ss.integrated.mg@reductions$pca@feature.loadings))),
             pearson.threshold = 0.15,partner.threshold = 5,
             bottom.dispersion.interval = 0.1)
  
  
  seurat = src.9ss.integrated.mg
  gene_list_sum = gsrc.9ss.integrated.mg.filtergene; p = 0.52
  #----------------
  pearson.pathway = WGCNA::cor(t(seurat@assays$RNA@data[gene_list_sum,]), method = "p")
  pearson.pathway = logistic(pearson.pathway, threshold = 0.25)
  pearson.pathway[is.na(pearson.pathway)]=0
  pearson.pathway[pearson.pathway<p]=0
  # pearson.pathway[pearson.pathway<0.58]=0  # Used in Background
  
  graph.pathway = graph.adjacency(pearson.pathway,mode = "undirected",weighted = T)
  graph.pathway = igraph::simplify(graph.pathway,remove.multiple = T,remove.loops = T)
  graph.pathway = delete.vertices(
    graph.pathway,
    v = V(graph.pathway)[clusters(graph.pathway)$membership%in%which(clusters(graph.pathway)$csize<=10)])
  degree.Foregut.pathway=igraph::degree(graph.pathway)
  layout.pathway = layout_nicely(graph.pathway)
  
  graph.cluster <- cluster_walktrap(graph.pathway,steps = 5)
  gene.cluster <- graph.cluster$membership
  names(gene.cluster) <- graph.cluster$names
  set.seed(3)
  
  pdf("gcn_9ss_mg.pdf",7,7)
  plot.igraph(graph.pathway,
              layout = layout.pathway,
              vertex.size=6,
              label.cex =6,
              vertex.label = NA, 
              vertex.label.color = "black",
              vertex.color = colors.num[gene.cluster]
              # vertex.color = cluster.endoderm.color.v5[gene_list_mg.re]
              )
  dev.off()
  #----------------
  
  gene_list_mg.re = c(
    gene.cluster[gene.cluster %in% c(1)],
    gene.cluster[gene.cluster %in% c(2)],
    gene.cluster[gene.cluster %in% c(3)])
  gene_list_mg.re[gene_list_mg.re%in%c(1)]="MG.1"
  gene_list_mg.re[gene_list_mg.re%in%c(2)]="MG.2"
  gene_list_mg.re[gene_list_mg.re%in%c(3)]="MG.3"
  
  gene_list_mg = names(gene_list_mg)
  names(gene_list_mg) = gene_list_mg
  
  pdf("heatmap_mg_9ss.pdf",7,7)
  #-------------
  data = src.9ss.integrated.mg@assays$RNA@data
  a = MyHeatmap(
    as.matrix(data[gene_list_mg,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    c.cov.method = "p",
    c.hc.method = "ward.D",
    RowSideColors = t(cbind(
      MyName2Col(
        names(gene_list_mg),
        cluster.endoderm.color.v5
      ))),
    # Rowv = "none", Colv = "none", 
    return.tree = "col",
    ColSideColorsSize = 2, RowSideColorsSize = 2,
    #labCol = name_list_all,  labRow = name_list_all,
    margins = c(8,8))
  #-------------
  dev.off()
  
  mg.coltree_1 = as.dendrogram(a) # cell
  cell_list_mg = c(
    as.dendrogram(labels(mg.coltree_1[[2]][[2]])),
    as.dendrogram(labels(mg.coltree_1[[2]][[1]])),
    as.dendrogram(labels(mg.coltree_1[[1]])))
  names(cell_list_mg) = c(
    rep("MG.1", length(c(as.dendrogram(labels(mg.coltree_1[[2]][[2]]))))),
    rep("MG.2", length(c(as.dendrogram(labels(mg.coltree_1[[2]][[1]]))))),
    rep("MG.3", length(c(as.dendrogram(labels(mg.coltree_1[[1]]))))))
  
  src.9ss.integrated_mg$cluster.v06.26.re..correct = NA
  src.9ss.integrated_mg@meta.data[cell_list_mg,]$cluster.v06.26.re..correct = names(cell_list_mg)
}

for(i.cell_type in c("HG")){
  src.9ss.integrated.hg = src.9ss.integrated[,src.9ss.integrated$cluster.v06.26.re%in%c("HG")]
  src.9ss.integrated.hg = FindVariableFeatures(src.9ss.integrated.hg,nfeatures = 2000)
  src.9ss.integrated.hg.filtergene =   
    Myfilter(as.matrix(src.9ss.integrated.hg@assays$RNA@data),
             gene = unique(c(src.9ss.integrated.hg @assays$RNA@var.features,
                             rownames(src.9ss.integrated.hg@reductions$pca@feature.loadings))),
             pearson.threshold = 0.15,partner.threshold = 5,
             bottom.dispersion.interval = 0.1)
  
  
  seurat = src.9ss.integrated.hg
  gene_list_sum = gsrc.9ss.integrated.hg.filtergene; p = 0.52
  #----------------
  pearson.pathway = WGCNA::cor(t(seurat@assays$RNA@data[gene_list_sum,]), method = "p")
  pearson.pathway = logistic(pearson.pathway, threshold = 0.25)
  pearson.pathway[is.na(pearson.pathway)]=0
  pearson.pathway[pearson.pathway<p]=0
  # pearson.pathway[pearson.pathway<0.58]=0  # Used in Background
  
  graph.pathway = graph.adjacency(pearson.pathway,mode = "undirected",weighted = T)
  graph.pathway = igraph::simplify(graph.pathway,remove.multiple = T,remove.loops = T)
  graph.pathway = delete.vertices(
    graph.pathway,
    v = V(graph.pathway)[clusters(graph.pathway)$membership%in%which(clusters(graph.pathway)$csize<=10)])
  degree.Foregut.pathway=igraph::degree(graph.pathway)
  layout.pathway = layout_nicely(graph.pathway)
  
  graph.cluster <- cluster_walktrap(graph.pathway,steps = 5)
  gene.cluster <- graph.cluster$membership
  names(gene.cluster) <- graph.cluster$names
  set.seed(3)
  
  pdf("gcn_9ss_hg.pdf",7,7)
  plot.igraph(graph.pathway,
              layout = layout.pathway,
              vertex.size=6,
              label.cex =6,
              vertex.label = NA, 
              vertex.label.color = "black",
              vertex.color = colors.num[gene.cluster]
              # vertex.color = cluster.endoderm.color.v5[gene_list_hg.re]
  )
  dev.off()
  #----------------
  
  gene_list_hg.re = c(
    gene.cluster[gene.cluster %in% c(2)],
    gene.cluster[gene.cluster %in% c(3)])
  gene_list_hg.re[gene_list_hg.re%in%c(2)]="HG.1"
  gene_list_hg.re[gene_list_hg.re%in%c(3)]="HG.2"
  
  gene_list_hg = names(gene_list_hg)
  names(gene_list_hg) = gene_list_hg
  
  pdf("heatmap_hg_9ss.pdf",7,7)
  #-------------
  data = src.9ss.integrated.hg@assays$RNA@data
  a = MyHeatmap(
    as.matrix(data[gene_list_hg,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    c.cov.method = "p",
    c.hc.method = "ward.D",
    RowSideColors = t(cbind(
      MyName2Col(
        names(gene_list_hg),
        cluster.endoderm.color.v5
      ))),
    # Rowv = "none", Colv = "none", 
    return.tree = "col",
    ColSideColorsSize = 2, RowSideColorsSize = 2,
    #labCol = name_list_all,  labRow = name_list_all,
    margins = c(8,8))
  #-------------
  dev.off()
  
  hg.coltree_1 = as.dendrogram(a) # cell
  cell_list_hg = c(
    as.dendrogram(labels(hg.coltree_1[[1]])),
    as.dendrogram(labels(hg.coltree_1[[2]])))
  names(cell_list_hg) = c(
    rep("HG.1", length(c(as.dendrogram(labels(hg.coltree_1[[1]]))))),
    rep("HG.2", length(c(as.dendrogram(labels(hg.coltree_1[[2]]))))))
  
  src.9ss.integrated_hg$cluster.v06.26.re..correct = NA
  src.9ss.integrated_hg@meta.data[cell_list_hg,]$cluster.v06.26.re..correct = names(cell_list_hg)
}

src.9ss.integrated$cluster.v06.26.re..correct = NA
src.9ss.integrated@meta.data[colnames(src.9ss.integrated_fg),]$cluster.v06.26.re..correct = src.9ss.integrated_fg$cluster.v06.26.re..correct
src.9ss.integrated@meta.data[colnames(src.9ss.integrated_al),]$cluster.v06.26.re..correct = src.9ss.integrated_al$cluster.v06.26.re..correct
src.9ss.integrated@meta.data[colnames(src.9ss.integrated_mg),]$cluster.v06.26.re..correct = src.9ss.integrated_mg$cluster.v06.26.re..correct
src.9ss.integrated@meta.data[colnames(src.9ss.integrated_hg),]$cluster.v06.26.re..correct = src.9ss.integrated_hg$cluster.v06.26.re..correct

src.9ss.integrated$cluster.endoderm = src.9ss.integrated$cluster.v06.26.re..correct
src.9ss.integrated$cluster.extract.v1.1 = src.9ss.integrated$cluster.v06.26.re..correct
#===============================================================================


#===============================================================================
#>> 2.2 27SS cell type
#===============================================================================
# identifying structures such as the pharynx (Pha), thyroid (Thy), esophagus 
# (Eso), lung (Lun), stomach (Sto), liver (Liv), extrahepatic biliary duct (EHBD), pancreas (Pan), 
# small intestine (SI), and large intestine (LI). 

src.27ss.integrated = FindNeighbors(src.27ss.integrated, reduction = "pca", dims = 1:30)
src.27ss.integrated = FindClusters(src.27ss.integrated, resolution = 3.5)
src.27ss.integrated = FindClusters(src.27ss.integrated, resolution = 8)
src.27ss.integrated = FindClusters(src.27ss.integrated, resolution = 10)

src.27ss.integrated$cluster.v06.26.re = NA
src.27ss.integrated@meta.data[src.27ss.integrated$RNA_snn_res.3.5%in%c(4,8,10,12,13,16,25,26),]$cluster.v06.26.re = "Pharynx.organ"
src.27ss.integrated@meta.data[src.27ss.integrated$RNA_snn_res.3.5%in%c(4),]$cluster.v06.26.re = "Thyroid"
src.27ss.integrated@meta.data[src.27ss.integrated$RNA_snn_res.3.5%in%c(9,10,23),]$cluster.v06.26.re = "Lung"
src.27ss.integrated@meta.data[src.27ss.integrated$RNA_snn_res.3.5%in%c(11),]$cluster.v06.26.re = "Esophagus"
src.27ss.integrated@meta.data[src.27ss.integrated$RNA_snn_res.3.5%in%c(5,15,24),]$cluster.v06.26.re = "Stomach"
src.27ss.integrated@meta.data[src.27ss.integrated$RNA_snn_res.3.5%in%c(18),]$cluster.v06.26.re = "EHBD" 
src.27ss.integrated@meta.data[src.27ss.integrated$RNA_snn_res.3.5%in%c(17,19,20,21),]$cluster.v06.26.re = "Liver"
src.27ss.integrated@meta.data[src.27ss.integrated$RNA_snn_res.3.5%in%c(18,22,27),]$cluster.v06.26.re = "Pancreas"
src.27ss.integrated@meta.data[src.27ss.integrated$RNA_snn_res.3.5%in%c(0,1,2,3,28),]$cluster.v06.26.re = "Small.intestine"
src.27ss.integrated@meta.data[src.27ss.integrated$RNA_snn_res.3.5%in%c(6,7,14),]$cluster.v06.26.re = "Large.intestine"

#>> correct res.5: 4 & 18
src.27ss.integrated@meta.data[src.27ss.integrated$RNA_snn_res.3.5%in%c(4)&
                                src.27ss.integrated$RNA_snn_res.10%in%c(0),]$cluster.v06.26.re = "Pharynx.organ"
src.27ss.integrated@meta.data[src.27ss.integrated$RNA_snn_res.3.5%in%c(4)&
                                !src.27ss.integrated$RNA_snn_res.10%in%c(0),]$cluster.v06.26.re = "Thyroid"
src.27ss.integrated@meta.data[src.27ss.integrated$RNA_snn_res.3.5%in%c(18)&
                                src.27ss.integrated$RNA_snn_res.10%in%c(30),]$cluster.v06.26.re = "EHBD" 
src.27ss.integrated@meta.data[src.27ss.integrated$RNA_snn_res.3.5%in%c(18)&
                                !src.27ss.integrated$RNA_snn_res.10%in%c(30),]$cluster.v06.26.re = "Pancreas" 

for(i.cell_type in c("Pharynx.organ")){
  src.27ss.integrated_pha = src.27ss.integrated[,src.27ss.integrated$cluster.v06.26.re%in%"Pharynx.organ"]
  src.27ss.integrated_pha = FindVariableFeatures(src.27ss.integrated_pha, nfeatures = 3000)
  src.27ss.integrated_pha.filtergene = 
    Myfilter(as.matrix(src.27ss.integrated_pha@assays$RNA@data),
             gene = src.27ss.integrated_pha@assays$RNA@var.features,
             pearson.threshold = 0.15,partner.threshold = 5,
             bottom.dispersion.interval = 0.1)
  
  gene_list_sum = unique(src.27ss.integrated_pha.filtergene); p = 0.6
  #----------------
  seurat = src.27ss.integrated_pha
  pearson.gene_pha = WGCNA::cor(t(seurat@assays$RNA@data[gene_list_sum,]), method = "p")
  pearson.gene_pha = logistic(pearson.gene_pha, threshold = 0.25)
  pearson.gene_pha[is.na(pearson.gene_pha)]=0
  pearson.gene_pha[pearson.gene_pha<p]=0
  # pearson.gene_pha[pearson.gene_pha<0.58]=0  # Used in Background
  
  graph.gene_pha = graph.adjacency(pearson.gene_pha,mode = "undirected",weighted = T)
  graph.gene_pha = igraph::simplify(graph.gene_pha,remove.multiple = T,remove.loops = T)
  graph.gene_pha = delete.vertices(
    graph.gene_pha,
    v = V(graph.gene_pha)[clusters(graph.gene_pha)$membership %in%
                            which(clusters(graph.gene_pha)$csize<=10)])
  degree.Foregut.gene_pha=igraph::degree(graph.gene_pha)
  layout.gene_pha = layout_nicely(graph.gene_pha)
  
  graph.cluster_pha <- cluster_walktrap(graph.gene_pha,steps = 5)
  gene.cluster_pha <- graph.cluster_pha$membership
  names(gene.cluster_pha) <- graph.cluster_pha$names
  set.seed(3)
  
  pdf("gcn_pha_27ss.pdf",7,7)
  plot.igraph(graph.gene_pha,
              layout = layout.gene_pha,
              vertex.size=6,
              label.cex =6,
              #vertex.label = NA, # gene.cluster_pha,
              vertex.label = gene.cluster_pha,
              vertex.label.color = "black",
              vertex.color = colors.num[gene.cluster_pha])
  #vertex.color = cluster.endoderm.color.v5[gene_list_pha])
  dev.off()
  #----------------
  
  gene_list_pha.re = c(
    gene.cluster_pha[gene.cluster_pha %in% c(3)],
    gene.cluster_pha[gene.cluster_pha %in% c(2)],
    gene.cluster_pha[gene.cluster_pha %in% c(1)],
    gene.cluster_pha[gene.cluster_pha %in% c(4)])
  gene_list_pha.re[gene_list_pha.re%in%c(3)]="Pharynx.organ.1"
  gene_list_pha.re[gene_list_pha.re%in%c(2)]="Pharynx.organ.2"
  gene_list_pha.re[gene_list_pha.re%in%c(1)]="Pharynx.organ.3"
  gene_list_pha.re[gene_list_pha.re%in%c(4)]="Pharynx.organ.4"
  
  gene_list_pha = names(gene_list_pha.re)
  names(gene_list_pha) = gene_list_pha.re
  
  pdf("heatmap_pha_27ss.pdf",7,7)
  #-------------
  data = src.27ss.integrated_pha@assays$RNA@data
  a = MyHeatmap(
    as.matrix(data[gene_list_pha,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    c.cov.method = "p",
    c.hc.method = "ward.D",
    RowSideColors = t(cbind(
      MyName2Col(
        names(gene_list_pha),
        cluster.endoderm.color.v5
      ))),
    # Rowv = "none", Colv = "none", 
    return.tree = "col",
    ColSideColorsSize = 2, RowSideColorsSize = 2,
    #labCol = name_list_all,  labRow = name_list_all,
    margins = c(8,8))
  #-------------
  dev.off()
  
  pha.coltree_1 = as.dendrogram(a) # cell
  cell_list_pha = c(
    as.dendrogram(labels(pha.coltree_1[[2]][[1]])),
    as.dendrogram(labels(pha.coltree_1[[1]][[1]])),
    as.dendrogram(labels(pha.coltree_1[[1]][[2]])),
    as.dendrogram(labels(pha.coltree_1[[2]][[2]][[1]])),
    as.dendrogram(labels(pha.coltree_1[[2]][[2]][[2]][[1]][[1]])),
    as.dendrogram(labels(pha.coltree_1[[2]][[2]][[2]][[1]][[2]])),
    as.dendrogram(labels(pha.coltree_1[[2]][[2]][[2]][[2]])))
  names(cell_list_pha) = c(
    rep("Pharynx.organ.2", length(c(as.dendrogram(labels(pha.coltree_1[[2]][[1]]))))),
    rep("Pharynx.organ.1", length(c(as.dendrogram(labels(pha.coltree_1[[1]][[1]]))))),
    rep("Pharynx.organ.4", length(c(as.dendrogram(labels(pha.coltree_1[[1]][[2]]))))),
    rep("Pharynx.organ.5", length(c(as.dendrogram(labels(pha.coltree_1[[2]][[2]][[1]])),
                                    as.dendrogram(labels(pha.coltree_1[[2]][[2]][[2]][[1]][[1]]))))),
    rep("Pharynx.organ.3", length(c(as.dendrogram(labels(pha.coltree_1[[2]][[2]][[2]][[1]][[2]])),
                                    as.dendrogram(labels(pha.coltree_1[[2]][[2]][[2]][[2]]))))))
  
  src.27ss.integrated_pha$cluster.v06.26.re..correct = NA
  src.27ss.integrated_pha@meta.data[cell_list_pha,]$cluster.v06.26.re..correct = names(cell_list_pha)
  
}

for(i.cell_type in c("Pancreas")){
  src.27ss.integrated_pan = src.27ss.integrated[,src.27ss.integrated$cluster.v06.26.re%in%c("Pancreas")]
  src.27ss.integrated_pan = src.27ss.integrated_pan[,src.27ss.integrated_pan$nFeature_RNA>2500]
  src.27ss.integrated_pan = FindVariableFeatures(src.27ss.integrated_pan, nfeatures = 3000)
  
  src.27ss.integrated_pan.filtergene = 
    Myfilter(as.matrix(src.27ss.integrated_pan@assays$RNA@data),
             gene = src.27ss.integrated_pan@assays$RNA@var.features,
             pearson.threshold = 0.15,partner.threshold = 5,
             bottom.dispersion.interval = 0.1)
  
  src.27ss.integrated_pan = SetIdent(src.27ss.integrated_pan,
                                     value = src.27ss.integrated_pan$cluster.v06.26.re..merge)
  
  gene_list_sum = src.27ss.integrated_pan.filtergene; p = 0.55
  #----------------
  seurat = src.27ss.integrated_pan
  pearson.gene_pan = WGCNA::cor(t(seurat@assays$RNA@data[gene_list_sum,]), method = "p")
  pearson.gene_pan = logistic(pearson.gene_pan, threshold = 0.25)
  pearson.gene_pan[is.na(pearson.gene_pan)]=0
  pearson.gene_pan[pearson.gene_pan<p]=0
  # pearson.gene_pan[pearson.gene_pan<0.58]=0  # Used in Background
  
  graph.gene_pan = graph.adjacency(pearson.gene_pan,mode = "undirected",weighted = T)
  graph.gene_pan = igraph::simplify(graph.gene_pan,remove.multiple = T,remove.loops = T)
  graph.gene_pan = delete.vertices(
    graph.gene_pan,
    v = V(graph.gene_pan)[clusters(graph.gene_pan)$membership %in%
                            which(clusters(graph.gene_pan)$csize<=10)])
  degree.Foregut.gene_pan=igraph::degree(graph.gene_pan)
  layout.gene_pan = layout_nicely(graph.gene_pan)
  
  graph.cluster_pan <- cluster_walktrap(graph.gene_pan,steps = 3)
  gene.cluster_pan <- graph.cluster_pan$membership
  names(gene.cluster_pan) <- graph.cluster_pan$names
  set.seed(1)
  
  pdf("gcn_pan_27ss.pdf",7,7)
  plot.igraph(graph.gene_pan,
              layout = layout.gene_pan,
              vertex.size=6,
              label.cex =6,
              vertex.label = NA, # gene.cluster_pan,
              vertex.label.color = "black",
              #vertex.color = colors.num[gene.cluster_pan])
              vertex.color = cluster.endoderm.color.v5[
                re_gene_list(gene.cluster_27ss_pan)[names(gene.cluster_pan)]
              ])
  dev.off()
  #----------------
  
  gene_list_pan.re = c(
    gene.cluster_pan[gene.cluster_pan %in% c(3)],
    gene.cluster_pan[gene.cluster_pan %in% c(1,2)],
    gene.cluster_pan[gene.cluster_pan %in% c(8)])
  gene_list_pan.re[gene_list_pan.re%in%c(3)]="DP"
  gene_list_pan.re[gene_list_pan.re%in%c(1,2)]="EP"
  gene_list_pan.re[gene_list_pan.re%in%c(8)]="VP"
  
  gene_list_pan = names(gene_list_pan.re)
  names(gene_list_pan) = gene_list_pan.re
  
  pdf("heatmap_pan_27ss.pdf",7,7)
  #-------------
  data = src.27ss.integrated_pan@assays$RNA@data
  a = MyHeatmap(
    as.matrix(data[gene_list_pan,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    c.cov.method = "p",
    c.hc.method = "ward.D",
    RowSideColors = t(cbind(
      MyName2Col(
        names(gene_list_pan),
        cluster.endoderm.color.v5
      ))),
    # Rowv = "none", Colv = "none", 
    return.tree = "col",
    ColSideColorsSize = 2, RowSideColorsSize = 2,
    #labCol = name_list_all,  labRow = name_list_all,
    margins = c(8,8))
  #-------------
  dev.off()
  
  pan.coltree_1 = as.dendrogram(a) # cell
  cell_list_pan = c(
    as.dendrogram(labels(pan.coltree_1[[2]][[2]])),
    as.dendrogram(labels(pan.coltree_1[[2]][[1]][[2]])),
    as.dendrogram(labels(pan.coltree_1[[2]][[1]][[1]])),
    as.dendrogram(labels(pan.coltree_1[[1]])))
  names(cell_list_pan) = c(
    rep("DP", length(c(as.dendrogram(labels(pan.coltree_1[[2]][[2]])),
                       as.dendrogram(labels(pan.coltree_1[[2]][[1]][[2]]))))),
    rep("EP", length(as.dendrogram(labels(pan.coltree_1[[2]][[1]][[1]])))),
    rep("VP", length(as.dendrogram(labels(pan.coltree_1[[1]])))))
  
  src.27ss.integrated_pan$cluster.v06.26.re..correct = NA
  src.27ss.integrated_pan@meta.data[cell_list_pan,]$cluster.v06.26.re..correct = names(cell_list_pan)
  
  
}

for(i.cell_type in c("EP")){
  src.27ss.integrated_ep = src.27ss.integrated_pan[,src.27ss.integrated_pan$cluster.v06.26.re..correct%in%c("EP")]
  src.27ss.integrated_ep = FindVariableFeatures(src.27ss.integrated_ep, nfeatures = 3000)
  
  src.27ss.integrated_ep.filtergene = 
    Myfilter(as.matrix(src.27ss.integrated_ep@assays$RNA@data),
             gene = src.27ss.integrated_ep@assays$RNA@var.features,
             pearson.threshold = 0.15,partner.threshold = 5,
             bottom.dispersion.interval = 0.1)
  
  gene_list_sum = src.27ss.integrated_ep.filtergene; p = 0.55
  #----------------
  seurat = src.27ss.integrated_ep
  pearson.gene_ep = WGCNA::cor(t(seurat@assays$RNA@data[gene_list_sum,]), method = "p")
  pearson.gene_ep = logistic(pearson.gene_ep, threshold = p)
  pearson.gene_ep[is.na(pearson.gene_ep)]=0
  pearson.gene_ep[pearson.gene_ep<p]=0
  # pearson.gene_ep[pearson.gene_ep<0.58]=0  # Used in Background
  
  graph.gene_ep = graph.adjacency(pearson.gene_ep,mode = "undirected",weighted = T)
  graph.gene_ep = igraph::simplify(graph.gene_ep,remove.multiple = T,remove.loops = T)
  graph.gene_ep = delete.vertices(
    graph.gene_ep,
    v = V(graph.gene_ep)[clusters(graph.gene_ep)$membership %in%
                           which(clusters(graph.gene_ep)$csize<=10)])
  degree.Foregut.gene_ep=igraph::degree(graph.gene_ep)
  layout.gene_ep = layout_nicely(graph.gene_ep)
  
  graph.cluster_ep <- cluster_walktrap(graph.gene_ep,steps = 2)
  gene.cluster_ep <- graph.cluster_ep$membership
  names(gene.cluster_ep) <- graph.cluster_ep$names
  set.seed(3)
  
  pdf("GCN_check/gcn_ep_27ss.pdf",7,7)
  plot.igraph(graph.gene_ep,
              layout = layout.gene_ep,
              vertex.size=6,
              label.cex =6,
              vertex.label = gene.cluster_ep,
              vertex.label.color = "black",
              vertex.color = colors.geneset[gene.cluster_ep ])
  dev.off()
  #----------------
  
  gene_list_ep.re = c(
    gene.cluster_ep[gene.cluster_ep %in% c(1)],
    rev(gene.cluster_ep[gene.cluster_ep %in% c(2,5,7)]))
  gene_list_ep.re[gene_list_ep.re%in%c(1)]="EP.1"
  gene_list_ep.re[gene_list_ep.re%in%c(2,5,7)]="EP.2"

  gene_list_ep = names(gene_list_ep.re)
  names(gene_list_ep.re) = gene_list_ep.re
  
  data = src.27ss.integrated_ep@assays$RNA@data
  pdf("GCN_check/heatmap_ep_27ss.pdf",7,7)
  #-------------
  a = MyHeatmap(as.matrix(data[gene_list_ep,]),
                type = "row.relat",
                hc.c.data.type = "row.relat",
                c.cov.method = "p",
                c.hc.method = "ward.D",
                RowSideColors = t(cbind(
                  MyName2Col(gene_list_ep, cluster.endoderm.color.v5))),
                # Rowv = "none", Colv = "none", 
                return.tree = "col",
                ColSideColorsSize = 2, RowSideColorsSize = 2,
                #labCol = name_list_all,  labRow = name_list_all,
                margins = c(8,8))
  #-------------
  dev.off()
  
  ep.coltree_1 = as.dendrogram(a) # cell
  cell_list_ep = c(
    as.dendrogram(labels(ep.coltree_1[[2]])),
    as.dendrogram(labels(ep.coltree_1[[1]])))
  names(cell_list_ep) = c(
    rep("EP.1", length(as.dendrogram(labels(ep.coltree_1[[2]])))),
    rep("EP.2", length(as.dendrogram(labels(ep.coltree_1[[1]])))))
  
  src.27ss.integrated_ep$cluster.v06.26.re..correct = NA
  src.27ss.integrated_ep@meta.data[cell_list_ep,]$cluster.v06.26.re..correct = names(cell_list_ep)
  
}

for(i.cell_type in c("Small.intestine")){
  src.27ss.integrated_sm = src.27ss.integrated[,src.27ss.integrated$cluster.v06.26.re%in%"Small.intestine"]
  src.27ss.integrated_sm = FindVariableFeatures(src.27ss.integrated_sm, nfeatures = 3000)
  src.27ss.integrated_sm.filtergene = 
    Myfilter(as.matrix(src.27ss.integrated_sm@assays$RNA@data),
             gene = src.27ss.integrated_sm@assays$RNA@var.features,
             pearson.threshold = 0.15,partner.threshold = 5,
             bottom.dispersion.interval = 0.1)
  
  gene_list_sum = src.27ss.integrated_sm.filtergene; p = 0.55
  gene_list_sum = setdiff(
    unique(src.27ss.integrated_sm.filtergene),
    src.27ss.integrated_sm.filtergene[grepl("mt-", src.27ss.integrated_sm.filtergene)])
  #----------------
  seurat = src.27ss.integrated_sm
  pearson.gene_sm = WGCNA::cor(t(seurat@assays$RNA@data[gene_list_sum,]), method = "p")
  pearson.gene_sm = logistic(pearson.gene_sm, threshold = 0.25)
  pearson.gene_sm[is.na(pearson.gene_sm)]=0
  pearson.gene_sm[pearson.gene_sm<p]=0
  # pearson.gene_sm[pearson.gene_sm<0.58]=0  # Used in Background
  
  graph.gene_sm = graph.adjacency(pearson.gene_sm,mode = "undirected",weighted = T)
  graph.gene_sm = igraph::simplify(graph.gene_sm,remove.multiple = T,remove.loops = T)
  graph.gene_sm = delete.vertices(
    graph.gene_sm,
    v = V(graph.gene_sm)[clusters(graph.gene_sm)$membership %in%
                           which(clusters(graph.gene_sm)$csize<=10)])
  degree.Foregut.gene_sm=igraph::degree(graph.gene_sm)
  layout.gene_sm = layout_nicely(graph.gene_sm)
  
  graph.cluster_sm <- cluster_walktrap(graph.gene_sm,steps = 3)
  gene.cluster_sm <- graph.cluster_sm$membership
  names(gene.cluster_sm) <- graph.cluster_sm$names
  set.seed(3)
  
  pdf("gcn_sm_27ss.update.pdf",7,7)
  plot.igraph(graph.gene_sm,
              layout = layout.gene_sm,
              vertex.size=6.5,
              label.cex =7,
              
              # vertex.size=5,
              # vertex.label.cex = 1.5, 
              # vertex.label.font= 4,
              # vertex.label.color = "black",
              # vertex.frame.color = "black", 
              # vertex.frame.width = 0.5,
              
              vertex.label = NA, #gene.cluster_sm,
              # vertex.label =  ifelse(names(gene.cluster_sm)%in%gi[gi$TF%in%TRUE,]$SymbolDedu,
              #                        names(gene.cluster_sm),NA),
              vertex.label.color = "black",
              # vertex.color = colors.num[gene.cluster_sm])
              vertex.color = cluster.endoderm.color.v5[
                gene_list_sm.re[names(gene.cluster_sm)]
              ])
  dev.off()
  #----------------
  
  gene_list_sm.re = c(
    gene.cluster_sm[gene.cluster_sm ==1],
    gene.cluster_sm[gene.cluster_sm %in% c(2,3)])
  gene_list_sm.re[gene_list_sm.re==1]="Small.intestine.1"
  gene_list_sm.re[gene_list_sm.re%in%c(2,3)]="Small.intestine.2"
  
  gene_list_sm = names(gene_list_sm.re)
  names(gene_list_sm) = gene_list_sm.re
  
  pdf("heatmap_sm_27ss.pdf",7,7)
  #-------------
  data = src.27ss.integrated_sm@assays$RNA@data
  a = MyHeatmap(as.matrix(data[
    gene_list_sm,]),
    type = "row.relat",
    hc.c.data.type = "row.relat",
    c.cov.method = "p",
    c.hc.method = "ward.D",
    RowSideColors = t(cbind(
      MyName2Col(
        names(gene_list_sm),
        cluster.endoderm.color.v5
      ))),
    # Rowv = "none", Colv = "none", 
    return.tree = "col",
    ColSideColorsSize = 2, RowSideColorsSize = 2,
    #labCol = name_list_all,  labRow = name_list_all,
    margins = c(8,8))
  #-------------
  dev.off()
  
  sm.coltree_1 = as.dendrogram(a) # cell
  cell_list_sm = c(
    as.dendrogram(labels(sm.coltree_1[[2]])),
    as.dendrogram(labels(sm.coltree_1[[1]])))
  names(cell_list_sm) = c(
    rep("Small.intestine.1", length(as.dendrogram(labels(sm.coltree_1[[2]])))),
    rep("Small.intestine.2", length(as.dendrogram(labels(sm.coltree_1[[1]])))))
  
  
  src.27ss.integrated_sm$cluster.v06.26.re..correct = NA
  src.27ss.integrated_sm@meta.data[cell_list_sm,]$cluster.v06.26.re..correct = names(cell_list_sm)
 
}

for(i.cell_type in c("Large.intestine")){
  src.27ss.integrated_lar = src.27ss.integrated[,src.27ss.integrated$cluster.v06.26.re%in%"Large.intestine"]
  src.27ss.integrated_lar = FindVariableFeatures(src.27ss.integrated_lar, nfeatures = 3000)
  src.27ss.integrated_lar.filtergene = 
    Myfilter(as.matrix(src.27ss.integrated_lar@assays$RNA@data),
             gene = src.27ss.integrated_lar@assays$RNA@var.features,
             pearson.threshold = 0.15,partner.threshold = 5,
             bottom.dispersion.interval = 0.1)
  
  gene_list_sum = src.27ss.integrated_lar.filtergene; p = 0.55
  gene_list_sum = unique(src.27ss.integrated_lar.filtergene)
  #----------------
  seurat = src.27ss.integrated_lar
  pearson.gene_lar = WGCNA::cor(t(seurat@assays$RNA@data[gene_list_sum,]), method = "p")
  pearson.gene_lar = logistic(pearson.gene_lar, threshold = 0.25)
  pearson.gene_lar[is.na(pearson.gene_lar)]=0
  pearson.gene_lar[pearson.gene_lar<p]=0
  
  graph.gene_lar = graph.adjacency(pearson.gene_lar,mode = "undirected",weighted = T)
  graph.gene_lar = igraph::simplify(graph.gene_lar,remove.multiple = T,remove.loops = T)
  graph.gene_lar = delete.vertices(
    graph.gene_lar,
    v = V(graph.gene_lar)[clusters(graph.gene_lar)$membership %in%
                            which(clusters(graph.gene_lar)$csize<=10)])
  degree.Foregut.gene_lar=igraph::degree(graph.gene_lar)
  layout.gene_lar = layout_nicely(graph.gene_lar)
  
  graph.cluster_lar <- cluster_walktrap(graph.gene_lar,steps = 3)
  gene.cluster_lar <- graph.cluster_lar$membership
  names(gene.cluster_lar) <- graph.cluster_lar$names
  set.seed(3)
  
  pdf("gcn_lar_27ss.pdf",7,7)
  plot.igraph(graph.gene_lar,
              layout = layout.gene_lar,
              vertex.size=6,
              label.cex =6,
              vertex.label = NA, #gene.cluster_lar,
              vertex.label.color = "black",
              # vertex.color = colors.num[gene.cluster_lar])
              vertex.color = cluster.endoderm.color.v5[gene_list_lar])
  dev.off()
  #----------------
  
  gene_list_lar.re = c(
    gene.cluster_lar[gene.cluster_lar ==1],
    gene.cluster_lar[gene.cluster_lar ==3],
    gene.cluster_lar[gene.cluster_lar ==2])
  
  gene_list_lar.re[gene_list_lar.re==1]="Large.intestine.1"
  gene_list_lar.re[gene_list_lar.re==3]="Large.intestine.2"
  gene_list_lar.re[gene_list_lar.re==2]="Large.intestine.3"
  
  gene_list_lar = names(gene_list_lar.re)
  names(gene_list_lar) = gene_list_lar.re
  
  data = src.27ss.integrated_lar@assays$RNA@data
  pdf("heatmap_lar_27ss.pdf",7,7)
  #-------------
  a = MyHeatmap(as.matrix(data[names(gene_list_lar),]),
                type = "row.relat",
                hc.c.data.type = "row.relat",
                c.cov.method = "p",
                c.hc.method = "ward.D",
                RowSideColors = t(cbind(
                  MyName2Col(gene_list_lar, cluster.endoderm.color.v5))),
                Rowv = "none",
                Colv = "none", 
                return.tree = "col",
                ColSideColorsSize = 2, RowSideColorsSize = 2,
                #labCol = name_list_all,  labRow = name_list_all,
                margins = c(8,8))
  #-------------
  dev.off()
  
  lar.coltree_1 = as.dendrogram(a) # cell
  cell_list_lar = c(
    as.dendrogram(labels(lar.coltree_1[[2]][[2]])),
    as.dendrogram(labels(lar.coltree_1[[2]][[1]])),
    as.dendrogram(labels(lar.coltree_1[[1]])))
  names(cell_list_lar) = c(
    rep("Large.intestine.1", length(as.dendrogram(labels(lar.coltree_1[[2]][[1]])))),
    rep("Large.intestine.2", length(as.dendrogram(labels(lar.coltree_1[[2]][[2]])))),
    rep("Large.intestine.3", length(as.dendrogram(labels(lar.coltree_1[[1]])))))
  
  
  src.27ss.integrated_lar$cluster.v06.26.re..correct = NA
  src.27ss.integrated_lar@meta.data[cell_list_lar,]$cluster.v06.26.re..correct = names(cell_list_lar)
}

src.27ss.integrated$cluster.v06.26.re..correct = NA
src.27ss.integrated@meta.data[colnames(src.27ss.integrated_pha),]$cluster.v06.26.re..correct = src.27ss.integrated_pha$cluster.v06.26.re..correct
src.27ss.integrated@meta.data[colnames(src.27ss.integrated_pan),]$cluster.v06.26.re..correct = src.27ss.integrated_pan$cluster.v06.26.re..correct
src.27ss.integrated@meta.data[colnames(src.27ss.integrated_ep),]$cluster.v06.26.re..correct = src.27ss.integrated_ep$cluster.v06.26.re..correct
src.27ss.integrated@meta.data[colnames(src.27ss.integrated_sm),]$cluster.v06.26.re..correct = src.27ss.integrated_sm$cluster.v06.26.re..correct
src.27ss.integrated@meta.data[colnames(src.27ss.integrated_lar),]$cluster.v06.26.re..correct = src.27ss.integrated_lar$cluster.v06.26.re..correct
#===============================================================================


#===============================================================================
#>> 2.3 27SS cell type for tracing cells & CCA
#===============================================================================

list.endoderm.name = c("Pharynx.organ.1","Pharynx.organ.2","Pharynx.organ.3",
                       "Pharynx.organ.4","Pharynx.organ.5",
                       "DP","VP","EP.1","EP.2",
                       "Small.intestine.1","Small.intestine.2",
                       "Large.intestine.1","Large.intestine.2","Large.intestine.3",
                       'Thyroid',"Lung","Stomach","Esophagus","Liver","EHBD")
names(list.endoderm.name) = c("pha.1","pha.2","pha.3","pha.4","pha.5",
                              "pan.1","pan.2","pan.4","pan.5",
                              "smai.1","smai.2","lari.1","lari.2","lari.3",
                              "other.1","other.2","other.3","other.4","other.5","other.6")

src.merge.tracing.27ss.all = src.sm3.merge[,colnames(src.sm3.merge)[src.sm3.merge$Time=="27ss"]]
#>> Need cell-type [Supplymental Table 1, Sheet4]
#>> Store them in: src.merge.tracing.27ss.all$cluster.v06.26.re_correct

#>> detail code of integration: 4.inference_on_LAI.R
src.merge.tracing.27ss.all$cluster.v06.26.re_correct = as.character(
  FNN::knn(
    src.27ss.integrated.merge@reductions$umap_integrated@cell.embeddings[
      colnames(src.27ss.integrated),],
    src.27ss.integrated.merge@reductions$umap_integrated@cell.embeddings[
      colnames(src.merge.tracing.27ss.all),],
    src.27ss.integrated$cluster.v06.26.re..correct, k = 10))

for(i.major_cluster in c("GCN for major cluster")){
  src.merge.tracing.27ss.all.pha = src.merge.tracing.27ss.all[,src.merge.tracing.27ss.all$cluster.v06.26.re_correct%in%c("Pharynx.organ.1","Pharynx.organ.2",'Pharynx.organ.3',"Pharynx.organ.4","Pharynx.organ.5")]
  src.merge.tracing.27ss.all.pan = src.merge.tracing.27ss.all[,src.merge.tracing.27ss.all$cluster.v06.26.re_correct%in%c("DP","VP","EP.1","EP.2")]
  src.merge.tracing.27ss.all.smai = src.merge.tracing.27ss.all[,src.merge.tracing.27ss.all$cluster.v06.26.re_correct%in%c("Small.intestine.1","Small.intestine.2")]
  src.merge.tracing.27ss.all.lari = src.merge.tracing.27ss.all[,src.merge.tracing.27ss.all$cluster.v06.26.re_correct%in%c("Large.intestine.1","Large.intestine.2","Large.intestine.3")]
  
  for(i.name in c("pha","pan","smai","lari")){
    src.temp = get(paste("src.merge.tracing.27ss.all.", i.name, sep=""))
    src.temp$cluster.temp = src.temp$cluster.v06.26.re_correct
    if("EP.1" %in% unique(src.temp$cluster.temp)){
      src.temp@meta.data[src.temp$cluster.v06.26.re_correct%in%c("EP.1","EP.2"),]$cluster.temp = "EP"
    }
    
    if(i.name == "pan"){
      src.temp = FindVariableFeatures(src.temp, nfeatures = 3000)
      src.temp = ScaleData(src.temp, features = rownames(src.temp))
      
      src.temp.selectgene =
        Myfilter(as.matrix(src.temp@assays$RNA@data),
                 gene = src.temp@assays$RNA@var.features,
                 pearson.threshold = 0.2, partner.threshold = 5)
    }else{
      src.temp = FindVariableFeatures(src.temp, nfeatures = 2000)
      src.temp = ScaleData(src.temp, features = rownames(src.temp))
      
      src.temp.selectgene =
        Myfilter(as.matrix(src.temp@assays$RNA@data),
                 gene = src.temp@assays$RNA@var.features,
                 pearson.threshold = 0.2, partner.threshold = 5)
    }
    
    for(i.tree in c("row")){
      pdf(paste("GCN_check/try.src.merge.tracing.27ss.all.", i.name, i.tree,".pdf", sep = ''), 9, 7)
      
      src = src.temp
      gene_src = src.temp.selectgene
      cellorder_src = colnames(src)
      
      src.tree = MyHeatmap(
        as.matrix(src@assays$RNA@data[gene_src, ]),
        type = "row.relat",
        hc.c.data.type = "row.relat",
        hc.r.data.type = "row.relat",
        c.cov.method = "s",
        r.cov.method = "s",
        ColSideColors = cbind(
          MyName2Col(src$Time[cellorder_src], colors.time.2),
          MyName2Col(src$lineage[cellorder_src], color.lineage),
          MyName2Col(src$batch[cellorder_src], colors.geneset),
          MyName2Col(src$cluster.v06.26.re_correct[cellorder_src], cluster.endoderm.color.v5)),
        c.hc.method = "ward.D",
        r.hc.method = "ward.D2",
        ColSideColorsSize = 3,
        return.tree = i.tree,
        graph = T)
      
      assign(paste("tree_src.merge.tracing.27ss.all.", i.name, ".", i.tree, "tree", sep = ''),
             as.dendrogram(src.tree))
      dev.off()
    }
  }
  
  gene.src.merge.tracing.27ss.all.list = list()
  gene.src.merge.tracing.27ss.all.list$pha = c(
    labels(tree_src.merge.tracing.27ss.all.pha.rowtree[[1]][[2]][[1]]),
    labels(tree_src.merge.tracing.27ss.all.pha.rowtree[[1]][[1]][[2]]),
    labels(tree_src.merge.tracing.27ss.all.pha.rowtree[[1]][[2]][[2]][[1]]),
    labels(tree_src.merge.tracing.27ss.all.pha.rowtree[[2]][[2]][[2]]),
    labels(tree_src.merge.tracing.27ss.all.pha.rowtree[[2]][[1]][[2]][[1]]),
    labels(tree_src.merge.tracing.27ss.all.pha.rowtree[[2]][[1]][[1]][[1]]))
  names(gene.src.merge.tracing.27ss.all.list$pha) = c(
    rep("Pharynx.organ.2", length(c(labels(tree_src.merge.tracing.27ss.all.pha.rowtree[[1]][[1]][[2]]),
                                    labels(tree_src.merge.tracing.27ss.all.pha.rowtree[[1]][[2]][[1]]),
                                    labels(tree_src.merge.tracing.27ss.all.pha.rowtree[[1]][[2]][[2]][[1]])))),
    rep("Pharynx.organ.1", length(c(labels(tree_src.merge.tracing.27ss.all.pha.rowtree[[2]][[2]][[2]])))),
    rep("Pharynx.organ.4", length(c(labels(tree_src.merge.tracing.27ss.all.pha.rowtree[[2]][[1]][[2]][[1]])))),
    rep("Pharynx.organ.3", length(c(labels(tree_src.merge.tracing.27ss.all.pha.rowtree[[2]][[1]][[1]][[1]])))))
  
  # -- set 3000 for more-genes
  gene.src.merge.tracing.27ss.all.list$pan = c(
    labels(tree_src.merge.tracing.27ss.all.pan.rowtree[[2]][[2]][[2]][[1]][[2]]),
    labels(tree_src.merge.tracing.27ss.all.pan.rowtree[[2]][[2]][[2]][[2]][[1]]),
    labels(tree_src.merge.tracing.27ss.all.pan.rowtree[[2]][[2]][[2]][[2]][[2]][[2]][[2]]),
    rev(c(
      labels(tree_src.merge.tracing.27ss.all.pan.rowtree[[2]][[1]][[2]][[2]][[2]][[2]]),
      rev(labels(tree_src.merge.tracing.27ss.all.pan.rowtree[[1]][[2]])),
      labels(tree_src.merge.tracing.27ss.all.pan.rowtree[[2]][[1]][[2]][[2]][[2]][[1]]),
      labels(tree_src.merge.tracing.27ss.all.pan.rowtree[[2]][[1]][[2]][[2]][[1]]),
      rev(labels(tree_src.merge.tracing.27ss.all.pan.rowtree[[1]][[1]]))
    )))
  names(gene.src.merge.tracing.27ss.all.list$pan) = c(
    rep("DP", length(c(labels(tree_src.merge.tracing.27ss.all.pan.rowtree[[2]][[2]][[2]][[2]][[1]]),
                       labels(tree_src.merge.tracing.27ss.all.pan.rowtree[[2]][[2]][[2]][[1]][[2]])))),
    rep("VP", length(c(labels(tree_src.merge.tracing.27ss.all.pan.rowtree[[2]][[2]][[2]][[2]][[2]][[2]][[2]])))),
    rep("EP", length(c(labels(tree_src.merge.tracing.27ss.all.pan.rowtree[[2]][[1]][[2]][[2]][[2]][[2]]),
                       rev(labels(tree_src.merge.tracing.27ss.all.pan.rowtree[[1]][[2]])),
                       labels(tree_src.merge.tracing.27ss.all.pan.rowtree[[2]][[1]][[2]][[2]][[2]][[1]]),
                       labels(tree_src.merge.tracing.27ss.all.pan.rowtree[[2]][[1]][[2]][[2]][[1]]),
                       rev(labels(tree_src.merge.tracing.27ss.all.pan.rowtree[[1]][[1]]))))))
  
  gene.src.merge.tracing.27ss.all.list$smai = c(
    setdiff(labels(tree_src.merge.tracing.27ss.all.smai.rowtree[[2]][[2]][[2]][[1]]),
            labels(tree_src.merge.tracing.27ss.all.smai.rowtree[[2]][[2]][[2]][[1]][[1]][[1]])),
    labels(tree_src.merge.tracing.27ss.all.smai.rowtree[[2]][[1]][[2]]))
  names(gene.src.merge.tracing.27ss.all.list$smai) = c(
    rep("Small.intestine.1", length(c(setdiff(labels(tree_src.merge.tracing.27ss.all.smai.rowtree[[2]][[2]][[2]][[1]]),
                                              labels(tree_src.merge.tracing.27ss.all.smai.rowtree[[2]][[2]][[2]][[1]][[1]][[1]]))))),
    rep("Small.intestine.2", length(c(labels(tree_src.merge.tracing.27ss.all.smai.rowtree[[2]][[1]][[2]])))))
  
  gene.src.merge.tracing.27ss.all.list$lari = c(
    labels(tree_src.merge.tracing.27ss.all.lari.rowtree[[2]][[1]][[2]][[2]]),
    labels(tree_src.merge.tracing.27ss.all.lari.rowtree[[2]][[1]][[1]][[2]]),
    setdiff(labels(tree_src.merge.tracing.27ss.all.lari.rowtree[[2]][[2]]),
            labels(tree_src.merge.tracing.27ss.all.lari.rowtree[[2]][[2]][[2]][[2]][[2]])),
    labels(tree_src.merge.tracing.27ss.all.lari.rowtree[[1]][[2]][[1]][[1]]),
    labels(tree_src.merge.tracing.27ss.all.lari.rowtree[[1]][[2]][[2]][[2]][[1]][[2]]))
  names(gene.src.merge.tracing.27ss.all.list$lari) = c(
    rep("Large.intestine.1", length(c(labels(tree_src.merge.tracing.27ss.all.lari.rowtree[[2]][[1]][[2]][[2]]),
                                      labels(tree_src.merge.tracing.27ss.all.lari.rowtree[[2]][[1]][[1]][[2]])))),
    rep("Large.intestine.2", length(c(setdiff(labels(tree_src.merge.tracing.27ss.all.lari.rowtree[[2]][[2]]),
                                              labels(tree_src.merge.tracing.27ss.all.lari.rowtree[[2]][[2]][[2]][[2]][[2]]))))),
    rep("Large.intestine.3", length(c(labels(tree_src.merge.tracing.27ss.all.lari.rowtree[[1]][[2]][[1]][[1]]),
                                      labels(tree_src.merge.tracing.27ss.all.lari.rowtree[[1]][[2]][[2]][[2]][[1]][[2]])))))
  
  for(i.name in c("pha","pan","smai","lari")){
    seurat = get(paste("src.merge.tracing.27ss.all.", i.name, sep=""))
    i.genelist = gene.src.merge.tracing.27ss.all.list[[i.name]]
    rev.i.genelist = names(i.genelist)
    names(rev.i.genelist) = i.genelist 
    
    pdf(paste("GCN_check/graph.src.merge.tracing.27ss.all.", i.name, ".pdf", sep = ""), 6, 6)
    temp = WGCNA::cor(t(as.matrix(seurat@assays$RNA@data[i.genelist,])))
    if(i.name %in% c("pha","lar1")){
      temp[temp<0.25]=0
    }else if(i.name %in% c("smai")){
      temp[temp<0.15]=0  
    }else{
      temp[temp<0.3]=0
    }
    
    temp = igraph::graph.adjacency(temp, mode = "undirected", weighted = T)
    temp = igraph::simplify(temp, remove.multiple = T, remove.loops = T)
    temp = igraph::delete.vertices(temp, v =  igraph::V(temp)[
      igraph::clusters(temp)$membership%in%which(igraph::clusters(temp)$csize<=10)])
    set.seed(10)
    graph.temp = temp
    layout.temp = igraph::layout_with_fr(graph.temp)
    igraph::plot.igraph(
      graph.temp,
      layout=layout.temp,
      edge.color="#CCCCCC66",
      vertex.label.color='#FFFFFF00',
      vertex.label.cex=1,
      vertex.frame.color="#FFFFFF00",
      vertex.frame.width=0.5,
      vertex.size=6,
      vertex.color=cluster.endoderm.color.v5[rev.i.genelist[names(igraph::V(graph.temp))]])
    
    dev.off()
  }
  
  cell.src.merge.tracing.27ss.all.list = list()
  for(i.name in c("pha","pan","smai","lari")){
    
    if(i.name=="pha"){type.list = c("Pharynx.organ.2","Pharynx.organ.1","Pharynx.organ.4","Pharynx.organ.5","Pharynx.organ.3")}
    if(i.name=="pan"){type.list = c("DP","VP","EP.1","EP.2")}
    if(i.name=="smai"){type.list = c("Small.intestine.1","Small.intestine.2")}
    if(i.name=="lari"){type.list = c("Large.intestine.1","Large.intestine.2","Large.intestine.3")}
    cell.temp = c()
    for(i.temp in type.list){
      if(i.temp == "EP.2"){next()}
      
      for(i.lineage in factor(list_type_27ss_lineage[[i.temp]], levels = list_type_lineage)){
        
        if(i.temp == "EP.1"){
          cell.list = rownames(src.merge.tracing.27ss.all@meta.data[
            src.merge.tracing.27ss.all$cluster.v06.26.re_correct%in%c("EP.1","EP.2")&
              src.merge.tracing.27ss.all$lineage%in%i.lineage,])
        }else{
          cell.list = rownames(src.merge.tracing.27ss.all@meta.data[
            src.merge.tracing.27ss.all$cluster.v06.26.re_correct%in%c(i.temp)&
              src.merge.tracing.27ss.all$lineage%in%i.lineage,])
        }
        
        cell.temp = c(cell.temp, sample(cell.list, length(cell.list)))
      }
    }
    cell.src.merge.tracing.27ss.all.list[[i.name]] = cell.temp
  }
  
  for(i.name in c("pha","pan","smai","lari")){
    for(i.tree in c("col")){
      pdf(paste("GCN_check/try.src.merge.tracing.27ss.all.", i.name, i.tree,".pdf", sep = ''), 7, 7)
      
      src = get(paste("src.merge.tracing.27ss.all.", i.name, sep=""))
      gene_src = gene.src.merge.tracing.27ss.all.list[[i.name]]
      cellorder_src = cell.src.merge.tracing.27ss.all.list[[i.name]]
      
      if(i.name == "pan"){
        src$cluster.v06.26.re_correct = gsub("EP.1","EP",gsub("EP.2","EP",(src$cluster.v06.26.re_correct)))
      }
      
      src.tree = MyHeatmap(
        as.matrix(src@assays$RNA@data[gene_src,cellorder_src]),
        type = "row.relat",
        hc.c.data.type = "row.relat",
        hc.r.data.type = "row.relat",
        c.cov.method = "s",
        r.cov.method = "s",
        ColSideColors = cbind(
          MyName2Col(src$lineage[cellorder_src], color.lineage.re),
          MyName2Col(src$cluster.v06.26.re_correct[cellorder_src], cluster.endoderm.color.v5)),
        RowSideColors = t(cbind(
          MyName2Col(names(gene_src), cluster.endoderm.color.v5))),
        c.hc.method = "ward.D",
        r.hc.method = "ward.D2",
        ColSideColorsSize = 3,
        Rowv = "none",
        Colv = "none",
        margins = c(10,10),
        # return.tree = i.tree,
        graph = T)
      dev.off()
    }
  }
  
  
}

#>>> Figure S4H
for(i.sub_cluster in c("GCN for sub-cluster")){
  gene.src.merge.tracing.27ss.all.list.sub = list()
  
  for(i.organ in c("Large intestine")){
    src.merge.tracing.27ss.all.lari.1 = src.merge.tracing.27ss.all.lari[,src.merge.tracing.27ss.all.lari$cluster.v06.26.re_correct%in%"Large.intestine.1"]
    src.merge.tracing.27ss.all.lari.2 = src.merge.tracing.27ss.all.lari[,src.merge.tracing.27ss.all.lari$cluster.v06.26.re_correct%in%"Large.intestine.2"]
    src.merge.tracing.27ss.all.lari.3 = src.merge.tracing.27ss.all.lari[,src.merge.tracing.27ss.all.lari$cluster.v06.26.re_correct%in%"Large.intestine.3"]
    
    for(i.name in c("lari.1","lari.2","lari.3")){
      src.temp = get(paste("src.merge.tracing.27ss.all.", i.name, sep=""))
      src.temp$cluster.temp = src.temp$cluster.v06.26.re_correct
      src.temp = FindVariableFeatures(src.temp, nfeatures = 2000)
      src.temp = ScaleData(src.temp, features = rownames(src.temp))
      src.temp.selectgene =
        Myfilter(as.matrix(src.temp@assays$RNA@data),
                 gene = src.temp@assays$RNA@var.features,
                 bottom.dispersion.interval = 0.1,
                 pearson.threshold = 0.2, partner.threshold = 5)
      src.temp.selectgene = intersect(src.temp.selectgene,
                                      gi[gi$GeneType%in%c("protein_coding"),]$Symbol2)
      for(i.tree in c("row")){
        pdf(paste("GCN_check/try.src.merge.tracing.27ss.all.", i.name, i.tree,".pdf", sep = ''), 9, 7)
        
        src = src.temp
        gene_src = src.temp.selectgene
        cellorder_src = colnames(src)
        
        src.tree = MyHeatmap(
          as.matrix(src@assays$RNA@data[gene_src, ]),
          type = "row.relat",
          hc.c.data.type = "row.relat",
          hc.r.data.type = "row.relat",
          c.cov.method = "s",
          r.cov.method = "s",
          ColSideColors = cbind(
            MyName2Col(src$Time[cellorder_src], colors.time.2),
            MyName2Col(src$lineage[cellorder_src], color.lineage),
            MyName2Col(src$batch[cellorder_src], colors.geneset),
            MyName2Col(src$cluster.v06.26.re_correct[cellorder_src], cluster.endoderm.color.v5)),
          c.hc.method = "ward.D",
          r.hc.method = "ward.D2",
          ColSideColorsSize = 3,
          return.tree = i.tree,
          graph = T)
        
        assign(paste("tree_src.merge.tracing.27ss.all.", i.name, ".", i.tree, "tree", sep = ''),
               as.dendrogram(src.tree))
        dev.off()
      }
      gene.src.merge.tracing.27ss.all.list.sub[[i.name]] = src.temp.selectgene
    }
    
    gene.src.merge.tracing.27ss.all.list.sub$lari.1 = setdiff(
      c(labels(tree_src.merge.tracing.27ss.all.lari.1.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.lari.1.rowtree[[2]])),
      c(labels(tree_src.merge.tracing.27ss.all.lari.1.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.lari.1.rowtree[[2]][[1]][[1]])))
    
    gene.src.merge.tracing.27ss.all.list.sub$lari.2 = setdiff(
      c(labels(tree_src.merge.tracing.27ss.all.lari.2.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.lari.2.rowtree[[2]])),
      c(labels(tree_src.merge.tracing.27ss.all.lari.2.rowtree[[1]])))
    
    gene.src.merge.tracing.27ss.all.list.sub$lari.3 = setdiff(
      c(labels(tree_src.merge.tracing.27ss.all.lari.3.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.lari.3.rowtree[[2]])),
      c(labels(tree_src.merge.tracing.27ss.all.lari.3.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.lari.3.rowtree[[2]][[1]][[1]]),
        labels(tree_src.merge.tracing.27ss.all.lari.3.rowtree[[2]][[2]][[2]][[1]])))
    
    for(i.name in c("lari.1","lari.2","lari.3")){
      seurat = get(paste("src.merge.tracing.27ss.all.", i.name, sep=""))
      i.genelist = gene.src.merge.tracing.27ss.all.list.sub[[i.name]]
      # rev.i.genelist = names(i.genelist)
      # names(rev.i.genelist) = i.genelist 
      
      pdf(paste("GCN_check/graph.src.merge.tracing.27ss.all.", i.name, ".pdf", sep = ""), 6, 6)
      
      temp = WGCNA::cor(t(as.matrix(seurat@assays$RNA@data[i.genelist,])))
      temp[temp<0.25]=0
      temp = igraph::graph.adjacency(temp, mode = "undirected", weighted = T)
      temp = igraph::simplify(temp, remove.multiple = T, remove.loops = T)
      temp = igraph::delete.vertices(temp, v =  igraph::V(temp)[
        igraph::clusters(temp)$membership%in%which(igraph::clusters(temp)$csize<=10)])
      set.seed(10)
      graph.temp = temp
      layout.temp = igraph::layout_with_fr(graph.temp)
      igraph::plot.igraph(
        graph.temp,
        layout=layout.temp,
        edge.color="#CCCCCC66",
        vertex.label.color='#FFFFFF00',
        vertex.label.cex=1,
        vertex.frame.color="#FFFFFF00",
        vertex.frame.width=0.5,
        vertex.size=4.5,
        vertex.alpha=0.5,
        vertex.color=cluster.endoderm.color.v5[list.endoderm.name[i.name]])
      dev.off()
    }
  }
  
  for(i.organ in c("Small intestine")){
    src.merge.tracing.27ss.all.smai.1 = src.merge.tracing.27ss.all.smai[,src.merge.tracing.27ss.all.smai$cluster.v06.26.re_correct%in%"Small.intestine.1"]
    src.merge.tracing.27ss.all.smai.2 = src.merge.tracing.27ss.all.smai[,src.merge.tracing.27ss.all.smai$cluster.v06.26.re_correct%in%"Small.intestine.2"]
    
    for(i.name in c("smai.1","smai.2")){
      src.temp = get(paste("src.merge.tracing.27ss.all.", i.name, sep=""))
      src.temp$cluster.temp = src.temp$cluster.v06.26.re_correct
      src.temp = FindVariableFeatures(src.temp, nfeatures = 2000)
      src.temp = ScaleData(src.temp, features = rownames(src.temp))
      src.temp.selectgene =
        Myfilter(as.matrix(src.temp@assays$RNA@data),
                 gene = src.temp@assays$RNA@var.features,
                 bottom.dispersion.interval = 0.1,
                 pearson.threshold = 0.2, partner.threshold = 5)
      src.temp.selectgene = intersect(src.temp.selectgene,
                                      gi[gi$GeneType%in%c("protein_coding"),]$Symbol2)
      for(i.tree in c("row")){
        pdf(paste("GCN_check/try.src.merge.tracing.27ss.all.", i.name, i.tree,".pdf", sep = ''), 9, 7)
        
        src = src.temp
        gene_src = src.temp.selectgene
        cellorder_src = colnames(src)
        
        src.tree = MyHeatmap(
          as.matrix(src@assays$RNA@data[gene_src, ]),
          type = "row.relat",
          hc.c.data.type = "row.relat",
          hc.r.data.type = "row.relat",
          c.cov.method = "s",
          r.cov.method = "s",
          ColSideColors = cbind(
            MyName2Col(src$Time[cellorder_src], colors.time.2),
            MyName2Col(src$lineage[cellorder_src], color.lineage),
            MyName2Col(src$batch[cellorder_src], colors.geneset),
            MyName2Col(src$cluster.v06.26.re_correct[cellorder_src], cluster.endoderm.color.v5)),
          c.hc.method = "ward.D",
          r.hc.method = "ward.D2",
          ColSideColorsSize = 3,
          return.tree = i.tree,
          graph = T)
        
        assign(paste("tree_src.merge.tracing.27ss.all.", i.name, ".", i.tree, "tree", sep = ''),
               as.dendrogram(src.tree))
        dev.off()
      }
      gene.src.merge.tracing.27ss.all.list.sub[[i.name]] = src.temp.selectgene
    }
    
    gene.src.merge.tracing.27ss.all.list.sub$smai.1 = setdiff(
      c(labels(tree_src.merge.tracing.27ss.all.smai.1.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.smai.1.rowtree[[2]])),
      c(labels(tree_src.merge.tracing.27ss.all.smai.1.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.smai.1.rowtree[[2]][[1]])))
    
    gene.src.merge.tracing.27ss.all.list.sub$smai.2 = setdiff(
      c(labels(tree_src.merge.tracing.27ss.all.smai.2.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.smai.2.rowtree[[2]])),
      c(labels(tree_src.merge.tracing.27ss.all.smai.2.rowtree[[1]])))
    
    for(i.name in c("smai.1","smai.2")){
      seurat = get(paste("src.merge.tracing.27ss.all.", i.name, sep=""))
      i.genelist = gene.src.merge.tracing.27ss.all.list.sub[[i.name]]
      # rev.i.genelist = names(i.genelist)
      # names(rev.i.genelist) = i.genelist 
      
      pdf(paste("GCN_check/graph.src.merge.tracing.27ss.all.", i.name, ".pdf", sep = ""), 6, 6)
      
      temp = WGCNA::cor(t(as.matrix(seurat@assays$RNA@data[i.genelist,])))
      temp[temp<0.25]=0
      temp = igraph::graph.adjacency(temp, mode = "undirected", weighted = T)
      temp = igraph::simplify(temp, remove.multiple = T, remove.loops = T)
      temp = igraph::delete.vertices(temp, v =  igraph::V(temp)[
        igraph::clusters(temp)$membership%in%which(igraph::clusters(temp)$csize<=10)])
      set.seed(10)
      graph.temp = temp
      layout.temp = igraph::layout_with_fr(graph.temp)
      igraph::plot.igraph(
        graph.temp,
        layout=layout.temp,
        edge.color="#CCCCCC66",
        vertex.label.color='#FFFFFF00',
        vertex.label.cex=1,
        vertex.frame.color="#FFFFFF00",
        vertex.frame.width=0.5,
        vertex.size=4.5,
        vertex.alpha=0.5,
        vertex.color=cluster.endoderm.color.v5[list.endoderm.name[i.name]])
      dev.off()
    }
  }
  
  for(i.organ in c("Pancreas")){
    src.merge.tracing.27ss.all.pan.1 = src.merge.tracing.27ss.all.pan[,src.merge.tracing.27ss.all.pan$cluster.v06.26.re_correct%in%"DP"]
    src.merge.tracing.27ss.all.pan.2 = src.merge.tracing.27ss.all.pan[,src.merge.tracing.27ss.all.pan$cluster.v06.26.re_correct%in%"VP"]
    src.merge.tracing.27ss.all.pan.3 = src.merge.tracing.27ss.all.pan[,src.merge.tracing.27ss.all.pan$cluster.v06.26.re_correct%in%c("EP.1","EP.2")]
    src.merge.tracing.27ss.all.pan.4 = src.merge.tracing.27ss.all.pan[,src.merge.tracing.27ss.all.pan$cluster.v06.26.re_correct%in%c("EP.1")]
    src.merge.tracing.27ss.all.pan.5 = src.merge.tracing.27ss.all.pan[,src.merge.tracing.27ss.all.pan$cluster.v06.26.re_correct%in%c("EP.2")]
    
    for(i.name in c("pan.1","pan.2","pan.3","pan.4","pan.5")){
      src.temp = get(paste("src.merge.tracing.27ss.all.", i.name, sep=""))
      src.temp$cluster.temp = src.temp$cluster.v06.26.re_correct
      src.temp = FindVariableFeatures(src.temp, nfeatures = 2000)
      src.temp = ScaleData(src.temp, features = rownames(src.temp))
      src.temp.selectgene =
        Myfilter(as.matrix(src.temp@assays$RNA@data),
                 gene = src.temp@assays$RNA@var.features,
                 bottom.dispersion.interval = 0.1,
                 pearson.threshold = 0.2, partner.threshold = 5)
      src.temp.selectgene = intersect(src.temp.selectgene,
                                      gi[gi$GeneType%in%c("protein_coding"),]$Symbol2)
      for(i.tree in c("row")){
        pdf(paste("GCN_check/try.src.merge.tracing.27ss.all.", i.name, i.tree,".pdf", sep = ''), 9, 7)
        
        src = src.temp
        gene_src = src.temp.selectgene
        cellorder_src = colnames(src)
        
        src.tree = MyHeatmap(
          as.matrix(src@assays$RNA@data[gene_src, ]),
          type = "row.relat",
          hc.c.data.type = "row.relat",
          hc.r.data.type = "row.relat",
          c.cov.method = "s",
          r.cov.method = "s",
          ColSideColors = cbind(
            MyName2Col(src$Time[cellorder_src], colors.time.2),
            MyName2Col(src$lineage[cellorder_src], color.lineage),
            MyName2Col(src$batch[cellorder_src], colors.geneset),
            MyName2Col(src$cluster.v06.26.re_correct[cellorder_src], cluster.endoderm.color.v5)),
          c.hc.method = "ward.D",
          r.hc.method = "ward.D2",
          ColSideColorsSize = 3,
          return.tree = i.tree,
          graph = T)
        
        assign(paste("tree_src.merge.tracing.27ss.all.", i.name, ".", i.tree, "tree", sep = ''),
               as.dendrogram(src.tree))
        dev.off()
      }
      gene.src.merge.tracing.27ss.all.list.sub[[i.name]] = src.temp.selectgene
    }
    
    gene.src.merge.tracing.27ss.all.list.sub$pan.1 = setdiff(
      c(labels(tree_src.merge.tracing.27ss.all.pan.1.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.pan.1.rowtree[[2]])),
      c(labels(tree_src.merge.tracing.27ss.all.pan.1.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.pan.1.rowtree[[2]][[1]][[1]])))
    
    gene.src.merge.tracing.27ss.all.list.sub$pan.2 = setdiff(
      c(labels(tree_src.merge.tracing.27ss.all.pan.2.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.pan.2.rowtree[[2]])),
      c())
    
    gene.src.merge.tracing.27ss.all.list.sub$pan.3 = 
      c(labels(tree_src.merge.tracing.27ss.all.pan.3.rowtree[[2]][[1]]),
        labels(tree_src.merge.tracing.27ss.all.pan.3.rowtree[[2]][[2]][[2]][[2]]),
        labels(tree_src.merge.tracing.27ss.all.pan.3.rowtree[[1]]))
    names(gene.src.merge.tracing.27ss.all.list.sub$pan.3) = c(
      rep("EP.1", length(c(labels(tree_src.merge.tracing.27ss.all.pan.3.rowtree[[2]][[1]])))),
      rep("EP.1", length(c(labels(tree_src.merge.tracing.27ss.all.pan.3.rowtree[[2]][[2]][[2]][[2]])))),
      rep("EP.2", length(c(labels(tree_src.merge.tracing.27ss.all.pan.3.rowtree[[1]])))))
    
    gene.src.merge.tracing.27ss.all.list.sub$pan.4 = setdiff(
      c(labels(tree_src.merge.tracing.27ss.all.pan.4.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.pan.4.rowtree[[2]])),
      c())
    
    gene.src.merge.tracing.27ss.all.list.sub$pan.5 = setdiff(
      c(labels(tree_src.merge.tracing.27ss.all.pan.5.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.pan.5.rowtree[[2]])),
      c())
    
    for(i.name in c("pan.1","pan.2","pan.4","pan.5")){
      seurat = get(paste("src.merge.tracing.27ss.all.", i.name, sep=""))
      i.genelist = gene.src.merge.tracing.27ss.all.list.sub[[i.name]]
      # rev.i.genelist = names(i.genelist)
      # names(rev.i.genelist) = i.genelist 
      
      pdf(paste("GCN_check/graph.src.merge.tracing.27ss.all.", i.name, ".pdf", sep = ""), 6, 6)
      
      temp = WGCNA::cor(t(as.matrix(seurat@assays$RNA@data[i.genelist,])))
      temp[temp<0.25]=0
      temp = igraph::graph.adjacency(temp, mode = "undirected", weighted = T)
      temp = igraph::simplify(temp, remove.multiple = T, remove.loops = T)
      temp = igraph::delete.vertices(temp, v =  igraph::V(temp)[
        igraph::clusters(temp)$membership%in%which(igraph::clusters(temp)$csize<=10)])
      set.seed(10)
      graph.temp = temp
      layout.temp = igraph::layout_with_fr(graph.temp)
      igraph::plot.igraph(
        graph.temp,
        layout=layout.temp,
        edge.color="#CCCCCC66",
        vertex.label.color='#FFFFFF00',
        vertex.label.cex=1,
        vertex.frame.color="#FFFFFF00",
        vertex.frame.width=0.5,
        vertex.size=4.5,
        vertex.alpha=0.5,
        vertex.color=cluster.endoderm.color.v5[list.endoderm.name[i.name]])
      dev.off()
    }
    
    for(i.name in c("pan.3")){
      seurat = get(paste("src.merge.tracing.27ss.all.", i.name, sep=""))
      i.genelist = gene.src.merge.tracing.27ss.all.list.sub[[i.name]]
      rev.i.genelist = names(i.genelist)
      names(rev.i.genelist) = i.genelist 
      
      pdf(paste("GCN_check/graph.src.merge.tracing.27ss.all.", i.name, ".pdf", sep = ""), 6, 6)
      
      temp = WGCNA::cor(t(as.matrix(seurat@assays$RNA@data[i.genelist,])))
      temp[temp<0.25]=0
      temp = igraph::graph.adjacency(temp, mode = "undirected", weighted = T)
      temp = igraph::simplify(temp, remove.multiple = T, remove.loops = T)
      temp = igraph::delete.vertices(temp, v =  igraph::V(temp)[
        igraph::clusters(temp)$membership%in%which(igraph::clusters(temp)$csize<=10)])
      set.seed(10)
      graph.temp = temp
      layout.temp = igraph::layout_with_fr(graph.temp)
      igraph::plot.igraph(
        graph.temp,
        layout=layout.temp,
        edge.color="#CCCCCC66",
        vertex.label.color='#FFFFFF00',
        vertex.label.cex=1,
        vertex.frame.color="#FFFFFF00",
        vertex.frame.width=0.5,
        vertex.size=4.5,
        vertex.alpha=0.5,
        vertex.color=cluster.endoderm.color.v5[rev.i.genelist[names(igraph::V(graph.temp))]])
      
      graph.temp.cluster <- cluster_walktrap(graph.temp,steps = 10)
      gene.cluster <- graph.temp.cluster$membership
      names(gene.cluster) <- graph.temp.cluster$names
      igraph::plot.igraph(
        graph.temp,
        layout=layout.temp,
        edge.color="#CCCCCC66",
        vertex.label = gene.cluster[names(igraph::V(graph.temp))],
        vertex.label.cex=0.5,
        vertex.frame.color="#FFFFFF00",
        vertex.frame.width=0.5,
        vertex.size=4.5,
        vertex.alpha=0.5,
        vertex.color=colors.geneset[gene.cluster[names(igraph::V(graph.temp))]])
      dev.off()
      
      gene.src.merge.tracing.27ss.all.list.sub$pan.3.gcn = c(
        names(gene.cluster[gene.cluster%in%c(3,4,7,8)]),
        names(gene.cluster[gene.cluster%in%c(1,2,6,9)]))
      names(gene.src.merge.tracing.27ss.all.list.sub$pan.3.gcn) = c(
        rep("EP.1", length(c(names(gene.cluster[gene.cluster%in%c(3,4,7,8)])))),
        rep("EP.2", length(c(names(gene.cluster[gene.cluster%in%c(1,2,6,9)])))))
      
      i.genelist = gene.src.merge.tracing.27ss.all.list.sub$pan.3.gcn
      rev.i.genelist = names(i.genelist)
      names(rev.i.genelist) = i.genelist 
      
      pdf(paste("GCN_check/graph.src.merge.tracing.27ss.all.", i.name, ".pdf", sep = ""), 6, 6)
      
      temp = WGCNA::cor(t(as.matrix(seurat@assays$RNA@data[i.genelist,])))
      temp[temp<0.25]=0
      temp = igraph::graph.adjacency(temp, mode = "undirected", weighted = T)
      temp = igraph::simplify(temp, remove.multiple = T, remove.loops = T)
      temp = igraph::delete.vertices(temp, v =  igraph::V(temp)[
        igraph::clusters(temp)$membership%in%which(igraph::clusters(temp)$csize<=10)])
      set.seed(10)
      graph.temp = temp
      layout.temp = igraph::layout_with_fr(graph.temp)
      igraph::plot.igraph(
        graph.temp,
        layout=layout.temp,
        edge.color="#CCCCCC66",
        vertex.label.color='#FFFFFF00',
        vertex.label.cex=1,
        vertex.frame.color="#FFFFFF00",
        vertex.frame.width=0.5,
        vertex.size=4.5,
        vertex.alpha=0.5,
        vertex.color=cluster.endoderm.color.v5[rev.i.genelist[names(igraph::V(graph.temp))]])
      dev.off()
    }  
  }
  
  for(i.organ in c("Pharynx.organ")){
    src.merge.tracing.27ss.all.pha.1 = src.merge.tracing.27ss.all.pha[,src.merge.tracing.27ss.all.pha$cluster.v06.26.re_correct%in%"Pharynx.organ.1"]
    src.merge.tracing.27ss.all.pha.2 = src.merge.tracing.27ss.all.pha[,src.merge.tracing.27ss.all.pha$cluster.v06.26.re_correct%in%"Pharynx.organ.2"]
    src.merge.tracing.27ss.all.pha.3 = src.merge.tracing.27ss.all.pha[,src.merge.tracing.27ss.all.pha$cluster.v06.26.re_correct%in%"Pharynx.organ.3"]
    src.merge.tracing.27ss.all.pha.4 = src.merge.tracing.27ss.all.pha[,src.merge.tracing.27ss.all.pha$cluster.v06.26.re_correct%in%"Pharynx.organ.4"]
    src.merge.tracing.27ss.all.pha.5 = src.merge.tracing.27ss.all.pha[,src.merge.tracing.27ss.all.pha$cluster.v06.26.re_correct%in%"Pharynx.organ.5"]
    
    for(i.name in c("pha.1","pha.2","pha.3","pha.4","pha.5")){
      src.temp = get(paste("src.merge.tracing.27ss.all.", i.name, sep=""))
      src.temp$cluster.temp = src.temp$cluster.v06.26.re_correct
      src.temp = FindVariableFeatures(src.temp, nfeatures = 2000)
      if(i.name == "pha.1"){
        src.temp = FindVariableFeatures(src.temp, nfeatures = 4000)
      }
      
      src.temp = ScaleData(src.temp, features = rownames(src.temp))
      src.temp = CellCycleScoring(src.temp, s.features = s_genes, g2m.features = g2m_genes)
      
      src.temp.selectgene =
        Myfilter(as.matrix(src.temp@assays$RNA@data),
                 gene = src.temp@assays$RNA@var.features,
                 bottom.dispersion.interval = 0.1,
                 pearson.threshold = 0.2, partner.threshold = 5)
      src.temp.selectgene = intersect(src.temp.selectgene,
                                      gi[gi$GeneType%in%c("protein_coding"),]$Symbol2)
      for(i.tree in c("row")){
        pdf(paste("GCN_check/try.src.merge.tracing.27ss.all.", i.name, i.tree,".pdf", sep = ''), 9, 7)
        
        src = src.temp
        gene_src = src.temp.selectgene
        cellorder_src = colnames(src)
        
        src.tree = MyHeatmap(
          as.matrix(src@assays$RNA@data[gene_src, ]),
          type = "row.relat",
          hc.c.data.type = "row.relat",
          hc.r.data.type = "row.relat",
          c.cov.method = "s",
          r.cov.method = "s",
          ColSideColors = cbind(
            MyName2Col(src$Time[cellorder_src], colors.time.2),
            MyName2Col(src$lineage[cellorder_src], color.lineage),
            MyName2Col(src$Phase[cellorder_src], as.character(colors.num)),
            MyName2Col(src$SeqDate[cellorder_src], as.character(colors.num)),
            MyName2Col(src$batch[cellorder_src], colors.geneset),
            MyName2Col(src$cluster.v06.26.re_correct[cellorder_src], cluster.endoderm.color.v5)),
          c.hc.method = "ward.D",
          r.hc.method = "ward.D2",
          ColSideColorsSize = 3,
          return.tree = i.tree,
          graph = T)
        
        assign(paste("tree_src.merge.tracing.27ss.all.", i.name, ".", i.tree, "tree", sep = ''),
               as.dendrogram(src.tree))
        dev.off()
      }
      gene.src.merge.tracing.27ss.all.list.sub[[i.name]] = src.temp.selectgene
    }
    
    gene.src.merge.tracing.27ss.all.list.sub$pha.1 = setdiff(
      c(labels(tree_src.merge.tracing.27ss.all.pha.1.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.pha.1.rowtree[[2]])),
      c(setdiff(labels(tree_src.merge.tracing.27ss.all.pha.1.rowtree[[2]]),
                labels(tree_src.merge.tracing.27ss.all.pha.1.rowtree[[2]][[2]][[2]][[2]][[2]][[2]]))))
    
    gene.src.merge.tracing.27ss.all.list.sub$pha.2 = setdiff(
      c(labels(tree_src.merge.tracing.27ss.all.pha.2.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.pha.2.rowtree[[2]])),
      c(labels(tree_src.merge.tracing.27ss.all.pha.2.rowtree[[1]])))
    
    gene.src.merge.tracing.27ss.all.list.sub$pha.3 = setdiff(
      c(labels(tree_src.merge.tracing.27ss.all.pha.3.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.pha.3.rowtree[[2]])),
      c(labels(tree_src.merge.tracing.27ss.all.pha.3.rowtree[[1]])))
    
    
    gene.src.merge.tracing.27ss.all.list.sub$pha.4 = setdiff(
      c(labels(tree_src.merge.tracing.27ss.all.pha.4.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.pha.4.rowtree[[2]])),
      c())
    
    gene.src.merge.tracing.27ss.all.list.sub$pha.5 = setdiff(
      c(labels(tree_src.merge.tracing.27ss.all.pha.5.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.pha.5.rowtree[[2]])),
      c())
    
    for(i.name in c("pha.2","pha.3","pha.4","pha.5")){
      seurat = get(paste("src.merge.tracing.27ss.all.", i.name, sep=""))
      i.genelist = gene.src.merge.tracing.27ss.all.list.sub[[i.name]]
      # rev.i.genelist = names(i.genelist)
      # names(rev.i.genelist) = i.genelist 
      
      pdf(paste("GCN_check/graph.src.merge.tracing.27ss.all.", i.name, ".pdf", sep = ""), 6, 6)
      
      temp = WGCNA::cor(t(as.matrix(seurat@assays$RNA@data[i.genelist,])))
      temp[temp<0.25]=0
      temp = igraph::graph.adjacency(temp, mode = "undirected", weighted = T)
      temp = igraph::simplify(temp, remove.multiple = T, remove.loops = T)
      temp = igraph::delete.vertices(temp, v =  igraph::V(temp)[
        igraph::clusters(temp)$membership%in%which(igraph::clusters(temp)$csize<=10)])
      set.seed(10)
      graph.temp = temp
      layout.temp = igraph::layout_with_fr(graph.temp)
      igraph::plot.igraph(
        graph.temp,
        layout=layout.temp,
        edge.color="#CCCCCC66",
        vertex.label.color='#FFFFFF00',
        vertex.label.cex=1,
        vertex.frame.color="#FFFFFF00",
        vertex.frame.width=0.5,
        vertex.size=4.5,
        vertex.alpha=0.5,
        vertex.color=cluster.endoderm.color.v5[list.endoderm.name[i.name]])
      dev.off()
    }
    
    for(i.name in c("pha.1")){
      seurat = get(paste("src.merge.tracing.27ss.all.", i.name, sep=""))
      i.genelist = gene.src.merge.tracing.27ss.all.list.sub[[i.name]]
      # rev.i.genelist = names(i.genelist)
      # names(rev.i.genelist) = i.genelist 
      
      pdf(paste("GCN_check/graph.src.merge.tracing.27ss.all.", i.name, ".pdf", sep = ""), 6, 6)
      
      temp = WGCNA::cor(t(as.matrix(seurat@assays$RNA@data[i.genelist,])))
      temp[temp<0.25]=0
      temp = igraph::graph.adjacency(temp, mode = "undirected", weighted = T)
      temp = igraph::simplify(temp, remove.multiple = T, remove.loops = T)
      temp = igraph::delete.vertices(temp, v =  igraph::V(temp)[
        igraph::clusters(temp)$membership%in%which(igraph::clusters(temp)$csize<=10)])
      set.seed(10)
      graph.temp = temp
      layout.temp = igraph::layout_with_fr(graph.temp)
      igraph::plot.igraph(
        graph.temp,
        layout=layout.temp,
        edge.color="#CCCCCC66",
        vertex.label.color='#FFFFFF00',
        vertex.label.cex=1,
        vertex.frame.color="#FFFFFF00",
        vertex.frame.width=0.5,
        vertex.size=3.5,
        vertex.alpha=0.5,
        vertex.color=cluster.endoderm.color.v5[list.endoderm.name[i.name]])
      dev.off()
    }
  }
  
  for(i.organ in c("Other")){
    src.merge.tracing.27ss.all.other = src.merge.tracing.27ss.all[,src.merge.tracing.27ss.all$cluster.v06.26.re_correct%in%c(
      "Thyroid","Lung","Stomach","Esophagus","Liver","EHBD")]
    src.merge.tracing.27ss.all.other.1 = src.merge.tracing.27ss.all.other[,src.merge.tracing.27ss.all.other$cluster.v06.26.re_correct%in%"Thyroid"]
    src.merge.tracing.27ss.all.other.2 = src.merge.tracing.27ss.all.other[,src.merge.tracing.27ss.all.other$cluster.v06.26.re_correct%in%"Lung"]
    src.merge.tracing.27ss.all.other.3 = src.merge.tracing.27ss.all.other[,src.merge.tracing.27ss.all.other$cluster.v06.26.re_correct%in%"Stomach"]
    src.merge.tracing.27ss.all.other.4 = src.merge.tracing.27ss.all.other[,src.merge.tracing.27ss.all.other$cluster.v06.26.re_correct%in%"Esophagus"]
    src.merge.tracing.27ss.all.other.5 = src.merge.tracing.27ss.all.other[,src.merge.tracing.27ss.all.other$cluster.v06.26.re_correct%in%"Liver"]
    src.merge.tracing.27ss.all.other.6 = src.merge.tracing.27ss.all.other[,src.merge.tracing.27ss.all.other$cluster.v06.26.re_correct%in%"EHBD"]
    
    for(i.name in c("other.1","other.2","other.3","other.4","other.5","other.6")){
      src.temp = get(paste("src.merge.tracing.27ss.all.", i.name, sep=""))
      src.temp$cluster.temp = src.temp$cluster.v06.26.re_correct
      src.temp = FindVariableFeatures(src.temp, nfeatures = 2000)
      if(i.name == "other.2"){
        src.temp = FindVariableFeatures(src.temp, nfeatures = 4000)
      }
      
      src.temp = ScaleData(src.temp, features = rownames(src.temp))
      src.temp = CellCycleScoring(src.temp, s.features = s_genes, g2m.features = g2m_genes)
      src.temp.selectgene =
        Myfilter(as.matrix(src.temp@assays$RNA@data),
                 gene = src.temp@assays$RNA@var.features,
                 bottom.dispersion.interval = 0.1,
                 pearson.threshold = 0.2, partner.threshold = 5)
      src.temp.selectgene = intersect(src.temp.selectgene,
                                      gi[gi$GeneType%in%c("protein_coding"),]$Symbol2)
      for(i.tree in c("row")){
        pdf(paste("GCN_check/try.src.merge.tracing.27ss.all.", i.name, i.tree,".pdf", sep = ''), 9, 7)
        
        src = src.temp
        gene_src = src.temp.selectgene
        cellorder_src = colnames(src)
        
        src.tree = MyHeatmap(
          as.matrix(src@assays$RNA@data[gene_src, ]),
          type = "row.relat",
          hc.c.data.type = "row.relat",
          hc.r.data.type = "row.relat",
          c.cov.method = "s",
          r.cov.method = "s",
          ColSideColors = cbind(
            MyName2Col(src$Time[cellorder_src], colors.time.2),
            MyName2Col(src$lineage[cellorder_src], color.lineage),
            MyName2Col(src$Phase[cellorder_src], as.character(colors.num)),
            MyName2Col(src$SeqDate[cellorder_src], as.character(colors.num)),
            MyName2Col(src$batch[cellorder_src], colors.geneset),
            MyName2Col(src$cluster.v06.26.re_correct[cellorder_src], cluster.endoderm.color.v5)),
          c.hc.method = "ward.D",
          r.hc.method = "ward.D2",
          ColSideColorsSize = 3,
          return.tree = i.tree,
          graph = T)
        
        assign(paste("tree_src.merge.tracing.27ss.all.", i.name, ".", i.tree, "tree", sep = ''),
               as.dendrogram(src.tree))
        dev.off()
      }
      gene.src.merge.tracing.27ss.all.list.sub[[i.name]] = src.temp.selectgene
    }
    
    gene.src.merge.tracing.27ss.all.list.sub$other.1 = setdiff(
      c(labels(tree_src.merge.tracing.27ss.all.other.1.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.other.1.rowtree[[2]])),
      c())
    
    gene.src.merge.tracing.27ss.all.list.sub$other.2 = setdiff(
      c(labels(tree_src.merge.tracing.27ss.all.other.2.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.other.2.rowtree[[2]])),
      c(labels(tree_src.merge.tracing.27ss.all.other.2.rowtree[[1]][[1]]),
        labels(tree_src.merge.tracing.27ss.all.other.2.rowtree[[2]][[1]]),
        labels(tree_src.merge.tracing.27ss.all.other.2.rowtree[[2]][[2]][[1]]),
        labels(tree_src.merge.tracing.27ss.all.other.2.rowtree[[2]][[2]][[2]][[2]])))
    
    gene.src.merge.tracing.27ss.all.list.sub$other.3 = setdiff(
      c(labels(tree_src.merge.tracing.27ss.all.other.3.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.other.3.rowtree[[2]])),
      c(labels(tree_src.merge.tracing.27ss.all.other.3.rowtree[[1]])))
    
    gene.src.merge.tracing.27ss.all.list.sub$other.4 = setdiff(
      c(labels(tree_src.merge.tracing.27ss.all.other.4.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.other.4.rowtree[[2]])),
      c(labels(tree_src.merge.tracing.27ss.all.other.4.rowtree[[2]][[1]]),
        labels(tree_src.merge.tracing.27ss.all.other.4.rowtree[[2]][[2]][[1]])))
    
    gene.src.merge.tracing.27ss.all.list.sub$other.5 = setdiff(
      c(labels(tree_src.merge.tracing.27ss.all.other.5.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.other.5.rowtree[[2]])),
      c())
    
    gene.src.merge.tracing.27ss.all.list.sub$other.6 = setdiff(
      c(labels(tree_src.merge.tracing.27ss.all.other.6.rowtree[[1]]),
        labels(tree_src.merge.tracing.27ss.all.other.6.rowtree[[2]])),
      c())
    
    for(i.name in c("other.1","other.2","other.3",
                    "other.4","other.5","other.6")){
      seurat = get(paste("src.merge.tracing.27ss.all.", i.name, sep=""))
      i.genelist = gene.src.merge.tracing.27ss.all.list.sub[[i.name]]
      # rev.i.genelist = names(i.genelist)
      # names(rev.i.genelist) = i.genelist 
      
      pdf(paste("GCN_check/graph.src.merge.tracing.27ss.all.", i.name, ".pdf", sep = ""), 6, 6)
      
      temp = WGCNA::cor(t(as.matrix(seurat@assays$RNA@data[i.genelist,])))
      temp[temp<0.25]=0
      temp = igraph::graph.adjacency(temp, mode = "undirected", weighted = T)
      temp = igraph::simplify(temp, remove.multiple = T, remove.loops = T)
      temp = igraph::delete.vertices(temp, v =  igraph::V(temp)[
        igraph::clusters(temp)$membership%in%which(igraph::clusters(temp)$csize<=10)])
      set.seed(10)
      graph.temp = temp
      layout.temp = igraph::layout_with_fr(graph.temp)
      igraph::plot.igraph(
        graph.temp,
        layout=layout.temp,
        edge.color="#CCCCCC66",
        vertex.label.color='#FFFFFF00',
        vertex.label.cex=1,
        vertex.frame.color="#FFFFFF00",
        vertex.frame.width=0.5,
        vertex.size=4.5,
        vertex.alpha=0.5,
        vertex.color=cluster.endoderm.color.v5[list.endoderm.name[i.name]])
      dev.off()
    }
  }
}

#>>> Figure S4I
for(i.cca in c("CCA")){
  src.merge.tracing.27ss.all = FindVariableFeatures(src.merge.tracing.27ss.all, nfeatures = 2000)
  src.merge.tracing.27ss.all = ScaleData(src.merge.tracing.27ss.all, features = rownames(src.merge.tracing.27ss.all))
  src.merge.tracing.27ss.all.selectgene =
    Myfilter(as.matrix(src.merge.tracing.27ss.all@assays$RNA@data),
             gene = src.merge.tracing.27ss.all@assays$RNA@var.features,
             pearson.threshold = 0.2, partner.threshold = 5)
  
  for(i.tree in c("col", "row")){
    pdf(paste("Milestones/try.src.merge.tracing.27ss.all." , i.tree,".pdf", sep = ''), 9, 7)
    
    src = src.merge.tracing.27ss.all
    gene_src = src.merge.tracing.27ss.all.selectgene
    cellorder_src = colnames(src)
    
    src.tree = MyHeatmap(
      as.matrix(src@assays$RNA@data[gene_src, ]),
      type = "row.relat",
      hc.c.data.type = "row.relat",
      hc.r.data.type = "row.relat",
      c.cov.method = "s",
      r.cov.method = "s",
      ColSideColors = cbind(
        MyName2Col(src$Time[cellorder_src], colors.time.2),
        MyName2Col(src$lineage[cellorder_src], color.lineage),
        MyName2Col(src$batch[cellorder_src], colors.geneset),
        MyName2Col(src$cluster.v06.26.re_correct[cellorder_src], cluster.endoderm.color.v5)),
      c.hc.method = "ward.D",
      r.hc.method = "ward.D2",
      ColSideColorsSize = 3,
      return.tree = i.tree,
      graph = T)
    
    assign(paste("tree_src.merge.tracing.27ss.all.", i.tree, "tree", sep = ''),
           as.dendrogram(src.tree))
    dev.off()
  }
  gene.src.merge.tracing.27ss.all.selectgene.rowtree = setdiff(
    src.merge.tracing.27ss.all.selectgene,
    c(labels(tree_src.merge.tracing.27ss.all.rowtree[[1]][[1]]),
      labels(tree_src.merge.tracing.27ss.all.rowtree[[2]][[1]][[2]][[1]]),
      c()))
  src.merge.tracing.27ss.all = RunPCA(src.merge.tracing.27ss.all, 
                                      features = gene.src.merge.tracing.27ss.all.selectgene.rowtree)
  
  #--- PCA-Cor
  cell_27ss_lineage_index = c()
  for(i.name in names(list_type_27ss_lineage)){
    cell_27ss_lineage_index = c(
      cell_27ss_lineage_index,
      rownames(src.merge.tracing.27ss.all@meta.data[
        src.merge.tracing.27ss.all$lineage%in%list_type_27ss_lineage[[i.name]] &
          src.merge.tracing.27ss.all$cluster.v06.26.re_correct%in%i.name,]))
  }
  
  src.merge.tracing.27ss.all$index = paste(src.merge.tracing.27ss.all$cluster.v06.26.re_correct,
                                           src.merge.tracing.27ss.all$lineage, sep = "_")
  
  cca.merge.27ss = by(as.data.frame(src.merge.tracing.27ss.all@reductions$pca@cell.embeddings[,1:20]),
                      INDICES = src.merge.tracing.27ss.all$index, FUN = colMeans)
  cca.merge.27ss = do.call(cbind, cca.merge.27ss)
  cor.cca.merge.27ss = cor(cca.merge.27ss)
  cor.cca.merge.27ss[cor.cca.merge.27ss<0]=0
  name_list1 = intersect(list_type_27ss_lineage_index, colnames(cor.cca.merge.27ss))
  name_list2 = intersect(list_type_27ss_lineage_index, colnames(cor.cca.merge.27ss))
  
  pdf("Milestones/cor.cca.merge.27ss.pdf", 5, 5)
  for(i.data in c("cor.cca.merge.27ss")){
    data = get(i.data)
    a = MyHeatmap(as.data.frame(data[name_list1, name_list2]),
                  type = "raw",
                  hc.c.data.type = "raw",
                  c.cov.method = "p",
                  c.hc.method = "ward.D",
                  RowSideColors = t(cbind(
                    MyName2Col(strsplit2(name_list1, split = "_")[,2], color.lineage.re),
                    MyName2Col(strsplit2(name_list1, split = "_")[,1], cluster.endoderm.color.v5))),
                  ColSideColors = cbind(
                    MyName2Col(strsplit2(name_list1, split = "_")[,1], cluster.endoderm.color.v5),
                    MyName2Col(strsplit2(name_list2, split = "_")[,2], color.lineage.re)),
                  color.palette = colorRampPalette(c(
                    # "#5aa7dd"
                    "#58a6dd","#EBEBEB","#fe0000"), space="Lab"),
                             Rowv = "none", Colv = "none", 
                  # return.tree = "col",
                  ColSideColorsSize = 1.5, RowSideColorsSize = 1.5,
                  labRow = name_list1,
                  labCol = name_list2, 
                  margins = c(8,8))
  }
  dev.off()
}

#===============================================================================


#===============================================================================
#>> 2.4 Endoderm intermediate state:: GCN
#> In [6.trajectory construction of endoderm]
#===============================================================================















