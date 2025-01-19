#----------------------------------
#>>> Plot-summary (for figure 2b)
#----------------------------------


# -- plot cell-type, time, lineage
#--------------------------------------------
for(i.plot in c("")){
  
  # -- FG.1
  #==========================
  DimPlot(src.fg1.integrated, dims = c(1,2), group.by = "cluster.v06.26.re", reduction = "umap")
  
  src.fg1.integrated.merge = merge(src.fg1.integrated, src.fg1.tracing.re)
  src.fg1.integrated.selectgene = src.endoderm.fg1.ext_v1.1.re.selectgene
  src.fg1.integrated.merge = batch_process_integration_Re(cluster_names = 'FG.1', set_src = T, 
                                                          seurat = src.fg1.integrated.merge,
                                                          red_refer = "fdl_pca")
  src.fg1.integrated.merge = integration_fta(
    seurat = src.fg1.integrated.merge, 
    cell_refer = rownames(src.fg1.integrated.merge@meta.data[src.fg1.integrated.merge$Source_tech%in%"refer",]),
    cell_query = rownames(src.fg1.integrated.merge@meta.data[src.fg1.integrated.merge$Source_tech%in%"query",]),
    seurat.selectgene = src.fg1.integrated.selectgene, 
    embeddings = src.fg1.integrated@reductions$umap@cell.embeddings)
  
  src.fg1.integrated.merge$cluster.v06.26.re_correct_mnn_umap_fta = 
    src.fg1.integrated.merge$cluster.v06.26.re
  src.fg1.integrated.merge@meta.data[
    colnames(src.fg1.tracing.re),]$cluster.v06.26.re_correct_mnn_umap_fta = 
    src.fg1.tracing.re$cluster.v06.26.re_correct
  
  save(src.fg1.integrated.merge, file = "figure.v08.07/organ_development_re_v240115/src.fg1.integrated.merge.Rdata")
  save(src.fg1.tracing, src.fg1.tracing.re, file = "figure.v08.07/organ_development_re_v240115/src.fg1.tracing.summary.Rdata")
  
  DimPlot(src.fg1.integrated.merge, dims = c(1,2),
          group.by = "cluster.v06.26.re_correct_mnn_umap_fta", reduction = "mnn_umap_fta")
  
  for(endoderm in endoderm_list){
    if(endoderm != "FG.1"){next()}
    
    endoderm = gsub("\\.","",tolower(endoderm))
    seurat = get(paste("src.",endoderm,".integrated.merge",sep=""))
    
    for(names in names(seurat@reductions)){
      if(names != "mnn_umap_fta"){next()}
      data_plot = reduct_to_meta(seurat, names)
      
      
      g1 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.v06.26.re_correct_mnn_umap_fta), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=cluster.v06.26.re_correct_mnn_umap_fta), 
                   # color = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g2 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2, color=lineage), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g3 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.extract.v1.1), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add + #ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g4 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=Time), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=Time), 
                   # color = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      
      for(k in c(1:4)){
        png(filename = paste("figure.v08.07/organ_development_re_v240115/png_re//",
                             "Pathway_",endoderm,"_pic",as.character(k),".png", sep=""),
            width = 1000,height = 1000,pointsize = 20)
        print(get(paste("g",as.character(k),sep="")))
        dev.off()
      }
    }
  }
  #==========================
  
  # -- FG.2
  #==========================
  DimPlot(src.fg2.integrated_re, dims = c(1,2),
          group.by = "Time", reduction = "umap")
  DimPlot(src.fg2.integrated_re, dims = c(1,2),
          group.by = "cluster.v06.26.re_hc", reduction = "umap")
  
  src.fg2.integrated.merge = integration_fta(
    seurat = src.fg2.integrated.merge,
    cell_refer = rownames(src.fg2.integrated.merge@meta.data[
      src.fg2.integrated.merge$Source_tech%in%'refer',]),
    cell_query = rownames(src.fg2.integrated.merge@meta.data[
      src.fg2.integrated.merge$Source_tech%in%'query',]),
    seurat.selectgene = src.endoderm.fg2.ext_v1.1.selectgene,
    embeddings = src.fg2.integrated_re[["umap"]]@cell.embeddings)
  
  src.fg2.integrated.merge = label_KNN_learn_Re(
    src.fg2.integrated.merge, reduction = "mnn_umap_fta", 
    label = "cluster.v06.26.re_hc", group = "Source_tech")
  
  for(endoderm in endoderm_list){
    if(endoderm != "FG.2"){next()}
    
    endoderm = gsub("\\.","",tolower(endoderm))
    seurat = get(paste("src.",endoderm,".integrated.merge",sep=""))
    
    for(names in names(seurat@reductions)){
      if(names != "mnn_umap_fta"){next()}
      data_plot = reduct_to_meta(seurat, names)
      
      
      g1 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.v06.26.re_hc), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=cluster.v06.26.re_hc_mnn_umap_fta), 
                   # color = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g2 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2, color=lineage), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g3 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.extract.v1.1), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add + #ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g4 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=Time), 
                   size=3)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=Time), 
                   # color = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      
      for(k in c(1:4)){
        png(filename = paste("figure.v08.07/organ_development_re_v240115/png_re//",
                             "Pathway_",endoderm,"_pic",as.character(k),".png", sep=""),
            width = 1000,height = 1000,pointsize = 20)
        print(get(paste("g",as.character(k),sep="")))
        dev.off()
      }
    }
  }
  #==========================
  
  # -- FG.3
  #==========================
  DimPlot(src.fg3.integrated_re, group.by =  "cluster.v06.26.re_correct_re", reduction = "umap_mon3_2d") 
  DimPlot(src.fg3.integrated_re, group.by =  "cluster.v06.26.re_correct_refine.re", reduction = "umap_mon3_2d")
  # DimPlot(src.fg3.integrated.merge, reduction = "mnn_umap_fta")
  
  src.fg3.integrated_re$cluster.v06.26.re_correct_refine.re = src.fg3.integrated_re$cluster.v06.26.re_correct_refine
  src.fg3.integrated_re@meta.data[
    src.fg3.integrated_re$cluster.v06.26.re_correct_refine%in%c('FG.3') &
      src.fg3.integrated_re$cluster.v06.26.re_correct_re%in%c('Pharynx.organ.4'),]$cluster.v06.26.re_correct_refine.re = "Pharynx.organ.4"
  
  src.fg3.integrated.merge$cluster.v06.26.re_correct_refine = NA
  src.fg3.integrated.merge@meta.data[colnames(src.fg3.integrated_re),]$cluster.v06.26.re_correct_refine = 
    src.fg3.integrated_re$cluster.v06.26.re_correct_refine.re
  
  src.fg3.integrated.merge = label_KNN_learn_Re(
    src.fg3.integrated.merge, reduction = "mnn_umap_fta", 
    label = "cluster.v06.26.re_correct_refine", group = "Source_tech")
  src.fg3.tracing$cluster.v06.26.re_correct_refine_mnn_umap_fta = 
    src.fg3.integrated.merge@meta.data[colnames(src.fg3.tracing),]$cluster.v06.26.re_correct_refine_mnn_umap_fta
  
  for(endoderm in endoderm_list){
    if(endoderm != "FG.3"){next()}
    
    endoderm = gsub("\\.","",tolower(endoderm))
    seurat = get(paste("src.",endoderm,".integrated.merge",sep=""))
    
    for(names in names(seurat@reductions)){
      if(names != "mnn_umap_fta"){next()}
      data_plot = reduct_to_meta(seurat, names)
      
      
      g1 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.v06.26.re_correct_refine_mnn_umap_fta), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=cluster.v06.26.re_correct_refine_mnn_umap_fta), 
                   # color = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g2 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2, color=lineage), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g3 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.extract.v1.1), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add + #ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g4 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=Time), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=Time), 
                   # color = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      
      for(k in c(1:4)){
        png(filename = paste("figure.v08.07/organ_development_re_v240115/png_re//",
                             "Pathway_",endoderm,"_pic",as.character(k),".png", sep=""),
            width = 1000,height = 1000,pointsize = 20)
        print(get(paste("g",as.character(k),sep="")))
        dev.off()
      }
    }
  }
  #==========================
  
  # -- FG.4
  #==========================
  DimPlot(src.fg4.integrated_re, reduction = "umap_mnn",
          group.by = "cluster.v06.26.re_correct")
  table(src.fg4.integrated_refine$cluster.v06.26z)
  
  src.fg4.tracing.re = src.fg4.tracing[,rownames(src.fg4.integrated.merge@meta.data[src.fg4.integrated.merge$Source_tech%in%"query",])]
  
  src.fg4.integrated.merge.re = merge(src.fg4.integrated_refine, src.fg4.tracing.re)
  src.fg4.integrated.selectgene = src.endoderm.fg4.ext.re.selectgene
  src.fg4.integrated.merge.re = batch_process_integration_Re(cluster_names = 'FG.4', set_src = T, 
                                                             seurat = src.fg4.integrated.merge.re,
                                                             red_refer = "umap_mnn")
  src.fg4.integrated.merge.re = label_KNN_learn_Re(
    src.fg4.integrated.merge.re, reduction = "mnn_umap_fta",
    label = "cluster.v06.26.re_correct_refine", group = "Source_tech")
  
  save(src.fg4.integrated.merge.re,
       file = "figure.v08.07/organ_development_re_v240115/src.fg4.integrated.merge.re.Rdata")
  DimPlot(src.fg4.integrated.merge.re, reduction = "mnn_umap_fta",
          group.by = "cluster.v06.26.re_correct_refine_mnn_umap_fta")
  
  
  src.fg4.integrated.merge.re@meta.data[colnames(src.fg4.tracing.re),]$cluster.v06.26.re_correct_refine_mnn_umap_fta = 
    src.fg4.tracing.re$cluster.v06.26.re_correct
  
  for(endoderm in endoderm_list){
    if(endoderm != "FG.4"){next()}
    
    endoderm = gsub("\\.","",tolower(endoderm))
    seurat = get(paste("src.",endoderm,".integrated.merge.re",sep=""))
    
    for(names in names(seurat@reductions)){
      if(names != "mnn_umap_fta"){next()}
      data_plot = reduct_to_meta(seurat, names)
      
      
      g1 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.v06.26.re_correct_refine_mnn_umap_fta), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=cluster.v06.26.re_correct_refine_mnn_umap_fta), 
                   # color = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g2 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2, color=lineage), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g3 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.extract.v1.1), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add + #ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g4 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=Time), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=Time), 
                   # color = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      
      for(k in c(1:4)){
        png(filename = paste("figure.v08.07/organ_development_re_v240115/png_re//",
                             "Pathway_",endoderm,"_pic",as.character(k),".png", sep=""),
            width = 1000,height = 1000,pointsize = 20)
        print(get(paste("g",as.character(k),sep="")))
        dev.off()
      }
    }
  }
  #==========================
  
  # -- FG.5
  #==========================
  DimPlot(src.fg5.integrated, reduction = 'umap_integrated',
          group.by =  "cluster.v06.26.re_correct")
  DimPlot(src.fg5.tracing, reduction = "mnn_umap_fta",
          group.by = "cluster.v06.26.re_mnn_umap_fta")
  
  save(src.fg5.integrated, file = "figure.v08.07/organ_development_re_v240115/src.fg5.integrated.Rdata")
  save(src.fg5.tracing, file = "figure.v08.07/organ_development_re_v240115/src.fg5.tracing.Rdata")
  
  src.fg5.integrated.merge = merge(src.fg5.integrated, src.fg5.tracing)
  src.fg5.integrated.selectgene = src.endoderm.fg5.ext.v1.1.selectgene
  src.fg5.integrated.merge = batch_process_integration_Re(cluster_names = 'FG.5', set_src = T, 
                                                          seurat = src.fg5.integrated.merge,
                                                          red_refer = "umap_integrated")
  
  src.fg5.integrated.merge$cluster.v06.26.re_correct_mnn_umap_fta = 
    src.fg5.integrated.merge$cluster.v06.26.re_mnn_umap_fta
  src.fg5.integrated.merge@meta.data[
    colnames(src.fg5.integrated),]$cluster.v06.26.re_correct_mnn_umap_fta = 
    src.fg5.integrated$cluster.v06.26.re_correct
  src.fg5.integrated.merge@meta.data[
    colnames(src.fg5.tracing),]$cluster.v06.26.re_correct_mnn_umap_fta = 
    src.fg5.tracing$cluster.v06.26.re_correct
  
  for(endoderm in endoderm_list){
    if(endoderm != "FG.5"){next()}
    
    endoderm = gsub("\\.","",tolower(endoderm))
    seurat = get(paste("src.",endoderm,".integrated.merge",sep=""))
    
    for(names in names(seurat@reductions)){
      if(names != "mnn_umap_fta"){next()}
      data_plot = reduct_to_meta(seurat, names)
      
      
      g1 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.v06.26.re_correct_mnn_umap_fta), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=cluster.v06.26.re_correct_mnn_umap_fta), 
                   # color = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g2 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2, color=lineage), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g3 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.extract.v1.1), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add + #ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g4 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=Time), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=Time), 
                   # color = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      
      for(k in c(1:4)){
        png(filename = paste("figure.v08.07/organ_development_re_v240115/png_re//",
                             "Pathway_",endoderm,"_pic",as.character(k),".png", sep=""),
            width = 1000,height = 1000,pointsize = 20)
        print(get(paste("g",as.character(k),sep="")))
        dev.off()
      }
    }
  }
  #==========================
  
  # -- FG.6
  #==========================
  DimPlot(src.fg6.integrated, reduction = "fdl_mnn")
  
  src.fg6.integrated.merge.re = src.fg6.integrated.merge[, colnames(src.fg6.integrated.merge)%in%c(
    colnames(src.fg6.integrated), colnames(src.fg6.tracing.re))]
  
  seurat = src.fg6.integrated.merge.re
  cell_refer = colnames(src.fg6.integrated)
  cell_query = setdiff(colnames(seurat), cell_refer)
  seurat.selectgene = src.fg6.integrated.filtergene
  umap_embedding = src.fg6.integrated@reductions$fdl_mnn@cell.embeddings
  
  
  for(assay in c("RNA", "mnnRNA")){
    anchor.integrated = 
      FindTransferAnchors(reference = seurat[,cell_refer],
                          query = seurat[,cell_query],
                          reference.assay =  assay,
                          query.assay =  assay, scale = T,
                          features = seurat.selectgene)
    
    umap.transfer = TransferData(anchor.integrated, t(umap_embedding))
    
    assay_name = gsub('RNA',"", gsub("mnn","mnn_",assay))
    
    seurat[[paste(assay_name, "umap_fta", sep="")]] =  seurat[["umap"]]
    
    seurat[[paste(assay_name, "umap_fta", sep="")]]@cell.embeddings[cell_refer,] = umap_embedding[,c(1:2)]
    
    seurat[[paste(assay_name, "umap_fta", sep="")]]@cell.embeddings[cell_query,] = as.matrix(t(umap.transfer@data))[,c(1:2)]
    
    print(paste(assay,"Find-Transfer-anchor processing done!"))
  }
  src.fg6.integrated.merge.re = seurat
  
  DimPlot(src.fg6.integrated.merge.re, reduction = 'mnn_umap_fta"',
          group.by =  "cluster.v06.26.re_correct_mnn_umap_fta")
  DimPlot(src.fg6.tracing.re, reduction = "mnn_umap_fta",
          group.by = "cluster.v06.26.re_correct")
  
  src.fg6.integrated.merge.re$cluster.v06.26.re_correct_mnn_umap_fta = 
    src.fg6.integrated.merge.re$cluster.v06.26.re
  src.fg6.integrated.merge.re@meta.data[
    colnames(src.fg6.tracing.re),]$cluster.v06.26.re_correct_mnn_umap_fta = 
    src.fg6.tracing.re$cluster.v06.26.re_correct
  
  src.fg6.integrated.merge.re = integration_fta(
    seurat = src.fg6.integrated.merge.re, 
    cell_refer = rownames(src.fg6.integrated.merge.re@meta.data[src.fg6.integrated.merge.re$Source_tech%in%"refer",]),
    cell_query = rownames(src.fg6.integrated.merge.re@meta.data[src.fg6.integrated.merge.re$Source_tech%in%"query",]),
    seurat.selectgene = src.fg6.integrated.filtergene,
    embeddings = src.fg6.integrated@reductions$umap@cell.embeddings)
  
  
  for(endoderm in endoderm_list){
    if(endoderm != "FG.6"){next()}
    
    endoderm = gsub("\\.","",tolower(endoderm))
    seurat = get(paste("src.",endoderm,".integrated.merge.re",sep=""))
    
    for(names in names(seurat@reductions)){
      if(names != "mnn_umap_fta"){next()}
      data_plot = reduct_to_meta(seurat, names)
      
      
      g1 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.v06.26.re_correct_mnn_umap_fta), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=cluster.v06.26.re_correct_mnn_umap_fta), 
                   # color = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g2 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2, color=lineage), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g3 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.extract.v1.1), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add + #ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g4 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=Time), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=Time), 
                   # color = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      
      for(k in c(1:4)){
        png(filename = paste("figure.v08.07/organ_development_re_v240115/png_re//",
                             "Pathway_",endoderm,"_pic",as.character(k),".png", sep=""),
            width = 1000,height = 1000,pointsize = 20)
        print(get(paste("g",as.character(k),sep="")))
        dev.off()
      }
    }
  }
  #==========================
  
  # -- AL.1
  #==========================
  DimPlot(src.al1.integrated.re, reduction = "umap_integrated")
  src.al1.integrated.merge.re = merge(src.al1.integrated.re, 
                                      src.al1.tracing_re)
  src.al1.integrated.selectgene = src.endoderm.al12.ext.v1.1.selectgene
  src.al1.integrated.merge.re = batch_process_integration_Re(cluster_names = 'AL.1', set_src = T, 
                                                             seurat = src.al1.integrated.merge.re,
                                                             red_refer = "umap_integrated")
  
  src.al1.integrated.merge.re$cluster.v06.26.re_correct_re_mnn_umap_fta = 
    src.al1.integrated.merge.re$cluster.v06.26.re_correct_re
  
  for(endoderm in endoderm_list){
    if(endoderm != "AL.1"){next()}
    
    endoderm = gsub("\\.","",tolower(endoderm))
    seurat = get(paste("src.",endoderm,".integrated.merge.re",sep=""))
    
    for(names in names(seurat@reductions)){
      if(names != "mnn_umap_fta"){next()}
      data_plot = reduct_to_meta(seurat, names)
      
      
      g1 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.v06.26.re_correct_re_mnn_umap_fta), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=cluster.v06.26.re_correct_re_mnn_umap_fta), 
                   # color = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g2 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2, color=lineage), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g3 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.extract.v1.1), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add + #ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g4 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=Time), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=Time), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      
      for(k in c(1:4)){
        png(filename = paste("figure.v08.07/organ_development_re_v240115/png_re//",
                             "Pathway_",endoderm,"_pic",as.character(k),".png", sep=""),
            width = 1000,height = 1000,pointsize = 20)
        print(get(paste("g",as.character(k),sep="")))
        dev.off()
      }
    }
  }
  #==========================
  
  # -- AL.2
  #==========================
  DimPlot(src.al2.integrated.re, reduction = "umap_integrated")
  src.al2.integrated.merge.re = merge(src.al2.integrated.re, src.al2.tracing)
  src.al2.integrated.selectgene = src.endoderm.al12.ext.v1.1.selectgene
  src.al2.integrated.merge.re = batch_process_integration_Re(cluster_names = 'AL.2', set_src = T, 
                                                             seurat = src.al2.integrated.merge.re,
                                                             red_refer = "umap_integrated")
  
  src.al2.integrated.merge.re$cluster.v06.26.re_correct_re_mnn_umap_fta = 
    src.al2.integrated.merge.re$cluster.v06.26.re_correct_re
  
  for(endoderm in endoderm_list){
    if(endoderm != "AL.2"){next()}
    
    endoderm = gsub("\\.","",tolower(endoderm))
    seurat = get(paste("src.",endoderm,".integrated.merge.re",sep=""))
    
    for(names in names(seurat@reductions)){
      if(names != "mnn_umap_fta"){next()}
      data_plot = reduct_to_meta(seurat, names)
      
      
      g1 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.v06.26.re_correct_re_mnn_umap_fta), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=cluster.v06.26.re_correct_re_mnn_umap_fta), 
                   # color = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g2 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2, color=lineage), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g3 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.extract.v1.1), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add + #ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g4 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=Time), 
                   size=3.3)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=Time), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      
      for(k in c(1:4)){
        png(filename = paste("figure.v08.07/organ_development_re_v240115/png_re//",
                             "Pathway_",endoderm,"_pic",as.character(k),".png", sep=""),
            width = 1000,height = 1000,pointsize = 20)
        print(get(paste("g",as.character(k),sep="")))
        dev.off()
      }
    }
  }
  #==========================
  
  # -- AL.3
  #==========================
  DimPlot(src.al3.integrated.re, reduction = "umap_integrated", group.by = "cluster.v06.26.re_correct")
  src.al3.tracing_re = src.al3.integrated.merge[,src.al3.integrated.merge$Source_tech%in%"query"]
  
  src.al3.integrated.merge.re = merge(src.al3.integrated.re, src.al3.tracing_re)
  src.al3.integrated.selectgene = src.endoderm.al3.ext.v1.1.selectgene
  src.al3.integrated.merge.re = batch_process_integration_Re(cluster_names = 'AL.3', set_src = T, 
                                                             seurat = src.al3.integrated.merge.re,
                                                             red_refer = "umap_integrated")
  
  src.al3.integrated.merge.re = label_KNN_learn_Re(
    src.al3.integrated.merge.re, reduction = "mnn_umap_fta",
    label = "cluster.v06.26.re_correct", group = "Source_tech")
  
  
  src.al3.integrated.merge.re@meta.data[
    colnames(src.al3.tracing_re),]$cluster.v06.26.re_correct_mnn_umap_fta =
    src.al3.tracing_re$cluster.v06.26.re_correct.refine
  
  for(endoderm in endoderm_list){
    if(endoderm != "AL.3"){next()}
    
    endoderm = gsub("\\.","",tolower(endoderm))
    seurat = get(paste("src.",endoderm,".integrated.merge.re",sep=""))
    
    for(names in names(seurat@reductions)){
      if(names != "mnn_umap_fta"){next()}
      data_plot = reduct_to_meta(seurat, names)
      
      
      g1 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.v06.26.re_correct_mnn_umap_fta), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=cluster.v06.26.re_correct_mnn_umap_fta), 
                   # color = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g2 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2, color=lineage), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g3 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.extract.v1.1), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add + #ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g4 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=Time), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=Time), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      
      for(k in c(1:4)){
        png(filename = paste("figure.v08.07/organ_development_re_v240115/png_re//",
                             "Pathway_",endoderm,"_pic",as.character(k),".png", sep=""),
            width = 1000,height = 1000,pointsize = 20)
        print(get(paste("g",as.character(k),sep="")))
        dev.off()
      }
    }
  }
  #==========================
  
  # -- MG.1
  #==========================
  DimPlot(src.mg1.integrated, reduction = 'umap_mnn',
          group.by =  "cluster.v06.26.re")
  DimPlot(src.mg1.tracing, reduction = "mnn_umap_fta",
          group.by = "lineage")
  # group.by = "cluster.v06.26.re_correct")
  
  src.mg1.integrated.merge$cluster.v06.26.re_correct_mnn_umap_fta =
    src.mg1.integrated.merge$cluster.v06.26.re_mnn_umap_fta
  src.mg1.integrated.merge@meta.data[src.mg1.integrated.merge$Source_tech%in%"query",]$cluster.v06.26.re_correct_mnn_umap_fta =
    src.mg1.tracing@meta.data[rownames(src.mg1.integrated.merge@meta.data[src.mg1.integrated.merge$Source_tech%in%"query",]),]$cluster.v06.26.re_correct
  
  src.mg1.integrated.merge = integration_fta(
    seurat = src.mg1.integrated.merge,
    cell_refer = rownames(src.mg1.integrated.merge@meta.data[
      src.mg1.integrated.merge$Source_tech%in%'refer',]),
    cell_query = rownames(src.mg1.integrated.merge@meta.data[
      src.mg1.integrated.merge$Source_tech%in%'query',]),
    seurat.selectgene = src.endoderm.mg1.ext.v1.1.selectgene,
    embeddings = src.mg1.integrated[["umap_mnn"]]@cell.embeddings)
  
  for(endoderm in endoderm_list){
    if(endoderm != "MG.1"){next()}
    
    endoderm = gsub("\\.","",tolower(endoderm))
    seurat = get(paste("src.",endoderm,".integrated.merge",sep=""))
    
    for(names in names(seurat@reductions)){
      if(names != "mnn_umap_fta"){next()}
      data_plot = reduct_to_meta(seurat, names)
      
      
      g1 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.v06.26.re_correct_mnn_umap_fta), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=cluster.v06.26.re_correct_mnn_umap_fta), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g2 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2, color=lineage), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g3 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.extract.v1.1), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add + #ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g4 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=Time), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=Time), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      
      for(k in c(1:4)){
        png(filename = paste("figure.v08.07/organ_development_re_v240115/png_re//",
                             "Pathway_",endoderm,"_pic",as.character(k),".png", sep=""),
            width = 1000,height = 1000,pointsize = 20)
        print(get(paste("g",as.character(k),sep="")))
        dev.off()
      }
    }
  }
  #==========================
  
  # -- MG.2
  #==========================
  DimPlot(src.mg2.integrated, reduction = 'umap_integrated',
          group.by =  "cluster.v06.26.re")
  DimPlot(src.mg2.tracing, reduction = "mnn_umap_fta",
          group.by = "cluster.v06.26.re_correct")
  
  src.mg2.integrated.merge$cluster.v06.26.re_correct_mnn_umap_fta = src.mg2.integrated.merge$cluster.v06.26.re_mnn_umap_fta
  src.mg2.integrated.merge$cluster.v06.26.re_correct_mnn_umap_fta[colnames(src.mg2.tracing)] = 
    src.mg2.tracing$cluster.v06.26.re_correct
  
  for(endoderm in endoderm_list){
    if(endoderm != "MG.2"){next()}
    
    endoderm = gsub("\\.","",tolower(endoderm))
    seurat = get(paste("src.",endoderm,".integrated.merge",sep=""))
    
    for(names in names(seurat@reductions)){
      if(names != "mnn_umap_fta"){next()}
      data_plot = reduct_to_meta(seurat, names)
      
      
      g1 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.v06.26.re_correct_mnn_umap_fta), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=cluster.v06.26.re_correct_mnn_umap_fta), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g2 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2, color=lineage), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g3 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.extract.v1.1), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add + #ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g4 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=Time), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=Time), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      
      for(k in c(1:4)){
        png(filename = paste("figure.v08.07/organ_development_re_v240115/png_re//",
                             "Pathway_",endoderm,"_pic",as.character(k),".png", sep=""),
            width = 1000,height = 1000,pointsize = 20)
        print(get(paste("g",as.character(k),sep="")))
        dev.off()
      }
    }
  }
  #==========================
  
  # -- MG.3
  #==========================
  DimPlot(src.mg3.integrated, reduction = 'umap_mnn',
          group.by =  "cluster.v06.26.re")
  DimPlot(src.mg3.integrated.merge, reduction = "mnn_umap_fta",
          group.by = "cluster.v06.26.re_mnn_umap_fta")
  
  src.mg3.integrated.merge = integration_fta(
    seurat = src.mg3.integrated.merge,
    cell_refer = rownames(src.mg3.integrated.merge@meta.data[
      src.mg3.integrated.merge$Source_tech%in%'refer',]),
    cell_query = rownames(src.mg3.integrated.merge@meta.data[
      src.mg3.integrated.merge$Source_tech%in%'query',]),
    seurat.selectgene = src.endoderm.mg3.ext.v1.1.selectgene,
    embeddings = src.mg3.integrated[["umap_mnn"]]@cell.embeddings)
  
  src.mg3.integrated.merge@meta.data[
    colnames(src.mg3.tracing),]$cluster.v06.26.re_mnn_umap_fta = 
    src.mg3.tracing$cluster.v06.26.re_correct
  
  
  for(endoderm in endoderm_list){
    if(endoderm != "MG.3"){next()}
    
    endoderm = gsub("\\.","",tolower(endoderm))
    seurat = get(paste("src.",endoderm,".integrated.merge",sep=""))
    
    for(names in names(seurat@reductions)){
      if(names != "mnn_umap_fta"){next()}
      data_plot = reduct_to_meta(seurat, names)
      
      
      g1 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.v06.26.re_mnn_umap_fta), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=cluster.v06.26.re_mnn_umap_fta), 
                   # color = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g2 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2, color=lineage), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g3 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.extract.v1.1), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add + #ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g4 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=Time), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=Time), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      
      for(k in c(1:4)){
        png(filename = paste("figure.v08.07/organ_development_re_v240115/png_re//",
                             "Pathway_",endoderm,"_pic",as.character(k),".png", sep=""),
            width = 1000,height = 1000,pointsize = 20)
        print(get(paste("g",as.character(k),sep="")))
        dev.off()
      }
    }
  }
  #==========================
  
  # -- HG.1
  #==========================
  DimPlot(src.hg1.integrated.re, reduction = 'umap_mnn_rot',
          group.by =  "cluster.v06.26.re_hc")
  DimPlot(src.hg1.integrated.merge.re, reduction = "mnn_umap_fta",
          group.by = "cluster.v06.26.re_mnn_umap_fta")
  
  src.hg1.integrated.merge = 
    label_KNN_learn_Re(seurat = src.hg1.integrated.merge.re,
                       reduction = "mnn_umap_fta", 
                       label = "cluster.v06.26.re_hc", 
                       group = "Source_tech")
  
  src.hg1.integrated.merge@meta.data[
    colnames(src.hg1.tracing),]$cluster.v06.26.re_hc_mnn_umap_fta =
    src.hg1.tracing$cluster.v06.26.re_correct.refine
  
  
  for(endoderm in endoderm_list){
    if(endoderm != "HG.1"){next()}
    
    lineage_list = t(hash::values(endoderm_lineage_raw[[endoderm]], keys="ss9"))
    endoderm = gsub("\\.","",tolower(endoderm))
    seurat = get(paste("src.",endoderm,".integrated.merge.re",sep=""))
    
    for(names in names(seurat@reductions)){
      if(names != "mnn_umap_fta"){next()}
      data_plot = reduct_to_meta(seurat, names)
      data_plot =  data_plot[data_plot$lineage%in%c(NA,lineage_list),]
      
      g1 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.v06.26.re_hc_mnn_umap_fta), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=cluster.v06.26.re_hc_mnn_umap_fta), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g2 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2, color=lineage), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g3 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.extract.v1.1), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add + #ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g4 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=Time), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=Time), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      
      for(k in c(1:4)){
        png(filename = paste("figure.v08.07/organ_development_re_v240115/png_re/",
                             "Pathway_",endoderm,"_pic",as.character(k),".re.png", sep=""),
            width = 1000,height = 1000,pointsize = 20)
        print(get(paste("g",as.character(k),sep="")))
        dev.off()
      }
    }
  }
  #==========================
  
  # -- HG.2
  #=========================
  DimPlot(src.hg2.integrated, reduction = 'umap_integrated',
          group.by =  "cluster.v06.26.re")
  DimPlot(src.hg2.integrated.merge, reduction = "mnn_umap_fta",
          group.by = "cluster.v06.26.re_mnn_umap_fta")
  
  src.hg2.integrated.merge = 
    label_KNN_learn_Re(seurat = src.hg2.integrated.merge,
                       reduction = "mnn_umap_fta", 
                       label = "cluster.v06.26.re", 
                       group = "Source_tech")
  
  src.hg2.integrated.merge@meta.data[
    colnames(src.hg2.tracing),]$cluster.v06.26.re_mnn_umap_fta =
    src.hg2.tracing$cluster.v06.26.re_correct
  
  
  for(endoderm in endoderm_list){
    if(endoderm != "HG.2"){next()}
    
    endoderm = gsub("\\.","",tolower(endoderm))
    seurat = get(paste("src.",endoderm,".integrated.merge",sep=""))
    
    for(names in names(seurat@reductions)){
      if(names != "mnn_umap_fta"){next()}
      data_plot = reduct_to_meta(seurat, names)
      
      
      g1 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.v06.26.re_mnn_umap_fta), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=cluster.v06.26.re_mnn_umap_fta), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g2 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2, color=lineage), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g3 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=cluster.extract.v1.1), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2), 
                   colour = "darkgray",
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
        theme_void() + p_add + #ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      g4 = ggplot()+ 
        geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                   mapping = aes(x=Coord_1, y =Coord_2, 
                                 color=Time), 
                   size=3.5)+
        geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                   mapping = aes(x=Coord_1, y =Coord_2,
                                 color=Time), 
                   size=5)+
        scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                      colors.time,colors.time.2))+
        scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                     colors.time,colors.time.2))+
        theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
        theme(legend.position = "none")
      
      
      for(k in c(1:4)){
        png(filename = paste("figure.v08.07/organ_development_re_v240115/png_re//",
                             "Pathway_",endoderm,"_pic",as.character(k),".png", sep=""),
            width = 1000,height = 1000,pointsize = 20)
        print(get(paste("g",as.character(k),sep="")))
        dev.off()
      }
    }
  }
  #=========================
}
#--------------------------------------------



# -- Summary-cell type For Tracing Data
#--------------------------------------------
src.sm3.merge$cluster.v06.26.re_mnn_umap_fta = NA
src.sm3.merge@meta.data[colnames(src.fg1.tracing.re),]$cluster.v06.26.re_mnn_umap_fta = 
  src.fg1.tracing.re$cluster.v06.26.re_correct
src.sm3.merge@meta.data[colnames(src.fg2.tracing),]$cluster.v06.26.re_mnn_umap_fta = 
  src.fg2.tracing$cluster.v06.26.re_correct
src.sm3.merge@meta.data[colnames(src.fg3.tracing),]$cluster.v06.26.re_mnn_umap_fta = 
  src.fg3.tracing$cluster.v06.26.re_correct
src.sm3.merge@meta.data[colnames(src.fg4.tracing),]$cluster.v06.26.re_mnn_umap_fta = 
  src.fg4.tracing$cluster.v06.26.re_mnn_umap_fta
src.sm3.merge@meta.data[colnames(src.fg5.tracing),]$cluster.v06.26.re_mnn_umap_fta = 
  src.fg5.tracing$cluster.v06.26.re_correct
src.sm3.merge@meta.data[colnames(src.fg6.tracing.re),]$cluster.v06.26.re_mnn_umap_fta = 
  src.fg6.tracing.re$cluster.v06.26.re_correct
src.sm3.merge@meta.data[colnames(src.al1.tracing_re),]$cluster.v06.26.re_mnn_umap_fta = 
  src.al1.tracing_re$cluster.v06.26.re_correct_re
src.sm3.merge@meta.data[colnames(src.al2.tracing),]$cluster.v06.26.re_mnn_umap_fta = 
  src.al2.tracing$cluster.v06.26.re_correct
src.sm3.merge@meta.data[colnames(src.al3.tracing),]$cluster.v06.26.re_mnn_umap_fta = 
  src.al3.tracing$cluster.v06.26.re_mnn_umap_fta
src.sm3.merge@meta.data[colnames(src.mg1.tracing),]$cluster.v06.26.re_mnn_umap_fta = 
  src.mg1.tracing$cluster.v06.26.re_mnn_umap_fta
src.sm3.merge@meta.data[colnames(src.mg2.tracing),]$cluster.v06.26.re_mnn_umap_fta = 
  src.mg2.tracing$cluster.v06.26.re_correct
src.sm3.merge@meta.data[colnames(src.mg3.tracing),]$cluster.v06.26.re_mnn_umap_fta = 
  src.mg3.tracing$cluster.v06.26.re_mnn_umap_fta
src.sm3.merge@meta.data[colnames(src.hg1.tracing),]$cluster.v06.26.re_mnn_umap_fta = 
  src.hg1.tracing$cluster.v06.26.re_mnn_umap_fta
src.sm3.merge@meta.data[colnames(src.hg2.tracing),]$cluster.v06.26.re_mnn_umap_fta = 
  src.hg2.tracing$cluster.v06.26.re_mnn_umap_fta

# table(src.sm3.merge@meta.data[!src.sm3.merge$Time%in%"E10.5",]$cluster.v06.26.re_mnn_umap_fta%in%NA) 

cell_tracing_refer = rownames(src.sm3.merge@meta.data[
  !src.sm3.merge$cluster.v06.26.re_mnn_umap_fta%in%NA,])
cell_tracing_query = rownames(src.sm3.merge@meta.data[
  src.sm3.merge$cluster.v06.26.re_mnn_umap_fta%in%NA & 
    !src.sm3.merge$Time%in%"E10.5",])

label_learning = FNN::knn(
    seurat@reductions$umap_fta_rpca@cell.embeddings[
      intersect(cell_sort_src.endoderm.merge, cell_tracing_refer), ],
    seurat@reductions$umap_fta_rpca@cell.embeddings[
      intersect(cell_sort_src.endoderm.merge, cell_tracing_query),],
    src.sm3.merge$cluster.v06.26.re_mnn_umap_fta[
      intersect(cell_sort_src.endoderm.merge, cell_tracing_refer)],
    k = 10)
src.sm3.merge$cluster.v06.26.re_mnn_umap_fta[
  intersect(cell_sort_src.endoderm.merge, cell_tracing_query)] = as.character(label_learning)
#--------------------------------------------

# -- Summary-cell type For Refernece data :: Cell type
#--------------------------------------------
src.endoderm$cluster.v06.26.re..correct = src.endoderm$cluster.v06.26.re..merge

src.endoderm@meta.data[colnames(src.fg1.integrated),]$cluster.v06.26.re..correct =
  src.fg1.integrated$cluster.v06.26.re_correct
src.endoderm@meta.data[colnames(src.fg2.integrated),]$cluster.v06.26.re..correct =
  src.fg2.integrated$cluster.v06.26.re_hc
src.endoderm@meta.data[colnames(src.fg3.integrated_re),]$cluster.v06.26.re..correct =
  src.fg3.integrated_re$cluster.v06.26.re_correct_re
src.endoderm@meta.data[colnames(src.fg4.integrated_re),]$cluster.v06.26.re..correct =
  src.fg4.integrated_re$cluster.v06.26.re_correct
src.endoderm@meta.data[colnames(src.fg5.integrated),]$cluster.v06.26.re..correct =
  src.fg5.integrated$cluster.v06.26.re_correct
src.endoderm@meta.data[colnames(src.fg6.integrated),]$cluster.v06.26.re..correct =
  src.fg6.integrated$cluster.v06.26.re

src.endoderm@meta.data[colnames(src.al1.integrated.re),]$cluster.v06.26.re..correct =
  src.al1.integrated.re$cluster.v06.26.re_correct
src.endoderm@meta.data[colnames(src.al2.integrated.re),]$cluster.v06.26.re..correct =
  src.al2.integrated.re$cluster.v06.26.re_correct
src.endoderm@meta.data[colnames(src.al3.integrated),]$cluster.v06.26.re..correct =
  src.al3.integrated$cluster.v06.26.re_correct

src.endoderm@meta.data[colnames(src.mg1.integrated),]$cluster.v06.26.re..correct =
  src.mg1.integrated$cluster.v06.26.re
src.endoderm@meta.data[colnames(src.mg2.integrated),]$cluster.v06.26.re..correct =
  src.mg2.integrated$cluster.v06.26.re
src.endoderm@meta.data[colnames(src.mg3.integrated),]$cluster.v06.26.re..correct =
  src.mg3.integrated$cluster.v06.26.re
src.endoderm@meta.data[colnames(src.hg1.integrated),]$cluster.v06.26.re..correct =
  src.hg1.integrated$cluster.v06.26.re_hc
src.endoderm@meta.data[colnames(src.hg2.integrated),]$cluster.v06.26.re..correct =
  src.hg2.integrated$cluster.v06.26.re
#--------------------------------------------


# -- Save :: Metadata_src.endoderm.merge
#-----------------------------------------
metadata_src.endoderm = cbind(src.endoderm@meta.data,
                              src.endoderm@reductions$umap@cell.embeddings)

metadata_src.sm3.merge = cbind.data.frame(
  src.sm3.merge@meta.data[
    intersect(cell_sort_src.endoderm.merge, 
              c(cell_tracing_refer, cell_tracing_query)),],
  seurat@reductions$umap_fta_pca@cell.embeddings[
    intersect(cell_sort_src.endoderm.merge, 
              c(cell_tracing_refer, cell_tracing_query)),],
  seurat@reductions$umap_fta_rpca@cell.embeddings[
    intersect(cell_sort_src.endoderm.merge, 
              c(cell_tracing_refer, cell_tracing_query)),])

save(metadata_src.endoderm,
     metadata_src.sm3.merge,
     file = "figure.v08.07/integrated_v240115/metadata_src.endoderm.merge.Rdata")
#------------------------------------------




















