

#--  RE-Name
#-----------------------
src.fg1.integrated = src.endoderm.fg1.ext_v1.1.re
src.fg2.integrated = src.endoderm.fg2.ext_v1.1
src.fg3.integrated = src.endoderm.fg3.ext_v1.1
# src.fg4.integrated = src.endoderm.fg4.ext_v1.1.re
src.fg4.integrated = src.fg4.integrated_re
src.fg5.integrated = src.endoderm.fg5.ext_v1.1
src.fg6.integrated = src.endoderm.fg6.ext_v1.1.re

src.al1.integrated = src.endoderm.al1.ext_v1.1
src.al2.integrated = src.endoderm.al2.ext_v1.1
src.al3.integrated = src.endoderm.al3.ext_v1.1
src.mg1.integrated = src.endoderm.mg1.ext_v1.1
src.mg2.integrated = src.endoderm.mg2.ext_v1.1
src.mg3.integrated = src.endoderm.mg3.ext_v1.1
src.hg1.integrated = src.endoderm.hg1.ext_v1.1
src.hg2.integrated = src.endoderm.hg2.ext_v1.1


src.fg1.integrated.selectgene = src.endoderm.fg1.ext_v1.1.re.filtergene
src.fg2.integrated.selectgene = src.endoderm.fg2.ext_v1.1.selectgene
src.fg3.integrated.selectgene = src.endoderm.fg3.ext_v1.1.selectgene
src.fg4.integrated.selectgene = src.endoderm.fg4.ext.v1.1.selectgene
src.fg5.integrated.selectgene = src.endoderm.fg5.ext_v1.1.selectgene
src.fg6.integrated.selectgene = src.endoderm.fg6.ext_v1.1.re.selectgene.fin
src.al1.integrated.selectgene = src.endoderm.al12.ext_v1.1.selectgene
src.al2.integrated.selectgene = src.endoderm.al12.ext_v1.1.selectgene
src.al3.integrated.selectgene = src.endoderm.al3.ext_v1.1.selectgene
src.mg1.integrated.selectgene = src.endoderm.mg1.ext.v1.1.selectgene
src.mg2.integrated.selectgene = src.endoderm.mg2.ext_v1.1.selectgene
src.mg3.integrated.selectgene = src.endoderm.mg3.ext_v1.1.selectgene
src.hg1.integrated.selectgene = src.endoderm.hg1.ext.v1.1.selectgene
src.hg2.integrated.selectgene = src.endoderm.hg2.ext_v1.1.selectgene

save(src.fg1.integrated, file = 'figure.v08.07/organ_development_re_v240115/src.fg1.integrated.Rdata')
save(src.fg2.integrated, file = 'figure.v08.07/organ_development_re_v240115/src.fg2.integrated.Rdata')
save(src.fg3.integrated, file = 'figure.v08.07/organ_development_re_v240115/src.fg3.integrated.Rdata')
save(src.fg4.integrated, file = 'figure.v08.07/organ_development_re_v240115/src.fg4.integrated.Rdata')
save(src.fg5.integrated, file = 'figure.v08.07/organ_development_re_v240115/src.fg5.integrated.Rdata')
save(src.fg6.integrated, file = 'figure.v08.07/organ_development_re_v240115/src.fg6.integrated.Rdata')
save(src.al1.integrated, file = 'figure.v08.07/organ_development_re_v240115/src.al1.integrated.Rdata')
save(src.al2.integrated, file = 'figure.v08.07/organ_development_re_v240115/src.al2.integrated.Rdata')
save(src.al3.integrated, file = 'figure.v08.07/organ_development_re_v240115/src.al3.integrated.Rdata')
save(src.mg1.integrated, file = 'figure.v08.07/organ_development_re_v240115/src.mg1.integrated.Rdata')
save(src.mg2.integrated, file = 'figure.v08.07/organ_development_re_v240115/src.mg2.integrated.Rdata')
save(src.mg3.integrated, file = 'figure.v08.07/organ_development_re_v240115/src.mg3.integrated.Rdata')
save(src.hg1.integrated, file = 'figure.v08.07/organ_development_re_v240115/src.hg1.integrated.Rdata')
save(src.hg2.integrated, file = 'figure.v08.07/organ_development_re_v240115/src.hg2.integrated.Rdata')
#-----------------------

label = "cluster.extract.v1.1"
reduction = c("umap_fta","umap_integrated","mnnumap_integrated")
for(names in reduction){
  src.sm3.merge@meta.data[,paste(label,"_",names,sep="")] = NA
}

for(time in c(9,12,15,18,21,24,27)){
  Time1 = paste(as.character(time),"ss",sep="")
  Time2 = paste("ss",as.character(time),sep="")
  
  src.Time1.integrated.merge = get(
    paste("src.",Time1,".integrated.merge",sep=""))
  
  cell_names = intersect(colnames(src.Time1.integrated.merge),
                         colnames(src.sm3.merge))
  
  
  if(time==9){
    reduction = c("umap_fta","umap_integrated")
  }else if(time==12){
    reduction = c("umap_fta","mnnumap_integrated")
  }else{
    reduction = c("umap_integrated","mnnumap_integrated")
  }
  
  for(names in reduction){
    src.sm3.merge@meta.data[cell_names, paste(label,"_",names,sep="")] = 
      src.Time1.integrated.merge@meta.data[cell_names, paste(label,"_",names,sep="")]
  }
}

src.sm3.merge$cluster.extract.v1.1_define = 
  src.sm3.merge$cluster.extract.v1.1_mnnumap_integrated
src.sm3.merge@meta.data[src.sm3.merge$cluster.extract.v1.1_mnnumap_integrated%in%NA,]$cluster.extract.v1.1_define = 
  src.sm3.merge@meta.data[src.sm3.merge$cluster.extract.v1.1_mnnumap_integrated%in%NA,]$cluster.extract.v1.1_umap_integrated

#-- Define Integrated Merge 
#-- Batch Process
#----------------------
for(names in endoderm_list){
  
  if(names %in% c(#"MG.2",
                  'FG.4')){}else{
    next() 
    }
  
  
  cluster =  gsub("\\.","",tolower(names))
  assign(paste("src.",cluster,".integrated.merge",sep=""),
         batch_process_integration_Re(cluster_names = names))
}

#--  Plot Integrated Merge
#----------------------
for(endoderm in endoderm_list){
  
  endoderm = gsub("\\.","",tolower(endoderm))
  seurat = get(paste("src.",endoderm,".integrated.merge",sep=""))
  
  pdf(paste("figure.v08.07/organ_development_re_v240115/pdf/",
             endoderm,"_embedding.pdf", sep=""))
  
  for(names in names(seurat@reductions)){
    
    data_plot = reduct_to_meta(seurat, names)
    
    g1 = ggplot()+ 
      geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                 mapping = aes(x=Coord_1, y =Coord_2), 
                 colour = "darkgray",size=1.5)+
      geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                 mapping = aes(x=Coord_1, y =Coord_2, color=lineage), 
                 size=1.5, shape=3, stroke=1)+
      scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
      theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
      theme(legend.position = "none")
    
    g2 = ggplot()+ 
      geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                 mapping = aes(x=Coord_1, y =Coord_2), 
                 colour = "darkgray", size=1.5)+
      geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                 mapping = aes(x=Coord_1, y =Coord_2, 
                               color=cluster.extract.v1.1_define), 
                 size=1.5, shape=3, stroke=1)+
      scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
      theme_void() + p_add + #ggtitle(paste(endoderm,names,sep="")) +
      theme(legend.position = "none")
    
    g3 = ggplot()+ 
      geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                 mapping = aes(x=Coord_1, y =Coord_2, 
                               color=cluster.extract.v1.1), 
                 size=1.5)+
      geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                 mapping = aes(x=Coord_1, y =Coord_2), 
                 colour = "darkgray", size=1.5, shape=3, stroke=1)+
      scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
      theme_void() + p_add + #ggtitle(paste(endoderm,names,sep="")) +
      theme(legend.position = "none")
    
    g4 = ggplot()+ 
      geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                 mapping = aes(x=Coord_1, y =Coord_2, 
                               color=Time), 
                 size=1.5)+
      geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                 mapping = aes(x=Coord_1, y =Coord_2), 
                 colour = "darkgray", size=1.5, shape=3, stroke=1)+
      scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,colors.time))+
      theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
      theme(legend.position = "none")
    
    print((g1+g2)/(g3+g4))
    rm(g1,g2,g3,g4)
  }
  dev.off()
}

#--  KNN Learning Annoation
#----------------------
for(endoderm in endoderm_list){
  
  if(endoderm=="FG.2"){
    reduction = c("mnn_umap_fta")
  }else{
    reduction = c("mnn_umap_fta")
  }
  
  endoderm = gsub("\\.","",tolower(endoderm))
  seurat = get(paste("src.",endoderm,".integrated.merge",sep=""))
  

  seurat = label_KNN_learn_Re(seurat = seurat,
                              reduction = reduction,
                              label = "cluster.v06.26.re",
                              group = "Source_tech")

  
  assign(paste("src.",endoderm,".integrated.merge",sep=""),
         seurat)
}

#--  RE-Merga Rdata
#-----------------------
save(src.fg1.integrated.merge, file = 'figure.v08.07/organ_development_re_v240115/src.fg1.integrated.merge.Rdata')
save(src.fg2.integrated.merge, file = 'figure.v08.07/organ_development_re_v240115/src.fg2.integrated.merge.Rdata')
save(src.fg3.integrated.merge, file = 'figure.v08.07/organ_development_re_v240115/src.fg3.integrated.merge.Rdata')
# save(src.fg4.integrated.merge, file = 'figure.v08.07/organ_development_re_v240115/src.fg4.integrated.merge.Rdata')
save(src.fg4.integrated.merge, file = 'figure.v08.07/organ_development_re_v240115/src.fg4.integrated.merge_re.Rdata')

save(src.fg5.integrated.merge, file = 'figure.v08.07/organ_development_re_v240115/src.fg5.integrated.merge.Rdata')
save(src.fg6.integrated.merge, file = 'figure.v08.07/organ_development_re_v240115/src.fg6.integrated.merge.Rdata')


save(src.al1.integrated.merge, file = 'figure.v08.07/organ_development_re_v240115/src.al1.integrated.merge.Rdata')
save(src.al2.integrated.merge, file = 'figure.v08.07/organ_development_re_v240115/src.al2.integrated.merge.Rdata')
save(src.al3.integrated.merge, file = 'figure.v08.07/organ_development_re_v240115/src.al3.integrated.merge.Rdata')
save(src.mg1.integrated.merge, file = 'figure.v08.07/organ_development_re_v240115/src.mg1.integrated.merge.Rdata')
save(src.mg2.integrated.merge, file = 'figure.v08.07/organ_development_re_v240115/src.mg2.integrated.merge.Rdata')
save(src.mg3.integrated.merge, file = 'figure.v08.07/organ_development_re_v240115/src.mg3.integrated.merge.Rdata')
save(src.hg1.integrated.merge, file = 'figure.v08.07/organ_development_re_v240115/src.hg1.integrated.merge.Rdata')
save(src.hg2.integrated.merge, file = 'figure.v08.07/organ_development_re_v240115/src.hg2.integrated.merge.Rdata')
#-----------------------


# Set Tracing Seurat
for(endoderm in endoderm_list){
  endoderm = gsub("\\.","",tolower(endoderm))
  seurat = get(paste("src.",endoderm,".integrated.merge",sep=""))
  seurat = seurat[,seurat$Source_tech%in%"query"]
  
  assign(paste("src.",endoderm,".tracing",sep=""),
         seurat)
}

#--  RE-Tracing Rdata
#-----------------------
save(src.fg1.tracing, file = 'figure.v08.07/organ_development_re_v240115/src.fg1.tracing.Rdata')
save(src.fg2.tracing, file = 'figure.v08.07/organ_development_re_v240115/src.fg2.tracing.Rdata')
save(src.fg3.tracing, file = 'figure.v08.07/organ_development_re_v240115/src.fg3.tracing.Rdata')
save(src.fg4.tracing, file = 'figure.v08.07/organ_development_re_v240115/src.fg4.tracing.Rdata')
save(src.fg5.tracing, file = 'figure.v08.07/organ_development_re_v240115/src.fg5.tracing.Rdata')
save(src.fg6.tracing, file = 'figure.v08.07/organ_development_re_v240115/src.fg6.tracing.Rdata')
save(src.al1.tracing, file = 'figure.v08.07/organ_development_re_v240115/src.al1.tracing.Rdata')
save(src.al2.tracing, file = 'figure.v08.07/organ_development_re_v240115/src.al2.tracing.Rdata')
save(src.al3.tracing, file = 'figure.v08.07/organ_development_re_v240115/src.al3.tracing.Rdata')
save(src.mg1.tracing, file = 'figure.v08.07/organ_development_re_v240115/src.mg1.tracing.Rdata')
save(src.mg2.tracing, file = 'figure.v08.07/organ_development_re_v240115/src.mg2.tracing.Rdata')
save(src.mg3.tracing, file = 'figure.v08.07/organ_development_re_v240115/src.mg3.tracing.Rdata')
save(src.hg1.tracing, file = 'figure.v08.07/organ_development_re_v240115/src.hg1.tracing.Rdata')
save(src.hg2.tracing, file = 'figure.v08.07/organ_development_re_v240115/src.hg2.tracing.Rdata')
#-----------------------

# Batch Tracing Process
for(endoderm in endoderm_list){
  seurat = batch_process_tracing(endoderm)
  
  endoderm_names = gsub("\\.","",tolower(endoderm))
  assign(paste("src.",endoderm_names,".tracing",sep=""),
         seurat)
}

# Batch Tracing Process
# Correction !!!
for(endoderm in  endoderm_list){
  
  endoderm_names = gsub("\\.","",tolower(endoderm))
  seurat = get(paste("src.",endoderm_names,".tracing",sep=""))
  
  lineage_list = t(values(endoderm_lineage_raw[[endoderm]], key="ss9"))
  
  
  cell_save = rownames(seurat@meta.data[seurat$lineage%in%lineage_list,])
  seurat = seurat[,colnames(seurat)%in%cell_save]
  assign(paste("src.",endoderm_names,".tracing",sep=""),
         seurat)
  
  
  seurat = batch_process_tracing(endoderm)
  assign(paste("src.",endoderm_names,".tracing",sep=""),
         seurat)
}

#--  Plot Tracing
for(endoderm in endoderm_list){
  endoderm_names = gsub("\\.","",tolower(endoderm))
  seurat = get(paste("src.",endoderm_names,".tracing",sep=""))
  
  pdf(paste("figure.v08.07/organ_development_re_v240115/pdf/",
            endoderm,"tracing_embedding.pdf", sep=""))
  
  for(reduction in c(names(src.fg1.tracing@reductions))){
    
    data_plot = reduct_to_meta(seurat, reduction = reduction)
    
    g1 = ggplot()+ 
      geom_point(data = data_plot,
                 mapping = aes(x= Coord_1, y = Coord_2,
                               color = Time), size=1.5)+
      scale_color_manual(values = c(color.lineage, colors.time.2,
                                    cluster.endoderm.color.v5))+
      theme_void() + p_add + 
      ggtitle(paste(endoderm,reduction,"Time",sep="-")) +
      theme(legend.position = "none")
    
    g2 = ggplot()+ 
      geom_point(data = data_plot,
                 mapping = aes(x= Coord_1, y = Coord_2,
                               color = cluster.v06.26.re_mnn_umap_fta), size=1.5)+
      scale_color_manual(values = c(color.lineage, colors.time.2,
                                    cluster.endoderm.color.v5))+
      theme_void() + p_add + 
      ggtitle(paste(endoderm,reduction,"Time",sep="-")) +
      theme(legend.position = "none")
    
    g3 = ggplot()+ 
      geom_point(data = data_plot,
                 mapping = aes(x= Coord_1, y = Coord_2,
                               color = cluster.extract.v1.1_define), size=1.5)+
      scale_color_manual(values = c(color.lineage, colors.time.2,
                                    cluster.endoderm.color.v5))+
      theme_void() + p_add + 
      ggtitle(paste(endoderm,reduction,"Time",sep="-")) +
      theme(legend.position = "none")
    
    g4 = ggplot()+ 
      geom_point(data = data_plot,
                 mapping = aes(x= Coord_1, y = Coord_2,
                               color = lineage), size=1.5)+
      scale_color_manual(values = c(color.lineage, colors.time.2,
                                    cluster.endoderm.color.v5))+
      theme_void() + p_add + 
      ggtitle(paste(endoderm,reduction,"Time",sep="-")) +
      theme(legend.position = "none")
    
    print(g1);print(g2);print(g3);print(g4)
    rm(g1,g2,g3,g4)
    
  }
  
  dev.off()
}






















