#===============================================================================
#>>> 9.signaling pathways analysis 
#===============================================================================

#===============================================================================
#>>> 9.2 Endoderm trajectory
#===============================================================================

#-- set :: cell-order
#--------------------------------------
for(i.cell.order in c("FG.1")){
  dim.reduct = "umap"
  type_define = "cluster.v06.26.re_correct"
  cellorder_src.fg1.integrated.fg1 = 
    rownames(src.fg1.integrated@meta.data[unlist(src.fg1.integrated[[type_define]]) %in%c("FG.1","Pharynx.organ.2"),])
  coord_src.fg1.integrated_fg1 = cbind(
    -src.fg1.integrated[[dim.reduct]]@cell.embeddings[cell_src.fg1.integrated_fg1,1],
    src.fg1.integrated[[dim.reduct]]@cell.embeddings[cell_src.fg1.integrated_fg1,2])
  pcurve_src.fg1.integrated_fg1 = princurve::principal_curve(x = coord_src.fg1.integrated_fg1, smoother = "smooth.spline")
  src.fg1.integrated$lambda_fg1 = pcurve_src.fg1.integrated_fg1$lambda[cell_src.fg1.integrated_fg1]
  src.fg1.integrated@meta.data[cell_src.fg1.integrated_fg1,]$lambda_fg1 = 
    pcurve_src.fg1.integrated_fg1$lambda[cell_src.fg1.integrated_fg1]
  src.fg1.integrated@meta.data[!src.fg1.integrated$lambda_fg1%in%NA,]$lambda_fg1 = 
    norm_range(src.fg1.integrated@meta.data[!src.fg1.integrated$lambda_fg1%in%NA,]$lambda_fg1)
  
  cellorder_src.fg1.integrated.fg1 = c(
    intersect(names(src.fg1.integrated$lambda_fg1[order(src.fg1.integrated$lambda_fg1)]),
              rownames(src.fg1.integrated@meta.data[src.fg1.integrated$cluster.v06.26.re_correct%in%c("FG.1","Pharynx.organ.2"),])))
}

for(i.cell.order in c("FG.2")){
  dim.reduct = "umap"
  type_define = "cluster.v06.26.re_hc"
  cell_src.fg2.integrated_re_fg2.pha.1 = 
    rownames(src.fg2.integrated_re@meta.data[unlist(src.fg2.integrated_re[[type_define]]) %in%c("FG.2","Pharynx.organ.1"),])
  coord_src.fg2.integrated_re_fg2.pha.1 = cbind(
    -src.fg2.integrated_re[[dim.reduct]]@cell.embeddings[cell_src.fg2.integrated_re_fg2.pha.1,1],
    src.fg2.integrated_re[[dim.reduct]]@cell.embeddings[cell_src.fg2.integrated_re_fg2.pha.1,2])
  pcurve_src.fg2.integrated_re_fg2.pha.1 = princurve::principal_curve(x = coord_src.fg2.integrated_re_fg2.pha.1, smoother = "smooth.spline")
  src.fg2.integrated_re$lambda_fg2.pha.1 = pcurve_src.fg2.integrated_re_fg2.pha.1$lambda[cell_src.fg2.integrated_re_fg2.pha.1]
  src.fg2.integrated_re@meta.data[cell_src.fg2.integrated_re_fg2.pha.1,]$lambda_fg2.pha.1 = 
    pcurve_src.fg2.integrated_re_fg2.pha.1$lambda[cell_src.fg2.integrated_re_fg2.pha.1]
  src.fg2.integrated_re@meta.data[!src.fg2.integrated_re$lambda_fg2.pha.1%in%NA,]$lambda_fg2.pha.1 = 
    norm_range(src.fg2.integrated_re@meta.data[!src.fg2.integrated_re$lambda_fg2.pha.1%in%NA,]$lambda_fg2.pha.1)
  
  cell_src.fg2.integrated_re_fg2.pha.1 = c(
    intersect(names(src.fg2.integrated_re$lambda_fg2.pha.1[order(src.fg2.integrated_re$lambda_fg2.pha.1)]),
              rownames(src.fg2.integrated_re@meta.data[src.fg2.integrated_re$cluster.v06.26.re_hc%in%c("FG.2","Pharynx.organ.1"),])))
}

for(i.cell.order in c("FG.3")){
  src.fg3.integrated_re$cluster.v06.26.re_correct_refine.re = 
    src.fg3.integrated_re$cluster.v06.26.re_correct_refine

  dim.reduct = "umap.rotated.3d"
  DimPlot(src.fg3.integrated_re, reduction = dim.reduct, pt.size = 1.5,
          group.by = "cluster.v06.26.re_correct_refine.re",
          cols = cluster.endoderm.color.v5)
  DimPlot(src.fg3.integrated_re, reduction = dim.reduct, pt.size = 1.5,
          group.by = "Time", cols = colors.time)
  
  type_define = "cluster.v06.26.re_correct_refine.re"
  #---------------------
  cell_src.fg3.integrated_re_fg3.lung = 
    rownames(src.fg3.integrated_re@meta.data[unlist(src.fg3.integrated_re[[type_define]]) %in%c("FG.3","Lung"),])
  coord_src.fg3.integrated_re_fg3.lung = src.fg3.integrated_re[["umap_mon3_2d"]]@cell.embeddings[cell_src.fg3.integrated_re_fg3.lung,c(1,2)]
  pcurve_src.fg3.integrated_re_fg3.lung = princurve::principal_curve(x = coord_src.fg3.integrated_re_fg3.lung, smoother = "smooth.spline")
  src.fg3.integrated_re$lambda_fg3.lung = pcurve_src.fg3.integrated_re_fg3.lung$lambda[cell_src.fg3.integrated_re_fg3.lung]
  src.fg3.integrated_re@meta.data[cell_src.fg3.integrated_re_fg3.lung,]$lambda_fg3.lung = 
    pcurve_src.fg3.integrated_re_fg3.lung$lambda[cell_src.fg3.integrated_re_fg3.lung]
  src.fg3.integrated_re@meta.data[!src.fg3.integrated_re$lambda_fg3.lung%in%NA,]$lambda_fg3.lung = 
    norm_range(src.fg3.integrated_re@meta.data[!src.fg3.integrated_re$lambda_fg3.lung%in%NA,]$lambda_fg3.lung)
  
  cell_src.fg3.integrated_re_fg3.pha4 = 
    rownames(src.fg3.integrated_re@meta.data[unlist(src.fg3.integrated_re[[type_define]]) %in%c("FG.3","Pharynx.organ.4"),])
  coord_src.fg3.integrated_re_fg3.pha4 = src.fg3.integrated_re[["umap_integrated"]]@cell.embeddings[cell_src.fg3.integrated_re_fg3.pha4,c(2,1)]
  pcurve_src.fg3.integrated_re_fg3.pha4 = princurve::principal_curve(x = coord_src.fg3.integrated_re_fg3.pha4, smoother = "smooth.spline")
  src.fg3.integrated_re$lambda_fg3.pha4 = pcurve_src.fg3.integrated_re_fg3.pha4$lambda[cell_src.fg3.integrated_re_fg3.pha4]
  src.fg3.integrated_re@meta.data[cell_src.fg3.integrated_re_fg3.pha4,]$lambda_fg3.pha4 = 
    pcurve_src.fg3.integrated_re_fg3.pha4$lambda[cell_src.fg3.integrated_re_fg3.pha4]
  src.fg3.integrated_re@meta.data[!src.fg3.integrated_re$lambda_fg3.pha4%in%NA,]$lambda_fg3.pha4 = 
    norm_range(src.fg3.integrated_re@meta.data[!src.fg3.integrated_re$lambda_fg3.pha4%in%NA,]$lambda_fg3.pha4)
  
  cell_src.fg3.integrated_re_fg3.pha5 = 
    rownames(src.fg3.integrated_re@meta.data[unlist(src.fg3.integrated_re[[type_define]]) %in%c("FG.3","Pharynx.organ.5"),])
  coord_src.fg3.integrated_re_fg3.pha5 = src.fg3.integrated_re[["umap"]]@cell.embeddings[cell_src.fg3.integrated_re_fg3.pha5,c(2,1)]
  pcurve_src.fg3.integrated_re_fg3.pha5 = princurve::principal_curve(x = coord_src.fg3.integrated_re_fg3.pha5, smoother = "smooth.spline")
  src.fg3.integrated_re$lambda_fg3.pha5 = pcurve_src.fg3.integrated_re_fg3.pha5$lambda[cell_src.fg3.integrated_re_fg3.pha5]
  src.fg3.integrated_re@meta.data[cell_src.fg3.integrated_re_fg3.pha5,]$lambda_fg3.pha5 = 
    pcurve_src.fg3.integrated_re_fg3.pha5$lambda[cell_src.fg3.integrated_re_fg3.pha5]
  src.fg3.integrated_re@meta.data[!src.fg3.integrated_re$lambda_fg3.pha5%in%NA,]$lambda_fg3.pha5 = 
    norm_range(src.fg3.integrated_re@meta.data[!src.fg3.integrated_re$lambda_fg3.pha5%in%NA,]$lambda_fg3.pha5)
  #---------------------
  
  cell_src.fg3.integrated_re_fg3 = intersect(
    names(src.fg3.integrated_re$lambda_fg3.lung[order(src.fg3.integrated_re$lambda_fg3.lung)]), 
    rownames(src.fg3.integrated_re@meta.data[unlist(src.fg3.integrated_re[[type_define]]) %in%"FG.3",]))
  cell_src.fg3.integrated_re_lung = intersect(
    names(src.fg3.integrated_re$lambda_fg3.lung[order(src.fg3.integrated_re$lambda_fg3.lung)]), 
    rownames(src.fg3.integrated_re@meta.data[unlist(src.fg3.integrated_re[[type_define]]) %in%"Lung",]))
  cell_src.fg3.integrated_re_pha4 = rev(intersect(
    names(src.fg3.integrated_re$lambda_fg3.pha4[order(src.fg3.integrated_re$lambda_fg3.pha4)]), 
    rownames(src.fg3.integrated_re@meta.data[unlist(src.fg3.integrated_re[[type_define]]) %in%"Pharynx.organ.4",])))
  cell_src.fg3.integrated_re_pha5 = intersect(
    names(src.fg3.integrated_re$lambda_fg3.pha5[order(src.fg3.integrated_re$lambda_fg3.pha5)]), 
    rownames(src.fg3.integrated_re@meta.data[unlist(src.fg3.integrated_re[[type_define]]) %in%"Pharynx.organ.5",]))
  
  cell_src.fg3.integrated_re_fg3.fin = c(cell_src.fg3.integrated_re_fg3, cell_src.fg3.integrated_re_pha4,
                                         cell_src.fg3.integrated_re_pha5, cell_src.fg3.integrated_re_lung) 
  cell_src.fg3.integrated_re_fg3.lung.fin = c(cell_src.fg3.integrated_re_fg3, cell_src.fg3.integrated_re_lung) 
  cell_src.fg3.integrated_re_fg3.pha5.fin = c(cell_src.fg3.integrated_re_fg3, cell_src.fg3.integrated_re_pha5) 
  
}

for(i.cell.order in c("FG.4")){
  dim.reduct = "umap_mnn"
  DimPlot(src.fg4.integrated_refine, reduction = dim.reduct, pt.size = 1.5,
          group.by = "cluster.v06.26.re_correct_refine",
          cols = cluster.endoderm.color.v5)
  DimPlot(src.fg4.integrated_refine, reduction = dim.reduct, pt.size = 1.5,
          group.by = "Time", cols = colors.time)
  
  type_define = "cluster.v06.26.re_correct_refine"
  #---------------------
  cell_src.fg4.integrated_refine_fg4.lung = 
    rownames(src.fg4.integrated_refine@meta.data[unlist(src.fg4.integrated_refine[[type_define]]) %in%c("FG.4","FG.4-Lung/Stomach","Lung"),])
  coord_src.fg4.integrated_refine_fg4.lung = src.fg4.integrated_refine[["umap_mnn"]]@cell.embeddings[cell_src.fg4.integrated_refine_fg4.lung,c(1,2)]
  pcurve_src.fg4.integrated_refine_fg4.lung = princurve::principal_curve(x = coord_src.fg4.integrated_refine_fg4.lung, smoother = "smooth.spline")
  src.fg4.integrated_refine$lambda_fg4.lung = pcurve_src.fg4.integrated_refine_fg4.lung$lambda[cell_src.fg4.integrated_refine_fg4.lung]
  src.fg4.integrated_refine@meta.data[cell_src.fg4.integrated_refine_fg4.lung,]$lambda_fg4.lung = 
    pcurve_src.fg4.integrated_refine_fg4.lung$lambda[cell_src.fg4.integrated_refine_fg4.lung]
  src.fg4.integrated_refine@meta.data[!src.fg4.integrated_refine$lambda_fg4.lung%in%NA,]$lambda_fg4.lung = 
    norm_range(src.fg4.integrated_refine@meta.data[!src.fg4.integrated_refine$lambda_fg4.lung%in%NA,]$lambda_fg4.lung)
  
  cell_src.fg4.integrated_refine_fg4.sto = 
    rownames(src.fg4.integrated_refine@meta.data[unlist(src.fg4.integrated_refine[[type_define]]) %in%c("FG.4","FG.4-Lung/Stomach","Stomach"),])
  coord_src.fg4.integrated_refine_fg4.sto = src.fg4.integrated_refine[["umap_mnn"]]@cell.embeddings[cell_src.fg4.integrated_refine_fg4.sto,c(1,2)]
  pcurve_src.fg4.integrated_refine_fg4.sto = princurve::principal_curve(x = coord_src.fg4.integrated_refine_fg4.sto, smoother = "smooth.spline")
  src.fg4.integrated_refine$lambda_fg4.sto = pcurve_src.fg4.integrated_refine_fg4.sto$lambda[cell_src.fg4.integrated_refine_fg4.sto]
  src.fg4.integrated_refine@meta.data[cell_src.fg4.integrated_refine_fg4.sto,]$lambda_fg4.sto = 
    pcurve_src.fg4.integrated_refine_fg4.sto$lambda[cell_src.fg4.integrated_refine_fg4.sto]
  src.fg4.integrated_refine@meta.data[!src.fg4.integrated_refine$lambda_fg4.sto%in%NA,]$lambda_fg4.sto = 
    norm_range(src.fg4.integrated_refine@meta.data[!src.fg4.integrated_refine$lambda_fg4.sto%in%NA,]$lambda_fg4.sto)
  
  cell_src.fg4.integrated_refine_fg4.liv = 
    rownames(src.fg4.integrated_refine@meta.data[unlist(src.fg4.integrated_refine[[type_define]]) %in%c("FG.4","FG.4-Liver","AL.1/2-Liver","Liver"),])
  coord_src.fg4.integrated_refine_fg4.liv = src.fg4.integrated_refine[["umap_mnn"]]@cell.embeddings[cell_src.fg4.integrated_refine_fg4.liv,c(1,2)]
  pcurve_src.fg4.integrated_refine_fg4.liv = princurve::principal_curve(x = coord_src.fg4.integrated_refine_fg4.liv, smoother = "smooth.spline")
  src.fg4.integrated_refine$lambda_fg4.liv = pcurve_src.fg4.integrated_refine_fg4.liv$lambda[cell_src.fg4.integrated_refine_fg4.liv]
  src.fg4.integrated_refine@meta.data[cell_src.fg4.integrated_refine_fg4.liv,]$lambda_fg4.liv = 
    pcurve_src.fg4.integrated_refine_fg4.liv$lambda[cell_src.fg4.integrated_refine_fg4.liv]
  src.fg4.integrated_refine@meta.data[!src.fg4.integrated_refine$lambda_fg4.liv%in%NA,]$lambda_fg4.liv = 
    norm_range(src.fg4.integrated_refine@meta.data[!src.fg4.integrated_refine$lambda_fg4.liv%in%NA,]$lambda_fg4.liv)
  
  cell_src.fg4.integrated_refine_fg4.pha5 = 
    rownames(src.fg4.integrated_refine@meta.data[unlist(src.fg4.integrated_refine[[type_define]]) %in%c("FG.4","Pharynx.organ.5"),])
  coord_src.fg4.integrated_refine_fg4.pha5 = src.fg4.integrated_refine[["umap_mnn"]]@cell.embeddings[cell_src.fg4.integrated_refine_fg4.pha5,c(1,2)]
  pcurve_src.fg4.integrated_refine_fg4.pha5 = princurve::principal_curve(x = coord_src.fg4.integrated_refine_fg4.pha5, smoother = "smooth.spline")
  src.fg4.integrated_refine$lambda_fg4.pha5 = pcurve_src.fg4.integrated_refine_fg4.pha5$lambda[cell_src.fg4.integrated_refine_fg4.pha5]
  src.fg4.integrated_refine@meta.data[cell_src.fg4.integrated_refine_fg4.pha5,]$lambda_fg4.pha5 = 
    pcurve_src.fg4.integrated_refine_fg4.pha5$lambda[cell_src.fg4.integrated_refine_fg4.pha5]
  src.fg4.integrated_refine@meta.data[!src.fg4.integrated_refine$lambda_fg4.pha5%in%NA,]$lambda_fg4.pha5 = 
    norm_range(src.fg4.integrated_refine@meta.data[!src.fg4.integrated_refine$lambda_fg4.pha5%in%NA,]$lambda_fg4.pha5)
  #---------------------
  
  cell.list = cell_src.fg4.integrated_refine_fg4.lung
  plot(x = c(1:length(cell.list), 1:length(cell.list)),
       y = c(1:length(cell.list), rep(1000, length(cell.list))),
       col = c(cluster.endoderm.color.v5[src.fg4.integrated_refine@meta.data[cell.list,]$cluster.v06.26.re_correct_refine],
               colors.time[src.fg4.integrated_refine@meta.data[cell.list,]$Time]))
  
  
  cell_src.fg4.integrated_refine_fg4 = intersect(
    names(src.fg4.integrated_refine$lambda_fg4.sto[order(src.fg4.integrated_refine$lambda_fg4.sto)]), 
    rownames(src.fg4.integrated_refine@meta.data[unlist(src.fg4.integrated_refine[[type_define]]) %in%"FG.4",]))
  cell_src.fg4.integrated_refine_fg4tolusto = intersect(
    names(src.fg4.integrated_refine$lambda_fg4.sto[order(src.fg4.integrated_refine$lambda_fg4.sto)]), 
    rownames(src.fg4.integrated_refine@meta.data[unlist(src.fg4.integrated_refine[[type_define]]) %in%"FG.4-Lung/Stomach",]))
  cell_src.fg4.integrated_refine_sto = intersect(
    names(src.fg4.integrated_refine$lambda_fg4.sto[order(src.fg4.integrated_refine$lambda_fg4.sto)]), 
    rownames(src.fg4.integrated_refine@meta.data[unlist(src.fg4.integrated_refine[[type_define]]) %in%"Stomach",]))
  cell_src.fg4.integrated_refine_lung = intersect(
    names(src.fg4.integrated_refine$lambda_fg4.lung[order(src.fg4.integrated_refine$lambda_fg4.lung)]), 
    rownames(src.fg4.integrated_refine@meta.data[unlist(src.fg4.integrated_refine[[type_define]]) %in%"Lung",]))
  cell_src.fg4.integrated_refine_pha5 = intersect(
    names(src.fg4.integrated_refine$lambda_fg4.pha5[order(src.fg4.integrated_refine$lambda_fg4.pha5)]), 
    rownames(src.fg4.integrated_refine@meta.data[unlist(src.fg4.integrated_refine[[type_define]]) %in%"Pharynx.organ.5",]))
  cell_src.fg4.integrated_refine_fg4toliv = intersect(
    names(src.fg4.integrated_refine$lambda_fg4.liv[order(src.fg4.integrated_refine$lambda_fg4.liv)]), 
    rownames(src.fg4.integrated_refine@meta.data[unlist(src.fg4.integrated_refine[[type_define]]) %in%"FG.4-Liver",]))
  cell_src.fg4.integrated_refine_al12toliv = intersect(
    names(src.fg4.integrated_refine$lambda_fg4.liv[order(src.fg4.integrated_refine$lambda_fg4.liv)]), 
    rownames(src.fg4.integrated_refine@meta.data[unlist(src.fg4.integrated_refine[[type_define]]) %in%"AL.1/2-Liver",]))
  cell_src.fg4.integrated_refine_liv = intersect(
    names(src.fg4.integrated_refine$lambda_fg4.liv[order(src.fg4.integrated_refine$lambda_fg4.liv)]), 
    rownames(src.fg4.integrated_refine@meta.data[unlist(src.fg4.integrated_refine[[type_define]]) %in%"Liver",]))
  
  cell_src.fg4.integrated_refine_fg4.fin = c(cell_src.fg4.integrated_refine_fg4, cell_src.fg4.integrated_refine_fg4tolusto, cell_src.fg4.integrated_refine_fg4toliv)
  cell_src.fg4.integrated_refine_fg4.pha5.fin = c(cell_src.fg4.integrated_refine_fg4, cell_src.fg4.integrated_refine_pha5) 
  cell_src.fg4.integrated_refine_fg4.lung.fin = c(cell_src.fg4.integrated_refine_fg4, cell_src.fg4.integrated_refine_fg4tolusto, cell_src.fg4.integrated_refine_lung)
  cell_src.fg4.integrated_refine_fg4.sto.fin = c(cell_src.fg4.integrated_refine_fg4, cell_src.fg4.integrated_refine_fg4tolusto, cell_src.fg4.integrated_refine_sto)
  cell_src.fg4.integrated_refine_fg4.liv.fin = c(cell_src.fg4.integrated_refine_fg4, cell_src.fg4.integrated_refine_fg4toliv, 
                                                 cell_src.fg4.integrated_refine_al12toliv, cell_src.fg4.integrated_refine_liv)
  cell_src.fg4.integrated_refine_fg4.lunsto.fin = c(cell_src.fg4.integrated_refine_fg4, cell_src.fg4.integrated_refine_fg4tolusto,
                                                    cell_src.fg4.integrated_refine_lung, cell_src.fg4.integrated_refine_sto)
}

for(i.cell.order in c("FG.5")){
  dim.reduct = "umap_integrated"
  DimPlot(src.fg5.integrated, reduction = dim.reduct, pt.size = 1.5,
          group.by = "cluster.v06.26.re_correct",
          cols = cluster.endoderm.color.v5)
  DimPlot(src.fg5.integrated, reduction = dim.reduct, pt.size = 1.5,
          group.by = "Time", cols = colors.time)
  DimPlot(src.fg5.integrated, reduction = dim.reduct, pt.size = 1.5,
          group.by = "batch")
  
  type_define = "cluster.v06.26.re_correct"
  #---------------------
  cell_src.fg5.integrated_fg5.th = 
    rownames(src.fg5.integrated@meta.data[unlist(src.fg5.integrated[[type_define]]) %in%c("FG.5","Thyroid"),])
  coord_src.fg5.integrated_fg5.th = cbind(src.fg5.integrated[[dim.reduct]]@cell.embeddings[cell_src.fg5.integrated_fg5.th,1],
                                          src.fg5.integrated[[dim.reduct]]@cell.embeddings[cell_src.fg5.integrated_fg5.th,2])
  pcurve_src.fg5.integrated_fg5.th = princurve::principal_curve(x = coord_src.fg5.integrated_fg5.th, smoother = "smooth.spline")
  src.fg5.integrated$lambda_fg5.th = pcurve_src.fg5.integrated_fg5.th$lambda[cell_src.fg5.integrated_fg5.th]
  src.fg5.integrated@meta.data[cell_src.fg5.integrated_fg5.th,]$lambda_fg5.th = 
    pcurve_src.fg5.integrated_fg5.th$lambda[cell_src.fg5.integrated_fg5.th]
  src.fg5.integrated@meta.data[!src.fg5.integrated$lambda_fg5.th%in%NA,]$lambda_fg5.th = 
    norm_range(src.fg5.integrated@meta.data[!src.fg5.integrated$lambda_fg5.th%in%NA,]$lambda_fg5.th)
  
  cell_src.fg5.integrated_fg5.pha.3 = 
    rownames(src.fg5.integrated@meta.data[unlist(src.fg5.integrated[[type_define]]) %in%c("FG.5","Pharynx.organ.3"),])
  coord_src.fg5.integrated_fg5.pha.3 = cbind(src.fg5.integrated[[dim.reduct]]@cell.embeddings[cell_src.fg5.integrated_fg5.pha.3,1],
                                             src.fg5.integrated[[dim.reduct]]@cell.embeddings[cell_src.fg5.integrated_fg5.pha.3,2])
  pcurve_src.fg5.integrated_fg5.pha.3 = princurve::principal_curve(x = coord_src.fg5.integrated_fg5.pha.3, smoother = "smooth.spline")
  src.fg5.integrated$lambda_fg5.pha.3 = pcurve_src.fg5.integrated_fg5.pha.3$lambda[cell_src.fg5.integrated_fg5.pha.3]
  src.fg5.integrated@meta.data[cell_src.fg5.integrated_fg5.pha.3,]$lambda_fg5.pha.3 = 
    pcurve_src.fg5.integrated_fg5.pha.3$lambda[cell_src.fg5.integrated_fg5.pha.3]
  src.fg5.integrated@meta.data[!src.fg5.integrated$lambda_fg5.pha.3%in%NA,]$lambda_fg5.pha.3 = 
    norm_range(src.fg5.integrated@meta.data[!src.fg5.integrated$lambda_fg5.pha.3%in%NA,]$lambda_fg5.pha.3)
  #---------------------
  
  cell_src.fg5.integrated_fg5.th = c(
    intersect(names(src.fg5.integrated$lambda_fg5.th[order(src.fg5.integrated$lambda_fg5.th)]),
              rownames(src.fg5.integrated@meta.data[src.fg5.integrated$cluster.v06.26.re_correct%in%c("FG.5"),])),
    intersect(names(src.fg5.integrated$lambda_fg5.th[order(src.fg5.integrated$lambda_fg5.th)]),
              rownames(src.fg5.integrated@meta.data[src.fg5.integrated$cluster.v06.26.re_correct%in%c("Thyroid"),])))
  
  cell_src.fg5.integrated_fg5.pha.3 = c(
    intersect(names(src.fg5.integrated$lambda_fg5.pha.3[order(src.fg5.integrated$lambda_fg5.pha.3)]),
              rownames(src.fg5.integrated@meta.data[src.fg5.integrated$cluster.v06.26.re_correct%in%c("FG.5"),])),
    intersect(names(src.fg5.integrated$lambda_fg5.pha.3[order(src.fg5.integrated$lambda_fg5.pha.3)]),
              rownames(src.fg5.integrated@meta.data[src.fg5.integrated$cluster.v06.26.re_correct%in%c("Pharynx.organ.3"),])))
}

for(i.cell.order in c("AL.1")){
  dim.reduct = "umap_integrated"
  DimPlot(src.al1.integrated.re, reduction = dim.reduct, pt.size = 1.5,
          group.by = "cluster.v06.26.re_correct_re",
          cols = cluster.endoderm.color.v5)
  DimPlot(src.al1.integrated.re, reduction = dim.reduct, pt.size = 1.5,
          group.by = "Time", cols = colors.time)
  DimPlot(src.al1.integrated.re, reduction = dim.reduct, pt.size = 1.5,
          group.by = "batch")
  
  type_define = "cluster.v06.26.re_correct_re"
  #---------------------
  cell_src.al1.integrated.re_al1.liv = 
    rownames(src.al1.integrated.re@meta.data[unlist(src.al1.integrated.re[[type_define]]) %in%c("AL.1","AL.1/2-Liver","Liver"),])
  coord_src.al1.integrated.re_al1.liv = cbind(-src.al1.integrated.re[[dim.reduct]]@cell.embeddings[cell_src.al1.integrated.re_al1.liv,1],
                                              src.al1.integrated.re[[dim.reduct]]@cell.embeddings[cell_src.al1.integrated.re_al1.liv,2])
  pcurve_src.al1.integrated.re_al1.liv = princurve::principal_curve(x = coord_src.al1.integrated.re_al1.liv, smoother = "smooth.spline")
  src.al1.integrated.re$lambda_al1.liv = pcurve_src.al1.integrated.re_al1.liv$lambda[cell_src.al1.integrated.re_al1.liv]
  src.al1.integrated.re@meta.data[cell_src.al1.integrated.re_al1.liv,]$lambda_al1.liv = 
    pcurve_src.al1.integrated.re_al1.liv$lambda[cell_src.al1.integrated.re_al1.liv]
  src.al1.integrated.re@meta.data[!src.al1.integrated.re$lambda_al1.liv%in%NA,]$lambda_al1.liv = 
    norm_range(src.al1.integrated.re@meta.data[!src.al1.integrated.re$lambda_al1.liv%in%NA,]$lambda_al1.liv)
  #---------------------
  
  cell_src.al1.integrated.re_al1.liv.step1 = c(
    intersect(names(src.al1.integrated.re$lambda_al1.liv[order(src.al1.integrated.re$lambda_al1.liv)]),
              rownames(src.al1.integrated.re@meta.data[src.al1.integrated.re$cluster.v06.26.re_correct_re%in%c("AL.1","AL.1/2-Liver"),])))
  cell_src.al1.integrated.re_al1.liv.step2 = c(
    intersect(names(src.al1.integrated.re$lambda_al1.liv[order(src.al1.integrated.re$lambda_al1.liv)]),
              rownames(src.al1.integrated.re@meta.data[src.al1.integrated.re$cluster.v06.26.re_correct_re%in%c("AL.1/2-Liver","Liver"),])))
}

for(i.cell.order in c("AL.2")){
  dim.reduct = "umap_integrated"
  DimPlot(src.al2.integrated.re, reduction = dim.reduct, pt.size = 1.5,
          group.by = "cluster.v06.26.re_correct_re",
          cols = cluster.endoderm.color.v5)
  DimPlot(src.al2.integrated.re, reduction = dim.reduct, pt.size = 1.5,
          group.by = "Time", cols = colors.time)
  DimPlot(src.al2.integrated.re, reduction = dim.reduct, pt.size = 1.5,
          group.by = "batch")
  
  type_define = "cluster.v06.26.re_correct_re"
  #---------------------
  cell_src.al2.integrated.re_al2.liv = 
    rownames(src.al2.integrated.re@meta.data[unlist(src.al2.integrated.re[[type_define]]) %in%c("AL.2","AL.1/2-Liver","Liver"),])
  coord_src.al2.integrated.re_al2.liv = cbind(-src.al2.integrated.re[[dim.reduct]]@cell.embeddings[cell_src.al2.integrated.re_al2.liv,1],
                                              src.al2.integrated.re[[dim.reduct]]@cell.embeddings[cell_src.al2.integrated.re_al2.liv,2])
  pcurve_src.al2.integrated.re_al2.liv = princurve::principal_curve(x = coord_src.al2.integrated.re_al2.liv, smoother = "smooth.spline")
  src.al2.integrated.re$lambda_al2.liv = pcurve_src.al2.integrated.re_al2.liv$lambda[cell_src.al2.integrated.re_al2.liv]
  src.al2.integrated.re@meta.data[cell_src.al2.integrated.re_al2.liv,]$lambda_al2.liv = 
    pcurve_src.al2.integrated.re_al2.liv$lambda[cell_src.al2.integrated.re_al2.liv]
  src.al2.integrated.re@meta.data[!src.al2.integrated.re$lambda_al2.liv%in%NA,]$lambda_al2.liv = 
    norm_range(src.al2.integrated.re@meta.data[!src.al2.integrated.re$lambda_al2.liv%in%NA,]$lambda_al2.liv)
  #---------------------
  cell_src.al2.integrated.re_al2.liv = c(
    intersect(names(src.al2.integrated.re$lambda_al2.liv[order(src.al2.integrated.re$lambda_al2.liv)]),
              rownames(src.al2.integrated.re@meta.data[src.al2.integrated.re$cluster.v06.26.re_correct_re%in%c("AL.2","AL.1/2-Liver","Liver"),])))
  
  cell_src.al2.integrated.re_al2.liv.step1 = c(
    intersect(names(src.al2.integrated.re$lambda_al2.liv[order(src.al2.integrated.re$lambda_al2.liv)]),
              rownames(src.al2.integrated.re@meta.data[src.al2.integrated.re$cluster.v06.26.re_correct_re%in%c("AL.2","AL.1/2-Liver"),])))
  cell_src.al2.integrated.re_al2.liv.step2 = c(
    intersect(names(src.al2.integrated.re$lambda_al2.liv[order(src.al2.integrated.re$lambda_al2.liv)]),
              rownames(src.al2.integrated.re@meta.data[src.al2.integrated.re$cluster.v06.26.re_correct_re%in%c("AL.1/2-Liver","Liver"),])))
  
}

for(i.cell.order in c("AL.3")){
  dim.reduct = "umap_integrated"
  DimPlot(src.al3.integrated.re, reduction = dim.reduct, pt.size = 1.5,
          group.by = "cluster.v06.26.re_correct",
          cols = cluster.endoderm.color.v5)
  DimPlot(src.al3.integrated.re, reduction = dim.reduct, pt.size = 1.5,
          group.by = "Time", cols = colors.time)
  DimPlot(src.al3.integrated.re, reduction = dim.reduct, pt.size = 1.5,
          group.by = "batch")
  
  type_define = "cluster.v06.26.re_correct"
  #---------------------
  cell_src.al3.integrated.re_al3.si1 = 
    rownames(src.al3.integrated.re@meta.data[unlist(src.al3.integrated.re[[type_define]]) %in%c("AL.3","AL.3-Small.intestine.1","Small.intestine.1"),])
  coord_src.al3.integrated.re_al3.si1 = cbind(src.al3.integrated.re[[dim.reduct]]@cell.embeddings[cell_src.al3.integrated.re_al3.si1,2],
                                              src.al3.integrated.re[[dim.reduct]]@cell.embeddings[cell_src.al3.integrated.re_al3.si1,1])
  pcurve_src.al3.integrated.re_al3.si1 = princurve::principal_curve(x = coord_src.al3.integrated.re_al3.si1, smoother = "smooth.spline")
  src.al3.integrated.re$lambda_al3.si1 = pcurve_src.al3.integrated.re_al3.si1$lambda[cell_src.al3.integrated.re_al3.si1]
  src.al3.integrated.re@meta.data[cell_src.al3.integrated.re_al3.si1,]$lambda_al3.si1 = 
    pcurve_src.al3.integrated.re_al3.si1$lambda[cell_src.al3.integrated.re_al3.si1]
  src.al3.integrated.re@meta.data[!src.al3.integrated.re$lambda_al3.si1%in%NA,]$lambda_al3.si1 = 
    norm_range(src.al3.integrated.re@meta.data[!src.al3.integrated.re$lambda_al3.si1%in%NA,]$lambda_al3.si1)
  
  cell_src.al3.integrated.re_al3.ehbd = 
    rownames(src.al3.integrated.re@meta.data[unlist(src.al3.integrated.re[[type_define]]) %in%c("AL.3","AL.3-EHBD/VP","EHBD"),])
  coord_src.al3.integrated.re_al3.ehbd = cbind(src.al3.integrated.re[[dim.reduct]]@cell.embeddings[cell_src.al3.integrated.re_al3.ehbd,2],
                                               src.al3.integrated.re[[dim.reduct]]@cell.embeddings[cell_src.al3.integrated.re_al3.ehbd,1])
  pcurve_src.al3.integrated.re_al3.ehbd = princurve::principal_curve(x = coord_src.al3.integrated.re_al3.ehbd, smoother = "smooth.spline")
  src.al3.integrated.re$lambda_al3.ehbd = pcurve_src.al3.integrated.re_al3.ehbd$lambda[cell_src.al3.integrated.re_al3.ehbd]
  src.al3.integrated.re@meta.data[cell_src.al3.integrated.re_al3.ehbd,]$lambda_al3.ehbd = 
    pcurve_src.al3.integrated.re_al3.ehbd$lambda[cell_src.al3.integrated.re_al3.ehbd]
  src.al3.integrated.re@meta.data[!src.al3.integrated.re$lambda_al3.ehbd%in%NA,]$lambda_al3.ehbd = 
    norm_range(src.al3.integrated.re@meta.data[!src.al3.integrated.re$lambda_al3.ehbd%in%NA,]$lambda_al3.ehbd)
  
  cell_src.al3.integrated.re_al3.vp = 
    rownames(src.al3.integrated.re@meta.data[unlist(src.al3.integrated.re[[type_define]]) %in%c("AL.3","AL.3-EHBD/VP","VP"),])
  coord_src.al3.integrated.re_al3.vp = cbind(src.al3.integrated.re[[dim.reduct]]@cell.embeddings[cell_src.al3.integrated.re_al3.vp,2],
                                             src.al3.integrated.re[[dim.reduct]]@cell.embeddings[cell_src.al3.integrated.re_al3.vp,1])
  pcurve_src.al3.integrated.re_al3.vp = princurve::principal_curve(x = coord_src.al3.integrated.re_al3.vp, smoother = "smooth.spline")
  src.al3.integrated.re$lambda_al3.vp = pcurve_src.al3.integrated.re_al3.vp$lambda[cell_src.al3.integrated.re_al3.vp]
  src.al3.integrated.re@meta.data[cell_src.al3.integrated.re_al3.vp,]$lambda_al3.vp = 
    pcurve_src.al3.integrated.re_al3.vp$lambda[cell_src.al3.integrated.re_al3.vp]
  src.al3.integrated.re@meta.data[!src.al3.integrated.re$lambda_al3.vp%in%NA,]$lambda_al3.vp = 
    norm_range(src.al3.integrated.re@meta.data[!src.al3.integrated.re$lambda_al3.vp%in%NA,]$lambda_al3.vp)
  
  cell_src.al3.integrated.re_al3.liv = 
    rownames(src.al3.integrated.re@meta.data[unlist(src.al3.integrated.re[[type_define]]) %in%c("AL.3","AL.3-Liver","AL.1/2-Liver","Liver"),])
  coord_src.al3.integrated.re_al3.liv = cbind(src.al3.integrated.re[[dim.reduct]]@cell.embeddings[cell_src.al3.integrated.re_al3.liv,2],
                                              src.al3.integrated.re[[dim.reduct]]@cell.embeddings[cell_src.al3.integrated.re_al3.liv,1])
  pcurve_src.al3.integrated.re_al3.liv = princurve::principal_curve(x = coord_src.al3.integrated.re_al3.liv, smoother = "smooth.spline")
  src.al3.integrated.re$lambda_al3.liv = pcurve_src.al3.integrated.re_al3.liv$lambda[cell_src.al3.integrated.re_al3.liv]
  src.al3.integrated.re@meta.data[cell_src.al3.integrated.re_al3.liv,]$lambda_al3.liv = 
    pcurve_src.al3.integrated.re_al3.liv$lambda[cell_src.al3.integrated.re_al3.liv]
  src.al3.integrated.re@meta.data[!src.al3.integrated.re$lambda_al3.liv%in%NA,]$lambda_al3.liv = 
    norm_range(src.al3.integrated.re@meta.data[!src.al3.integrated.re$lambda_al3.liv%in%NA,]$lambda_al3.liv)
  #---------------------
  cell_src.al3.integrated.re_al3.si1.step1 = c(
    intersect(names(src.al3.integrated.re$lambda_al3.si1[order(src.al3.integrated.re$lambda_al3.si1)]),
              rownames(src.al3.integrated.re@meta.data[src.al3.integrated.re$cluster.v06.26.re_correct%in%c("AL.3","AL.3-Small.intestine.1"),])))
  cell_src.al3.integrated.re_al3.si1.step2 = c(
    intersect(names(src.al3.integrated.re$lambda_al3.si1[order(src.al3.integrated.re$lambda_al3.si1)]),
              rownames(src.al3.integrated.re@meta.data[src.al3.integrated.re$cluster.v06.26.re_correct%in%c("AL.3-Small.intestine.1","Small.intestine.1"),])))
  
  cell_src.al3.integrated.re_al3.ehbd.step1 = c(
    intersect(names(src.al3.integrated.re$lambda_al3.ehbd[order(src.al3.integrated.re$lambda_al3.ehbd)]),
              rownames(src.al3.integrated.re@meta.data[src.al3.integrated.re$cluster.v06.26.re_correct%in%c("AL.3","AL.3-EHBD/VP"),])))
  cell_src.al3.integrated.re_al3.ehbd.step2 = c(
    intersect(names(src.al3.integrated.re$lambda_al3.ehbd[order(src.al3.integrated.re$lambda_al3.ehbd)]),
              rownames(src.al3.integrated.re@meta.data[src.al3.integrated.re$cluster.v06.26.re_correct%in%c("AL.3-EHBD/VP"),])),
    intersect(names(src.al3.integrated.re$lambda_al3.ehbd[order(src.al3.integrated.re$lambda_al3.ehbd)]),
              rownames(src.al3.integrated.re@meta.data[src.al3.integrated.re$cluster.v06.26.re_correct%in%c("EHBD"),])))
  
  cell_src.al3.integrated.re_al3.vp.step1 = c(
    intersect(names(src.al3.integrated.re$lambda_al3.vp[order(src.al3.integrated.re$lambda_al3.vp)]),
              rownames(src.al3.integrated.re@meta.data[src.al3.integrated.re$cluster.v06.26.re_correct%in%c("AL.3","AL.3-EHBD/VP"),])))
  cell_src.al3.integrated.re_al3.vp.step2 = c(
    intersect(names(src.al3.integrated.re$lambda_al3.vp[order(src.al3.integrated.re$lambda_al3.vp)]),
              rownames(src.al3.integrated.re@meta.data[src.al3.integrated.re$cluster.v06.26.re_correct%in%c("AL.3-EHBD/VP"),])),
    intersect(names(src.al3.integrated.re$lambda_al3.vp[order(src.al3.integrated.re$lambda_al3.vp)]),
              rownames(src.al3.integrated.re@meta.data[src.al3.integrated.re$cluster.v06.26.re_correct%in%c("VP"),])))
  
  cell_src.al3.integrated.re_al3.liv.step1 = c(
    intersect(names(src.al3.integrated.re$lambda_al3.liv[order(src.al3.integrated.re$lambda_al3.liv)]),
              rownames(src.al3.integrated.re@meta.data[src.al3.integrated.re$cluster.v06.26.re_correct%in%c("AL.3","AL.3-Liver"),])))
  cell_src.al3.integrated.re_al3.liv.step2 = c(
    intersect(names(src.al3.integrated.re$lambda_al3.liv[order(src.al3.integrated.re$lambda_al3.liv)]),
              rownames(src.al3.integrated.re@meta.data[src.al3.integrated.re$cluster.v06.26.re_correct%in%c("AL.3-Liver"),])),
    intersect(names(src.al3.integrated.re$lambda_al3.liv[order(src.al3.integrated.re$lambda_al3.liv)]),
              rownames(src.al3.integrated.re@meta.data[src.al3.integrated.re$cluster.v06.26.re_correct%in%c("AL.1/2-Liver"),])))
  cell_src.al3.integrated.re_al3.liv.step3 = c(
    intersect(names(src.al3.integrated.re$lambda_al3.liv[order(src.al3.integrated.re$lambda_al3.liv)]),
              rownames(src.al3.integrated.re@meta.data[src.al3.integrated.re$cluster.v06.26.re_correct%in%c("AL.1/2-Liver","Liver"),])))
}

for(i.cell.order in c("MG.1")){
  dim.reduct = "umap_mnn"
  DimPlot(src.mg1.integrated, reduction = dim.reduct, pt.size = 1.5,
          group.by = "cluster.v06.26.re",
          cols = cluster.endoderm.color.v5)
  DimPlot(src.mg1.integrated, reduction = dim.reduct, pt.size = 1.5,
          group.by = "Time", cols = colors.time)
  DimPlot(src.mg1.integrated, reduction = dim.reduct, pt.size = 1.5,
          group.by = "batch")
  
  type_define = "cluster.v06.26.re"
  #---------------------
  cell_src.mg1.integrated_mg1.sto = 
    rownames(src.mg1.integrated@meta.data[unlist(src.mg1.integrated[[type_define]]) %in%c("MG.1","Stomach"),])
  coord_src.mg1.integrated_mg1.sto = cbind(-src.mg1.integrated[[dim.reduct]]@cell.embeddings[cell_src.mg1.integrated_mg1.sto,2],
                                           src.mg1.integrated[[dim.reduct]]@cell.embeddings[cell_src.mg1.integrated_mg1.sto,1])
  pcurve_src.mg1.integrated_mg1.sto = princurve::principal_curve(x = coord_src.mg1.integrated_mg1.sto, smoother = "smooth.spline")
  src.mg1.integrated$lambda_mg1.sto = pcurve_src.mg1.integrated_mg1.sto$lambda[cell_src.mg1.integrated_mg1.sto]
  src.mg1.integrated@meta.data[cell_src.mg1.integrated_mg1.sto,]$lambda_mg1.sto = 
    pcurve_src.mg1.integrated_mg1.sto$lambda[cell_src.mg1.integrated_mg1.sto]
  src.mg1.integrated@meta.data[!src.mg1.integrated$lambda_mg1.sto%in%NA,]$lambda_mg1.sto = 
    norm_range(src.mg1.integrated@meta.data[!src.mg1.integrated$lambda_mg1.sto%in%NA,]$lambda_mg1.sto)
  
  cell_src.mg1.integrated_mg1.si1 = 
    rownames(src.mg1.integrated@meta.data[unlist(src.mg1.integrated[[type_define]]) %in%c("MG.1","Small.intestine.1"),])
  coord_src.mg1.integrated_mg1.si1 = cbind(-src.mg1.integrated[[dim.reduct]]@cell.embeddings[cell_src.mg1.integrated_mg1.si1,2],
                                           src.mg1.integrated[[dim.reduct]]@cell.embeddings[cell_src.mg1.integrated_mg1.si1,1])
  pcurve_src.mg1.integrated_mg1.si1 = princurve::principal_curve(x = coord_src.mg1.integrated_mg1.si1, smoother = "smooth.spline")
  src.mg1.integrated$lambda_mg1.si1 = pcurve_src.mg1.integrated_mg1.si1$lambda[cell_src.mg1.integrated_mg1.si1]
  src.mg1.integrated@meta.data[cell_src.mg1.integrated_mg1.si1,]$lambda_mg1.si1 = 
    pcurve_src.mg1.integrated_mg1.si1$lambda[cell_src.mg1.integrated_mg1.si1]
  src.mg1.integrated@meta.data[!src.mg1.integrated$lambda_mg1.si1%in%NA,]$lambda_mg1.si1 = 
    norm_range(src.mg1.integrated@meta.data[!src.mg1.integrated$lambda_mg1.si1%in%NA,]$lambda_mg1.si1)
  #---------------------
  
  cell_src.mg1.integrated_mg1.sto = c(
    intersect(names(src.mg1.integrated$lambda_mg1.sto[order(src.mg1.integrated$lambda_mg1.sto)]),
              rownames(src.mg1.integrated@meta.data[src.mg1.integrated$cluster.v06.26.re%in%c("MG.1"),])),
    intersect(names(src.mg1.integrated$lambda_mg1.sto[order(src.mg1.integrated$lambda_mg1.sto)]),
              rownames(src.mg1.integrated@meta.data[src.mg1.integrated$cluster.v06.26.re%in%c("Stomach"),])))
  
  cell_src.mg1.integrated_mg1.si1 = c(
    intersect(names(src.mg1.integrated$lambda_mg1.si1[order(src.mg1.integrated$lambda_mg1.si1)]),
              rownames(src.mg1.integrated@meta.data[src.mg1.integrated$cluster.v06.26.re%in%c("MG.1"),])),
    intersect(names(src.mg1.integrated$lambda_mg1.si1[order(src.mg1.integrated$lambda_mg1.si1)]),
              rownames(src.mg1.integrated@meta.data[src.mg1.integrated$cluster.v06.26.re%in%c("Small.intestine.1"),])))
  
}

for(i.cell.order in c("MG.2")){
  dim.reduct = "umap_integrated"
  DimPlot(src.mg2.integrated, reduction = dim.reduct, pt.size = 1.5,
          group.by = "cluster.v06.26.re",
          cols = cluster.endoderm.color.v5)
  DimPlot(src.mg2.integrated, reduction = dim.reduct, pt.size = 1.5,
          group.by = "Time", cols = colors.time)
  DimPlot(src.mg2.integrated, reduction = dim.reduct, pt.size = 1.5,
          group.by = "batch")
  
  type_define = "cluster.v06.26.re"
  #---------------------
  cell_src.mg2.integrated_mg2.si1 = 
    rownames(src.mg2.integrated@meta.data[unlist(src.mg2.integrated[[type_define]]) %in%c("MG.2","Small.intestine.1"),])
  coord_src.mg2.integrated_mg2.si1 = cbind(-src.mg2.integrated[[dim.reduct]]@cell.embeddings[cell_src.mg2.integrated_mg2.si1,1],
                                           src.mg2.integrated[[dim.reduct]]@cell.embeddings[cell_src.mg2.integrated_mg2.si1,2])
  pcurve_src.mg2.integrated_mg2.si1 = princurve::principal_curve(x = coord_src.mg2.integrated_mg2.si1, smoother = "smooth.spline")
  src.mg2.integrated$lambda_mg2.si1 = pcurve_src.mg2.integrated_mg2.si1$lambda[cell_src.mg2.integrated_mg2.si1]
  src.mg2.integrated@meta.data[cell_src.mg2.integrated_mg2.si1,]$lambda_mg2.si1 = 
    pcurve_src.mg2.integrated_mg2.si1$lambda[cell_src.mg2.integrated_mg2.si1]
  src.mg2.integrated@meta.data[!src.mg2.integrated$lambda_mg2.si1%in%NA,]$lambda_mg2.si1 = 
    norm_range(src.mg2.integrated@meta.data[!src.mg2.integrated$lambda_mg2.si1%in%NA,]$lambda_mg2.si1)
  
  cell_src.mg2.integrated_mg2.si2 = 
    rownames(src.mg2.integrated@meta.data[unlist(src.mg2.integrated[[type_define]]) %in%c("MG.2","Small.intestine.2"),])
  coord_src.mg2.integrated_mg2.si2 = cbind(-src.mg2.integrated[[dim.reduct]]@cell.embeddings[cell_src.mg2.integrated_mg2.si2,1],
                                           src.mg2.integrated[[dim.reduct]]@cell.embeddings[cell_src.mg2.integrated_mg2.si2,2])
  pcurve_src.mg2.integrated_mg2.si2 = princurve::principal_curve(x = coord_src.mg2.integrated_mg2.si2, smoother = "smooth.spline")
  src.mg2.integrated$lambda_mg2.si2 = pcurve_src.mg2.integrated_mg2.si2$lambda[cell_src.mg2.integrated_mg2.si2]
  src.mg2.integrated@meta.data[cell_src.mg2.integrated_mg2.si2,]$lambda_mg2.si2 = 
    pcurve_src.mg2.integrated_mg2.si2$lambda[cell_src.mg2.integrated_mg2.si2]
  src.mg2.integrated@meta.data[!src.mg2.integrated$lambda_mg2.si2%in%NA,]$lambda_mg2.si2 = 
    norm_range(src.mg2.integrated@meta.data[!src.mg2.integrated$lambda_mg2.si2%in%NA,]$lambda_mg2.si2)
  #---------------------
  
  cell_src.mg2.integrated_mg2.si1 = c(
    intersect(names(src.mg2.integrated$lambda_mg2.si1[order(src.mg2.integrated$lambda_mg2.si1)]),
              rownames(src.mg2.integrated@meta.data[src.mg2.integrated$cluster.v06.26.re%in%c("MG.2"),])),
    intersect(names(src.mg2.integrated$lambda_mg2.si1[order(src.mg2.integrated$lambda_mg2.si1)]),
              rownames(src.mg2.integrated@meta.data[src.mg2.integrated$cluster.v06.26.re%in%c("Small.intestine.1"),])))
  
  cell_src.mg2.integrated_mg2.si2 = c(
    intersect(names(src.mg2.integrated$lambda_mg2.si2[order(src.mg2.integrated$lambda_mg2.si2)]),
              rownames(src.mg2.integrated@meta.data[src.mg2.integrated$cluster.v06.26.re%in%c("MG.2"),])),
    intersect(names(src.mg2.integrated$lambda_mg2.si2[order(src.mg2.integrated$lambda_mg2.si2)]),
              rownames(src.mg2.integrated@meta.data[src.mg2.integrated$cluster.v06.26.re%in%c("Small.intestine.2"),])))
  
}

for(i.cell.order in c("MG.3")){
  dim.reduct = "umap_mnn"
  DimPlot(src.mg3.integrated, reduction = dim.reduct, pt.size = 1.5,
          group.by = "cluster.v06.26.re",
          cols = cluster.endoderm.color.v5)
  DimPlot(src.mg3.integrated, reduction = dim.reduct, pt.size = 1.5,
          group.by = "Time", cols = colors.time)
  DimPlot(src.mg3.integrated, reduction = dim.reduct, pt.size = 1.5,
          group.by = "batch")
  
  type_define = "cluster.v06.26.re"
  #---------------------
  cell_src.mg3.integrated_mg3.sto = 
    rownames(src.mg3.integrated@meta.data[unlist(src.mg3.integrated[[type_define]]) %in%c("MG.3","MG.3.A","MG.3.M","Stomach"),])
  coord_src.mg3.integrated_mg3.sto = cbind(src.mg3.integrated[[dim.reduct]]@cell.embeddings[cell_src.mg3.integrated_mg3.sto,2],
                                           src.mg3.integrated[[dim.reduct]]@cell.embeddings[cell_src.mg3.integrated_mg3.sto,1])
  pcurve_src.mg3.integrated_mg3.sto = princurve::principal_curve(x = coord_src.mg3.integrated_mg3.sto, smoother = "smooth.spline")
  src.mg3.integrated$lambda_mg3.sto = pcurve_src.mg3.integrated_mg3.sto$lambda[cell_src.mg3.integrated_mg3.sto]
  src.mg3.integrated@meta.data[cell_src.mg3.integrated_mg3.sto,]$lambda_mg3.sto = 
    pcurve_src.mg3.integrated_mg3.sto$lambda[cell_src.mg3.integrated_mg3.sto]
  src.mg3.integrated@meta.data[!src.mg3.integrated$lambda_mg3.sto%in%NA,]$lambda_mg3.sto = 
    norm_range(src.mg3.integrated@meta.data[!src.mg3.integrated$lambda_mg3.sto%in%NA,]$lambda_mg3.sto)
  
  cell_src.mg3.integrated_mg3.dp = 
    rownames(src.mg3.integrated@meta.data[unlist(src.mg3.integrated[[type_define]]) %in%c("MG.3","MG.3.A","MG.3.M","DP"),])
  coord_src.mg3.integrated_mg3.dp = cbind(src.mg3.integrated[[dim.reduct]]@cell.embeddings[cell_src.mg3.integrated_mg3.dp,2],
                                          src.mg3.integrated[[dim.reduct]]@cell.embeddings[cell_src.mg3.integrated_mg3.dp,1])
  pcurve_src.mg3.integrated_mg3.dp = princurve::principal_curve(x = coord_src.mg3.integrated_mg3.dp, smoother = "smooth.spline")
  src.mg3.integrated$lambda_mg3.dp = pcurve_src.mg3.integrated_mg3.dp$lambda[cell_src.mg3.integrated_mg3.dp]
  src.mg3.integrated@meta.data[cell_src.mg3.integrated_mg3.dp,]$lambda_mg3.dp = 
    pcurve_src.mg3.integrated_mg3.dp$lambda[cell_src.mg3.integrated_mg3.dp]
  src.mg3.integrated@meta.data[!src.mg3.integrated$lambda_mg3.dp%in%NA,]$lambda_mg3.dp = 
    norm_range(src.mg3.integrated@meta.data[!src.mg3.integrated$lambda_mg3.dp%in%NA,]$lambda_mg3.dp)
  
  cell_src.mg3.integrated_mg3.ep = 
    rownames(src.mg3.integrated@meta.data[unlist(src.mg3.integrated[[type_define]]) %in%c("MG.3","MG.3.A","MG.3.M","EP.1","EP.2"),])
  coord_src.mg3.integrated_mg3.ep = cbind(src.mg3.integrated[[dim.reduct]]@cell.embeddings[cell_src.mg3.integrated_mg3.ep,2],
                                          src.mg3.integrated[[dim.reduct]]@cell.embeddings[cell_src.mg3.integrated_mg3.ep,1])
  pcurve_src.mg3.integrated_mg3.ep = princurve::principal_curve(x = coord_src.mg3.integrated_mg3.ep, smoother = "smooth.spline")
  src.mg3.integrated$lambda_mg3.ep = pcurve_src.mg3.integrated_mg3.ep$lambda[cell_src.mg3.integrated_mg3.ep]
  src.mg3.integrated@meta.data[cell_src.mg3.integrated_mg3.ep,]$lambda_mg3.ep = 
    pcurve_src.mg3.integrated_mg3.ep$lambda[cell_src.mg3.integrated_mg3.ep]
  src.mg3.integrated@meta.data[!src.mg3.integrated$lambda_mg3.ep%in%NA,]$lambda_mg3.ep = 
    norm_range(src.mg3.integrated@meta.data[!src.mg3.integrated$lambda_mg3.ep%in%NA,]$lambda_mg3.ep)
  
  cell_src.mg3.integrated_mg3.si1 = 
    rownames(src.mg3.integrated@meta.data[unlist(src.mg3.integrated[[type_define]]) %in%c("MG.3","MG.3.P","Small.intestine.1"),])
  coord_src.mg3.integrated_mg3.si1 = cbind(src.mg3.integrated[[dim.reduct]]@cell.embeddings[cell_src.mg3.integrated_mg3.si1,2],
                                           src.mg3.integrated[[dim.reduct]]@cell.embeddings[cell_src.mg3.integrated_mg3.si1,1])
  pcurve_src.mg3.integrated_mg3.si1 = princurve::principal_curve(x = coord_src.mg3.integrated_mg3.si1, smoother = "smooth.spline")
  src.mg3.integrated$lambda_mg3.si1 = pcurve_src.mg3.integrated_mg3.si1$lambda[cell_src.mg3.integrated_mg3.si1]
  src.mg3.integrated@meta.data[cell_src.mg3.integrated_mg3.si1,]$lambda_mg3.si1 = 
    pcurve_src.mg3.integrated_mg3.si1$lambda[cell_src.mg3.integrated_mg3.si1]
  src.mg3.integrated@meta.data[!src.mg3.integrated$lambda_mg3.si1%in%NA,]$lambda_mg3.si1 = 
    norm_range(src.mg3.integrated@meta.data[!src.mg3.integrated$lambda_mg3.si1%in%NA,]$lambda_mg3.si1)
  
  cell_src.mg3.integrated_mg3.si2 = 
    rownames(src.mg3.integrated@meta.data[unlist(src.mg3.integrated[[type_define]]) %in%c("MG.3","MG.3.P","Small.intestine.2"),])
  coord_src.mg3.integrated_mg3.si2 = cbind(src.mg3.integrated[[dim.reduct]]@cell.embeddings[cell_src.mg3.integrated_mg3.si2,2],
                                           src.mg3.integrated[[dim.reduct]]@cell.embeddings[cell_src.mg3.integrated_mg3.si2,1])
  pcurve_src.mg3.integrated_mg3.si2 = princurve::principal_curve(x = coord_src.mg3.integrated_mg3.si2, smoother = "smooth.spline")
  src.mg3.integrated$lambda_mg3.si2 = pcurve_src.mg3.integrated_mg3.si2$lambda[cell_src.mg3.integrated_mg3.si2]
  src.mg3.integrated@meta.data[cell_src.mg3.integrated_mg3.si2,]$lambda_mg3.si2 = 
    pcurve_src.mg3.integrated_mg3.si2$lambda[cell_src.mg3.integrated_mg3.si2]
  src.mg3.integrated@meta.data[!src.mg3.integrated$lambda_mg3.si2%in%NA,]$lambda_mg3.si2 = 
    norm_range(src.mg3.integrated@meta.data[!src.mg3.integrated$lambda_mg3.si2%in%NA,]$lambda_mg3.si2)
  #---------------------
  
  cell_src.mg3.integrated_mg3.sto.step1 = intersect(names(src.mg3.integrated$lambda_mg3.sto[order(src.mg3.integrated$lambda_mg3.sto)]),
                                                    rownames(src.mg3.integrated@meta.data[src.mg3.integrated$cluster.v06.26.re%in%c("MG.3","MG.3.A","MG.3.M"),]))
  cell_src.mg3.integrated_mg3.sto.step2 = intersect(names(src.mg3.integrated$lambda_mg3.sto[order(src.mg3.integrated$lambda_mg3.sto)]),
                                                    rownames(src.mg3.integrated@meta.data[src.mg3.integrated$cluster.v06.26.re%in%c("MG.3.A","MG.3.M","Stomach"),]))
  
  cell_src.mg3.integrated_mg3.dp.step1 = intersect(names(src.mg3.integrated$lambda_mg3.dp[order(src.mg3.integrated$lambda_mg3.dp)]),
                                                   rownames(src.mg3.integrated@meta.data[src.mg3.integrated$cluster.v06.26.re%in%c("MG.3","MG.3.A","MG.3.M"),]))
  cell_src.mg3.integrated_mg3.dp.step2 = intersect(names(src.mg3.integrated$lambda_mg3.dp[order(src.mg3.integrated$lambda_mg3.dp)]),
                                                   rownames(src.mg3.integrated@meta.data[src.mg3.integrated$cluster.v06.26.re%in%c("MG.3.A","MG.3.M","DP"),]))
  
  cell_src.mg3.integrated_mg3.ep.step1 = c(
    intersect(names(src.mg3.integrated$lambda_mg3.dp[order(src.mg3.integrated$lambda_mg3.ep)]),
              rownames(src.mg3.integrated@meta.data[src.mg3.integrated$cluster.v06.26.re%in%c("MG.3","MG.3.A","MG.3.M"),])))
  cell_src.mg3.integrated_mg3.ep.step2 = c(
    intersect(names(src.mg3.integrated$lambda_mg3.dp[order(src.mg3.integrated$lambda_mg3.ep)]),
              rownames(src.mg3.integrated@meta.data[src.mg3.integrated$cluster.v06.26.re%in%c("MG.3.A","MG.3.M"),])),
    intersect(names(src.mg3.integrated$lambda_mg3.dp[order(src.mg3.integrated$lambda_mg3.ep)]),
              rownames(src.mg3.integrated@meta.data[src.mg3.integrated$cluster.v06.26.re%in%c("EP.1"),])))
  cell_src.mg3.integrated_mg3.ep.step3 = c(
    intersect(names(src.mg3.integrated$lambda_mg3.dp[order(src.mg3.integrated$lambda_mg3.ep)]),
              rownames(src.mg3.integrated@meta.data[src.mg3.integrated$cluster.v06.26.re%in%c("EP.1"),])),
    intersect(names(src.mg3.integrated$lambda_mg3.dp[order(src.mg3.integrated$lambda_mg3.ep)]),
              rownames(src.mg3.integrated@meta.data[src.mg3.integrated$cluster.v06.26.re%in%c("EP.2"),])))
  
  cell_src.mg3.integrated_mg3.si1.step1 = c(
    intersect(names(src.mg3.integrated$lambda_mg3.si2[order(src.mg3.integrated$lambda_mg3.si2)]),
              rownames(src.mg3.integrated@meta.data[src.mg3.integrated$cluster.v06.26.re%in%c("MG.3"),])),
    intersect(names(src.mg3.integrated$lambda_mg3.si2[order(src.mg3.integrated$lambda_mg3.si2)]),
              rownames(src.mg3.integrated@meta.data[src.mg3.integrated$cluster.v06.26.re%in%c("MG.3.P"),])))
  cell_src.mg3.integrated_mg3.si1.step2 = c(
    intersect(names(src.mg3.integrated$lambda_mg3.si1[order(src.mg3.integrated$lambda_mg3.si1)]),
              rownames(src.mg3.integrated@meta.data[src.mg3.integrated$cluster.v06.26.re%in%c("MG.3.P"),])),
    intersect(names(src.mg3.integrated$lambda_mg3.si1[order(src.mg3.integrated$lambda_mg3.si1)]),
              rownames(src.mg3.integrated@meta.data[src.mg3.integrated$cluster.v06.26.re%in%c("Small.intestine.1"),])))
  
  cell_src.mg3.integrated_mg3.si2.step1 = c(
    intersect(names(src.mg3.integrated$lambda_mg3.si2[order(src.mg3.integrated$lambda_mg3.si2)]),
              rownames(src.mg3.integrated@meta.data[src.mg3.integrated$cluster.v06.26.re%in%c("MG.3"),])),
    intersect(names(src.mg3.integrated$lambda_mg3.si2[order(src.mg3.integrated$lambda_mg3.si2)]),
              rownames(src.mg3.integrated@meta.data[src.mg3.integrated$cluster.v06.26.re%in%c("MG.3.P"),])))
  cell_src.mg3.integrated_mg3.si2.step2 = c(
    intersect(names(src.mg3.integrated$lambda_mg3.si1[order(src.mg3.integrated$lambda_mg3.si1)]),
              rownames(src.mg3.integrated@meta.data[src.mg3.integrated$cluster.v06.26.re%in%c("MG.3.P"),])),
    intersect(names(src.mg3.integrated$lambda_mg3.si1[order(src.mg3.integrated$lambda_mg3.si1)]),
              rownames(src.mg3.integrated@meta.data[src.mg3.integrated$cluster.v06.26.re%in%c("Small.intestine.2"),])))
}

for(i.cell.order in c("HG.1")){
  dim.reduct = "umap_mnn_rot"
  DimPlot(src.hg1.integrated.re, reduction = dim.reduct, pt.size = 1.5,
          group.by = "cluster.v06.26.re_hc",
          cols = cluster.endoderm.color.v5)
  DimPlot(src.hg1.integrated.re, reduction = dim.reduct, pt.size = 1.5,
          group.by = "Time", cols = colors.time)
  DimPlot(src.hg1.integrated.re, reduction = dim.reduct, pt.size = 1.5,
          group.by = "batch")
  
  type_define = "cluster.v06.26.re_hc"
  #---------------------
  cell_src.hg1.integrated.re_hg1.si2 = 
    rownames(src.hg1.integrated.re@meta.data[unlist(src.hg1.integrated.re[[type_define]]) %in%c("HG.1","Small.intestine.2"),])
  coord_src.hg1.integrated.re_hg1.si2 = cbind(-src.hg1.integrated.re[[dim.reduct]]@cell.embeddings[cell_src.hg1.integrated.re_hg1.si2,1],
                                              src.hg1.integrated.re[[dim.reduct]]@cell.embeddings[cell_src.hg1.integrated.re_hg1.si2,2])
  pcurve_src.hg1.integrated.re_hg1.si2 = princurve::principal_curve(x = coord_src.hg1.integrated.re_hg1.si2, smoother = "smooth.spline")
  src.hg1.integrated.re$lambda_hg1.si2 = pcurve_src.hg1.integrated.re_hg1.si2$lambda[cell_src.hg1.integrated.re_hg1.si2]
  src.hg1.integrated.re@meta.data[cell_src.hg1.integrated.re_hg1.si2,]$lambda_hg1.si2 = 
    pcurve_src.hg1.integrated.re_hg1.si2$lambda[cell_src.hg1.integrated.re_hg1.si2]
  src.hg1.integrated.re@meta.data[!src.hg1.integrated.re$lambda_hg1.si2%in%NA,]$lambda_hg1.si2 = 
    norm_range(src.hg1.integrated.re@meta.data[!src.hg1.integrated.re$lambda_hg1.si2%in%NA,]$lambda_hg1.si2)
  
  cell_src.hg1.integrated.re_hg1.li1 = 
    rownames(src.hg1.integrated.re@meta.data[unlist(src.hg1.integrated.re[[type_define]]) %in%c("HG.1","Large.intestine.1"),])
  coord_src.hg1.integrated.re_hg1.li1 = cbind(-src.hg1.integrated.re[[dim.reduct]]@cell.embeddings[cell_src.hg1.integrated.re_hg1.li1,1],
                                              src.hg1.integrated.re[[dim.reduct]]@cell.embeddings[cell_src.hg1.integrated.re_hg1.li1,2])
  pcurve_src.hg1.integrated.re_hg1.li1 = princurve::principal_curve(x = coord_src.hg1.integrated.re_hg1.li1, smoother = "smooth.spline")
  src.hg1.integrated.re$lambda_hg1.li1 = pcurve_src.hg1.integrated.re_hg1.li1$lambda[cell_src.hg1.integrated.re_hg1.li1]
  src.hg1.integrated.re@meta.data[cell_src.hg1.integrated.re_hg1.li1,]$lambda_hg1.li1 = 
    pcurve_src.hg1.integrated.re_hg1.li1$lambda[cell_src.hg1.integrated.re_hg1.li1]
  src.hg1.integrated.re@meta.data[!src.hg1.integrated.re$lambda_hg1.li1%in%NA,]$lambda_hg1.li1 = 
    norm_range(src.hg1.integrated.re@meta.data[!src.hg1.integrated.re$lambda_hg1.li1%in%NA,]$lambda_hg1.li1)
  
  cell_src.hg1.integrated.re_hg1.li2 = 
    rownames(src.hg1.integrated.re@meta.data[unlist(src.hg1.integrated.re[[type_define]]) %in%c("HG.1-Large.intestine.2","Large.intestine.2"),])
  coord_src.hg1.integrated.re_hg1.li2 = cbind(src.hg1.integrated.re[[dim.reduct]]@cell.embeddings[cell_src.hg1.integrated.re_hg1.li2,1],
                                              src.hg1.integrated.re[[dim.reduct]]@cell.embeddings[cell_src.hg1.integrated.re_hg1.li2,2])
  pcurve_src.hg1.integrated.re_hg1.li2 = princurve::principal_curve(x = coord_src.hg1.integrated.re_hg1.li2, smoother = "smooth.spline")
  src.hg1.integrated.re$lambda_hg1.li2 = pcurve_src.hg1.integrated.re_hg1.li2$lambda[cell_src.hg1.integrated.re_hg1.li2]
  src.hg1.integrated.re@meta.data[cell_src.hg1.integrated.re_hg1.li2,]$lambda_hg1.li2 = 
    pcurve_src.hg1.integrated.re_hg1.li2$lambda[cell_src.hg1.integrated.re_hg1.li2]
  src.hg1.integrated.re@meta.data[!src.hg1.integrated.re$lambda_hg1.li2%in%NA,]$lambda_hg1.li2 = 
    norm_range(src.hg1.integrated.re@meta.data[!src.hg1.integrated.re$lambda_hg1.li2%in%NA,]$lambda_hg1.li2)
  #---------------------
  
  cell_src.hg1.integrated.re_hg1.si2 = c(
    intersect(names(src.hg1.integrated.re$lambda_hg1.si2[order(src.hg1.integrated.re$lambda_hg1.si2)]),
              rownames(src.hg1.integrated.re@meta.data[src.hg1.integrated.re$cluster.v06.26.re_hc%in%c("HG.1"),])),
    intersect(names(src.hg1.integrated.re$lambda_hg1.si2[order(src.hg1.integrated.re$lambda_hg1.si2)]),
              rownames(src.hg1.integrated.re@meta.data[src.hg1.integrated.re$cluster.v06.26.re_hc%in%c("Small.intestine.2"),]))) 
  
  cell_src.hg1.integrated.re_hg1.li1 = c(
    intersect(names(src.hg1.integrated.re$lambda_hg1.li1[order(src.hg1.integrated.re$lambda_hg1.li1)]),
              rownames(src.hg1.integrated.re@meta.data[src.hg1.integrated.re$cluster.v06.26.re_hc%in%c("HG.1"),])),
    intersect(names(src.hg1.integrated.re$lambda_hg1.li1[order(src.hg1.integrated.re$lambda_hg1.li1)]),
              rownames(src.hg1.integrated.re@meta.data[src.hg1.integrated.re$cluster.v06.26.re_hc%in%c("Large.intestine.1"),]))) 
  
  
  type_define = "cluster.v06.26.re_hc.refine"
  #---------------------
  dim.reduct = "umap_mnn"
  cell.step_src.hg1.integrated.re_hg1.hg1 = 
    rownames(src.hg1.integrated.re@meta.data[unlist(src.hg1.integrated.re[[type_define]]) %in%c("HG.1","HG.1-Large.intestine.2"),])
  coord_src.hg1.integrated.re_hg1.hg1 = cbind(src.hg1.integrated.re[[dim.reduct]]@cell.embeddings[cell.step_src.hg1.integrated.re_hg1.hg1,2],
                                              src.hg1.integrated.re[[dim.reduct]]@cell.embeddings[cell.step_src.hg1.integrated.re_hg1.hg1,1])
  pcurve_src.hg1.integrated.re_hg1.hg1 = princurve::principal_curve(x = coord_src.hg1.integrated.re_hg1.hg1, smoother = "smooth.spline")
  src.hg1.integrated.re$lambda_hg1.hg1 = pcurve_src.hg1.integrated.re_hg1.hg1$lambda[cell.step_src.hg1.integrated.re_hg1.hg1]
  src.hg1.integrated.re@meta.data[cell.step_src.hg1.integrated.re_hg1.hg1,]$lambda_hg1.hg1 = 
    pcurve_src.hg1.integrated.re_hg1.hg1$lambda[cell.step_src.hg1.integrated.re_hg1.hg1]
  src.hg1.integrated.re@meta.data[!src.hg1.integrated.re$lambda_hg1.hg1%in%NA,]$lambda_hg1.hg1 = 
    norm_range(src.hg1.integrated.re@meta.data[!src.hg1.integrated.re$lambda_hg1.hg1%in%NA,]$lambda_hg1.hg1) 
  
  dim.reduct = "umap_mnn_rot"
  cell.step_src.hg1.integrated.re_hg1.li2 = 
    rownames(src.hg1.integrated.re@meta.data[unlist(src.hg1.integrated.re[[type_define]]) %in%c("HG.1-Large.intestine.2","Large.intestine.2"),])
  coord_src.hg1.integrated.re_hg1.li2 = cbind(src.hg1.integrated.re[[dim.reduct]]@cell.embeddings[cell.step_src.hg1.integrated.re_hg1.li2,2],
                                              src.hg1.integrated.re[[dim.reduct]]@cell.embeddings[cell.step_src.hg1.integrated.re_hg1.li2,1])
  pcurve_src.hg1.integrated.re_hg1.li2 = princurve::principal_curve(x = coord_src.hg1.integrated.re_hg1.li2, smoother = "smooth.spline")
  src.hg1.integrated.re$lambda_hg1.li2 = pcurve_src.hg1.integrated.re_hg1.li2$lambda[cell.step_src.hg1.integrated.re_hg1.li2]
  src.hg1.integrated.re@meta.data[cell.step_src.hg1.integrated.re_hg1.li2,]$lambda_hg1.li2 = 
    pcurve_src.hg1.integrated.re_hg1.li2$lambda[cell.step_src.hg1.integrated.re_hg1.li2]
  src.hg1.integrated.re@meta.data[!src.hg1.integrated.re$lambda_hg1.li2%in%NA,]$lambda_hg1.li2 = 
    norm_range(src.hg1.integrated.re@meta.data[!src.hg1.integrated.re$lambda_hg1.li2%in%NA,]$lambda_hg1.li2)
  #---------------------
  
  cell.step_src.hg1.integrated.re_hg1.li2.step1 = c(
    intersect(names(src.hg1.integrated.re$lambda_hg1.hg1[order(src.hg1.integrated.re$lambda_hg1.hg1)]),
              c(rownames(src.hg1.integrated.re@meta.data[src.hg1.integrated.re$cluster.v06.26.re_hc.refine%in%c("HG.1") &
                                                           !src.hg1.integrated.re$RNA_snn_res.0.8%in%c(3,5,8,11,12,19),]))),
    intersect(rev(names(src.hg1.integrated.re$lambda_hg1.li2[order(src.hg1.integrated.re$lambda_hg1.li2)])),
              rownames(src.hg1.integrated.re@meta.data[src.hg1.integrated.re$cluster.v06.26.re_hc.refine%in%c("HG.1-Large.intestine.2"),]))) 
  cell.step_src.hg1.integrated.re_hg1.li2.step2 = c(
    intersect(rev(names(src.hg1.integrated.re$lambda_hg1.li2[order(src.hg1.integrated.re$lambda_hg1.li2)])),
              rownames(src.hg1.integrated.re@meta.data[src.hg1.integrated.re$cluster.v06.26.re_hc.refine%in%c("HG.1-Large.intestine.2"),])),
    intersect(rev(names(src.hg1.integrated.re$lambda_hg1.li2[order(src.hg1.integrated.re$lambda_hg1.li2)])),
              rownames(src.hg1.integrated.re@meta.data[src.hg1.integrated.re$cluster.v06.26.re_hc.refine%in%c("Large.intestine.2"),]))) 
}

for(i.cell.order in c("HG.2")){
  dim.reduct = "umap_integrated"
  DimPlot(src.hg2.integrated, reduction = dim.reduct, pt.size = 1.5,
          group.by = "cluster.v06.26.re",
          cols = cluster.endoderm.color.v5)
  DimPlot(src.hg2.integrated, reduction = dim.reduct, pt.size = 1.5,
          group.by = "Time", cols = colors.time)
  DimPlot(src.hg2.integrated, reduction = dim.reduct, pt.size = 1.5,
          group.by = "batch")
  
  type_define = "cluster.v06.26.re"
  #---------------------
  cell_src.hg2.integrated_hg2.li1 = 
    rownames(src.hg2.integrated@meta.data[unlist(src.hg2.integrated[[type_define]]) %in%c("HG.2","Large.intestine.1"),])
  coord_src.hg2.integrated_hg2.li1 = cbind(-src.hg2.integrated[[dim.reduct]]@cell.embeddings[cell_src.hg2.integrated_hg2.li1,2],
                                           src.hg2.integrated[[dim.reduct]]@cell.embeddings[cell_src.hg2.integrated_hg2.li1,1])
  pcurve_src.hg2.integrated_hg2.li1 = princurve::principal_curve(x = coord_src.hg2.integrated_hg2.li1, smoother = "smooth.spline")
  src.hg2.integrated$lambda_hg2.li1 = pcurve_src.hg2.integrated_hg2.li1$lambda[cell_src.hg2.integrated_hg2.li1]
  src.hg2.integrated@meta.data[cell_src.hg2.integrated_hg2.li1,]$lambda_hg2.li1 = 
    pcurve_src.hg2.integrated_hg2.li1$lambda[cell_src.hg2.integrated_hg2.li1]
  src.hg2.integrated@meta.data[!src.hg2.integrated$lambda_hg2.li1%in%NA,]$lambda_hg2.li1 = 
    norm_range(src.hg2.integrated@meta.data[!src.hg2.integrated$lambda_hg2.li1%in%NA,]$lambda_hg2.li1)
  
  cell_src.hg2.integrated_hg2.li3 = 
    rownames(src.hg2.integrated@meta.data[unlist(src.hg2.integrated[[type_define]]) %in%c("HG.2","Large.intestine.3"),])
  coord_src.hg2.integrated_hg2.li3 = cbind(-src.hg2.integrated[[dim.reduct]]@cell.embeddings[cell_src.hg2.integrated_hg2.li3,2],
                                           src.hg2.integrated[[dim.reduct]]@cell.embeddings[cell_src.hg2.integrated_hg2.li3,1])
  pcurve_src.hg2.integrated_hg2.li3 = princurve::principal_curve(x = coord_src.hg2.integrated_hg2.li3, smoother = "smooth.spline")
  src.hg2.integrated$lambda_hg2.li3 = pcurve_src.hg2.integrated_hg2.li3$lambda[cell_src.hg2.integrated_hg2.li3]
  src.hg2.integrated@meta.data[cell_src.hg2.integrated_hg2.li3,]$lambda_hg2.li3 = 
    pcurve_src.hg2.integrated_hg2.li3$lambda[cell_src.hg2.integrated_hg2.li3]
  src.hg2.integrated@meta.data[!src.hg2.integrated$lambda_hg2.li3%in%NA,]$lambda_hg2.li3 = 
    norm_range(src.hg2.integrated@meta.data[!src.hg2.integrated$lambda_hg2.li3%in%NA,]$lambda_hg2.li3)
  #---------------------
  
  cell_src.hg2.integrated_hg2.li1 = c(
    intersect(names(src.hg2.integrated$lambda_hg2.li1[order(src.hg2.integrated$lambda_hg2.li1)]),
              rownames(src.hg2.integrated@meta.data[src.hg2.integrated$cluster.v06.26.re%in%c("HG.2"),])),
    intersect(names(src.hg2.integrated$lambda_hg2.li1[order(src.hg2.integrated$lambda_hg2.li1)]),
              rownames(src.hg2.integrated@meta.data[src.hg2.integrated$cluster.v06.26.re%in%c("Large.intestine.1"),])))
  
  cell_src.hg2.integrated_hg2.li3 = c(
    rev(intersect(names(src.hg2.integrated$lambda_hg2.li3[order(src.hg2.integrated$lambda_hg2.li3)]),
                  rownames(src.hg2.integrated@meta.data[src.hg2.integrated$cluster.v06.26.re%in%c("HG.2"),]))),
    rev(intersect(names(src.hg2.integrated$lambda_hg2.li3[order(src.hg2.integrated$lambda_hg2.li3)]),
                  rownames(src.hg2.integrated@meta.data[src.hg2.integrated$cluster.v06.26.re%in%c("Large.intestine.3"),]))))
  
}
#--------------------------------------

list_cellorder_endoderm.merge.step = list()
for(i.type.plot in c("list_cellorder_endoderm.merge.step")){
  #--------------------------------------
  list_cellorder_endoderm.merge.step[["FG.1_Pha.2_pathway"]] = cellorder_src.fg1.integrated.fg1
  list_cellorder_endoderm.merge.step[["FG.2_Pha.1_pathway"]] = rev(cell_src.fg2.integrated_re_fg2.pha.1)
  
  list_cellorder_endoderm.merge.step[["FG.3_Lung_pathway"]] = c(cell_src.fg3.integrated_re_fg3, cell_src.fg3.integrated_re_lung) 
  list_cellorder_endoderm.merge.step[["FG.3_Pha.4_pathway"]] = c(cell_src.fg3.integrated_re_fg3,cell_src.fg3.integrated_re_pha4)
  list_cellorder_endoderm.merge.step[["FG.3_Pha.5_pathway"]] = c(cell_src.fg3.integrated_re_fg3, cell_src.fg3.integrated_re_pha5) 
  
  list_cellorder_endoderm.merge.step[["FG.4_Pha.5_pathway"]] = c(cell_src.fg4.integrated_refine_fg4, cell_src.fg4.integrated_refine_pha5) 
  list_cellorder_endoderm.merge.step[["FG.4_Lung_pathway.step1"]]  = c(cell_src.fg4.integrated_refine_fg4, cell_src.fg4.integrated_refine_fg4tolusto)
  list_cellorder_endoderm.merge.step[["FG.4_Lung_pathway.step2"]]  = c(cell_src.fg4.integrated_refine_fg4tolusto, cell_src.fg4.integrated_refine_lung)
  list_cellorder_endoderm.merge.step[["FG.4_Sto_pathway.step1"]]  = c(cell_src.fg4.integrated_refine_fg4, cell_src.fg4.integrated_refine_fg4tolusto)
  list_cellorder_endoderm.merge.step[["FG.4_Sto_pathway.step2"]]  = c(cell_src.fg4.integrated_refine_fg4tolusto, cell_src.fg4.integrated_refine_sto)
  list_cellorder_endoderm.merge.step[["FG.4_Liv_pathway.step1"]]  = c(cell_src.fg4.integrated_refine_fg4, cell_src.fg4.integrated_refine_fg4toliv)
  list_cellorder_endoderm.merge.step[["FG.4_Liv_pathway.step2"]]  = c(cell_src.fg4.integrated_refine_fg4toliv, cell_src.fg4.integrated_refine_al12toliv)
  list_cellorder_endoderm.merge.step[["FG.4_Liv_pathway.step3"]]  = c(cell_src.fg4.integrated_refine_al12toliv, cell_src.fg4.integrated_refine_liv)
  
  list_cellorder_endoderm.merge.step[["FG.5_Pha.3_pathway"]] = cell_src.fg5.integrated_fg5.pha.3
  list_cellorder_endoderm.merge.step[["FG.5_Th_pathway"]] = cell_src.fg5.integrated_fg5.th
  list_cellorder_endoderm.merge.step[["FG.6_Eso_pathway"]] = cellorder_src.fg6.integrated.fg6
  
  list_cellorder_endoderm.merge.step[["AL.1_Liv_pathway.step1"]] = cell_src.al1.integrated.re_al1.liv.step1 
  list_cellorder_endoderm.merge.step[["AL.1_Liv_pathway.step2"]] = cell_src.al1.integrated.re_al1.liv.step2 
  list_cellorder_endoderm.merge.step[["AL.2_Liv_pathway.step1"]] = cell_src.al2.integrated.re_al2.liv.step1
  list_cellorder_endoderm.merge.step[["AL.2_Liv_pathway.step2"]] = cell_src.al2.integrated.re_al2.liv.step2
  
  list_cellorder_endoderm.merge.step[["AL.3_SI.1_pathway.step1"]] = cell_src.al3.integrated.re_al3.si1.step1
  list_cellorder_endoderm.merge.step[["AL.3_SI.1_pathway.step2"]] = cell_src.al3.integrated.re_al3.si1.step2
  list_cellorder_endoderm.merge.step[["AL.3_EHBD_pathway.step1"]] = cell_src.al3.integrated.re_al3.ehbd.step1
  list_cellorder_endoderm.merge.step[["AL.3_EHBD_pathway.step2"]] = cell_src.al3.integrated.re_al3.ehbd.step2
  list_cellorder_endoderm.merge.step[["AL.3_VP_pathway.step1"]] = cell_src.al3.integrated.re_al3.vp.step1
  list_cellorder_endoderm.merge.step[["AL.3_VP_pathway.step2"]] = cell_src.al3.integrated.re_al3.vp.step2
  list_cellorder_endoderm.merge.step[["AL.3_Liv_pathway.step1"]] = cell_src.al3.integrated.re_al3.liv.step1
  list_cellorder_endoderm.merge.step[["AL.3_Liv_pathway.step2"]] = cell_src.al3.integrated.re_al3.liv.step2
  list_cellorder_endoderm.merge.step[["AL.3_Liv_pathway.step3"]] = cell_src.al3.integrated.re_al3.liv.step3
  
  list_cellorder_endoderm.merge.step[["MG.1_Sto_pathway"]] = cell_src.mg1.integrated_mg1.sto
  list_cellorder_endoderm.merge.step[["MG.1_SI.1_pathway"]] = cell_src.mg1.integrated_mg1.si1
  list_cellorder_endoderm.merge.step[["MG.2_SI.1_pathway"]] = cell_src.mg2.integrated_mg2.si1
  list_cellorder_endoderm.merge.step[["MG.2_SI.2_pathway"]] = cell_src.mg2.integrated_mg2.si2
  
  list_cellorder_endoderm.merge.step[["MG.3_Sto_pathway.step1"]] = cell_src.mg3.integrated_mg3.sto.step1
  list_cellorder_endoderm.merge.step[["MG.3_Sto_pathway.step2"]] = cell_src.mg3.integrated_mg3.sto.step2
  list_cellorder_endoderm.merge.step[["MG.3_DP_pathway.step1"]] = cell_src.mg3.integrated_mg3.dp.step1
  list_cellorder_endoderm.merge.step[["MG.3_DP_pathway.step1"]] = cell_src.mg3.integrated_mg3.dp.step2
  list_cellorder_endoderm.merge.step[["MG.3_EP_pathway.step1"]] = cell_src.mg3.integrated_mg3.ep.step1
  list_cellorder_endoderm.merge.step[["MG.3_EP_pathway.step2"]] = cell_src.mg3.integrated_mg3.ep.step2
  list_cellorder_endoderm.merge.step[["MG.3_EP_pathway.step3"]] = cell_src.mg3.integrated_mg3.ep.step3
  list_cellorder_endoderm.merge.step[["MG.3_SI.1_pathway.step1"]] = cell_src.mg3.integrated_mg3.si1.step1
  list_cellorder_endoderm.merge.step[["MG.3_SI.1_pathway.step2"]] = cell_src.mg3.integrated_mg3.si1.step2
  list_cellorder_endoderm.merge.step[["MG.3_SI.2_pathway.step1"]] = cell_src.mg3.integrated_mg3.si2.step1
  list_cellorder_endoderm.merge.step[["MG.3_SI.2_pathway.step2"]] = cell_src.mg3.integrated_mg3.si2.step2
  
  list_cellorder_endoderm.merge.step[["HG.1_SI.2_pathway"]] = cell_src.hg1.integrated.re_hg1.si2
  list_cellorder_endoderm.merge.step[["HG.1_LI.1_pathway"]] = cell_src.hg1.integrated.re_hg1.li1
  list_cellorder_endoderm.merge.step[["HG.1_LI.2_pathway.step1"]] = cell.step_src.hg1.integrated.re_hg1.li2.step1
  list_cellorder_endoderm.merge.step[["HG.1_LI.2_pathway.step2"]] = cell.step_src.hg1.integrated.re_hg1.li2.step2
  
  list_cellorder_endoderm.merge.step[["HG.2_LI.1_pathway"]] = cell_src.hg2.integrated_hg2.li1
  list_cellorder_endoderm.merge.step[["HG.2_LI.3_pathway"]] = cell_src.hg2.integrated_hg2.li3
  #--------------------------------------
}

#-- set :: seurat.plsda.pathway.endoderm.list.step
#--------------------------------------
seurat.auc.pathway.endoderm.list.step = list()
seurat.plsda.pathway.endoderm.list.step = list()
for(i.type in c("endoderm")){
  for(i.name in names(seurat.endoderm.list.pathway[[i.type]])){
    # if(i.name=="FG.5"){}else{next()}
    seurat = seurat.endoderm.list.pathway[[i.type]][[i.name]]
    seurat$cluster.temp = type.endoderm.list[[i.type]][[i.name]]$type
    
    seurat.auc.pathway.endoderm.list.step[[i.type]][[i.name]] = list()
    seurat.plsda.pathway.endoderm.list.step[[i.type]][[i.name]] = list()
    
    list_cellorder = names(list_cellorder_endoderm.merge.step)
    list_cellorder = list_cellorder[grepl(i.name, list_cellorder)]
    
    for(j.name in list_cellorder){
      cellorder = intersect(list_cellorder_endoderm.merge.step[[j.name]], colnames(seurat))
      
      if(length(cellorder)==0){
        print(paste(i.type, i.name, j.name, "Error for Selection!"))
        next()
      }
      
      seurat.temp = seurat[,cellorder]
      
      if(length(unique(seurat.temp$cluster.temp))<=1){
        print(paste(i.type, i.name, j.name, 'Error for Types!'))
        next()
      }else{
        print(paste(i.type, i.name, j.name))
      }
      
      seurat.temp = seurat_plsda_pathway(seurat.temp, type = seurat.temp$cluster.temp)
      seurat.auc.pathway.endoderm.list.step[[i.type]][[i.name]][[j.name]] = seurat.temp$auc
      seurat.plsda.pathway.endoderm.list.step[[i.type]][[i.name]][[j.name]] = seurat.temp$seurat
    }
  }
}

pathway.plsda.select.endoderm.list.step = list()
for(i.type in c("endoderm")){
  for(i.name in names(seurat.plsda.pathway.endoderm.list.step[[i.type]])){
    pathway.plsda.select.endoderm.list.step[[i.type]][[i.name]] = list()
    for(j.name in names(seurat.plsda.pathway.endoderm.list.step[[i.type]][[i.name]])){
      pathway.plsda.select.endoderm.list.step[[i.type]][[i.name]][[j.name]] =
        plsda_select(seurat.plsda.pathway.endoderm.list.step[[i.type]][[i.name]][[j.name]],
                     comp = 3, threshold = 25)
    }
  }
}

tree.plsda.pathway.endoderm.list.step = list()
for(i.type in names(seurat.plsda.pathway.endoderm.list.step)){
  for(i.name in names(seurat.plsda.pathway.endoderm.list.step[[i.type]])){
    # if(i.name=="MG.2"){}else{next()}
    tree.plsda.pathway.endoderm.list.step[[i.type]][[i.name]] = list()
    
    pdf(paste("Milestones/try.step.", i.type, "_", i.name, '.pdf', sep = ""), 10, 10)
    for(j.name in names(seurat.plsda.pathway.endoderm.list.step[[i.type]][[i.name]])){
      
      seurat = seurat.plsda.pathway.endoderm.list.step[[i.type]][[i.name]][[j.name]]
      gene_list = pathway.plsda.select.endoderm.list.step[[i.type]][[i.name]][[j.name]]
      cell_name = intersect(list_cellorder_endoderm.merge.step[[j.name]], colnames(seurat))
      
      seurat = ScaleData(seurat, features = rownames(seurat))
      data = seurat@assays$Pathway@scale.data[gene_list, cell_name]
      data_re = t(apply(data,1,kernelsmooth))
      rownames(data_re) = rownames(data); colnames(data_re) = colnames(data)
      
      plot(c(1:10), c(1:10))
      text(x=5,y=5,labels=paste(i.type, i.name, j.name, sep = " "), cex = 5)
      
      pathway.heatmap  =
        MyHeatmap(data_re,
                  type = "raw",
                  ColSideColors = cbind(
                    MyName2Col(seurat@meta.data[cell_name,]$batch, colors.type),
                    MyName2Col(seurat@meta.data[cell_name,]$Time, c(colors.time, colors.time.2)),
                    MyName2Col(seurat@meta.data[cell_name,]$cluster.temp, cluster.endoderm.color.v5)),
                  color.palette = colorRampPalette(c("#5aa7dd","#EBEBEB","red"), space="Lab"),
                  ColSideColorsSize = 3,
                  RowSideColorsSize = 2,
                  #Rowv = "none", Colv = "none",
                  return.tree = "row",
                  labCol="none", graph = T)
      
      tree.plsda.pathway.endoderm.list.step[[i.type]][[i.name]][[j.name]] = as.dendrogram(pathway.heatmap) 
    }
    dev.off()
  }
}
#--------------------------------------

gene.tree.plsda.pathway.endoderm.list.step = list()
for(i.tree.plsda in c("gene.tree.plsda.pathway.endoderm.list.step")){
  #-- Add-tree
  for(i.type in names(tree.plsda.pathway.endoderm.list.step)){
    gene.tree.plsda.pathway.endoderm.list.step[[i.type]] = list()
    for(j.type in names(tree.plsda.pathway.endoderm.list.step[[i.type]])){
      gene.tree.plsda.pathway.endoderm.list.step[[i.type]][[j.type]] = list()
      for(k.type in names(tree.plsda.pathway.endoderm.list.step[[i.type]][[j.type]])){
        gene.tree.plsda.pathway.endoderm.list.step[[i.type]][[j.type]][[k.type]] = NA
      }
    }
  }
  
  #-- FG.1::Pha.2
  #----------------------------------
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.1$FG.1_Pha.2_pathway = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.1$FG.1_Pha.2_pathway[[1]][[1]]),
    # labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.1$FG.1_Pha.2_pathway[[1]][[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.1$FG.1_Pha.2_pathway[[1]][[2]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.1$FG.1_Pha.2_pathway[[1]][[2]][[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.1$FG.1_Pha.2_pathway[[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.1$FG.1_Pha.2_pathway[[2]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.1$FG.1_Pha.2_pathway[[2]][[1]])
  )
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.1$FG.1_Pha.2_pathway) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.1$FG.1_Pha.2_pathway[[1]][[1]]),
                    # labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.1$FG.1_Pha.2_pathway[[1]][[2]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.1$FG.1_Pha.2_pathway[[1]][[2]][[2]][[2]])))),
    rep(3, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.1$FG.1_Pha.2_pathway[[1]][[2]][[1]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.1$FG.1_Pha.2_pathway[[2]][[2]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.1$FG.1_Pha.2_pathway[[2]][[2]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.1$FG.1_Pha.2_pathway[[2]][[1]])))))
  #----------------------------------
  
  #-- FG.2::Pha.1
  #----------------------------------
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.2$FG.2_Pha.1_pathway = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.2$FG.2_Pha.1_pathway[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.2$FG.2_Pha.1_pathway[[1]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.2$FG.2_Pha.1_pathway[[1]][[2]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.2$FG.2_Pha.1_pathway[[1]][[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.2$FG.2_Pha.1_pathway[[2]][[2]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.2$FG.2_Pha.1_pathway) = c(
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.2$FG.2_Pha.1_pathway[[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.2$FG.2_Pha.1_pathway[[1]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.2$FG.2_Pha.1_pathway[[1]][[2]][[2]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.2$FG.2_Pha.1_pathway[[1]][[2]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.2$FG.2_Pha.1_pathway[[2]][[2]])))))
  #----------------------------------
  
  #-- FG.3::Lung, Pha.4, Pha.5
  #----------------------------------
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Lung_pathway = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Lung_pathway[[1]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Lung_pathway[[1]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Lung_pathway[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Lung_pathway[[2]][[2]][[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Lung_pathway[[2]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Lung_pathway) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Lung_pathway[[1]][[2]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Lung_pathway[[1]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Lung_pathway[[1]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Lung_pathway[[2]][[2]][[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Lung_pathway[[2]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Pha.4_pathway = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Pha.4_pathway[[2]][[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Pha.4_pathway[[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Pha.4_pathway[[2]][[2]][[2]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Pha.4_pathway) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Pha.4_pathway[[2]][[1]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Pha.4_pathway[[2]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Pha.4_pathway[[2]][[2]][[2]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Pha.5_pathway = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Pha.5_pathway[[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Pha.5_pathway[[2]][[2]][[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Pha.5_pathway[[2]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Pha.5_pathway) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Pha.5_pathway[[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Pha.5_pathway[[2]][[2]][[2]][[2]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.3$FG.3_Pha.5_pathway[[2]][[1]])))))
  #----------------------------------
  
  #-- FG.4::Pha.5, Lung, Sto, Live
  #----------------------------------
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Pha.5_pathway = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Pha.5_pathway[[1]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Pha.5_pathway[[1]][[2]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Pha.5_pathway[[1]][[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Pha.5_pathway[[2]][[2]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Pha.5_pathway) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Pha.5_pathway[[1]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Pha.5_pathway[[1]][[2]][[2]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Pha.5_pathway[[1]][[2]][[2]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Pha.5_pathway[[2]][[2]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Lung_pathway.step1 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Lung_pathway.step1[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Lung_pathway.step1[[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Lung_pathway.step1[[2]][[2]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Lung_pathway.step1) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Lung_pathway.step1[[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Lung_pathway.step1[[1]][[2]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Lung_pathway.step1[[2]][[2]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Lung_pathway.step2 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Lung_pathway.step2[[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Lung_pathway.step2[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Lung_pathway.step2[[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Lung_pathway.step2[[2]][[2]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Lung_pathway.step2) = c(
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Lung_pathway.step2[[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Lung_pathway.step2[[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Lung_pathway.step2[[1]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Lung_pathway.step2[[2]][[2]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Sto_pathway.step1 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Sto_pathway.step1[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Sto_pathway.step1[[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Sto_pathway.step1[[2]][[2]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Sto_pathway.step1) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Sto_pathway.step1[[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Sto_pathway.step1[[1]][[2]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Sto_pathway.step1[[2]][[2]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Sto_pathway.step2 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Sto_pathway.step2[[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Sto_pathway.step2[[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Sto_pathway.step2[[1]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Sto_pathway.step2) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Sto_pathway.step2[[2]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Sto_pathway.step2[[1]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Sto_pathway.step2[[1]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step1 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step1[[2]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step1[[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step1[[1]][[2]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step1) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step1[[2]][[2]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step1[[2]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step1[[1]][[2]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step2 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step2[[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step2[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step2[[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step2[[2]][[2]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step2) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step2[[1]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step2[[1]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step2[[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step2[[2]][[2]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step3 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step3[[1]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step3[[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step3[[2]][[2]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step3) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step3[[1]][[2]][[2]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step3[[2]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.4$FG.4_Liv_pathway.step3[[2]][[2]])))))
  #----------------------------------
  
  #-- FG.5::Pha.3, Th
  #----------------------------------
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.5$FG.5_Pha.3_pathway = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.5$FG.5_Pha.3_pathway[[1]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.5$FG.5_Pha.3_pathway[[2]][[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.5$FG.5_Pha.3_pathway[[2]][[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.5$FG.5_Pha.3_pathway[[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.5$FG.5_Pha.3_pathway[[2]][[2]][[2]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.5$FG.5_Pha.3_pathway) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.5$FG.5_Pha.3_pathway[[1]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.5$FG.5_Pha.3_pathway[[2]][[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.5$FG.5_Pha.3_pathway[[2]][[1]][[2]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.5$FG.5_Pha.3_pathway[[2]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.5$FG.5_Pha.3_pathway[[2]][[2]][[2]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.5$FG.5_Th_pathway = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.5$FG.5_Th_pathway[[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.5$FG.5_Th_pathway[[2]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.5$FG.5_Th_pathway[[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.5$FG.5_Th_pathway[[1]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.5$FG.5_Th_pathway) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.5$FG.5_Th_pathway[[2]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.5$FG.5_Th_pathway[[2]][[2]][[2]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.5$FG.5_Th_pathway[[1]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.5$FG.5_Th_pathway[[1]][[1]])))))
  #----------------------------------
  
  #-- FG.6::Eso
  #----------------------------------
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.6$FG.6_Eso_pathway = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.6$FG.6_Eso_pathway[[2]][[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.6$FG.6_Eso_pathway[[2]][[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.6$FG.6_Eso_pathway[[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.6$FG.6_Eso_pathway[[1]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$FG.6$FG.6_Eso_pathway) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.6$FG.6_Eso_pathway[[2]][[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.6$FG.6_Eso_pathway[[2]][[1]][[2]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.6$FG.6_Eso_pathway[[1]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$FG.6$FG.6_Eso_pathway[[1]][[1]])))))
  #----------------------------------
  
  #-- AL.1::Liv
  #----------------------------------
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.1$AL.1_Liv_pathway.step1 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.1$AL.1_Liv_pathway.step1[[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.1$AL.1_Liv_pathway.step1[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.1$AL.1_Liv_pathway.step1[[2]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.1$AL.1_Liv_pathway.step1[[2]][[2]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.1$AL.1_Liv_pathway.step1) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.1$AL.1_Liv_pathway.step1[[1]][[2]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.1$AL.1_Liv_pathway.step1[[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.1$AL.1_Liv_pathway.step1[[2]][[2]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.1$AL.1_Liv_pathway.step1[[2]][[2]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.1$AL.1_Liv_pathway.step2 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.1$AL.1_Liv_pathway.step2[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.1$AL.1_Liv_pathway.step2[[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.1$AL.1_Liv_pathway.step2[[2]][[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.1$AL.1_Liv_pathway.step2[[2]][[2]][[2]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.1$AL.1_Liv_pathway.step2) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.1$AL.1_Liv_pathway.step2[[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.1$AL.1_Liv_pathway.step2[[1]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.1$AL.1_Liv_pathway.step2[[2]][[1]][[2]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.1$AL.1_Liv_pathway.step2[[2]][[2]][[2]])))))
  #----------------------------------
  
  #-- AL.2::Liv
  #----------------------------------
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.2$AL.2_Liv_pathway.step1 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.2$AL.2_Liv_pathway.step1[[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.2$AL.2_Liv_pathway.step1[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.2$AL.2_Liv_pathway.step1[[2]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.2$AL.2_Liv_pathway.step1[[2]][[2]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.2$AL.2_Liv_pathway.step1) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.2$AL.2_Liv_pathway.step1[[1]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.2$AL.2_Liv_pathway.step1[[1]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.2$AL.2_Liv_pathway.step1[[2]][[2]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.2$AL.2_Liv_pathway.step1[[2]][[2]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.2$AL.2_Liv_pathway.step2 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.2$AL.2_Liv_pathway.step2[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.2$AL.2_Liv_pathway.step2[[1]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.2$AL.2_Liv_pathway.step2[[1]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.2$AL.2_Liv_pathway.step2[[2]][[2]][[2]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.2$AL.2_Liv_pathway.step2) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.2$AL.2_Liv_pathway.step2[[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.2$AL.2_Liv_pathway.step2[[1]][[2]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.2$AL.2_Liv_pathway.step2[[1]][[2]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.2$AL.2_Liv_pathway.step2[[2]][[2]][[2]])))))
  
  #----------------------------------
  
  #-- AL.3::Liv, VP, EHBD, SI.1
  #----------------------------------
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step1 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step1[[1]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step1[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step1[[1]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step1[[2]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step1) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step1[[1]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step1[[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step1[[1]][[2]][[2]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step1[[2]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step2 = c(
    rev(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step2[[2]][[2]])),
    rev(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step2[[2]][[1]])),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step2[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step2[[1]][[2]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step2) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step2[[2]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step2[[2]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step2[[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step2[[1]][[2]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step3 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step3[[1]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step3[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step3[[2]][[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step3[[2]][[2]][[2]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step3) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step3[[1]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step3[[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step3[[2]][[2]][[2]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_Liv_pathway.step3[[2]][[2]][[2]])))))
  
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_VP_pathway.step1 = c(
    rev(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_VP_pathway.step1[[2]][[2]])),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_VP_pathway.step1[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_VP_pathway.step1[[1]][[2]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_VP_pathway.step1) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_VP_pathway.step1[[2]][[2]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_VP_pathway.step1[[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_VP_pathway.step1[[1]][[2]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_VP_pathway.step2 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_VP_pathway.step2[[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_VP_pathway.step2[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_VP_pathway.step2[[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_VP_pathway.step2[[2]][[2]][[2]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_VP_pathway.step2) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_VP_pathway.step2[[1]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_VP_pathway.step2[[1]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_VP_pathway.step2[[2]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_VP_pathway.step2[[2]][[2]][[2]])))))
  
  
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_EHBD_pathway.step1 = c(
    rev(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_EHBD_pathway.step1[[2]][[2]])),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_EHBD_pathway.step1[[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_EHBD_pathway.step1[[1]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_EHBD_pathway.step1) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_EHBD_pathway.step1[[2]][[2]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_EHBD_pathway.step1[[1]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_EHBD_pathway.step1[[1]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_EHBD_pathway.step2 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_EHBD_pathway.step2[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_EHBD_pathway.step2[[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_EHBD_pathway.step2[[2]][[2]][[2]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_EHBD_pathway.step2) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_EHBD_pathway.step2[[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_EHBD_pathway.step2[[2]][[2]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_EHBD_pathway.step2[[2]][[2]][[2]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_SI.1_pathway.step1 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_SI.1_pathway.step1[[2]][[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_SI.1_pathway.step1[[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_SI.1_pathway.step1[[2]][[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_SI.1_pathway.step1[[2]][[2]][[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_SI.1_pathway.step1[[2]][[2]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_SI.1_pathway.step1) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_SI.1_pathway.step1[[2]][[1]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_SI.1_pathway.step1[[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_SI.1_pathway.step1[[2]][[2]][[2]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_SI.1_pathway.step1[[2]][[2]][[2]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_SI.1_pathway.step1[[2]][[2]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_SI.1_pathway.step2 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_SI.1_pathway.step2[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_SI.1_pathway.step2[[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_SI.1_pathway.step2[[2]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_SI.1_pathway.step2) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_SI.1_pathway.step2[[1]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_SI.1_pathway.step2[[2]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$AL.3$AL.3_SI.1_pathway.step2[[2]][[1]])))))
  #----------------------------------
  
  #-- MG.1::Sto, SI.1 
  #----------------------------------
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.1$MG.1_Sto_pathway = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.1$MG.1_Sto_pathway[[2]][[2]][[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.1$MG.1_Sto_pathway[[2]][[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.1$MG.1_Sto_pathway[[2]][[2]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.1$MG.1_Sto_pathway[[2]][[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.1$MG.1_Sto_pathway[[2]][[1]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.1$MG.1_Sto_pathway) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.1$MG.1_Sto_pathway[[2]][[2]][[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.1$MG.1_Sto_pathway[[2]][[2]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.1$MG.1_Sto_pathway[[2]][[2]][[2]][[2]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.1$MG.1_Sto_pathway[[2]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.1$MG.1_SI.1_pathway = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.1$MG.1_SI.1_pathway[[1]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.1$MG.1_SI.1_pathway[[1]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.1$MG.1_SI.1_pathway[[2]][[2]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.1$MG.1_SI.1_pathway) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.1$MG.1_SI.1_pathway[[1]][[2]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.1$MG.1_SI.1_pathway[[1]][[2]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.1$MG.1_SI.1_pathway[[2]][[2]][[1]])))))
  #----------------------------------
  
  #-- MG.2::SI.1, SI.2 
  #----------------------------------
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.2$MG.2_SI.1_pathway = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.2$MG.2_SI.1_pathway[[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.2$MG.2_SI.1_pathway[[2]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.2$MG.2_SI.1_pathway[[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.2$MG.2_SI.1_pathway[[1]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.2$MG.2_SI.1_pathway) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.2$MG.2_SI.1_pathway[[2]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.2$MG.2_SI.1_pathway[[2]][[2]][[2]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.2$MG.2_SI.1_pathway[[1]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.2$MG.2_SI.1_pathway[[1]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.2$MG.2_SI.2_pathway = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.2$MG.2_SI.2_pathway[[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.2$MG.2_SI.2_pathway[[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.2$MG.2_SI.2_pathway[[2]][[2]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.2$MG.2_SI.2_pathway) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.2$MG.2_SI.2_pathway[[1]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.2$MG.2_SI.2_pathway[[2]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.2$MG.2_SI.2_pathway[[2]][[2]][[1]])))))
  #----------------------------------
  
  #-- MG.3::Sto, DP, EP, SI.1, SI.2
  #----------------------------------
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_Sto_pathway.step1 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_Sto_pathway.step1[[2]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_Sto_pathway.step1[[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_Sto_pathway.step1[[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_Sto_pathway.step1[[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_Sto_pathway.step1[[1]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_Sto_pathway.step1) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_Sto_pathway.step1[[2]][[2]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_Sto_pathway.step1[[2]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_Sto_pathway.step1[[2]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_Sto_pathway.step1[[1]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_Sto_pathway.step1[[1]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_Sto_pathway.step2 = c(
    rev(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_Sto_pathway.step2[[1]][[1]])),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_Sto_pathway.step2[[2]][[2]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_Sto_pathway.step2) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_Sto_pathway.step2[[1]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_Sto_pathway.step2[[2]][[2]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_DP_pathway.step1 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_DP_pathway.step1[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_DP_pathway.step1[[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_DP_pathway.step1[[2]][[2]][[2]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_DP_pathway.step1) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_DP_pathway.step1[[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_DP_pathway.step1[[1]][[2]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_DP_pathway.step1[[2]][[2]][[2]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step1 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step1[[2]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step1[[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step1[[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step1[[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step1[[1]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step1) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step1[[2]][[2]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step1[[2]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step1[[2]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step1[[1]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step1[[1]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step2 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step2[[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step2[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step2[[2]][[2]][[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step2[[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step2[[2]][[2]][[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step2[[2]][[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step2[[2]][[2]][[2]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step2[[2]][[2]][[2]][[2]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step2) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step2[[1]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step2[[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step2[[2]][[2]][[1]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step2[[2]][[2]][[1]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step2[[2]][[2]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step2[[2]][[2]][[2]][[2]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step2[[2]][[2]][[2]][[2]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step3 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step3[[1]][[1]]),
    rev(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step3[[2]][[2]])),
    rev(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step3[[2]][[1]])))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step3) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step3[[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step3[[2]][[2]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_EP_pathway.step3[[2]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.1_pathway.step1 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.1_pathway.step1[[2]][[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.1_pathway.step1[[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.1_pathway.step1[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.1_pathway.step1[[1]][[2]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.1_pathway.step1) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.1_pathway.step1[[2]][[2]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.1_pathway.step1[[2]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.1_pathway.step1[[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.1_pathway.step1[[1]][[2]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.1_pathway.step2 = c(
    rev(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.1_pathway.step2[[1]])),
    rev(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.1_pathway.step2[[2]][[1]])),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.1_pathway.step2[[2]][[2]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.1_pathway.step2) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.1_pathway.step2[[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.1_pathway.step2[[2]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.1_pathway.step2[[2]][[2]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.2_pathway.step1 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.2_pathway.step1[[2]][[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.2_pathway.step1[[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.2_pathway.step1[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.2_pathway.step1[[1]][[2]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.2_pathway.step1) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.2_pathway.step1[[2]][[2]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.2_pathway.step1[[2]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.2_pathway.step1[[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.2_pathway.step1[[1]][[2]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.2_pathway.step2 = c(
    rev(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.2_pathway.step2[[2]][[2]][[2]])),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.2_pathway.step2[[2]][[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.2_pathway.step2[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.2_pathway.step2[[1]][[2]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.2_pathway.step2) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.2_pathway.step2[[2]][[2]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.2_pathway.step2[[2]][[1]][[2]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.2_pathway.step2[[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$MG.3$MG.3_SI.2_pathway.step2[[1]][[2]])))))
  
  
  #----------------------------------
  
  #-- HG.1::SI.2, LI.1, LI.2
  #---------------------------------- 
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_SI.2_pathway = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_SI.2_pathway[[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_SI.2_pathway[[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_SI.2_pathway[[2]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_SI.2_pathway[[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_SI.2_pathway[[1]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_SI.2_pathway) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_SI.2_pathway[[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_SI.2_pathway[[2]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_SI.2_pathway[[2]][[2]][[2]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_SI.2_pathway[[1]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_SI.2_pathway[[1]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.1_pathway = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.1_pathway[[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.1_pathway[[2]][[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.1_pathway[[2]][[2]][[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.1_pathway[[1]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.1_pathway[[1]][[2]][[2]][[2]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.1_pathway) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.1_pathway[[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.1_pathway[[2]][[2]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.1_pathway[[2]][[2]][[2]][[2]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.1_pathway[[1]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.1_pathway[[1]][[2]][[2]][[2]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.2_pathway.step1 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.2_pathway.step1[[2]][[1]]),
    rev(labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.2_pathway.step1[[2]][[2]])),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.2_pathway.step1[[1]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.2_pathway.step1) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.2_pathway.step1[[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.2_pathway.step1[[2]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.2_pathway.step1[[1]][[1]])))))
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.2_pathway.step2 = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.2_pathway.step2[[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.2_pathway.step2[[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.2_pathway.step2[[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.2_pathway.step2[[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.2_pathway.step2[[2]][[2]][[2]][[2]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.2_pathway.step2) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.2_pathway.step2[[2]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.2_pathway.step2[[2]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.2_pathway.step2[[1]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.2_pathway.step2[[1]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.1$HG.1_LI.2_pathway.step2[[2]][[2]][[2]][[2]])))))
  #----------------------------------
  
  #-- HG.2::LI.1, LI.3
  #----------------------------------
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.1_pathway = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.1_pathway[[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.1_pathway[[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.1_pathway[[1]][[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.1_pathway[[1]][[2]][[2]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.1_pathway[[1]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.1_pathway[[1]][[2]][[2]][[1]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.1_pathway) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.1_pathway[[2]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.1_pathway[[2]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.1_pathway[[1]][[1]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.1_pathway[[1]][[2]][[2]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.1_pathway[[1]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.1_pathway[[1]][[2]][[2]][[1]])))))
  
  
  gene.tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.3_pathway = c(
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.3_pathway[[2]][[2]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.3_pathway[[2]][[2]][[2]][[1]][[2]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.3_pathway[[2]][[2]][[2]][[1]][[1]]),
    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.3_pathway[[2]][[1]][[2]]))
  names(gene.tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.3_pathway) = c(
    rep(4, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.3_pathway[[2]][[2]][[1]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.3_pathway[[2]][[2]][[2]][[1]][[2]]),
                    labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.3_pathway[[2]][[2]][[2]][[1]][[1]])))),
    rep(7, length(c(labels(tree.plsda.pathway.endoderm.list.step$endoderm$HG.2$HG.2_LI.3_pathway[[2]][[1]][[2]])))))
  #----------------------------------
  
} 

#-- Filter by anova
#--------------------------------------
gene.tree.plsda.pathway.endoderm.list.step.filter = list()
for(i.type in names(gene.tree.plsda.pathway.endoderm.list.step)){
  gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]] = list()
  
  for(j.type in names(gene.tree.plsda.pathway.endoderm.list.step[[i.type]])){
    gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]][[j.type]] = list()
    
    for(k.type in names(gene.tree.plsda.pathway.endoderm.list.step[[i.type]][[j.type]])){
      print(paste(i.type,j.type,k.type))
      
      src.rna = seurat.endoderm.list[[i.type]][[j.type]]
      type.src = type.endoderm.list[[i.type]][[j.type]]$type
      
      src.plsda = seurat.plsda.pathway.endoderm.list.step[[i.type]][[j.type]][[k.type]]
      gene.plsda = gene.tree.plsda.pathway.endoderm.list.step[[i.type]][[j.type]][[k.type]]
      
      gene.anova = anova.test(tpm = src.plsda@assays$Pathway@data[gene.plsda, colnames(src.plsda)],
                              variable = type.src[colnames(src.plsda)])
      gene.anova.scale = anova.test(tpm = src.plsda@assays$Pathway@scale.data[gene.plsda, colnames(src.plsda)],
                                    variable = type.src[colnames(src.plsda)])
      
      gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]][[j.type]][[k.type]][["gene"]] = gene.plsda
      gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]][[j.type]][[k.type]][["p.val.data"]] = gene.anova
      gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]][[j.type]][[k.type]][["p.val.scaledata"]] = gene.anova.scale
    }
  }
}

anova.curve.gene.tree.plsda.pathway.endoderm.list.step.filter = list()
anova.df.gene.tree.plsda.pathway.endoderm.list.step.filter = cbind(c(1:50))
for(i.type in names(gene.tree.plsda.pathway.endoderm.list.step.filter)){
  for(j.type in names(gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]])){
    for(k.type in names(gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]][[j.type]])){
      gene = gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]][[j.type]][[k.type]][["gene"]]
      pvalue.data = gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]][[j.type]][[k.type]][["p.val.data"]]
      num = c()
      for(i in c(1:50)){
        num = c(num, length(gene[pvalue.data<10^(-i)]))
        if(i==1 & length(gene[pvalue.data<10^(-i)])==0){print(paste(i.type,j.type,k.type))}
      }
      names(num) = c(1:50)
      anova.curve.gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]][[j.type]][[k.type]] = num 
      anova.df.gene.tree.plsda.pathway.endoderm.list.step.filter = cbind(
        anova.df.gene.tree.plsda.pathway.endoderm.list.step.filter, num)
    }
  }
}

for(i.plot in c("gene.tree.plsda.pathway.endoderm.list.step.filter")){
  df.plot = melt(anova.df.gene.tree.plsda.pathway.endoderm.list.step.filter[,2:ncol(anova.df.gene.tree.plsda.pathway.endoderm.list.step.filter)])
  list.index = c()
  for(i in c(1:(ncol(anova.df.gene.tree.plsda.pathway.endoderm.list.step.filter)-1))){list.index = c(list.index, rep(i,50))}
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

for(i.type in names(gene.tree.plsda.pathway.endoderm.list.step.filter)){
  for(j.type in names(gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]])){
    for(k.type in names(gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]][[j.type]])){
      gene = gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]][[j.type]][[k.type]][["gene"]]
      pvalue.data = gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]][[j.type]][[k.type]][["p.val.data"]]
      gene.rev = names(gene)
      names(gene.rev) = gene
      
      #-- p10 & p20 (Value range)
      #---------------------------- 
      gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]][[j.type]][[k.type]][["gene.p10"]] = gene[pvalue.data<10^(-10)]
      gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]][[j.type]][[k.type]][["gene.p20"]] = gene[pvalue.data<10^(-20)]
      
      if(length(gene)==0){print(paste(i.type, j.type, k.type));next()}
      if(length(gene[pvalue.data<10^(-10)])/length(gene) < 0.5){
        gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]][[j.type]][[k.type]][["gene.set"]] = "gene.p10"
      }else{
        gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]][[j.type]][[k.type]][["gene.set"]] = "gene.p20"
      }
      #---------------------------- 
      
      #-- p10 with top30 & top50
      #---------------------------- 
      gene.p10.top30 = intersect(gene, names(sort(pvalue.data)[1:min(30, length(gene[pvalue.data<10^(-10)]))]))
      names(gene.p10.top30) = gene.rev[gene.p10.top30]
      
      gene.p10.top50 = intersect(gene, names(sort(pvalue.data)[1:min(50, length(gene[pvalue.data<10^(-10)]))]))
      names(gene.p10.top50) = gene.rev[gene.p10.top50]
      
      gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]][[j.type]][[k.type]][["gene.p10.top30"]] = gene.p10.top30
      gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]][[j.type]][[k.type]][["gene.p10.top50"]] = gene.p10.top50
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
      
      
      gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]][[j.type]][[k.type]][["gene.p10.top30.stratified"]] = gene.p10.top30.stratified
      gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]][[j.type]][[k.type]][["gene.p10.top50.stratified"]] = gene.p10.top50.stratified
      #---------------------------- 
    }
  }
}
#--------------------------------------

save(list_cellorder_endoderm.merge.step,
     seurat.plsda.pathway.endoderm.list.step,
     pathway.plsda.select.endoderm.list.step,
     tree.plsda.pathway.endoderm.list.step,
     gene.tree.plsda.pathway.endoderm.list.step,
     gene.tree.plsda.pathway.endoderm.list.step.filter,
     file = "trajectory_signal/seurat.plsda.pathway.endoderm.list.step.parameter.Rdata")

#-- fin-plot
for(i.type in c("endoderm")){
  for(i.name in names(seurat.plsda.pathway.endoderm.list.step[[i.type]])){
    # if(i.name != "HG.1"){next()}
    pdf(paste("trajectory_signal/try.step.", i.type, "_", i.name, '.pdf', sep = ""), 8, 8)
    for(j.name in names(seurat.plsda.pathway.endoderm.list.step[[i.type]][[i.name]])){
      
      seurat = seurat.plsda.pathway.endoderm.list.step[[i.type]][[i.name]][[j.name]]
      seurat = ScaleData(seurat, features = rownames(seurat))
      
      seurat$cluster.temp = type.endoderm.list$endoderm[[i.name]]$type[colnames(seurat)]
      cell_name = intersect(list_cellorder_endoderm.merge.step[[j.name]], colnames(seurat))
      
      if(j.name %in% names(gene.tree.plsda.pathway.endoderm.list.step[[i.type]][[i.name]])){
        gene_list = gene.tree.plsda.pathway.endoderm.list.step.filter[[i.type]][[i.name]][[j.name]][["gene.p10.top30.stratified"]]
        gene_list = gene_list[is.na(gene_list)==F]
      }else{next()}
      
      print(paste(i.type, i.name, j.name))
      
      data = seurat@assays$Pathway@scale.data[gene_list, cell_name]
      data_re = t(apply(data,1,kernelsmooth))
      rownames(data_re) = rownames(data); colnames(data_re) = colnames(data)
      
      pathway.heatmap  =
        MyHeatmap(data_re,
                  type = "raw",
                  ColSideColors = cbind(
                    MyName2Col(seurat@meta.data[cell_name,]$Time, c(colors.time, colors.time.2)),
                    MyName2Col(seurat@meta.data[cell_name,]$cluster.temp, cluster.endoderm.color.v5)
                  ),
                  RowSideColors = t(cbind(
                    MyName2Col(names(gene_list), colors.geneset)
                  )),
                  color.palette = colorRampPalette(c("#5aa7dd","#EBEBEB","red"), space="Lab"),
                  ColSideColorsSize = 3,
                  RowSideColorsSize = 2,
                  Rowv = "none", 
                  Colv = "none",
                  # return.tree = "row",
                  margins = c(10,10),
                  labCol="none", graph = T)
    }
    dev.off()
  }
}

#===============================================================================




