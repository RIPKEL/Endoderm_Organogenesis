
# Set UMAP from View
#----------------------------------
#-- lung & pha5
b = t(view_src.lung.integrated.pha5[1:3,1:3] %*% 
  t(as.matrix(src.lung.integrated.pha5@reductions$umap@cell.embeddings)))
colnames(b) = colnames(src.lung.integrated.pha5@reductions$umap@cell.embeddings)
src.lung.integrated.pha5@reductions$umap.rot = 
  src.lung.integrated.pha5@reductions$umap
src.lung.integrated.pha5@reductions$umap.rot@cell.embeddings = b
DimPlot(src.lung.integrated.pha5, reduction = 'umap.rot',
        group.by = "cluster.v06.26.re_correct")

#-- lung
src.lung.integrated@reductions$umap.rot
DimPlot(src.lung.integrated, group.by = "Time", reduction = 'umap.rot')

#-- pha5
src.pha5.integrated.re@reductions$umap.rot = 
  src.pha5.integrated.re@reductions$umap
src.pha5.integrated.re@reductions$umap.rot@cell.embeddings =
  src.pha5.integrated@reductions$umap.rot@cell.embeddings[colnames(src.pha5.integrated.re),]
DimPlot(src.pha5.integrated.re, group.by = "Time", reduction = 'umap.rot')
#----------------------------------


src.fg3.tracing$cluster.v06.26.re_correct_refine_mnn_umap_fta =
  src.fg3.integrated.merge@meta.data[colnames(src.fg3.tracing),]$cluster.v06.26.re_correct_refine_mnn_umap_fta
src.fg4.tracing$cluster.v06.26.re_correct_refine_mnn_umap_fta =
  src.fg4.integrated.merge.re@meta.data[colnames(src.fg4.tracing),]$cluster.v06.26.re_correct_refine_mnn_umap_fta


# Set Organ-merging seurat
#-------------------------------------------------------------------------------------------------------------------------------------------

# Fig 4
#=============================================================================
src.lung.integrated@meta.data = src.lung.integrated@meta.data[!src.lung.integrated$Time%in%NA,]
src.pha5.integrated.re@meta.data = src.pha5.integrated.re@meta.data[!src.pha5.integrated.re$Time%in%NA,]
src.lung.integrated.pha5@meta.data = src.lung.integrated.pha5@meta.data[!src.lung.integrated.pha5$Time%in%NA,]

#-- pha5 & lung
cell.lung_pha5.tracing.fg3 = rownames(
  src.fg3.tracing@meta.data[
    src.fg3.tracing$cluster.v06.26.re_correct_refine_mnn_umap_fta%in%c("FG.3","Pharynx.organ.5","Lung"),])
cell.lung_pha5.tracing.fg4 = rownames(
  src.fg4.tracing@meta.data[
    src.fg4.tracing$cluster.v06.26.re_correct_refine_mnn_umap_fta%in%c("FG.4","Pharynx.organ.5","FG.4-Lung/Stomach","Lung"),])

src.lung_pha5.tracing = merge(
  src.fg3.tracing[,cell.lung_pha5.tracing.fg3],
  src.fg4.tracing[,setdiff(cell.lung_pha5.tracing.fg4, cell.lung_pha5.tracing.fg3)])
src.lung.integrated.pha5 = NormalizeData(src.lung.integrated.pha5, assay = "RNA", scale.factor = 10^5)
src.lung_pha5.integrated.merge = merge(src.lung.integrated.pha5, src.lung_pha5.tracing)


#-- pha5
cell.pha5.tracing.fg3 = rownames(src.fg3.tracing@meta.data[src.fg3.tracing$cluster.v06.26.re_correct_refine_mnn_umap_fta%in%c("FG.3","Pharynx.organ.5"),])
cell.pha5.tracing.fg4 = rownames(src.fg4.tracing@meta.data[src.fg4.tracing$cluster.v06.26.re_correct_refine_mnn_umap_fta%in%c("FG.4","Pharynx.organ.5"),])
src.pha5.tracing = merge(
  src.fg3.tracing[,cell.pha5.tracing.fg3],
  src.fg4.tracing[,setdiff(cell.pha5.tracing.fg4, cell.pha5.tracing.fg3)])
src.pha5.integrated.re = NormalizeData(src.pha5.integrated.re, assay = "RNA", scale.factor = 10^5)
src.pha5.integrated.merge = merge(src.pha5.integrated.re, src.pha5.tracing)


#-- lung
cell.lung.tracing.fg3 = rownames(src.fg3.tracing@meta.data[src.fg3.tracing$cluster.v06.26.re_correct_refine_mnn_umap_fta%in%c("FG.3","Lung"),])
cell.lung.tracing.fg4 = rownames(src.fg4.tracing@meta.data[src.fg4.tracing$cluster.v06.26.re_correct_refine_mnn_umap_fta%in%c("FG.4",'FG.4-Lung/Stomach',"Lung"),])
src.lung.tracing = merge(
  src.fg3.tracing[,cell.lung.tracing.fg3],
  src.fg4.tracing[,setdiff(cell.lung.tracing.fg4, cell.lung.tracing.fg3)])
src.lung.integrated = NormalizeData(src.lung.integrated, assay = "RNA", scale.factor = 10^5)
src.lung.integrated.merge = merge(src.lung.integrated, src.lung.tracing)




# Set cell type
# ————————————————————————————————————————————————————————————
#-- lung
cell_refer_fg3_lung = rownames(src.fg3.integrated_re@meta.data[
  src.fg3.integrated_re$cluster.v06.26.re_correct_refine.re%in%c("FG.3","Lung"),])
cell_refer_fg4_lung = rownames(src.fg4.integrated_refine@meta.data[
  src.fg4.integrated_refine$cluster.v06.26.re_correct_refine%in%c("FG.4","FG.4-Lung/Stomach","Lung"),])

src.lung.integrated$cluster.v06.26.re_correct_refine = NA
src.lung.integrated@meta.data[cell_refer_fg3_lung,]$cluster.v06.26.re_correct_refine = 
  src.fg3.integrated_re@meta.data[cell_refer_fg3_lung,]$cluster.v06.26.re_correct_refine 
src.lung.integrated@meta.data[cell_refer_fg4_lung,]$cluster.v06.26.re_correct_refine = 
  src.fg4.integrated_refine@meta.data[cell_refer_fg4_lung,]$cluster.v06.26.re_correct_refine 

cell_src.lung.integrated = rownames(src.lung.integrated@meta.data[
  !src.lung.integrated$cluster.v06.26.re_correct_refine%in%NA,])

src.lung.integrated.merge$cluster.v06.26.re_mnn_umap_fta = 
  src.lung.integrated.merge$cluster.v06.26.re_correct_refine_mnn_umap_fta
src.lung.integrated.merge@meta.data[cell_src.lung.integrated,]$cluster.v06.26.re_mnn_umap_fta = 
  src.lung.integrated@meta.data[cell_src.lung.integrated,]$cluster.v06.26.re_correct_refine


#-- pha5
cell_refer_fg3_pha5 = rownames(src.fg3.integrated_re@meta.data[
  src.fg3.integrated_re$cluster.v06.26.re_correct_refine.re%in%c("FG.3","Pharynx.organ.5"),])
cell_refer_fg4_pha5 = rownames(src.fg4.integrated_refine@meta.data[
  src.fg4.integrated_refine$cluster.v06.26.re_correct_refine%in%c("FG.4","Pharynx.organ.5"),])

src.pha5.integrated.re$cluster.v06.26.re_correct_refine = NA
src.pha5.integrated.re@meta.data[cell_refer_fg3_pha5,]$cluster.v06.26.re_correct_refine = 
  src.fg3.integrated_re@meta.data[cell_refer_fg3_pha5,]$cluster.v06.26.re_correct_refine 
src.pha5.integrated.re@meta.data[cell_refer_fg4_pha5,]$cluster.v06.26.re_correct_refine = 
  src.fg4.integrated_refine@meta.data[cell_refer_fg4_pha5,]$cluster.v06.26.re_correct_refine 
cell_src.pha5.integrated.re = rownames(src.pha5.integrated.re@meta.data[
  !src.pha5.integrated.re$cluster.v06.26.re_correct_refine%in%NA,])

src.pha5.integrated.merge$cluster.v06.26.re_mnn_umap_fta = 
  src.pha5.integrated.merge$cluster.v06.26.re_correct_refine_mnn_umap_fta
src.pha5.integrated.merge@meta.data[cell_src.pha5.integrated.re, ]$cluster.v06.26.re_mnn_umap_fta = 
  src.pha5.integrated.re@meta.data[cell_src.pha5.integrated.re, ]$cluster.v06.26.re_correct_refine

src.pha5.integrated.merge@meta.data[colnames(src.pha5.integrated.re),]$cluster.extract.v1.1 = 
  src.pha5.integrated.re$cluster.extract.v1.1.re

#-- lung_pha5
cell_refer_fg3_lung_pha5 = rownames(src.fg3.integrated_re@meta.data[
  src.fg3.integrated_re$cluster.v06.26.re_correct_refine.re%in%c("FG.3","Pharynx.organ.5","Lung"),])
cell_refer_fg4_lung_pha5 = rownames(src.fg4.integrated_refine@meta.data[
  src.fg4.integrated_refine$cluster.v06.26.re_correct_refine%in%c("FG.4","Pharynx.organ.5",
                                                                  "FG.4-Lung/Stomach","Lung"),])
src.lung.integrated.pha5_fg5$cluster.v06.26.re_correct_refine = NA
src.lung.integrated.pha5_fg5@meta.data[cell_refer_fg3_lung_pha5,]$cluster.v06.26.re_correct_refine = 
  src.fg3.integrated_re@meta.data[cell_refer_fg3_lung_pha5,]$cluster.v06.26.re_correct_refine 
src.lung.integrated.pha5_fg5@meta.data[cell_refer_fg4_lung_pha5,]$cluster.v06.26.re_correct_refine = 
  src.fg4.integrated_refine@meta.data[cell_refer_fg4_lung_pha5,]$cluster.v06.26.re_correct_refine 
cell_src.lung.integrated.pha5_fg5 = rownames(src.lung.integrated.pha5_fg5@meta.data[
  !src.lung.integrated.pha5_fg5$cluster.v06.26.re_correct_refine%in%NA,])

src.lung_pha5.integrated.merge$cluster.v06.26.re_mnn_umap_fta = 
  src.lung_pha5.integrated.merge$cluster.v06.26.re_correct_refine_mnn_umap_fta
src.lung_pha5.integrated.merge@meta.data[cell_src.lung.integrated.pha5_fg5,]$cluster.v06.26.re_mnn_umap_fta = 
  src.lung.integrated.pha5_fg5@meta.data[cell_src.lung.integrated.pha5_fg5,]$cluster.v06.26.re_correct_refine

src.lung_pha5.integrated.merge@meta.data[colnames(src.pha5.integrated.re),]$cluster.extract.v1.1 =
  src.pha5.integrated.re$cluster.extract.v1.1.re
# ————————————————————————————————————————————————————————————



#=============================================================================

# Fig 5
#=============================================================================
src.al3.tracing_re$cluster.v06.26.re_correct_mnn_umap_fta = 
  src.al3.integrated.merge.re@meta.data[colnames(src.al3.tracing_re),]$cluster.v06.26.re_correct_mnn_umap_fta


#-- liver
cell.liver.tracing.fg4 = rownames(src.fg4.tracing@meta.data[
  src.fg4.tracing$cluster.v06.26.re_correct_refine_mnn_umap_fta%in%c("FG.4","FG.4-Liver","AL.1/2-Liver","Liver"),])
cell.liver.tracing.al1 = rownames(src.al1.tracing_re@meta.data[
  src.al1.tracing_re$cluster.v06.26.re_correct_re%in%c("AL.1","AL.1/2-Liver","Liver"),])
cell.liver.tracing.al2 = rownames(src.al2.tracing@meta.data[
  src.al2.tracing$cluster.v06.26.re_correct %in%c("AL.2","AL.1/2-Liver","Liver"),])
cell.liver.tracing.al3 = rownames(src.al3.tracing_re@meta.data[
  src.al3.tracing_re$cluster.v06.26.re_correct_mnn_umap_fta%in%c("AL.3","AL.3-Liver","AL.1/2-Liver","Liver"),])

src.liver.tracing = merge(
  src.fg4.tracing[,cell.liver.tracing.fg4], 
  c(src.al1.tracing_re[, setdiff(cell.liver.tracing.al1, 
                                 cell.liver.tracing.fg4)],
    src.al2.tracing[, setdiff(cell.liver.tracing.al2, 
                              unique(c(cell.liver.tracing.fg4, cell.liver.tracing.al1)))],
    src.al3.tracing[, setdiff(cell.liver.tracing.al3, 
                              unique(c(cell.liver.tracing.fg4, cell.liver.tracing.al1, cell.liver.tracing.al2)))]))

src.endoderm.liver.ext_v1.1.re = NormalizeData(src.endoderm.liver.ext_v1.1.re, assay = "RNA", scale.factor = 10^5)
cell_src.fg4.integrated.re_lv = rownames(src.fg4.integrated_refine@meta.data[
  src.fg4.integrated_refine$cluster.v06.26.re_correct_refine%in%c('FG.4',"AL.1/2-Liver","Liver"),])
cell_src.al1.integrated.re_lv = colnames(src.al1.integrated.re)
cell_src.al2.integrated.re_lv = colnames(src.al2.integrated.re)
cell_src.al3.integrated.re_lv = rownames(src.al3.integrated.re@meta.data[
  src.al3.integrated.re$cluster.v06.26.re_correct%in%c('AL.3',"AL.3-Liver","AL.1/2-Liver","Liver"),])

cell_src.liver.integrated.merge = unique(c(
  cell_src.fg4.integrated.re_lv,
  cell_src.al1.integrated.re_lv,
  cell_src.al2.integrated.re_lv,
  cell_src.al3.integrated.re_lv))

src.liver.integrated.merge = merge(
  src.endoderm.liver.ext_v1.1.re[,cell_src.liver.integrated.merge], 
  src.liver.tracing)

# -- cluster.v06.26.re_correct_refine
src.liver.integrated.merge$cluster.v06.26.re_correct_refine = NA
src.liver.integrated.merge@meta.data[cell_src.fg4.integrated.re_lv,]$cluster.v06.26.re_correct_refine =
  src.fg4.integrated_refine@meta.data[cell_src.fg4.integrated.re_lv,]$cluster.v06.26.re_correct_refine
src.liver.integrated.merge@meta.data[cell_src.al1.integrated.re_lv,]$cluster.v06.26.re_correct_refine =
  src.al1.integrated.re@meta.data[cell_src.al1.integrated.re_lv,]$cluster.v06.26.re_correct_re
src.liver.integrated.merge@meta.data[cell_src.al2.integrated.re_lv,]$cluster.v06.26.re_correct_refine =
  src.al2.integrated.re@meta.data[cell_src.al2.integrated.re_lv,]$cluster.v06.26.re_correct_re
src.liver.integrated.merge@meta.data[cell_src.al3.integrated.re_lv,]$cluster.v06.26.re_correct_refine =
  src.al3.integrated.re@meta.data[cell_src.al3.integrated.re_lv,]$cluster.v06.26.re_correct

# -- cluster.v06.26.re_mnn_umap_fta
src.liver.integrated.merge$cluster.v06.26.re_mnn_umap_fta = 
  src.liver.integrated.merge$cluster.v06.26.re_correct_refine 
src.liver.integrated.merge@meta.data[cell.liver.tracing.fg4,]$cluster.v06.26.re_mnn_umap_fta =
  src.fg4.tracing@meta.data[cell.liver.tracing.fg4,]$cluster.v06.26.re_correct_refine_mnn_umap_fta
src.liver.integrated.merge@meta.data[cell.liver.tracing.al1,]$cluster.v06.26.re_mnn_umap_fta =
  src.al1.tracing_re@meta.data[cell.liver.tracing.al1,]$cluster.v06.26.re_correct_re
src.liver.integrated.merge@meta.data[cell.liver.tracing.al2,]$cluster.v06.26.re_mnn_umap_fta =
  src.al2.tracing@meta.data[cell.liver.tracing.al2,]$cluster.v06.26.re_correct
src.liver.integrated.merge@meta.data[cell.liver.tracing.al3,]$cluster.v06.26.re_mnn_umap_fta =
  src.al3.tracing_re@meta.data[cell.liver.tracing.al3,]$cluster.v06.26.re_correct_mnn_umap_fta

save(src.liver.integrated.merge, file = "figure.v08.07/organ_development_re_v240115/src.liver.integrated.merge.Rdata")
#=============================================================================

# Fig 6
#=============================================================================

#-- stomach
cell.sto.tracing.fg4 = rownames(src.fg4.tracing@meta.data[
  src.fg4.tracing$cluster.v06.26.re_correct_refine_mnn_umap_fta%in%c("FG.4","FG.4-Lung/Stomach","Stomach"),])
cell.sto.tracing.mg1 = rownames(src.mg1.tracing.re@meta.data[
  src.mg1.tracing.re$cluster.v06.26.re_correct%in%c("MG.1","Stomach"),])
cell.sto.tracing.mg3 = rownames(src.mg3.tracing@meta.data[
  src.mg3.tracing$cluster.v06.26.re_mnn_umap_fta%in%c("MG.3","MG.3.A","MG.3.M","Stomach"),])

src.sto.tracing = merge(src.fg4.tracing[,cell.sto.tracing.fg4], 
  c(src.mg1.tracing[, setdiff(cell.sto.tracing.mg1, cell.sto.tracing.fg4)],
    src.mg3.tracing[, setdiff(cell.sto.tracing.mg3, unique(c(cell.sto.tracing.fg4, cell.sto.tracing.mg1)))]))

cell.sto.integrated.fg4 = rownames(src.fg4.integrated_refine@meta.data[
  src.fg4.integrated_refine$cluster.v06.26.re_correct_refine%in%c('FG.4',"FG.4-Lung/Stomach","Stomach"),])
cell.sto.integrated.mg1 = rownames(src.mg1.integrated@meta.data[
  src.mg1.integrated$cluster.v06.26.re%in%c("MG.1", "Stomach"),])
cell.sto.integrated.mg3 = rownames(src.mg3.integrated@meta.data[
  src.mg3.integrated$cluster.v06.26.re%in%c("MG.3","MG.3.A","MG.3.M","Stomach"),])

src.endoderm.sto.ext_v1.1.refine = NormalizeData(src.endoderm.sto.ext_v1.1.refine,
                                                 assay = "RNA", scale.factor = 10^5)
src.sto.integrated.merge = merge(
  src.endoderm.sto.ext_v1.1.refine[, intersect(colnames(src.endoderm.sto.ext_v1.1.refine),
                                               unique(c(cell.sto.integrated.fg4,
                                                        cell.sto.integrated.mg1,
                                                        cell.sto.integrated.mg3)))], 
  src.sto.tracing)

src.sto.integrated.merge$cluster.v06.26.re_correct_refine = NA
src.sto.integrated.merge@meta.data[cell.sto.integrated.fg4,]$cluster.v06.26.re_correct_refine = 
  src.fg4.integrated_refine@meta.data[cell.sto.integrated.fg4,]$cluster.v06.26.re_correct_refine
src.sto.integrated.merge@meta.data[cell.sto.integrated.mg1,]$cluster.v06.26.re_correct_refine = 
  src.mg1.integrated@meta.data[cell.sto.integrated.mg1,]$cluster.v06.26.re
src.sto.integrated.merge@meta.data[cell.sto.integrated.mg3,]$cluster.v06.26.re_correct_refine = 
  src.mg3.integrated@meta.data[cell.sto.integrated.mg3,]$cluster.v06.26.re
src.sto.integrated.merge@meta.data = src.sto.integrated.merge@meta.data[
  !src.sto.integrated.merge$Time%in%NA,]

src.sto.integrated.merge$cluster.v06.26.re_mnn_umap_fta =
  src.sto.integrated.merge$cluster.v06.26.re_correct_refine 
src.sto.integrated.merge@meta.data[cell.sto.tracing.fg4,]$cluster.v06.26.re_mnn_umap_fta = 
  src.fg4.tracing@meta.data[cell.sto.tracing.fg4,]$cluster.v06.26.re_correct_refine_mnn_umap_fta
src.sto.integrated.merge@meta.data[cell.sto.tracing.mg1,]$cluster.v06.26.re_mnn_umap_fta = 
  src.mg1.tracing.re@meta.data[cell.sto.tracing.mg1,]$cluster.v06.26.re_correct
src.sto.integrated.merge@meta.data[cell.sto.tracing.mg3,]$cluster.v06.26.re_mnn_umap_fta = 
  src.mg3.tracing@meta.data[cell.sto.tracing.mg3,]$cluster.v06.26.re_mnn_umap_fta



#-- small intestine 1
cell.sm1.tracing.mg1 = rownames(src.mg1.tracing@meta.data[
  src.mg1.tracing.re$cluster.v06.26.re_correct%in%c("MG.1","Small.intestine.1"),])
cell.sm1.tracing.mg2 = rownames(src.mg2.tracing@meta.data[
  src.mg2.tracing$cluster.v06.26.re_correct%in%c("MG.2","Small.intestine.1"),])
cell.sm1.tracing.mg3 = rownames(src.mg3.tracing@meta.data[
  src.mg3.tracing$cluster.v06.26.re_mnn_umap_fta%in%c("MG.3","MG.3.P","Small.intestine.1"),])
cell.sm1.tracing.al3 = rownames(src.al3.tracing@meta.data[
  src.al3.tracing$cluster.v06.26.re_mnn_umap_fta%in%c("AL.3","AL.3-Small.intestine.1","Small.intestine.1"),])

src.sm1.tracing = merge(
  src.mg1.tracing[,cell.sm1.tracing.mg1], 
  c(src.mg2.tracing[, setdiff(cell.sm1.tracing.mg2, cell.sm1.tracing.mg1)],
    src.mg3.tracing[, setdiff(cell.sm1.tracing.mg3, unique(c(cell.sm1.tracing.mg1, cell.sm1.tracing.mg2)))],
    src.al3.tracing[, setdiff(cell.sm1.tracing.al3, unique(c(cell.sm1.tracing.mg1, cell.sm1.tracing.mg2, cell.sm1.tracing.mg3)))]))


cell.sm1.integrated.al3 = rownames(src.al3.integrated.re@meta.data[
  src.al3.integrated.re$cluster.v06.26.re_correct%in%c('AL.3',"AL.3-Small.intestine.1","Small.intestine.1"),])
cell.sm1.integrated.mg1 = rownames(src.mg1.integrated@meta.data[
  src.mg1.integrated$cluster.v06.26.re%in%c("MG.1", "Small.intestine.1"),])
cell.sm1.integrated.mg2 = rownames(src.mg2.integrated@meta.data[
  src.mg2.integrated$cluster.v06.26.re%in%c("MG.2", "Small.intestine.1"),])
cell.sm1.integrated.mg3 = rownames(src.mg3.integrated@meta.data[
  src.mg3.integrated$cluster.v06.26.re%in%c("MG.3","MG.3.P","Small.intestine.1"),])


src.endoderm.sm1.ext_v1.1.re = NormalizeData(src.endoderm.sm1.ext_v1.1.re, assay = "RNA", scale.factor = 10^5)
src.sm1.integrated.merge = merge(
  src.endoderm.sm1.ext_v1.1.re[, intersect(colnames(src.endoderm.sm1.ext_v1.1.re),
                                           unique(c(cell.sm1.integrated.al3,
                                                    cell.sm1.integrated.mg1,
                                                    cell.sm1.integrated.mg2,
                                                    cell.sm1.integrated.mg3)))], 
  src.sm1.tracing)


src.sm1.integrated.merge$cluster.v06.26.re_correct_refine = NA
src.sm1.integrated.merge@meta.data[cell.sm1.integrated.al3,]$cluster.v06.26.re_correct_refine = 
  src.al3.integrated.re@meta.data[cell.sm1.integrated.al3,]$cluster.v06.26.re_correct
src.sm1.integrated.merge@meta.data[cell.sm1.integrated.mg1,]$cluster.v06.26.re_correct_refine = 
  src.mg1.integrated@meta.data[cell.sm1.integrated.mg1,]$cluster.v06.26.re
src.sm1.integrated.merge@meta.data[cell.sm1.integrated.mg2,]$cluster.v06.26.re_correct_refine = 
  src.mg2.integrated@meta.data[cell.sm1.integrated.mg2,]$cluster.v06.26.re
src.sm1.integrated.merge@meta.data[cell.sm1.integrated.mg3,]$cluster.v06.26.re_correct_refine = 
  src.mg3.integrated@meta.data[cell.sm1.integrated.mg3,]$cluster.v06.26.re

src.sm1.integrated.merge$cluster.v06.26.re_mnn_umap_fta =
  src.sm1.integrated.merge$cluster.v06.26.re_correct_refine 
src.sm1.integrated.merge@meta.data[cell.sm1.tracing.al3,]$cluster.v06.26.re_mnn_umap_fta = 
  src.al3.tracing@meta.data[cell.sm1.tracing.al3,]$cluster.v06.26.re_mnn_umap_fta
src.sm1.integrated.merge@meta.data[cell.sm1.tracing.mg1,]$cluster.v06.26.re_mnn_umap_fta = 
  src.mg1.tracing.re@meta.data[cell.sm1.tracing.mg1,]$cluster.v06.26.re_correct
src.sm1.integrated.merge@meta.data[cell.sm1.tracing.mg2,]$cluster.v06.26.re_mnn_umap_fta = 
  src.mg2.tracing@meta.data[cell.sm1.tracing.mg2,]$cluster.v06.26.re_correct
src.sm1.integrated.merge@meta.data[cell.sm1.tracing.mg3,]$cluster.v06.26.re_mnn_umap_fta = 
  src.mg3.tracing@meta.data[cell.sm1.tracing.mg3,]$cluster.v06.26.re_mnn_umap_fta

src.sm1.integrated.merge@meta.data = 
  src.sm1.integrated.merge@meta.data[!src.sm1.integrated.merge$Time%in%NA,]


#-- small intestine 2
src.hg1.tracing$cluster.v06.26.re_correct_mnn_umap_fta = 
  src.hg1.integrated.merge.re@meta.data[colnames(src.hg1.tracing),]$cluster.v06.26.re_hc_mnn_umap_fta

cell.sm2.tracing.hg1 = rownames(src.hg1.tracing@meta.data[
  src.hg1.tracing$cluster.v06.26.re_correct_mnn_umap_fta%in%c("HG.1","Small.intestine.2"),])
cell.sm2.tracing.mg2 = rownames(src.mg2.tracing@meta.data[
  src.mg2.tracing$cluster.v06.26.re_correct%in%c("MG.2","Small.intestine.2"),])
cell.sm2.tracing.mg3 = rownames(src.mg3.tracing@meta.data[
  src.mg3.tracing$cluster.v06.26.re_mnn_umap_fta%in%c("MG.3","MG.3.P","Small.intestine.2"),])

src.sm2.tracing = merge(
  src.hg1.tracing[,cell.sm2.tracing.hg1], 
  c(src.mg2.tracing[, setdiff(cell.sm2.tracing.mg2, cell.sm2.tracing.hg1)],
    src.mg3.tracing[, setdiff(cell.sm2.tracing.mg3, unique(c(cell.sm2.tracing.hg1, cell.sm2.tracing.mg2)))]))

cell.sm2.integrated.hg1 = rownames(src.hg1.integrated.re@meta.data[
  src.hg1.integrated.re$cluster.v06.26.re_hc%in%c('HG.1',"Small.intestine.2"),])
cell.sm2.integrated.mg2 = rownames(src.mg2.integrated@meta.data[
  src.mg2.integrated$cluster.v06.26.re%in%c("MG.2", "Small.intestine.2"),])
cell.sm2.integrated.mg3 = rownames(src.mg3.integrated@meta.data[
  src.mg3.integrated$cluster.v06.26.re%in%c("MG.3","MG.3.P","Small.intestine.2"),])


src.endoderm.sm2.ext_v1.1.re = NormalizeData(src.endoderm.sm2.ext_v1.1.re, assay = "RNA", scale.factor = 10^5)
src.sm2.integrated.merge = merge(
  src.endoderm.sm2.ext_v1.1.re[, intersect(colnames(src.endoderm.sm2.ext_v1.1.re),
                                           unique(c(cell.sm2.integrated.hg1,
                                                    cell.sm2.integrated.mg2,
                                                    cell.sm2.integrated.mg3)))], 
  src.sm2.tracing)


src.sm2.integrated.merge$cluster.v06.26.re_correct_refine = NA
src.sm2.integrated.merge@meta.data[cell.sm2.integrated.hg1,]$cluster.v06.26.re_correct_refine = 
  src.hg1.integrated.re@meta.data[cell.sm2.integrated.hg1,]$cluster.v06.26.re_hc
src.sm2.integrated.merge@meta.data[cell.sm2.integrated.mg2,]$cluster.v06.26.re_correct_refine = 
  src.mg2.integrated@meta.data[cell.sm2.integrated.mg2,]$cluster.v06.26.re
src.sm2.integrated.merge@meta.data[cell.sm2.integrated.mg3,]$cluster.v06.26.re_correct_refine = 
  src.mg3.integrated@meta.data[cell.sm2.integrated.mg3,]$cluster.v06.26.re

src.sm2.integrated.merge$cluster.v06.26.re_mnn_umap_fta =
  src.sm2.integrated.merge$cluster.v06.26.re_correct_refine 
src.sm2.integrated.merge@meta.data[cell.sm2.tracing.hg1,]$cluster.v06.26.re_mnn_umap_fta = 
  src.hg1.tracing@meta.data[cell.sm2.tracing.hg1,]$cluster.v06.26.re_correct_mnn_umap_fta 
src.sm2.integrated.merge@meta.data[cell.sm2.tracing.mg2,]$cluster.v06.26.re_mnn_umap_fta = 
  src.mg2.tracing@meta.data[cell.sm2.tracing.mg2,]$cluster.v06.26.re_correct
src.sm2.integrated.merge@meta.data[cell.sm2.tracing.mg3,]$cluster.v06.26.re_mnn_umap_fta = 
  src.mg3.tracing@meta.data[cell.sm2.tracing.mg3,]$cluster.v06.26.re_mnn_umap_fta

src.sm2.integrated.merge@meta.data = 
  src.sm2.integrated.merge@meta.data[!src.sm2.integrated.merge$Time%in%NA,]

#-- large intestine 1
cell.lar1.tracing.hg1 = rownames(src.hg1.tracing@meta.data[
  src.hg1.tracing$cluster.v06.26.re_correct_mnn_umap_fta%in%c("HG.1","Large.intestine.1"),])
cell.lar1.tracing.hg2 = rownames(src.hg2.tracing@meta.data[
  src.hg2.tracing$cluster.v06.26.re_mnn_umap_fta%in%c("HG.2","Large.intestine.1"),])

src.lar1.tracing = merge(
  src.hg1.tracing[,cell.lar1.tracing.hg1], 
  c(src.hg2.tracing[, setdiff(cell.lar1.tracing.hg2, cell.lar1.tracing.hg1)]))

cell.lar1.integrated.hg1 = rownames(src.hg1.integrated@meta.data[
  src.hg1.integrated.re$cluster.v06.26.re_hc%in%c("HG.1", "Large.intestine.1"),])
cell.lar1.integrated.hg2 = rownames(src.hg2.integrated@meta.data[
  src.hg2.integrated$cluster.v06.26.re%in%c("HG.2", "Large.intestine.1"),])

src.endoderm.lar1.ext_v1.1.re = NormalizeData(src.endoderm.lar1.ext_v1.1.re, assay = "RNA", scale.factor = 10^5)
src.lar1.integrated.merge = merge(
  src.endoderm.lar1.ext_v1.1.re[, intersect(colnames(src.endoderm.lar1.ext_v1.1.re),
                                            unique(c(cell.lar1.integrated.hg1,
                                                     cell.lar1.integrated.hg2)))], 
  src.lar1.tracing)


src.lar1.integrated.merge$cluster.v06.26.re_correct_refine = NA
src.lar1.integrated.merge@meta.data[cell.lar1.integrated.hg1,]$cluster.v06.26.re_correct_refine = 
  src.hg1.integrated@meta.data[cell.lar1.integrated.hg1,]$cluster.v06.26.re_hc
src.lar1.integrated.merge@meta.data[cell.lar1.integrated.hg2,]$cluster.v06.26.re_correct_refine = 
  src.hg2.integrated@meta.data[cell.lar1.integrated.hg2,]$cluster.v06.26.re

src.lar1.integrated.merge$cluster.v06.26.re_mnn_umap_fta =
  src.lar1.integrated.merge$cluster.v06.26.re_correct_refine 
src.lar1.integrated.merge@meta.data[cell.lar1.tracing.hg1,]$cluster.v06.26.re_mnn_umap_fta = 
  src.hg1.tracing@meta.data[cell.lar1.tracing.hg1,]$cluster.v06.26.re_correct_mnn_umap_fta
src.lar1.integrated.merge@meta.data[cell.lar1.tracing.hg2,]$cluster.v06.26.re_mnn_umap_fta = 
  src.hg2.tracing@meta.data[cell.lar1.tracing.hg2,]$cluster.v06.26.re_mnn_umap_fta

src.lar1.integrated.merge@meta.data = 
  src.lar1.integrated.merge@meta.data[!src.lar1.integrated.merge$Time%in%NA,]

#=============================================================================

# Fig 7
#=============================================================================
#-- pancreas
cell.pan.tracing.mg3 = rownames(
  src.mg3.tracing@meta.data[
    src.mg3.tracing$cluster.v06.26.re_mnn_umap_fta%in%c("MG.3","MG.3.M","EP.1","EP.2","DP"),])
cell.pan.tracing.al3 = rownames(
  src.al3.tracing@meta.data[
    src.al3.tracing$cluster.v06.26.re_mnn_umap_fta%in%c("AL.3","AL.3-EHBD/VP","VP"),])
cell.ep.tracing.ngn3cre = rownames(
  src.sm3.endoderm_Ngn3Cre@meta.data[
    src.sm3.endoderm_Ngn3Cre$cluster.refine%in%c("EP.1",'EP.2',"DP"),])
cell.ep.tracing.ngn3gfp = rownames(
  src.sm3.endoderm_Ngn3GFP@meta.data[
    src.sm3.endoderm_Ngn3GFP$cluster.predict.fin.umap.v06.26.re%in%c("MG.3.M","DP","EP.1","EP.2"),])
src.ep.tracing = merge(src.sm3.endoderm_Ngn3GFP[,cell.ep.tracing.ngn3gfp],
                       src.sm3.endoderm_Ngn3Cre[,cell.ep.tracing.ngn3cre])
src.pan.tracing = merge(
  src.mg3.tracing[,cell.pan.tracing.mg3], 
  c(src.al3.tracing[, setdiff(cell.pan.tracing.al3, cell.pan.tracing.mg3)]))

src.endoderm.pan.ext_v1.1.refine = NormalizeData(src.endoderm.pan.ext_v1.1.refine,
                                                 assay = "RNA", scale.factor = 10^4)
src.pan.integrated.merge = merge(src.endoderm.pan.ext_v1.1.refine,
                                 c(src.pan.tracing,
                                   src.ep.tracing))

src.pan.integrated.merge$cluster.v06.26.re_mnn_umap_fta = 
  src.pan.integrated.merge$cluster.v06.26.re..merge
src.pan.integrated.merge@meta.data[cell.ep.tracing.ngn3cre,]$cluster.v06.26.re_mnn_umap_fta =
  src.sm3.endoderm_Ngn3Cre$cluster.refine[cell.ep.tracing.ngn3cre]
src.pan.integrated.merge@meta.data[cell.ep.tracing.ngn3gfp,]$cluster.v06.26.re_mnn_umap_fta =
  src.sm3.endoderm_Ngn3GFP$cluster.refine[cell.ep.tracing.ngn3gfp]
src.pan.integrated.merge@meta.data[cell.pan.tracing.mg3,]$cluster.v06.26.re_mnn_umap_fta =
  src.mg3.tracing$cluster.v06.26.re_correct[cell.pan.tracing.mg3]
src.pan.integrated.merge@meta.data[cell.pan.tracing.al3,]$cluster.v06.26.re_mnn_umap_fta =
  src.al3.tracing_re$cluster.v06.26.re_correct.refine[cell.pan.tracing.al3]

src.pan.integrated.merge@meta.data[
  src.pan.integrated.merge$cluster.v06.26.re_mnn_umap_fta%in%c("AL.3-Small.intestine.1", "EHBD",
                                                               "Small.intestine.1", "Small.intestine.2"),]$cluster.v06.26.re_mnn_umap_fta = NA

#=============================================================================


#-- selected gene
#=================================================
src.lung_pha5.integrated.merge.selectgene = rownames(src.lung_pha5.integrated.merge@assays$mnnRNA)
src.pha5.integrated.merge.selectgene = src.pha5.integrated_filtergene
src.lung.integrated.merge.selectgene = rownames(src.lung.integrated@assays$mnnRNA)
src.liver.integrated.merge.selectgene = rownames(src.endoderm.liver.ext_v1.1.re@assays$mnnRNA)

src.sto.integrated.merge.selectgene = src.endoderm.sto.ext_v1.1.refine.gene.fin
src.sm1.integrated.merge.selectgene = src.endoderm.sm1.ext_v1.1.gene.fin
src.sm2.integrated.merge.selectgene = src.endoderm.sm2.ext_v1.1.gene.fin
src.lar1.integrated.merge.selectgene = src.endoderm.lar1.ext_v1.1.gene.fin

src.pan.integrated.merge.selectgene = src.endoderm.pan.ext.re.selectgene.fin
src.ep.integrated.merge.selectgene = src.endoderm.pan.ext.re.selectgene.fin

#-- reduction
red_organ = hash::hash()
red_organ["src.lung_pha5.integrated.merge"] = "umap.rot"
red_organ["src.pha5.integrated.merge"] = "umap.rot"
red_organ["src.lung.integrated.merge"] = "umap.rot"

red_organ["src.liver.integrated.merge"] = "umap_mnn"

red_organ["src.sto.integrated.merge"] = "umap_integrated"
red_organ["src.sm1.integrated.merge"] = "umap_integrated"
red_organ["src.sm2.integrated.merge"] = "umap_mnn"
red_organ["src.lar1.integrated.merge"] = "umap_mnn.rotated" # "umap_mnn"

red_organ["src.pan.integrated.merge"] = "umap_integrated"
red_organ["src.ep.integrated.merge"] = "umap_integrated"

red_organ_raw = hash::hash()
red_organ_raw["src.lung_pha5.integrated.merge"] = "src.lung.integrated.pha5"
red_organ_raw["src.pha5.integrated.merge"] = "src.pha5.integrated.re"
red_organ_raw["src.lung.integrated.merge"] = "src.lung.integrated"
red_organ_raw["src.liver.integrated.merge"] = "src.endoderm.liver.ext_v1.1.re"
red_organ_raw["src.sto.integrated.merge"] = "src.endoderm.sto.ext_v1.1.refine" # "src.endoderm.sto.ext_v1.1.re"
red_organ_raw["src.sm1.integrated.merge"] = "src.endoderm.sm1.ext_v1.1.re"
red_organ_raw["src.sm2.integrated.merge"] = "src.endoderm.sm2.ext_v1.1.re"
red_organ_raw["src.lar1.integrated.merge"] = "src.endoderm.lar1.ext_v1.1.re"
red_organ_raw["src.pan.integrated.merge"] = "src.endoderm.pan.ext_v1.1.refine"
red_organ_raw["src.ep.integrated.merge"] = "src.endoderm.pan.ext_v1.1.re"

DimPlot(src.endoderm.sto.ext_v1.1.refine, reduction = "umap_integrated")

embedding_rotated = function(seurat, reduction = "umap", theta= pi/2){
  ncol.red = ncol(seurat[[reduction]]@cell.embeddings)
  
  rotation_matrix = matrix(c(cos(theta), sin(theta), 
                             -sin(theta), cos(theta)), ncol = 2, byrow = TRUE)
  reduction_rotated = paste(reduction, ".rotated", sep = "")
  
  seurat[[reduction_rotated]] = seurat[[reduction]]
  if(ncol.red ==2 ){
    
    seurat[[reduction_rotated]]@cell.embeddings =
      seurat[[reduction]]@cell.embeddings %*% rotation_matrix
  
    }else{
    
    seurat[[reduction_rotated]]@cell.embeddings =
      cbind(seurat[[reduction]]@cell.embeddings[,c(1:2)] %*% rotation_matrix,
            seurat[[reduction]]@cell.embeddings[,c(3:ncol.red)])
  }
  
  # colnames(seurat[[reduction_rotated]]@cell.embeddings) = paste("UMAP_", c(1:ncol.red), sep = "")
  # seurat[[reduction_rotated]]@key = "UMAP_"
  return(seurat)
}

src.endoderm.lar1.ext_v1.1.re = embedding_rotated(src.endoderm.lar1.ext_v1.1.re, 
                                                  reduction = "umap_mnn", 
                                                  theta = pi/2.75)
colnames(src.endoderm.lar1.ext_v1.1.re[["umap_mnn.rotated"]]@cell.embeddings) = paste('Coord_', c(1:3), sep = "")
src.endoderm.lar1.ext_v1.1.re[["umap_mnn.rotated"]]@key = "Coord_"
DimPlot(src.endoderm.lar1.ext_v1.1.re, reduction = "umap_mnn.rotated",
        group.by = "cluster.v06.26.re..merge")
save(src.endoderm.lar1.ext_v1.1.re, file = "figure.v08.07/organ_development_re_v240115/src.endoderm.lar1.ext_v1.1.re.Rdata")



#=================================================

for(name_seurat in c(
  # "src.lung_pha5.integrated.merge",
  # "src.pha5.integrated.merge", "src.lung.integrated.merge" #,
  # "src.liver.integrated.merge"
  # "src.sto.integrated.merge"
  # "src.sm1.integrated.merge", 
  # "src.sm2.integrated.merge",
  # "src.lar1.integrated.merge"# , 
  "src.pan.integrated.merge"
  # "src.ep.integrated.merge"
)){
  
  seurat = get(name_seurat)
  DefaultAssay(seurat) = "RNA"
  seurat.selectgene = get(paste(name_seurat, "selectgene", sep="."))
  # seurat = CellCycleScoring(seurat, s.features = s_genes, g2m.features = g2m_genes)
  
  seurat$Source_tech = 'refer'
  seurat@meta.data[!seurat$lineage%in%NA,]$Source_tech = "query"
  seurat$index = paste(seurat$Source_tech, seurat$Phase, sep="_")
  
  seurat = ScaleData(seurat, split.by = "index")
  seurat = RunPCA(seurat, features = seurat.selectgene)
  seurat = RunUMAP(seurat, dims = 1:30, reduction = "pca")
  
  MNN.res =  mnnCorrect(
      as.matrix(seurat@assays$RNA@data[seurat.selectgene, 
                                       seurat$Source_tech%in%"refer"&
                                         seurat$batch%in%1]),
      as.matrix(seurat@assays$RNA@data[seurat.selectgene, 
                                       seurat$Source_tech%in%"refer"&
                                         seurat$batch%in%2]),
      as.matrix(seurat@assays$RNA@data[seurat.selectgene,
                                       seurat$Source_tech%in%"query"]),
      k = 20, cos.norm.out=F)
  
  
  seurat@assays$mnnRNA = CreateAssayObject(data = MNN.res@assays@data$corrected[,colnames(seurat)])
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
  
  assign(name_seurat, seurat)
}
  

for(name_seurat in c(
  # "src.lung_pha5.integrated.merge",
  # "src.pha5.integrated.merge", "src.lung.integrated.merge" #,
  # "src.liver.integrated.merge"
  # "src.sto.integrated.merge" #,
  # "src.sm1.integrated.merge", 
  # "src.sm2.integrated.merge",
  # "src.lar1.integrated.merge"# , 
  "src.pan.integrated.merge"
  # "src.ep.integrated.merge"
  )){
  
  seurat = get(name_seurat)
  seurat.selectgene = get(paste(name_seurat, "selectgene", sep="."))
  
  if(name_seurat %in% c("src.lung_pha5.integrated.merge",
                        "src.pha5.integrated.merge", 
                        "src.lung.integrated.merge")){  
    reference.assay =  "RNA" 
    
  }else if(name_seurat %in% c(#"src.ep.integrated.merge", 
                              "src.liver.integrated.merge", 
                              "src.sm2.integrated.merge", 
                              "src.lar1.integrated.merge")){
    reference.assay =  "mnnRNA" 
    
  }else{
    reference.assay =  "integrated" 
    }
  
  anchor.integrated = 
    FindTransferAnchors(reference = seurat[, seurat$lineage%in%NA],
                        query = seurat[, !seurat$lineage%in%NA],
                        reference.assay = "RNA",
                        query.assay =  "RNA", scale = T, 
                        features = seurat.selectgene, 
                        reduction = 'rpca', 
                        npcs = 30) # 50
  
  seurat_test = get(t(hash::values(red_organ_raw, key=name_seurat))[1])
  umap_embedding = seurat_test[[t(hash::values(red_organ, key=name_seurat))[1]]]@cell.embeddings[
    rownames(seurat@meta.data[seurat$lineage%in%NA,]),]
  umap.transfer = TransferData(anchor.integrated,  t(umap_embedding))
  seurat[["umap_fta"]] = seurat[["umap"]]
  seurat[["umap_fta"]]@cell.embeddings[rownames(seurat@meta.data[seurat$lineage%in%NA,]),] = umap_embedding[,c(1:2)]
  seurat[["umap_fta"]]@cell.embeddings[rownames(seurat@meta.data[!seurat$lineage%in%NA,]),] = 
    as.matrix(t(umap.transfer@data))[,c(1:2)]
  rm(seurat_test)
  
  assign(name_seurat, seurat)
}


# Set Plot parameter
#----------------------------------
#-- pancreas
src.pan.integrated.merge = label_KNN_learn_Re(
  src.pan.integrated.merge, reduction = "umap_fta", 
  label = "cluster.v06.26.re..merge", group = "Source_tech")
src.pan.integrated.merge$cluster.v06.26.re_mnn_umap_fta = 
  src.pan.integrated.merge$cluster.v06.26.re..merge_umap_fta

#-- ep
src.ep.integrated.merge = label_KNN_learn_Re(
  src.ep.integrated.merge, reduction = "umap_fta", 
  label = "cluster.v06.26.re..merge", group = "Source_tech")

src.ep.integrated.merge$cluster.v06.26.re_mnn_umap_fta = 
  src.ep.integrated.merge$cluster.v06.26.re..merge_umap_fta

src.ep.integrated.merge@meta.data[
  src.ep.integrated.merge$Source_tech%in%"query",]$cluster.v06.26.re_mnn_umap_fta = 
  src.ep.integrated.merge@meta.data[
    src.ep.integrated.merge$Source_tech%in%"query",]$cluster.predict.fin.umap.v06.26.re
#----------------------------------

reduct_to_meta = function(seurat,reduction){
  select_col = colnames(seurat@meta.data)[!grepl("Coord_",colnames(seurat@meta.data))]
  data_plot = cbind(
    seurat@meta.data[!seurat$Time%in%NA,select_col],
    seurat[[reduction]]@cell.embeddings[,c(1:2)])
  colnames(data_plot) = c(select_col, "Coord_1","Coord_2")
  return(data_plot)
}


#--  Plot Tracing
#----------------------------------------------------------------------------------
for(name_seurat in c(
  "src.lung_pha5.integrated.merge",
  "src.pha5.integrated.merge", "src.lung.integrated.merge" #,
  # "src.liver.integrated.merge"
  # "src.sto.integrated.merge" #,
  #"src.sm1.integrated.merge", 
  #"src.sm2.integrated.merge",
  #"src.lar1.integrated.merge"# , 
  # "src.pan.integrated.merge"
  # "src.ep.integrated.merge"
  )){
  
  # endoderm_names = gsub("\\.","",tolower(endoderm))
  seurat = get(name_seurat)
  
  for(names in c("umap_fta"
                 # "umap_integrated"
                 )){
    
    data_plot = reduct_to_meta(seurat, names)
    data_plot = data_plot[(!data_plot$cluster.v06.26.re_mnn_umap_fta%in%c(NA, "Small.intestine.1", "Small.intestine.2"))&
                            (!data_plot$lineage%in%c("Ngn3GFPHM", "Ngn3GFPHZ",
                                                    "Ngn3Cre"
                                                    )|
                               data_plot$lineage%in%NA),]
    
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
    
    g5 = ggplot()+ 
      geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                 mapping = aes(x=Coord_1, y =Coord_2),
                               # color=cluster.v06.26.re_mnn_umap_fta),
                 colour = "darkgray",
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
    
    g6 = ggplot()+ 
      geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                 mapping = aes(x=Coord_1, y =Coord_2),
                 # color=cluster.v06.26.re_mnn_umap_fta),
                 colour = "darkgray",
                 size=3.5)+
      # geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
      #            mapping = aes(x=Coord_1, y =Coord_2,
      #                          color=cluster.v06.26.re_mnn_umap_fta), 
      #            # color = "darkgray",
      #            size=2.55, shape=24, stroke=1.5)+
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
                 mapping = aes(x=Coord_1, y =Coord_2#, 
                               #color=Time
                               ), 
                 colour = "darkgray",
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
    
    for(k in c(1:6)){
      png(filename = paste("figure.v08.07/organ_development_re_v240115/organ/tracing_embedding.",
                           "Pathway_",gsub("src.",'', gsub("integrated.merge", "", name_seurat)),
                           "reduce_",names,"_pic",as.character(k),".png", sep=""),
          width = 1000,height = 1000,pointsize = 20)
      print(get(paste("g",as.character(k),sep="")))
      dev.off()
    }
    
  }
}

table(src.sm3.endoderm_Ngn3GFP$lineage,
      src.sm3.endoderm_Ngn3GFP$cluster.refine)

#----------------------------------------------------------------------------------



#-- pancreas
#------------------------------
seurat.selectgene = src.ep.integrated.merge.selectgene
pancreas.integrated = RunUMAP(src.endoderm.pan.ext_v1.1.re, 
                              dims = 1:30, 
                              reduction = "pca_integrated",
                              n.components = 3, 
                              return.model = TRUE)
DimPlot(pancreas.integrated, group.by = "Time")
# pancreas.integrated@reductions$umap@cell.embeddings = 
#  pancreas.integrated@reductions$umap_integrated@cell.embeddings

anchor.integrated = 
  FindTransferAnchors(reference = pancreas.integrated, 
                      query = src.ep.tracing,
                      reference.assay = "integrated", # "integrated",
                      query.assay =  "RNA", scale = T, # reduction = "cca", 
                      features = seurat.selectgene, 
                      reduction = "cca",
                      # reference.reduction = "umap_integrated", 
                      npcs = 50)

tracing.query = TransferData(anchorset = anchor.integrated, 
                             refdata = t(pancreas.integrated@reductions$umap_integrated@cell.embeddings))

pancreas.query = IntegrateEmbeddings(anchorset = anchor.integrated,
                                     reference = pancreas.integrated,
                                     query = src.ep.tracing, 
                                     new.reduction.name = "ref.pca")
pancreas.query = ProjectUMAP(query = pancreas.query, 
                             query.reduction = "ref.pca", 
                             reference = pancreas.integrated,
                             reference.reduction = "pca_integrated", 
                             reduction.model = "umap")
DimPlot(pancreas.query, reduction = "ref.umap",group.by = "lineage")

src.ep.integrated.merge@reductions$umap_fta@cell.embeddings[colnames(pancreas.integrated),] = 
  pancreas.integrated@reductions$umap_integrated@cell.embeddings[,c(1,2)]
src.ep.integrated.merge@reductions$umap_fta@cell.embeddings[colnames(pancreas.query),] = 
  as.matrix(t(tracing.query@data))[,c(1,2)]
  # pancreas.query@reductions$ref.umap@cell.embeddings[,c(1,3)]
  # as.matrix(t(epmach.query.tran@data))[colnames(epmach.query),c(1,2)]
DimPlot(src.ep.integrated.merge, reduction = "umap_fta",
        # group.by = "Time",
        group.by = "lineage", na.value = "#eeeeee", pt.size = 1.3,
        cols = c(color.lineage, cluster.endoderm.color.v5, colors.time.2, colors.time))

#------------------------------
save(src.ep.integrated.merge, file = "organ_development_re_v240115/organ/src.ep.integrated.merge.Rdata")

src.pan.integrated = src.endoderm.pan.ext_v1.1.re
# Check for pancreas and ep
#------------------------------
seurat = src.ep.integrated.merge
# seurat_q2 = src.pan.integrated.merge
cell_refer1 = rownames(seurat@meta.data[seurat$Source_tech%in%"refer" & seurat$batch%in%1,])
cell_refer2 = rownames(seurat@meta.data[seurat$Source_tech%in%"refer" & seurat$batch%in%2,])
cell_query_1 = rownames(seurat@meta.data[seurat$Source_tech%in%"query"&seurat$lineage%in%"Ngn3Cre",])
cell_query_2 = rownames(seurat@meta.data[seurat$Source_tech%in%"query"&seurat$lineage%in%"Ngn3GFP",])
# cell_query2 = rownames(seurat_q2@meta.data[seurat_q2$Source_tech%in%"query",])
seurat.selectgene = src.endoderm.pan.ext.re.selectgene.fin

for(assay in c("RNA")){
  
  anchor.integrated.merge = 
    FindIntegrationAnchors(c(seurat[,cell_refer1],
                             seurat[,cell_refer2],
                             seurat[,cell_query_1],
                             seurat[,cell_query_2]# ,
                             # seurat_q2[,cell_query2]
                             ),
                           assay = c(assay, assay, assay,
                                     assay, assay), 
                           reduction = "rpca", 
                           anchor.features = seurat.selectgene ,
                           dims = 1:30)
  
  seurat.re = IntegrateData(anchorset = anchor.integrated.merge, dims = 1:30)
  DefaultAssay(seurat.re) = "integrated"
  seurat.re = ScaleData(seurat.re, features = rownames(seurat.re))
  seurat.re = RunPCA(seurat.re, features = seurat.selectgene)
  seurat.re = RunUMAP(seurat.re, dims = 1:30, n.neighbors = 100,
                      n.components=2, umap.method = "uwot-learn")
  
  # seurat.re@reductions$umap
  DimPlot(seurat.re, reduction = 'umap',
          #group.by = "Treatment",
          #group.by = "Time",
          group.by = "lineage",
          pt.size = 1.5,
          # group.by = "cluster.v06.26.re_mnn_umap_fta",
          na.value = "#eeeeee",
          cols = c(# colors.time,
            colors.time.2, color.lineage, cluster.endoderm.color.v5,colors.type))
  
  
  assay_name_re = gsub('RNA',"", gsub("mnn_","mnn",assay))
  seurat.re[[paste(assay_name_re,"pca_integrated",sep="")]] = seurat.re[["pca"]]
  seurat.re[[paste(assay_name_re,"umap_integrated",sep="")]] = seurat.re[["umap"]]
  
  seurat[[paste(assay_name_re,"pca_integrated",sep="")]] = seurat[["pca"]]
  seurat[[paste(assay_name_re,"pca_integrated",sep="")]]@cell.embeddings[colnames(seurat),] = 
    seurat.re[[paste(assay_name_re,"pca_integrated",sep="")]]@cell.embeddings[colnames(seurat),] 
  seurat[[paste(assay_name_re,"umap_integrated",sep="")]] = seurat[["umap"]]
  seurat[[paste(assay_name_re,"umap_integrated",sep="")]]@cell.embeddings[colnames(seurat),] = 
    seurat.re[[paste(assay_name_re,"umap_integrated",sep="")]]@cell.embeddings[colnames(seurat),] 
  
  # seurat_q2[[paste(assay_name_re,"pca_integrated",sep="")]] = seurat_q2[["pca"]]
  # seurat_q2[[paste(assay_name_re,"pca_integrated",sep="")]]@cell.embeddings[colnames(seurat_q2),] = 
  #   seurat.re[[paste(assay_name_re,"pca_integrated",sep="")]]@cell.embeddings[colnames(seurat_q2),] 
  # seurat_q2[[paste(assay_name_re,"umap_integrated",sep="")]] = seurat_q2[["umap"]]
  # seurat_q2[[paste(assay_name_re,"umap_integrated",sep="")]]@cell.embeddings[colnames(seurat_q2),] = 
  #   seurat.re[[paste(assay_name_re,"umap_integrated",sep="")]]@cell.embeddings[colnames(seurat_q2),] 
  
  # rm(seurat.re)
}

# src.pan.integrated.merge = seurat_q2
src.ep.integrated.merge = seurat

src.endoderm.pan.ext_v1.1.re@reductions$umap_integrated_merge = 
  src.endoderm.pan.ext_v1.1.re@reductions$umap
src.endoderm.pan.ext_v1.1.re@reductions$umap_integrated_merge@cell.embeddings =
  src.pan.integrated.merge@reductions$umap_integrated@cell.embeddings[colnames(src.endoderm.pan.ext_v1.1.re),]
src.endoderm.pan.ext_v1.1.re@reductions$umap@cell.embeddings = 
  src.endoderm.pan.ext_v1.1.re@reductions$umap_integrated_merge@cell.embeddings

DimPlot(src.endoderm.pan.ext_v1.1.re, reduction = "umap_integrated_merge")
DimPlot(src.ep.integrated.merge, reduction = "mnnumap_integrated")

src.ep.integrated.merge@meta.data[
  src.ep.integrated.merge$lineage%in%"Ngn3GFP", ]$Treatment = 
  src.sm3.endoderm_Ngn3GFP$index[
    rownames(src.ep.integrated.merge@meta.data[
      src.ep.integrated.merge$lineage%in%"Ngn3GFP",])]

for(name_seurat in c("Ngn3Cre","Ngn3GFP")){

  seurat = src.ep.integrated.merge[,!src.ep.integrated.merge$lineage%in%name_seurat]
  
  pdf(paste("organ_development_re_v240115/organ/tracing_embedding.",
            gsub("src.",'', gsub("integrated.merge", "", name_seurat)),
            "sum.pdf", sep=""),12,10)
  
  for(names in c(#"umap_fta"
    "umap_integrated"
  )){
    data_plot = reduct_to_meta(seurat, names)
    
    g1 = ggplot()+ 
      geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                 mapping = aes(x=Coord_1, y =Coord_2, 
                               color=cluster.v06.26.re_mnn_umap_fta), 
                 size=2)+
      geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                 mapping = aes(x=Coord_1, y =Coord_2,
                               color=cluster.v06.26.re_mnn_umap_fta), 
                 # color = "darkgray",
                 size=2.55, shape=24, stroke=1.5)+
      scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                    colors.time,colors.time.2))+
      scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                   colors.time,colors.time.2))+
      theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
      theme(legend.position = "none")
    
    g5 = ggplot()+ 
      geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                 mapping = aes(x=Coord_1, y =Coord_2),
                 # color=cluster.v06.26.re_mnn_umap_fta),
                 colour = "#eeeeee",
                 size=2)+
      geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                 mapping = aes(x=Coord_1, y =Coord_2,
                               color=cluster.v06.26.re_mnn_umap_fta), 
                 # color = "darkgray",
                 size=2.55, shape=24, stroke=1.5)+
      scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                    colors.time,colors.time.2))+
      scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                   colors.time,colors.time.2))+
      theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
      theme(legend.position = "none")
    
    g2 = ggplot()+ 
      geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                 mapping = aes(x=Coord_1, y =Coord_2), 
                 colour = "#eeeeee",size=3)+
      geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                 mapping = aes(x=Coord_1, y =Coord_2, color=lineage), 
                 size=2.55, shape=2, stroke=1.2)+
      scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
      theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
      theme(legend.position = "none")
    
    g6 = ggplot()+ 
      geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                 mapping = aes(x=Coord_1, y =Coord_2), 
                 colour = "#eeeeee",size=3)+
      geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                 mapping = aes(x=Coord_1, y =Coord_2, color=Treatment), 
                 size=2.55, shape=2, stroke=1.2)+
      scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,colors.type))+
      theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
      theme(legend.position = "none")
    
    g3 = ggplot()+ 
      geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                 mapping = aes(x=Coord_1, y =Coord_2, 
                               color=cluster.extract.v1.1), 
                 size=2)+
      geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                 mapping = aes(x=Coord_1, y =Coord_2), 
                 colour = "#eeeeee",
                 size=2.55, shape=2, stroke=1.2)+
      scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5))+
      theme_void() + p_add + #ggtitle(paste(endoderm,names,sep="")) +
      theme(legend.position = "none")
    
    g4 = ggplot()+ 
      geom_point(data = data_plot[data_plot$Source_tech%in%"refer",],
                 mapping = aes(x=Coord_1, y =Coord_2), 
                 colour = "#eeeeee",
                 size=2)+
      geom_point(data = data_plot[data_plot$Source_tech%in%"query",],
                 mapping = aes(x=Coord_1, y =Coord_2,
                               color=Time), 
                 # color = "darkgray",
                 size=2.55, shape=24, stroke=1.5)+
      scale_color_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                    colors.time,colors.time.2))+
      scale_fill_manual(values = c(color.lineage, cluster.endoderm.color.v5,
                                   colors.time,colors.time.2))+
      theme_void() + p_add +# ggtitle(paste(endoderm,names,sep="")) +
      theme(legend.position = "none")
    
    for(k in c(1:6)){
      print(get(paste("g",as.character(k),sep="")))
    }
  }
  dev.off()
}
#------------------------------



