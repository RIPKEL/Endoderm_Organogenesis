#===============================================================================
#>>> 9.signaling pathways analysis 
#===============================================================================

#>>>> The entire section summarizes the reference pathways, 
#> tracking group pathways, and organ pathways for the endoderm in Seurat, 
#> as well as the cell classification results. 
#> 
#>>>> The cell classification part can be found in the supplementary table 1, Sheets 3/4. 
#>
#>>>> The names of the variables below are the names used in 
#> my analysis and do not have any actual significance. 
#> 
#>>>> Here, only the `type.endoderm.list` and `seurat.endoderm.list` 
#> need to be organized for subsequent analysis.

#===============================================================================
#>>> 9.0 Summary of Seurta, Cell orders and Types : Tracing, Endoderm, Endoderm.organ
#===============================================================================
list_cluster_order.9ss = c("FG.1","FG.2","FG.3","FG.4","FG.5","FG.6",
                           "AL.1","AL.2","AL.3","MG.1","MG.2","MG.3","HG.1","HG.2")

list_cluster_order.15ss = c("Pharynx.organ.2","FG.1","FG.6","Esophagus","MG.3.A","DP","EP.1","MG.3.M","MG.3","MG.3.P",
                            'FG.2',"Pharynx.organ.1","FG.3","Pharynx.organ.4","Stomach","MG.1",
                            "Small.intestine.1","Small.intestine.2","MG.2","HG.1","Large.intestine.1",
                            "HG.1-Large.intestine.2","Large.intestine.2","HG.2","Large.intestine.3",
                            "FG.5","Thyroid","Pharynx.organ.3","Pharynx.organ.5","Lung","FG.4","FG.4-Lung/Stomach",
                            "FG.4-Liver","AL.1/2-Liver","Liver","EHBD","AL.3-Liver","VP","AL.3-EHBD/VP",'AL.3',"AL.3-Small.intestine.1")

list_cluster_order.27ss = c("Pharynx.organ.2","Pharynx.organ.1","Pharynx.organ.4","Pharynx.organ.5","Pharynx.organ.3",
                            "Thyroid","Lung","Esophagus","Stomach","Liver","EHBD","VP","DP","EP.1", "EP.2",
                            "Small.intestine.1","Small.intestine.2","Large.intestine.1","Large.intestine.2","Large.intestine.3")

#-- seurat:
list_seurat_tracing =
  c("src.fg1.tracing","src.fg2.tracing","src.fg3.tracing",
    "src.fg4.tracing","src.fg5.tracing","src.fg6.tracing",
    "src.al1.tracing","src.al2.tracing","src.al3.tracing",
    "src.mg1.tracing","src.mg2.tracing","src.mg3.tracing",
    "src.hg1.tracing","src.hg2.tracing")
names(list_seurat_tracing) = list_type_9ss

list_seurat_endoderm =
  c("src.fg1.integrated","src.fg2.integrated","src.fg3.integrated",
    "src.fg4.integrated","src.fg5.integrated","src.fg6.integrated",
    "src.al1.integrated","src.al2.integrated","src.al3.integrated",
    "src.mg1.integrated","src.mg2.integrated","src.mg3.integrated",
    "src.hg1.integrated","src.hg2.integrated")
names(list_seurat_endoderm) = list_type_9ss

list_seurat_endoderm.organ =
  c("src.endoderm.lung.ext_v1.1.re","src.endoderm.pha5.ext_v1.1.re",
    "src.endoderm.liver.ext_v1.1.re","src.endoderm.pan.ext_v1.1.re",
    "src.endoderm.sto.ext_v1.1.re","src.endoderm.sm1.ext_v1.1.re",
    "src.endoderm.sm2.ext_v1.1.re","src.endoderm.lar1.ext_v1.1.re")
names(list_seurat_endoderm.organ) =
  c("Lung","Pharynx.organ.5",'Liver',"Pancreas","Stomach",
    "Small.intestine.1","Small.intestine.2","Large.intestine.1")

seurat.endoderm.list = list()
for(i.name in c("tracing", "endoderm", "endoderm.organ")){
  seurat.endoderm.list[[i.name]] = list()
  list.try = get(paste("list_seurat_", i.name, sep = ""))
  for(j.name in names(list.try)){
    seurat.endoderm.list[[i.name]][[j.name]] = list.try[j.name]
  }
}

save(seurat.endoderm.list, file = "Milestones/seurat.endoderm.list.Rdata")


#-- cell types:
seurat.tracing.type.list = list()
for(i.type in list_type_9ss){
  seurat.tracing.type.list[[i.type]] = list()}
seurat.tracing.type.list$FG.1$type = src.fg1.tracing$cluster.v06.26.re_correct
seurat.tracing.type.list$FG.2$type = src.fg2.tracing$cluster.v06.26.re_umap_fta
seurat.tracing.type.list$FG.3$type = src.fg3.tracing$cluster.v06.26.re_correct_refine_umap_fta
seurat.tracing.type.list$FG.4$type = src.fg4.tracing$cluster.v06.26.re_correct
seurat.tracing.type.list$FG.5$type = src.fg5.tracing$cluster.v06.26.re_correct
seurat.tracing.type.list$FG.6$type = src.fg6.tracing$cluster.v06.26.re_correct
seurat.tracing.type.list$AL.1$type = src.al1.tracing$cluster.v06.26.re_correct_re
seurat.tracing.type.list$AL.2$type = src.al2.tracing$cluster.v06.26.re_correct
seurat.tracing.type.list$AL.3$type = src.al3.tracing$cluster.v06.26.re_correct.refine
seurat.tracing.type.list$MG.1$type = src.mg1.tracing$cluster.v06.26.re_correct
seurat.tracing.type.list$MG.2$type = src.mg2.tracing$cluster.v06.26.re_correct
seurat.tracing.type.list$MG.3$type = src.mg3.tracing$cluster.v06.26.re_correct
seurat.tracing.type.list$HG.1$type = src.hg1.tracing$cluster.v06.26.re_correct.refine
seurat.tracing.type.list$HG.2$type = src.hg2.tracing$cluster.v06.26.re_correct

seurat.endoderm.type.list = list()
for(i.type in list_type_9ss){
  seurat.endoderm.type.list[[i.type]] = list()}
seurat.endoderm.type.list$FG.1$type = src.fg1.integrated$cluster.v06.26.re_correct
seurat.endoderm.type.list$FG.2$type = src.fg2.integrated_re$cluster.v06.26.re_hc
seurat.endoderm.type.list$FG.3$type = src.fg3.integrated_re$cluster.v06.26.re_correct_refine.re
seurat.endoderm.type.list$FG.4$type = src.fg4.integrated_refine$cluster.v06.26.re_correct_refine
seurat.endoderm.type.list$FG.5$type = src.fg5.integrated$cluster.v06.26.re_correct
seurat.endoderm.type.list$FG.6$type = src.fg6.integrated$cluster.v06.26.re
seurat.endoderm.type.list$AL.1$type = src.al1.integrated.re$cluster.v06.26.re_correct_re
seurat.endoderm.type.list$AL.2$type = src.al2.integrated.re$cluster.v06.26.re_correct_re
seurat.endoderm.type.list$AL.3$type = src.al3.integrated.re$cluster.v06.26.re_correct
seurat.endoderm.type.list$MG.1$type = src.mg1.integrated$cluster.v06.26.re
seurat.endoderm.type.list$MG.2$type = src.mg2.integrated$cluster.v06.26.re
seurat.endoderm.type.list$MG.3$type = src.mg3.integrated$cluster.v06.26.re
seurat.endoderm.type.list$HG.1$type = src.hg1.integrated.re$cluster.v06.26.re_hc
seurat.endoderm.type.list$HG.2$type = src.hg2.integrated$cluster.v06.26.re

seurat.endoderm.organ.type.list = list()
for(i.type in c("Lung","Pharynx.organ.5",'Liver',"Pancreas","Stomach",
                "Small.intestine.1","Small.intestine.2","Large.intestine.1")){
  seurat.endoderm.organ.type.list[[i.type]] = list()}
seurat.endoderm.organ.type.list$Lung$type = src.endoderm.lung.ext_v1.1.re$cluster.v06.26.re_correct
seurat.endoderm.organ.type.list$Pharynx.organ.5$type = src.endoderm.pha5.ext_v1.1.re$cluster.v06.26.re_correct
seurat.endoderm.organ.type.list$Liver$type = src.endoderm.liver.ext_v1.1.re$cluster.v06.26.re_correct_refine
seurat.endoderm.organ.type.list$Pancreas$type = src.endoderm.pan.ext_v1.1.re$cluster.v06.26.re..merge
seurat.endoderm.organ.type.list$Stomach$type = src.endoderm.sto.ext_v1.1.re$cluster.v06.26.re..merge
seurat.endoderm.organ.type.list$Small.intestine.1$type = src.endoderm.sm1.ext_v1.1.re$cluster.v06.26.re..merge
seurat.endoderm.organ.type.list$Small.intestine.2$type = src.endoderm.sm2.ext_v1.1.re$cluster.v06.26.re..merge
seurat.endoderm.organ.type.list$Large.intestine.1$type = src.endoderm.lar1.ext_v1.1.re$cluster.v06.26.re..merge

type.endoderm.list = list()
for(i.name in c("tracing", "endoderm","endoderm.organ")){
  type.endoderm.list[[i.name]] = list()
  list.try = get(paste("seurat.", i.name, ".type.list", sep = ""))
  for(j.name in names(list.try)){
    type.endoderm.list[[i.name]][[j.name]] = list.try[[j.name]]
  }
}

save(type.endoderm.list, file = "Milestones/type.endoderm.list.Rdata")
#===============================================================================




