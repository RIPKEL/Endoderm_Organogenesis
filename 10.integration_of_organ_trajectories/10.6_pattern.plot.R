
src.endoderm.pan.ext_v1.1.re@reductions$umap@cell.embeddings = 
  src.endoderm.pan.ext_v1.1.re@reductions$umap_integrated@cell.embeddings[,c(1,2)]
pdf("organ_development_re_v240115/try.pattern.pan.pdf",9,7)
for(gene in c("Pdx1","Ptf1a","Mnx1",
              "Neurog3", "Sim1","Neurod1","Sox17",
              "Gcg","Nkx6-2","Nkx6-1","Arx")){
  print(FeaturePlot(src.endoderm.pan.ext_v1.1.re, 
                    features = gene, pt.size = 1.25) + 
          theme(aspect.ratio = 1))
}
dev.off()


pdf(paste("organ_development_re_v240115/organ/pattern_of_endoderm_15ss.pdf",sep=''),9,7)
for(gene in c("Nepn","Pdx1","Sox17","Prox1","Mnx1","Sim1","Irs4","Cckar")){
    print(
      FeaturePlot(src.15ss.integrated[,src.15ss.integrated$nFeature_RNA>2500], reduction = "umap",
                  features = gene, pt.size = 1.25) + 
        theme(aspect.ratio = 1))
  }
dev.off()




# Set Marker plot
#-----------------------------------
# Pharynx organ 5
src.pha5.integrated = SetIdent(src.pha5.integrated, 
                               value = src.pha5.integrated$cluster.v06.26.re_correct)
marker.pha5.integrated = FindAllMarkers(src.pha5.integrated)
marker.pha5.integrated$pct.ratio =
  marker.pha5.integrated$pct.1 / marker.pha5.integrated$pct.2 * 
  - log(marker.pha5.integrated$p_val_adj)

src.pha5.integrated.re@reductions$umap@cell.embeddings = 
  src.pha5.integrated.re@reductions$umap.rot@cell.embeddings
DimPlot(src.pha5.integrated.re, 
        reduction = "umap.rot", group.by = "Time")

marker_test_fg3 = marker.pha5.integrated[marker.pha5.integrated$cluster%in%'FG.3',]
gene_test_fg3 = marker_test_fg3[
  marker_test_fg3$pct.ratio > quantile(marker_test_fg3$pct.ratio, probs = 0.994),]$gene

marker_test_fg4 = marker.pha5.integrated[marker.pha5.integrated$cluster%in%'FG.4',]
gene_test_fg4 = marker_test_fg4[
  marker_test_fg4$pct.ratio > quantile(marker_test_fg4$pct.ratio, probs = 0.994),]$gene

marker_test_pha5 = marker.pha5.integrated[marker.pha5.integrated$cluster%in%'Pharynx.organ.5',]
gene_test_pha5 = marker_test_pha5[
  marker_test_pha5$pct.ratio > quantile(marker_test_pha5$pct.ratio, probs = 0.994),]$gene

pdf("organ_development_re_v240115/organ/pattern_of_pathway.pha5.pdf",9,7)
for(type in c("fg3","fg4","pha5")){
  gene_list = get(paste("gene_test_",type,sep=""))
  for(gene in gene_list){
    print(
      FeaturePlot(src.pha5.integrated.re, reduction = "umap.rot",
                  features = gene, pt.size = 1.25) + 
        theme(aspect.ratio = 1))
  }
}
dev.off()


# Lung
src.lung.integrated = SetIdent(src.lung.integrated, 
                               value = src.lung.integrated$cluster.v06.26.re_correct)
marker.lung.integrated = FindAllMarkers(src.lung.integrated)
marker.lung.integrated$pct.ratio =
  marker.lung.integrated$pct.1 / marker.lung.integrated$pct.2 * 
  - log(marker.lung.integrated$p_val_adj)

src.lung.integrated@reductions$umap@cell.embeddings = 
  src.lung.integrated@reductions$umap.rot@cell.embeddings
DimPlot(src.lung.integrated, 
        reduction = "umap", group.by = "Time")

marker_test_fg3 = marker.lung.integrated[marker.lung.integrated$cluster%in%'FG.3',]
gene_test_fg3 = marker_test_fg3[
  marker_test_fg3$pct.ratio > quantile(marker_test_fg3$pct.ratio, probs = 0.994),]$gene

marker_test_fg4 = marker.lung.integrated[marker.lung.integrated$cluster%in%'FG.4',]
gene_test_fg4 = marker_test_fg4[
  marker_test_fg4$pct.ratio > quantile(marker_test_fg4$pct.ratio, probs = 0.994),]$gene

marker_test_fg4lung = marker.lung.integrated[marker.lung.integrated$cluster%in%'FG.4-Lung/Stomach',]
gene_test_fg4lung = marker_test_fg4lung[
  marker_test_fg4lung$pct.ratio > quantile(marker_test_fg4lung$pct.ratio, probs = 0.994),]$gene

marker_test_lung = marker.lung.integrated[marker.lung.integrated$cluster%in%'Lung',]
gene_test_lung = marker_test_lung[
  marker_test_lung$pct.ratio > quantile(marker_test_lung$pct.ratio, probs = 0.994),]$gene

pdf("organ_development_re_v240115/organ/pattern_of_pathway.lung.pdf",9,7)
for(type in c("fg3","fg4","fg4lung","lung")){
  gene_list = get(paste("gene_test_",type,sep=""))
  for(gene in gene_list){
    print(
      FeaturePlot(src.lung.integrated, reduction = "umap.rot",
                  features = gene, pt.size = 1.25) + 
        theme(aspect.ratio = 1))
  }
}
dev.off()

# Liver
src.endoderm.liver.ext_v1.1.re = SetIdent(src.endoderm.liver.ext_v1.1.re, 
                                          value = src.endoderm.liver.ext_v1.1.re$cluster.v06.26.re..correct)
marker.liver.integrated = FindAllMarkers(src.endoderm.liver.ext_v1.1.re)
marker.liver.integrated$pct.ratio =
  marker.liver.integrated$pct.1 / marker.liver.integrated$pct.2 * 
  - log(marker.liver.integrated$p_val_adj)

exact_cell = intersect(
  colnames(src.endoderm.liver.ext_v1.1.re),
  rownames(src.fg4.integrated_re@meta.data[src.fg4.integrated_re$cluster.v06.26.re_correct_refine%in%c("FG.4-Lung/Stomach"),]))

src.endoderm.liver.ext_v1.1.re@reductions$umap@cell.embeddings = 
  src.endoderm.liver.ext_v1.1.re@reductions$umap_mnn@cell.embeddings[,c(1:2)]
DimPlot(src.endoderm.liver.ext_v1.1.re, 
        reduction = "umap_mnn", group.by = "cluster.v06.26.re..correct")
src.endoderm.liver.ext_v1.1.re.ext = 
  src.endoderm.liver.ext_v1.1.re[, setdiff(colnames(src.endoderm.liver.ext_v1.1.re),exact_cell)]

marker_test_al1 = marker.liver.integrated[marker.liver.integrated$cluster%in%'AL.1',]
gene_test_al1 = marker_test_al1[
  marker_test_al1$pct.ratio > quantile(marker_test_al1$pct.ratio, probs = 0.99),]$gene

marker_test_al2 = marker.liver.integrated[marker.liver.integrated$cluster%in%'AL.2',]
gene_test_al2 = marker_test_al2[
  marker_test_al2$pct.ratio > quantile(marker_test_al2$pct.ratio, probs = 0.99),]$gene

marker_test_al3 = marker.liver.integrated[marker.liver.integrated$cluster%in%'AL.3',]
gene_test_al3 = marker_test_al3[
  marker_test_al3$pct.ratio > quantile(marker_test_al3$pct.ratio, probs = 0.99),]$gene

marker_test_al13lv = marker.liver.integrated[marker.liver.integrated$cluster%in%'AL.3-Liver',]
gene_test_al13lv = marker_test_al13lv[
  marker_test_al13lv$pct.ratio > quantile(marker_test_al13lv$pct.ratio, probs = 0.99),]$gene

marker_test_fg4 = marker.liver.integrated[marker.liver.integrated$cluster%in%'FG.4',]
gene_test_fg4 = marker_test_fg4[
  marker_test_fg4$pct.ratio > quantile(marker_test_fg4$pct.ratio, probs = 0.985),]$gene

marker_test_fg4liver = marker.liver.integrated[marker.liver.integrated$cluster%in%'FG.4-Liver',]
gene_test_fg4liver = marker_test_fg4liver[
  marker_test_fg4liver$pct.ratio > quantile(marker_test_fg4liver$pct.ratio, probs = 0.99),]$gene

marker_test_al12lv = marker.liver.integrated[marker.liver.integrated$cluster%in%'AL.1/2-Liver',]
gene_test_al12lv = marker_test_al12lv[
  marker_test_al12lv$pct.ratio > quantile(marker_test_al12lv$pct.ratio, probs = 0.97),]$gene

marker_test_liver = marker.liver.integrated[marker.liver.integrated$cluster%in%'Liver',]
gene_test_liver = marker_test_liver[
  marker_test_liver$pct.ratio == Inf,]$gene

for(type in c(#"al1",
              #"al2",
              #'al3',
              #'fg4', 
              "al12lv",
              #"al13lv",
              #"fg4liver",
              #"liver",
              c())){
  pdf(paste("organ_development_re_v240115/organ/pattern_of_pathway.liver_",type,".pdf",sep=''),9,7)
  gene_list = get(paste("gene_test_",type,sep=""))
  for(gene in gene_list){
    print(
      FeaturePlot(src.endoderm.liver.ext_v1.1.re.ext, reduction = "umap",
                  features = gene, pt.size = 1.25) + 
        theme(aspect.ratio = 1))
  }
  dev.off()
}



# Stomach
src.endoderm.sto.ext_v1.1.re = SetIdent(src.endoderm.sto.ext_v1.1.re, 
                                        value = src.endoderm.sto.ext_v1.1.re$cluster.v06.26.re..merge)
marker.sto.integrated = FindAllMarkers(src.endoderm.sto.ext_v1.1.re)
marker.sto.integrated$pct.ratio =
  marker.sto.integrated$pct.1 / marker.sto.integrated$pct.2 * 
  - log(marker.sto.integrated$p_val_adj)

src.endoderm.sto.ext_v1.1.re@reductions$umap@cell.embeddings = 
  src.endoderm.sto.ext_v1.1.re@reductions$umap_integrated@cell.embeddings[,c(1:2)]
DimPlot(src.endoderm.sto.ext_v1.1.re,
        reduction = "umap_integrated", group.by = "Time")

marker_test_mg1 = marker.sto.integrated[marker.sto.integrated$cluster%in%'MG.1',]
gene_test_mg1 = marker_test_mg1[
  marker_test_mg1$pct.ratio > quantile(marker_test_mg1$pct.ratio, probs = 0.994),]$gene

marker_test_mg3 = marker.sto.integrated[marker.sto.integrated$cluster%in%'MG.3',]
gene_test_mg3 = marker_test_mg3[
  marker_test_mg3$pct.ratio > quantile(marker_test_mg3$pct.ratio, probs = 0.994),]$gene

marker_test_mg3sto = marker.sto.integrated[marker.sto.integrated$cluster%in%'MG.3.A',]
gene_test_mg3sto = marker_test_mg3sto[
  marker_test_mg3sto$pct.ratio > quantile(marker_test_mg3sto$pct.ratio, probs = 0.994),]$gene

marker_test_fg4 = marker.sto.integrated[marker.sto.integrated$cluster%in%'FG.4',]
gene_test_fg4 = marker_test_fg4[
  marker_test_fg4$pct.ratio > quantile(marker_test_fg4$pct.ratio, probs = 0.994),]$gene

marker_test_fg4sto = marker.sto.integrated[marker.sto.integrated$cluster%in%'FG.4-Lung/Stomach',]
gene_test_fg4sto = marker_test_fg4sto[
  marker_test_fg4sto$pct.ratio > quantile(marker_test_fg4sto$pct.ratio, probs = 0.98),]$gene

marker_test_sto = marker.sto.integrated[marker.sto.integrated$cluster%in%'Stomach',]
gene_test_sto = marker_test_sto[
  marker_test_sto$pct.ratio > quantile(marker_test_sto$pct.ratio, probs = 0.98),]$gene

pdf("organ_development_re_v240115/organ/pattern_of_pathway.sto.re.pdf",9,7)
for(type in c("fg4sto","sto")){
  gene_list = get(paste("gene_test_",type,sep=""))
  for(gene in gene_list){
    print(
      FeaturePlot(src.endoderm.sto.ext_v1.1.re, reduction = "umap",
                  features = gene, pt.size = 1.25) + 
        theme(aspect.ratio = 1))
  }
}
dev.off()


# small intestine 1
src.endoderm.sm1.ext_v1.1.re = SetIdent(src.endoderm.sm1.ext_v1.1.re, 
                                        value = src.endoderm.sm1.ext_v1.1.re$cluster.v06.26.re..merge)
marker.sm1.integrated = FindAllMarkers(src.endoderm.sm1.ext_v1.1.re)
marker.sm1.integrated$pct.ratio =
  marker.sm1.integrated$pct.1 / marker.sm1.integrated$pct.2 * 
  - log(marker.sm1.integrated$p_val_adj)

src.endoderm.sm1.ext_v1.1.re@reductions$umap@cell.embeddings = 
  src.endoderm.sm1.ext_v1.1.re@reductions$umap_integrated@cell.embeddings
DimPlot(src.endoderm.sm1.ext_v1.1.re, 
        reduction = "umap_integrated", group.by = "Time")

marker_test_mg1 = marker.sm1.integrated[marker.sm1.integrated$cluster%in%'MG.1',]
gene_test_mg1 = intersect(
  marker_test_mg1[
    marker_test_mg1$pct.ratio > quantile(marker_test_mg1$pct.ratio, probs = 0.91),]$gene,
  gi[gi$TF%in%T,]$Symbol2)

marker_test_mg2 = marker.sm1.integrated[marker.sm1.integrated$cluster%in%'MG.2',]
gene_test_mg2 = intersect(
  marker_test_mg2[
    marker_test_mg2$pct.ratio > quantile(marker_test_mg2$pct.ratio, probs = 0.91),]$gene,
  gi[gi$TF%in%T,]$Symbol2)

marker_test_mg3 = marker.sm1.integrated[marker.sm1.integrated$cluster%in%'MG.3',]
gene_test_mg3 = marker_test_mg3[
  marker_test_mg3$pct.ratio > quantile(marker_test_mg3$pct.ratio, probs = 0.992),]$gene

marker_test_mg3sm1 = marker.sm1.integrated[marker.sm1.integrated$cluster%in%'MG.3.P',]
gene_test_mg3sm1 = marker_test_mg3sm1[
  marker_test_mg3sm1$pct.ratio > quantile(marker_test_mg3sm1$pct.ratio, probs = 0.992),]$gene

marker_test_al3 = marker.sm1.integrated[marker.sm1.integrated$cluster%in%'AL.3',]
gene_test_al3 = marker_test_al3[
  marker_test_al3$pct.ratio > quantile(marker_test_al3$pct.ratio, probs = 0.994),]$gene

marker_test_al3sm1 = marker.sm1.integrated[marker.sm1.integrated$cluster%in%'AL.3-Small.intestine.1',]
gene_test_al3sm1 = marker_test_al3sm1[
  marker_test_al3sm1$pct.ratio > quantile(marker_test_al3sm1$pct.ratio, probs = 0.994),]$gene

marker_test_sm1 = marker.sm1.integrated[marker.sm1.integrated$cluster%in%'Small.intestine.1',]
gene_test_sm1 = intersect(
  marker_test_sm1[
    marker_test_sm1$pct.ratio > quantile(marker_test_sm1$pct.ratio, probs = 0.9),]$gene,
  gi[gi$TF%in%T,]$Symbol2)

for(type in c("mg1",
              "mg2",
              #"mg3",
              #"al3",
              #"mg3sm1",
              #"al3sm1",
              #"sm1",
              c())){
  pdf(paste("organ_development_re_v240115/organ/pattern_of_pathway.sm1_",type,".pdf",sep=""),9,7)
  gene_list = get(paste("gene_test_",type,sep=""))
  for(gene in gene_list){
    print(
      FeaturePlot(src.endoderm.sm1.ext_v1.1.re, reduction = "umap",
                  features = gene, pt.size = 1.25) + 
        theme(aspect.ratio = 1))
  }
  dev.off()
}


# small intestine 2
src.endoderm.sm2.ext_v1.1.re = SetIdent(src.endoderm.sm2.ext_v1.1.re, 
                                        value = src.endoderm.sm2.ext_v1.1.re$cluster.v06.26.re..merge)
marker.sm2.integrated = FindAllMarkers(src.endoderm.sm2.ext_v1.1.re)
marker.sm2.integrated$pct.ratio =
  marker.sm2.integrated$pct.1 / marker.sm2.integrated$pct.2 * 
  - log(marker.sm2.integrated$p_val_adj)

src.endoderm.sm2.ext_v1.1.re@reductions$umap@cell.embeddings = 
  src.endoderm.sm2.ext_v1.1.re@reductions$umap_mnn@cell.embeddings[,c(1:ncol(src.endoderm.sm2.ext_v1.1.re[["umap"]]@cell.embeddings))]
DimPlot(src.endoderm.sm2.ext_v1.1.re,
        reduction = "umap_mnn", group.by = "Time")

marker_test_mg2 = marker.sm2.integrated[marker.sm2.integrated$cluster%in%'MG.2',]
gene_test_mg2 = intersect(
  marker_test_mg2[
    marker_test_mg2$pct.ratio > quantile(marker_test_mg2$pct.ratio, probs = 0.9),]$gene,
  gi[gi$TF%in%T,]$Symbol2)

marker_test_mg3 = marker.sm2.integrated[marker.sm2.integrated$cluster%in%'MG.3',]
gene_test_mg3 = marker_test_mg3[
  marker_test_mg3$pct.ratio > quantile(marker_test_mg3$pct.ratio, probs = 0.988),]$gene

marker_test_mg3sm2 = marker.sm2.integrated[marker.sm2.integrated$cluster%in%'MG.3.P',]
gene_test_mg3sm2 = intersect(
  marker_test_mg3sm2[
    marker_test_mg3sm2$pct.ratio > quantile(marker_test_mg3sm2$pct.ratio, probs = 0.96),]$gene,
  gi[gi$TF%in%T,]$Symbol2)

marker_test_hg1 = marker.sm2.integrated[marker.sm2.integrated$cluster%in%'HG.1',]
gene_test_hg1 = marker_test_hg1[
  marker_test_hg1$pct.ratio > quantile(marker_test_hg1$pct.ratio, probs = 0.994),]$gene

marker_test_sm2 = marker.sm2.integrated[marker.sm2.integrated$cluster%in%'Small.intestine.2',]
gene_test_sm2 = marker_test_sm2[
    marker_test_sm2$pct.ratio > quantile(marker_test_sm2$pct.ratio, probs = 0.966),]$gene

for(type in c(#"hg1",
              "mg2"
              #,"mg3","mg3sm2","sm2"
              )){
  pdf(paste("organ_development_re_v240115/organ/pattern_of_pathway.sm2_",type,".pdf",sep=""),9,7)
  gene_list = get(paste("gene_test_",type,sep=""))
  for(gene in gene_list){
    print(
      FeaturePlot(src.endoderm.sm2.ext_v1.1.re, reduction = "umap",
                  features = gene, pt.size = 1.25) + 
        theme(aspect.ratio = 1))
  }
  dev.off()
}



# large intestine 1
src.endoderm.lar1.ext_v1.1.re = SetIdent(src.endoderm.lar1.ext_v1.1.re, 
                                         value = src.endoderm.lar1.ext_v1.1.re$cluster.v06.26.re..merge)
marker.lar1.integrated = FindAllMarkers(src.endoderm.lar1.ext_v1.1.re)
marker.lar1.integrated$pct.ratio =
  marker.lar1.integrated$pct.1 / marker.lar1.integrated$pct.2 * 
  - log(marker.lar1.integrated$p_val_adj)

src.endoderm.lar1.ext_v1.1.re@reductions$umap@cell.embeddings = 
  src.endoderm.lar1.ext_v1.1.re@reductions$umap_mnn@cell.embeddings[,1:3]
DimPlot(src.endoderm.lar1.ext_v1.1.re,
        reduction = "umap_mnn", group.by = "Time")


marker_test_hg1 = marker.lar1.integrated[marker.lar1.integrated$cluster%in%'HG.1',]
gene_test_hg1 = marker_test_hg1[
  marker_test_hg1$pct.ratio > quantile(marker_test_hg1$pct.ratio, probs = 0.99),]$gene

marker_test_hg2 = marker.lar1.integrated[marker.lar1.integrated$cluster%in%'HG.2',]
gene_test_hg2 = marker_test_hg2[
  marker_test_hg2$pct.ratio > quantile(marker_test_hg2$pct.ratio, probs = 0.99),]$gene

marker_test_lar1 = marker.lar1.integrated[marker.lar1.integrated$cluster%in%'Large.intestine.1',]
gene_test_lar1 = marker_test_lar1[
  marker_test_lar1$pct.ratio > quantile(marker_test_lar1$pct.ratio, probs = 0.99),]$gene

pdf("organ_development_re_v240115/organ/pattern_of_pathway.lar1.pdf",9,7)
for(type in c("hg1","hg2","lar1")){
  gene_list = get(paste("gene_test_",type,sep=""))
  for(gene in gene_list){
    print(
      FeaturePlot(src.endoderm.lar1.ext_v1.1.re, reduction = "umap",
                  features = gene, pt.size = 1.25) + 
        theme(aspect.ratio = 1))
  }
}
dev.off()


pdf("organ_development_re_v240115/organ/pattern_of_fg4.pdf",9,7)
DimPlot(src.fg4.integrated_re, reduction = "umap_mnn")
src.fg4.integrated_re@reductions$umap@cell.embeddings = 
  src.fg4.integrated_re@reductions$umap_mnn@cell.embeddings 
for(gene in c("Prox1","Nkx2-3","Pax9","Mnx1","Hhex","Irx1","Pyy","Gata3","Shh")){
  print(
    FeaturePlot(src.fg4.integrated_re, 
                features = gene, reduction ="umap", pt.size = 1.25) + 
      theme(aspect.ratio = 1))
}
dev.off()


