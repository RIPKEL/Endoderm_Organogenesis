setwd("/data/yupeilong/NewBioinformatic/endoderm_organogenesis/")

#============================================================================================================
#--project: Spatiotemporal and genetic cell lineage tracing of endodermal organogenesis at single-cell resolution 
#============================================================================================================

#============================================================================================================
#>>>> Outline of project
#============================================================================================================
#> 
#> 0.my function source and color set
#> 
#> 1.pre-processing
#>  
#> 2.gene co-expression network
#> 
#> 3.integration of 10x with Smart2
#> 
#> 4.inference on lineage affinity index 
#> 
#> 5.simulation of lineage affinity index 
#> 
#> 6.trajectory construction of endoderm
#> 
#> 7.gene dynamic heatmap of endoderm trajectories
#> 
#> 8.optimal transport (Python+R)
#> 
#> 9.signaling pathways analysis 
#> 
#> 10.integration of organ trajectories
#> 
#> 11.Tracing_by_Ngn3Cre_Ngn3GFP
#============================================================================================================


#============================================================================================================
#-->>>> Note !!!:
#============================================================================================================

#> The names of cell types in the original analysis code are different from those in the final paper, 
#> and the corresponding relationships are as follows. 

#> It should be noted that these differences are only in the naming of cells 
#> and do not affect the actual properties of the cells or the final conclusions.
#> 
#> Original code ::  Paper
#> 
#> (1) Intermeidate state
#> (1.1) FG.4-Lung/Stomach :: FG.4.1
#> (1.2) FG.4-Liver :: FG.4.2
#> (1.3) AL.3-EHBD.VP :: AL3.2
#> (1.4) AL.3-Small.intestine.1 :: AL3.3
#> (1.5) AL.3-Liver :: AL3.1
#> (1.6) AL.1/2-Liver :: pre-Liver
#> (1.7) MG.3.A/M / MG.3.A / MG.3.M :: MG.3.1
#> (1.8) MG.3.P :: MG.3.2
#> (1.9) HG.1-Large.intestine.2 :: pre-Lagre.intestine.2
#> 
#> (2) Organ state
#> (2.1) Phayrnx.organ.2 :: Pharynx.organ.1
#> (2.2) Phayrnx.organ.1 :: Pharynx.organ.2
#> (2.3) Phayrnx.organ.4 :: Pharynx.organ.3
#> (2.4) Phayrnx.organ.5 :: Pharynx.organ.4
#> (2.5) Phayrnx.organ.3 :: Pharynx.organ.5
#> (2.6) EP.1 :: PEP (pancreatic endocrine progentior)
#> (2.7) EP.2 :: alpha cell (1st wave)

#============================================================================================================


#============================================================================================================
#-->>>> Used R-package 
#============================================================================================================
#>> (1) R-version :R 4.1.2
#>> (2) package version
# AnnotationDbi	1.56.2
# base	4.1.2
# batchelor	1.10.0
# Biobase	2.54.0
# BiocGenerics	0.40.0
# Category	2.60.0
# ChIPseeker	1.30.3
# clusterProfiler	4.2.2
# cmdstanr	0.8.1
# data.table	1.14.8
# datasets	4.1.2
# DESeq2	1.34.0
# doParallel	1.0.17
# dplyr	1.1.2
# FactoMineR	2.8
# foreach	1.5.2
# furrr	0.3.1
# futile.logger	1.4.3
# future	1.32.0
# genomation	1.26.0
# GenomeInfoDb	1.30.1
# GenomicRanges	1.46.1
# ggplot2	3.4.2
# ggrepel	0.9.3
# ggsci	3.0.0
# ggsignif	0.6.4
# ggvenn	0.1.10
# GO.db	3.14.0
# GOstats	2.60.0
# gplots	3.1.3
# graph	1.78.0
# graphics	4.1.2
# grDevices	4.1.2
# grid	4.1.2
# hash	2.2.6.3
# igraph	1.5.0
# IRanges	2.28.0
# ISLR2	1.3-2
# iterators	1.0.14
# KEGGREST	1.34.0
# lattice	0.21-8
# limma	3.50.3
# MASS	7.3-60
# Matrix	1.5-4.1
# MatrixGenerics	1.6.0
# matrixStats	1.0.0
# methods	4.1.2
# mixOmics	6.18.1
# nucleR	2.26.0
# org.Mm.eg.db	3.14.0
# parallel	4.1.2
# patchwork	1.1.2
# plotly	4.10.2
# princurve	2.1.6
# processx	3.8.1
# progress	1.2.2
# progressr	0.13.0
# RColorBrewer	1.1-3
# readr	2.1.4
# readxl	1.4.2
# reshape2	1.4.4
# reticulate	1.3
# rgl	1.1.3
# rlang	1.1.1
# rmspc	1.0.0
# S4Vectors	0.32.4
# scales	1.2.1
# scatterplot3d	0.3-44
# scran	1.22.1
# scuttle	1.4.0
# Seurat	4.3.0.1
# SeuratObject	4.1.3
# SeuratWrappers	0.3.0
# Signac	1.10.0
# SingleCellExperiment	1.16.0
# stats	4.1.2
# stats4	4.1.2
# stringr	1.5.0
# SummarizedExperiment	1.24.0
# tools	4.1.2
# tradeSeq	1.8.0
# utils	4.1.2
# VennDiagram	1.7.3
# Vennerable	3.1.0.9000
# viridis	0.6.3
# viridisLite	0.4.2
#============================================================================================================













