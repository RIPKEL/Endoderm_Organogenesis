#===============================================================================
#>> 8.optimal transport
#>>>> Python-part
#>
#>>>> (8.2) calculate OT transition matirx using wot in [Python] 
#===============================================================================

import wot
import os
import pandas as pd
import numpy as np
# import pegasus as pg
import igraph as ig
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
import anndata

sns.set_style("white")
s_genes=["Mcm5","Pcna","Tyms","Fen1","Mcm7","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl",
         "Prim1","Uhrf1","Cenpu","Hells","Rfc2","Polr1b","Nasp","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2",
         "Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm",
         "Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Mrpl36","E2f8"]
g2m_genes=["Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2","Cks1b","Mki67","Tmpo",
           "Cenpf","Tacc3","Pimreg","Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1",
           "Kif20b","Hjurp","Cdca3","Jpt1","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8",
           "Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2e3","Gas2l3",
           "Cbx5","Cenpa"]
os.chdir("D:/project180311/endoderm 10x/OT")
############################################
#              load
#############################################
adata = wot.io.read_dataset("../cluster/endoderm.loom")
adata.shape
select_gene = pd.read_table("../cluster/endoderm.merge.selectgene.tab",header=None)[0].values.tolist()
adata.var["highly_variable"]=adata.var["Accession"].isin(select_gene)
proliferation = pd.read_table("../cluster/endoderm.merge.cellcycle.score.tab",header=0)['Cell.cycle']
sns.histplot(proliferation)
plt.show(block=True)
# adata=anndata.read_h5ad("./wot.h5ad")
############################################
#               DR
#############################################
sc.tl.pca(adata)
g1 = sns.scatterplot(x=adata.obsm["X_pca"][:,0],
                     y=adata.obsm["X_pca"][:,1],
                     hue=adata.obs["Time"],
                     edgecolor="none")
g1.legend(bbox_to_anchor=(0.8,0.7), ncol=1)
plt.show(block=True)
adata.obs["day"]=adata.obs["Time"].str.replace("ss","").astype("float")/3-2
#################################################
#            cell cycle regression
##################################################
def logistic(x, L, k, x0=0):
    f = L / (1 + np.exp(-k * (x - x0)))
    return f
def gen_logistic(p, beta_max, beta_min, center, width):
    return beta_min + logistic(p, L=beta_max - beta_min, k=4 / width, x0=center)

def beta(p, beta_max=0.04, beta_min=0.016, center=-0.4):
    return gen_logistic(p, beta_max, beta_min, center, width=0.5)

sns.scatterplot(x=np.linspace(-2,2,1000),
                y=beta(np.linspace(-2,2,1000)))
plt.show(block=True)

birth = beta(proliferation)

nbrs = NearestNeighbors(n_neighbors=50, algorithm='ball_tree').fit(adata.obsm["X_pca"])
nbrs_matrix=nbrs.kneighbors_graph(adata.obsm["X_pca"]).toarray()

adata.obs["cell_growth_rate"]=[sum(birth[nbrs_matrix[i]==1])/50 for i in range(nbrs_matrix.shape[0])]

g2 = sns.scatterplot(x=adata.obsm["X_pca"][:,0],
                     y=adata.obsm["X_pca"][:,1],
                     hue=adata.obs["cell_growth_rate"],
                     edgecolor="none")
g2.legend(bbox_to_anchor=(0.8,0.7), ncol=1)
plt.show(block=True)
##############################################################3
#                    OT
###############################################################
ot_model = wot.ot.OTModel(adata,epsilon = 0.05)
trans_9_12=ot_model.compute_transport_map(1,2)
trans_12_15=ot_model.compute_transport_map(2,3)
trans_15_18=ot_model.compute_transport_map(3,4)
trans_18_21=ot_model.compute_transport_map(4,5)
trans_21_24=ot_model.compute_transport_map(5,6)
trans_24_27=ot_model.compute_transport_map(6,7)
adata.write("wot.h5ad")

trans_matrix_9_12=trans_9_12.X
trans_matrix_12_15=trans_12_15.X
trans_matrix_15_18=trans_15_18.X
trans_matrix_18_21=trans_18_21.X
trans_matrix_21_24=trans_21_24.X
trans_matrix_24_27=trans_24_27.X

np.savetxt("trans_matrix_9_12.csv", trans_matrix_9_12, delimiter=",")
np.savetxt("trans_matrix_12_15.csv", trans_matrix_12_15, delimiter=",")
np.savetxt("trans_matrix_15_18.csv", trans_matrix_15_18, delimiter=",")
np.savetxt("trans_matrix_18_21.csv", trans_matrix_18_21, delimiter=",")
np.savetxt("trans_matrix_21_24.csv", trans_matrix_21_24, delimiter=",")
np.savetxt("trans_matrix_24_27.csv", trans_matrix_24_27, delimiter=",")
