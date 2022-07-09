###
#   File name  : ReadH5AD_WriteObjects.py
#   Author     : Hyunjin Kim
#   Date       : Jul 7, 2022
#   Email      : firadazer@gmail.com
#   Purpose    : Read scRNA-Seq H5AD file and save it as raw counts/scaled counts/meta data/PCA/UMAP
#
#   Instruction
#               1. Run ReadH5AD_WriteObjects.py
#               3. The results will be generated under the specified directory
###

### import modules
import scanpy as sc
import pandas as pd
import anndata as ad

### read h5ad
adata = sc.read_h5ad("/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_gex.h5ad")

### write raw counts
mat = adata.raw.X.toarray()
pd.DataFrame(data=mat, index=adata.obs_names, columns=adata.raw.var_names).to_csv("/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_raw_counts.csv")

### write scaled counts
pd.DataFrame(data=adata.X, index=adata.obs_names, columns=adata.var_names).to_csv("/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_scaled_counts.csv")

### write meta data
pd.DataFrame(data=adata.obs, index=adata.obs.index, columns=adata.obs.columns).to_csv("/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_metadata.csv")

### write PCA/UMAP
dim_methods = list(adata.obsm.keys())

for method in dim_methods:
    print(method)
    mat = adata.obsm[method]
    pd.DataFrame(data=mat, index=adata.obsm.dim_names, columns=[(method.split("X_")[1] + str(x)) for x in range(1,(mat.shape[1]+1))]).to_csv(
        ("/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_" + method.split("X_")[1] + ".csv"))

### create new anndata with necessary parts only
adata2 = ad.AnnData(adata.raw.X)
adata2.raw = adata2
adata2.X = adata.X
# adata2.obs = adata.obs
# adata2.var = adata.var
adata2.obsm = adata.obsm

### write out the new anndata as h5ad
adata2.write("/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_gex_new.h5ad")
