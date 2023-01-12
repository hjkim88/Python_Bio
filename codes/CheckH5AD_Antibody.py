###
#   File name  : CheckH5AD_Antibody.py
#   Author     : Hyunjin Kim
#   Date       : Jul 12, 2023
#   Email      : firadazer@gmail.com
#   Purpose    : Read scRNA-Seq H5AD file and check the antibody data
#
#   Instruction
#               1. Run CheckH5AD_Antibody.py
#               3. The results will be generated on the console
###

### import modules
import scanpy as sc
import pandas as pd
import anndata as ad

### read h5ad
adata = sc.read_h5ad("C:/Users/hyunjin.kim2/Documents/RProjects/SimpleTasks/data/PID4835/PID4835_concat_data_010423_gex.h5ad",
                     gex_only=False)
