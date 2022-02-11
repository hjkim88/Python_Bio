###
#   File name  : ReadCR_WriteH5AD.py
#   Author     : Hyunjin Kim
#   Date       : Feb 10, 2022
#   Email      : firadazer@gmail.com
#   Purpose    : Read scRNA-Seq Cell Ranger output and save it as h5ad file
#                This is because the file size is too big to use in R (unable to load directly)
#
#   Instruction
#               1. import ReadCR_WriteH5AD.py
#               2. Run the function ReadFastq_RmReps.start()
#               3. The results will be generated under the same directory as the input
###

### import modules
import numpy as np
import pandas as pd
import scanpy as sc

### load the data
adata = sc.read_10x_mtx(path="Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/JCC212_SJCAR19_AggregOct2020ns/filtered_feature_bc_matrix/")

### save the data as h5ad
adata.write("Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/JCC212_SJCAR19_AggregOct2020ns/filtered_feature_bc_matrix/mat.h5ad")