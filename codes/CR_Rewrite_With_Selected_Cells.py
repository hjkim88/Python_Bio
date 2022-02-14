###
#   File name  : CR_Rewrite_With_Selected_Cells.py
#   Author     : Hyunjin Kim
#   Date       : Feb 12, 2022
#   Email      : firadazer@gmail.com
#   Purpose    : Read scRNA-Seq Cell Ranger output and save it with selected cells
#                This is because the file size is too big to use in R (unable to load directly)
#
#   Instruction
#               1. Run ReadCR_WriteH5AD.py
#               3. The results will be generated under the spcified directory
###

### import modules
import numpy as np
import pandas as pd
import os
import re
import gzip
import scanpy as sc

### load barcodes
target_path = "Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/JCC212_SJCAR19_AggregOct2020ns/"
f = [(target_path + x) for x in os.listdir(target_path)]
barcode_path = []
for x in f:
    if re.search(r'barcodes.tsv.gz$', x):
        barcode_path.append(x)

total_barcodes=[]
for x in barcode_path:
    barcodes = gzip.open(x, "rb")
    total_barcodes.append(barcodes.read().decode("utf-8").splitlines())

### load the data
adata = sc.read_10x_mtx(path="Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/JCC212_SJCAR19_AggregOct2020ns/filtered_feature_bc_matrix/")
adata_barcodes = list(adata.obs_names)

### find indices in a list
def indices(lst, element):
    result = []
    offset = -1
    while True:
        try:
            offset = lst.index(element, offset+1)
        except ValueError:
            return result
        result.append(offset)

### barcodes are different
### total_barcodes are from indivisual samples (all -1)
### and adata_barcodes has different '-*' after the actual barcodes (ex. -1, -2, -3, -4, etc.)
### so subset the adata_barcodes based on their suffix and find which are the barcodes of our interest
actual_bcd = []
suffices = []
for x in adata_barcodes:
    c = x.split('-')
    actual_bcd.append(c[0])
    suffices.append(c[1])

### save the indices of the adata suffices
unique_suffix = sorted(set(suffices))
unique_suffix = sorted([int(x) for x in unique_suffix])
unique_suffix = [str(x) for x in unique_suffix]

suffix_idx = []
for x in range(len(unique_suffix)):
    suffix_idx.append(indices(suffices, unique_suffix[x]))

### find barcodes that are matched
### there are 6 samples matched out of 8 - adata has only 6 healthy donors while we have barcodes of 8
total_barcodes_idx = []
for x in range(len(total_barcodes)):
    for y in range(len(unique_suffix)):
        if len(total_barcodes[x]) == len(suffix_idx[y]):
            total_barcodes_idx.append(y)
            break

### get the combined indices of the adata_barcodes
target_barcodes = []
for x in total_barcodes_idx:
    temp = suffix_idx[x]
    for y in temp:
        target_barcodes.append(adata_barcodes[y])

### subset the data with selected cells
adata2 = adata[target_barcodes].copy()

### save the data as h5ad
adata2.write("Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/JCC212_SJCAR19_AggregOct2020ns/filtered_feature_bc_matrix/mat_healthy_donors.h5ad")
