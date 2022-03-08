###
#   File name  : RunScenic.py
#   Author     : Hyunjin Kim
#   Date       : Feb 19, 2022
#   Email      : firadazer@gmail.com
#   Purpose    : Run Scenic example, load SJCAR19 data and run Scenic
#
#   https://pyscenic.readthedocs.io/en/latest/tutorial.html
#
#   Instruction
#               1. Run RunScenic.py
#               3. The results will be generated under the specified directory
###

import os
import glob
import pickle
import pandas as pd
import numpy as np
import scanpy
import scanpy as sc
import loompy as lp

from dask.diagnostics import ProgressBar

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

import seaborn as sns


if __name__ == '__main__':
    ### set parameters
    DATABASES_GLOB = os.path.join("Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/hg19-*.mc9nr.feather")

    db_fnames = glob.glob(DATABASES_GLOB)

    def name(fname):
        return os.path.splitext(os.path.basename(fname))[0]

    dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]

    MOTIF_ANNOTATIONS_FNAME = os.path.join(r"Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/motifs-v9-nr.hgnc-m0.001-o0.0.tbl")

    ### load the data
    adata = sc.read_10x_h5("C:/Users/hkim8/Documents/SJCAR/data/pbmc_10k_v3_filtered_feature_bc_matrix.h5")
    tf_names = load_tf_names("C:/Users/hkim8/Documents/SJCAR/data/tf_list/hs_hgnc_tfs.txt")

    ### preprocess the anndata
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.log1p(adata)
    adata.to_df().to_csv("C:/Users/hkim8/Documents/SJCAR/data/pbmc_10k_v3_filtered_feature_bc_matrix.csv")

    ### load the data again
    ex_matrix = pd.read_csv("C:/Users/hkim8/Documents/SJCAR/data/pbmc_10k_v3_filtered_feature_bc_matrix.csv", header=0, index_col=0)
    ex_matrix.shape

    ### run scenic
    ### since this grnboost2 makes an error, we used arboreto_with_multiprocessing.py instead
    adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)

    ### load adj matrices
    precursor_adjacencies = pd.read_csv("Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/precursor_adj.tsv", sep="\t")
    precursor_mat = sc.read_loom("Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/precursor_seurat.loom")
    precursor_mat = precursor_mat.to_df()

    ### set modules from grnboost2
    modules = list(modules_from_adjacencies(precursor_adjacencies, precursor_mat))

    ### Calculate a list of enriched motifs and the corresponding target genes for all modules.
    with ProgressBar():
        df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)

    # Create regulons from this table of enriched motifs.
    regulons = df2regulons(df)

    # Save the enriched motifs and the discovered regulons to disk.
    df.to_csv("Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/precursor_motifs.csv")
    with open("Z:/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/precursor_regulons.p", "wb") as f:
        pickle.dump(regulons, f)

