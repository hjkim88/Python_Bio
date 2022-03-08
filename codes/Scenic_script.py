sudo mkdir /mnt/z
sudo mount -t drvfs Z: /mnt/z
# sudo umount /mnt/d

python /mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/arboreto_with_multiprocessing.py \
  /mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/precursor_seurat.loom \
  /mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/hs_hgnc_tfs.txt \
  --method grnboost2 \
  --output /mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/precursor_adj.tsv \
  --num_workers 20 \
  --seed 777

python /mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/arboreto_with_multiprocessing.py \
  /mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/non_precursor_seurat_subset.loom \
  /mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/hs_hgnc_tfs.txt \
  --method grnboost2 \
  --output /mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/non_precursor_subset_adj.tsv \
  --num_workers 20 \
  --seed 777

python /mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/arboreto_with_multiprocessing.py \
  /mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/combined_seurat.loom \
  /mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/hs_hgnc_tfs.txt \
  --method grnboost2 \
  --output /mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/combined_adj.tsv \
  --num_workers 20 \
  --seed 777

python /mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/arboreto_with_multiprocessing.py \
  /mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/cluster3_8_seurat.loom \
  /mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/hs_hgnc_tfs.txt \
  --method grnboost2 \
  --output /mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/cluster3_8_adj.tsv \
  --num_workers 20 \
  --seed 777

### now use those adj.tsv, run the following code in python - linux environment
sudo mkdir /mnt/z
sudo mount -t drvfs Z: /mnt/z

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

### set parameters
DATABASES_GLOB = os.path.join("/mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/hg19-*.mc9nr.feather")

db_fnames = glob.glob(DATABASES_GLOB)

def name(fname):
	return os.path.splitext(os.path.basename(fname))[0]

dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]

MOTIF_ANNOTATIONS_FNAME = os.path.join("/mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/motifs-v9-nr.hgnc-m0.001-o0.0.tbl")

### load adj matrices
precursor_adjacencies = pd.read_csv("/mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/precursor_adj.tsv", sep="\t")
precursor_mat = sc.read_loom("/mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/precursor_seurat.loom")
precursor_mat = precursor_mat.to_df()

### set modules from grnboost2
modules = list(modules_from_adjacencies(precursor_adjacencies, precursor_mat))

### Calculate a list of enriched motifs and the corresponding target genes for all modules.
with ProgressBar():
	df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)

### Create regulons from this table of enriched motifs.
regulons = df2regulons(df)

### Save the enriched motifs and the discovered regulons to disk.
df.to_csv("/mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/precursor_motifs.csv")
with open("/mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/precursor_regulons.p", "wb") as f:
	pickle.dump(regulons, f)

### characterize the different cells in a single-cell transcriptomics experiment
### via the enrichment of the previously discovered regulons.
### Enrichment of a regulon is measured as the Area Under the recovery Curve (AUC)
### of the genes that define this regulon.
auc_mtx = aucell(precursor_mat, regulons, num_workers=4)
auc_mtx.to_csv("/mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/precursor_aucell.csv")
sns.clustermap(auc_mtx, figsize=(8,8))


### Now it's non-precursor's turn

### load adj matrices
non_precursor_adjacencies = pd.read_csv("/mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/non_precursor_subset_adj.tsv", sep="\t")
non_precursor_mat = sc.read_loom("/mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/non_precursor_seurat_subset.loom")
non_precursor_mat = non_precursor_mat.to_df()

### set modules from grnboost2
modules2 = list(modules_from_adjacencies(non_precursor_adjacencies, non_precursor_mat))

### Calculate a list of enriched motifs and the corresponding target genes for all modules.
with ProgressBar():
	df2 = prune2df(dbs, modules2, MOTIF_ANNOTATIONS_FNAME)

### Create regulons from this table of enriched motifs.
regulons2 = df2regulons(df2)

### Save the enriched motifs and the discovered regulons to disk.
df2.to_csv("/mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/non_precursor_subset_motifs.csv")
with open("/mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/non_precursor_subset_regulons.p", "wb") as f:
	pickle.dump(regulons2, f)

### characterize the different cells in a single-cell transcriptomics experiment
### via the enrichment of the previously discovered regulons.
### Enrichment of a regulon is measured as the Area Under the recovery Curve (AUC)
### of the genes that define this regulon.
auc_mtx2 = aucell(non_precursor_mat, regulons2, num_workers=4)
auc_mtx2.to_csv("/mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/non_precursor_subset_aucell.csv")


### Now it's combined data's turn

### load adj matrices
combined_adjacencies = pd.read_csv("/mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/combined_adj.tsv", sep="\t")
combined_mat = sc.read_loom("/mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/combined_seurat.loom")
combined_mat = combined_mat.to_df()

### set modules from grnboost2
modules3 = list(modules_from_adjacencies(combined_adjacencies, combined_mat))

### Calculate a list of enriched motifs and the corresponding target genes for all modules.
with ProgressBar():
	df3 = prune2df(dbs, modules3, MOTIF_ANNOTATIONS_FNAME)

### Create regulons from this table of enriched motifs.
regulons3 = df2regulons(df3)

### Save the enriched motifs and the discovered regulons to disk.
df3.to_csv("/mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/combined_subset_motifs.csv")
with open("/mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/combined_regulons.p", "wb") as f:
	pickle.dump(regulons3, f)

### characterize the different cells in a single-cell transcriptomics experiment
### via the enrichment of the previously discovered regulons.
### Enrichment of a regulon is measured as the Area Under the recovery Curve (AUC)
### of the genes that define this regulon.
auc_mtx3 = aucell(combined_mat, regulons3, num_workers=4)
auc_mtx3.to_csv("/mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/combined_aucell.csv")



### Now it's cluster3+8's turn

### load adj matrices
cluster3_8_adjacencies = pd.read_csv("/mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/cluster3_8_adj.tsv", sep="\t")
cluster3_8_mat = sc.read_loom("/mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/cluster3_8_seurat.loom")
cluster3_8_mat = cluster3_8_mat.to_df()

### set modules from grnboost2
modules4 = list(modules_from_adjacencies(cluster3_8_adjacencies, cluster3_8_mat))

### Calculate a list of enriched motifs and the corresponding target genes for all modules.
with ProgressBar():
	df4 = prune2df(dbs, modules4, MOTIF_ANNOTATIONS_FNAME)

### Create regulons from this table of enriched motifs.
regulons4 = df2regulons(df4)

### Save the enriched motifs and the discovered regulons to disk.
df4.to_csv("/mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/cluster3_8_motifs.csv")
with open("/mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/cluster3_8_regulons.p", "wb") as f:
	pickle.dump(regulons4, f)

### characterize the different cells in a single-cell transcriptomics experiment
### via the enrichment of the previously discovered regulons.
### Enrichment of a regulon is measured as the Area Under the recovery Curve (AUC)
### of the genes that define this regulon.
auc_mtx4 = aucell(cluster3_8_mat, regulons4, num_workers=4)
auc_mtx4.to_csv("/mnt/z/ResearchHome/SharedResources/Immunoinformatics/hkim8/Scenic/cluster3_8_aucell.csv")
