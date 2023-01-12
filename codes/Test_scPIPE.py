
### load modules
import os
import sys
import scanpy as sc
import pandas as pd
# sys.path.append('C:/Users/hyunjin.kim2/Downloads/scpipe-master/packages/')
# import scPIPE

### load the data
data_dir="C:/Users/hyunjin.kim2/Documents/RProjects/SimpleTasks/data/PID4835/"

try:
    features = pd.read_csv(os.path.join(data_dir, 'features.tsv.gz'), sep='\t', header=None, index_col=0,
                           names=['gene_symbols', 'feature_types'])
except:
    try:
        features = pd.read_csv(os.path.join(data_dir, 'genes.tsv.gz'), sep='\t', header=None, index_col=0,
                               names=['gene_symbols'])
    except:
        features = pd.read_csv(os.path.join(data_dir, 'genes.tsv'), sep='\t', header=None, index_col=0,
                               names=['gene_symbols'])

try:
    barcodes = pd.read_csv(os.path.join(data_dir, 'barcodes.tsv.gz'), sep='\t', header=None, index_col=0, names=[''])
except:
    barcodes = pd.read_csv(os.path.join(data_dir, 'barcodes.tsv'), sep='\t', header=None, index_col=0, names=[''])

try:
    data = sc.read(os.path.join(data_dir, 'matrix.mtx.gz'), first_column_names=False).transpose()
except:
    data = sc.read(os.path.join(data_dir, 'matrix.mtx'), first_column_names=False).transpose()

# add feature and barcode information to anndata object
data.var = features
data.obs = barcodes


adata = sc.read_h5ad("C:/Users/hyunjin.kim2/Documents/RProjects/SimpleTasks/data/PID4835/PID4835_concat_data3_gex.h5ad")




