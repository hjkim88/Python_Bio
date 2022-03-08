###
#   File name  : Run_hscscore.py
#   Author     : Hyunjin Kim
#   Date       : Mar 7, 2022
#   Email      : firadazer@gmail.com
#   Purpose    : Run hscscore on the GSM4645153 data
#
#   https://nbviewer.org/github/fionahamey/hscScore/blob/master/hscScore_demonstration_notebook.ipynb
#   https://zenodo.org/record/3332151
#
#   Instruction
#               1. Run Run_hscscore.py
#               3. The results will be generated under the specified directory
###

import pandas as pd
import numpy as np
import pickle
import sklearn
import platform

data_dir = 'C:/Users/hkim8/Documents/Hematopoiesis/data/GSM4645153/Zenodo_data/'


def total_count_normalise(count_matrix):
    """Normalise count matrix for input into hscScore model.
    Performs read depth normalisation normalising each cell so that normalised
    counts sum to the same value.

    Parameters
    ----------
    count_matrix : pandas dataframe
        Gene count matrix of dimension cells x genes with column names as genes
        and index as cell names

    Returns
    -------
    **norm_matrix** : pandas dataframe
        Normalised count matrix of dimension cells x genes
    """

    # Set the value normalised counts will sum to for each cell
    wilson_molo_genes_median_counts = 18704.5

    # Scale rows
    count_matrix_expression = np.array(count_matrix, dtype='float')
    counts_per_cell = np.sum(count_matrix_expression, axis=1)
    counts_per_cell += (counts_per_cell == 0)
    counts_per_cell /= wilson_molo_genes_median_counts
    norm_matrix_expression = count_matrix_expression / counts_per_cell[:, None]
    norm_matrix = pd.DataFrame(norm_matrix_expression, index=count_matrix.index,
                               columns=count_matrix.columns)
    # log + 1 transform the data
    norm_matrix = np.log(norm_matrix + 1)

    return norm_matrix

### load the trained model
hsc_score = pickle.load(open(data_dir + 'hscScore_model3.pkl', 'rb'))

### load the target count matrix
count_data = pd.read_csv(data_dir + 'GSM4645153.tsv', sep=' ')
count_data.head()

### load the genes of the trained model
model_genes = np.genfromtxt(data_dir + 'model_molo_genes.txt', dtype='str')
model_genes = model_genes.tolist()

### Subset to the same genes in the same order as used for training the model
### for hscsocre pipeline all the genes should be in the count data
### we are missing 4 genes, so we just put all zeros for those 4 genes
count_data_genes = count_data.columns.values.tolist()
common_genes = list(set(model_genes).intersection(count_data_genes))
common_genes.sort()
diff_genes = set(model_genes).difference(count_data_genes)
### add missing 4 genes to the count data with all zeros
for g in diff_genes:
    count_data[g] = 0
count_data_molo = count_data[model_genes]
count_data_molo.head()

### normalize the data
normalised_data_molo = total_count_normalise(count_data_molo)
normalised_data_molo.head()

### calculate the HSCscore
predicted_hsc_scores = hsc_score.predict(np.array(normalised_data_molo))
predicted_hsc_scores

### create a data frame with the hscscore and the cell barcodes
df = pd.DataFrame()
df['Barcode'] = normalised_data_molo.index.to_list()
df['hscScore'] = predicted_hsc_scores.tolist()

### save the result
df.to_csv(data_dir + 'GSM4645153_hscScore.csv', index=False)
