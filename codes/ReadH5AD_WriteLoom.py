###
#   File name  : ReadH5AD_WriteLoom.py
#   Author     : Hyunjin Kim
#   Date       : Jul 7, 2022
#   Email      : firadazer@gmail.com
#   Purpose    : Read scRNA-Seq H5AD file and save it as loom file
#
#   Instruction
#               1. Run ReadH5AD_WriteLoom.py
#               3. The results will be generated under the specified directory
###

### import modules
import scanpy as sc

### read h5ad
adata = sc.read_h5ad("/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_gex.h5ad")
adata = sc.read_h5ad("/Users/hyunjin.kim2/Downloads/GSE193816_all_data_raw_counts.h5ad")
adata = sc.read_h5ad("/Users/hyunjin.kim2/Downloads/GSE193816_all_cells_data.h5ad")

### save as loom data
adata.write_loom("/Users/hyunjin.kim2/Documents/SimpleTasks/results/PID5202/MSR1/concat_gex.loom",
                 write_obsm_varm=True)
adata.write_loom("/Users/hyunjin.kim2/Downloads/GSE193816_all_data_raw_counts.loom",
                 write_obsm_varm=True)
adata.write_loom("/Users/hyunjin.kim2/Downloads/GSE193816_all_cells_data.loom",
                 write_obsm_varm=True)