import scanpy as sc
import os,sys,glob
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt

dataset='dataset'

## QC
adata = sc.read_h5ad(file)

# mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-') 
# ribosomal genes
adata.var['ribo'] = adata.var_names.str.startswith(("RPS","RPL"))
# hemoglobin genes.
adata.var['hb'] = adata.var_names.str.contains(("^HB[^(P)]"))

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

q1=np.quantile(adata.obs['n_genes_by_counts'], .25)
q3=np.quantile(adata.obs['n_genes_by_counts'], .75)
print(q1, q3)
upper_gene=q3+3*(q3-q1)
under_gene=q1-3*(q3-q1)
print(under_gene, upper_gene)

q1=np.quantile(adata.obs['total_counts'], .25)
q3=np.quantile(adata.obs['total_counts'], .75)
print(q1, q3)
upper_total=q3+3*(q3-q1)
under_total=q1-3*(q3-q1)
print(under_total, upper_total)

q1=np.quantile(adata.obs['pct_counts_mt'], .25)
q3=np.quantile(adata.obs['pct_counts_mt'], .75)
print(q1, q3)
upper_mt=q3+3*(q3-q1)
under_mt=q1-3*(q3-q1)
print(under_mt, upper_mt)

adata = adata[adata.obs.n_genes_by_counts > under_gene, :]
adata = adata[adata.obs.n_genes_by_counts < upper_gene, :]
adata = adata[adata.obs.total_counts > under_total, :]
adata = adata[adata.obs.total_counts < upper_total, :]
adata = adata[adata.obs.pct_counts_mt > under_mt, :]
adata = adata[adata.obs.pct_counts_mt < upper_mt, :]

sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_genes(adata, min_cells=10)
adata = adata[adata.obs.pct_counts_mt < 5, :]

## Doublet
# split per batch into new objects.
batches = adata.obs['batch'].cat.categories.tolist()
alldata = {}
for batch in batches:
    tmp = adata[adata.obs['batch'] == batch,]
    print(batch, ":", tmp.shape[0], " cells")
    scrub = scr.Scrublet(tmp.X)
    out = scrub.scrub_doublets(verbose=False, n_prin_comps = 20)
    alldata[batch] = pd.DataFrame({'doublet_score':out[0],'predicted_doublets':out[1]},index = tmp.obs.index)
    print(alldata[batch].predicted_doublets.sum(), " predicted_doublets")
    
# add predictions to the adata object.
scrub_pred = pd.concat(adata.values())
adata.obs['doublet_scores'] = scrub_pred['doublet_score'] 
adata.obs['predicted_doublets'] = scrub_pred['predicted_doublets'] 

# also revert back to the raw counts as the main matrix in adata
print(adata.shape)
adata = adata[adata.obs['doublet_info'] == 'False',:]
print(adata.shape)

## Convert h5ad to Seurat

class_list=adata.obs.mainclass.unique()

for cls in class_list :
    print(cls)
    
    try: 
        os.mkdir(path_dir + cls)
    except OSError as error: 
        print(error) 
        
    adata1=adata[adata.obs.mainclass==cls]
    
    with open(path_dir + cls + '/barcodes.tsv', 'w') as f:
        for item in adata1.obs_names:
            f.write(item + '\n')
            
    with open(path_dir + cls + '/features.tsv', 'w') as f:
        for item in ['\t'.join([x,x,'Gene Expression']) for x in adata1.var_names]:
            f.write(item + '\n')
    
    io.mmwrite(path_dir + cls + '/matrix', adata1.X.T)
    
    adata1.obs.to_csv( path_dir + cls + '/metadata.csv')






