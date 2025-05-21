import scanpy as sc
import os,sys,glob
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
import pyreadr
from scipy import sparse, io
import matplotlib

### Load reference file
adata_sel = sc.read_h5ad(reference_file)

adata_sel.layers["counts"] = adata_sel.X.copy()
sc.pp.normalize_total(adata_sel, target_sum=1e4)
sc.pp.log1p(adata_sel)

sc.pp.highly_variable_genes(adata_sel, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=5000)
adata_sel = adata_sel[:, adata_sel.var.highly_variable]
sc.pp.scale(adata_sel, max_value=10)

sc.tl.pca(adata_sel, svd_solver='arpack')
sc.pp.neighbors(adata_sel, n_neighbors=100, n_pcs=20)

adata_sel.obs['batch'] = adata_sel.obs['batch'].astype('category')
adata_sel.obsm['X_pca'] = adata_sel.obsm['X_pca_harmony']

sc.pp.neighbors(adata_sel, n_neighbors=100, n_pcs=20)
sc.tl.umap(adata_sel)

### Load the target file
adata_new = sc.read_h5ad(singlebrain_file)
sc.pp.scale(adata_new, max_value=10)
sc.pp.pca(adata_new, n_comps=50)

adata_ref = adata_sel.copy()

var_names = adata_ref.var_names.intersection(adata_new.var_names)
print(len(var_names))

adata_ref = adata_ref[:, var_names]
adata_new = adata_new[:, var_names]

sc.pp.neighbors(adata_ref, n_neighbors=10, n_pcs=40)

# Map the new data onto the reference data
sc.tl.ingest(adata_new, adata_ref, obs='paper')

adata_new.obs.to_csv('~/projection_result.csv')



