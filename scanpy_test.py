import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pylab as pl
import anndata as ad

GENE_facs = 'data/facs/tabula-muris-senis-facs-processed-official-annotations.h5ad'
adata = sc.read_h5ad(GENE_facs)  # read the data

sc.pl.highest_expr_genes(adata, n_top=30)

