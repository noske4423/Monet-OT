import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pylab as pl
import anndata as ad
import math
import matplotlib.colors as mcolors
import ot
from matplotlib import ticker, cm, colors
import pickle
import random
from statistics import stdev, variance, median
from PIL import Image
import scipy.stats
import define_function
import csv
from datetime import datetime
import os

# Set up the folder to save the images
now = datetime.now()
date_str = now.strftime("%Y-%m-%d")
folder_name='facs_female'
folder_path = f'./image/{date_str}/{folder_name}'
if not os.path.exists(folder_path):
    os.makedirs(folder_path)
original_dir = os.getcwd()  # Save the original directory

print(1)
GENE_facs = 'data/facs/tabula-muris-senis-facs-processed-official-annotations.h5ad'
adata = sc.read_h5ad(GENE_facs)  # read the data
print(2)
print(adata)
print(adata.obs)
print(adata.var)

# adata.X: the data matrix X stores the feature expression values in the rows
# (observations/cells) and the features in the columns.
# adata.obs: stores observation annotation (usually one row per cell) as data frame.
# adata.var: stores variable annotation (usually one column per gene)

adata = adata[adata.obs.sex == 'female']
print(3)
mean_expression = adata.X.mean(axis=0)  # mean expression of each gene
fig = pl.figure(figsize=(3, 3))
pl.plot(mean_expression.tolist()[0])
pl.xlabel('Genes')
pl.ylabel('Mean expression')
pl.title('Mean expression of each gene')
pl.show()
fig.savefig(f'{folder_path}/{folder_name}_mean_expression.png')  # save the figure

adata = define_function.sort_adata_by_attribute(adata, 'tissue')  # sort by tissue
age_dict = define_function.split_adata_by_attribute(adata, 'age')  # split by age
ages = list(age_dict.keys())  # get the list of ages
tissue_dict = {}
for age in ages:
    tissue_dict[age] = define_function.split_adata_by_attribute(age_dict[age], 'tissue')  # split by tissue

os.chdir(folder_path)  # Change the current directory to the desired folder

# save the PCA and umap plot for each tissue
for tissue in list(tissue_dict[ages[0]].keys()):
