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

# explain the code
#  1. Set up the folder to save the images
#  2. Read the data
#  3. Preprocessing
#  4. Plot the data
#  5. Save the image
#  6. Save the data

# Set up the folder to save the images
now = datetime.now()
date_str = now.strftime("%Y-%m-%d")
folder_name = 'facs_male'
folder_path = f'./image/{date_str}/{folder_name}'
if not os.path.exists(folder_path):
    os.makedirs(folder_path)

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
# adata.var: stores variable annotation (usually one column per gene) as data frame.

adata = adata[adata.obs.sex == 'male']  # select only male cells

# pca
first_component_pca = adata.obsm['X_pca'][:, 0]
second_component_pca = adata.obsm['X_pca'][:, 1]

max_value_pca_1 = first_component_pca.max()
min_value_pca_1 = first_component_pca.min()
max_value_pca_2 = second_component_pca.max()
min_value_pca_2 = second_component_pca.min()

define_function.adjust_and_save_plot(adata, folder_path, f'pca_{folder_name}_age,', 'pca', min_value_pca_1,
                                     max_value_pca_1,
                                     min_value_pca_2, max_value_pca_2, 'PC1', 'PC2', 'age')
define_function.adjust_and_save_plot(adata, folder_path, f'pca_{folder_name}_tissue,', 'pca', min_value_pca_1,
                                     max_value_pca_1,
                                     min_value_pca_2, max_value_pca_2, 'PC1', 'PC2', 'tissue')

# umap
first_component_umap = adata.obsm['X_umap'][:, 0]
second_component_umap = adata.obsm['X_umap'][:, 1]

max_value_umap_1 = first_component_umap.max()
min_value_umap_1 = first_component_umap.min()
max_value_umap_2 = second_component_umap.max()
min_value_umap_2 = second_component_umap.min()

define_function.adjust_and_save_plot(adata, folder_path, f'umap_{folder_name}_age,', 'umap', min_value_umap_1,
                                     max_value_umap_1, min_value_umap_2, max_value_umap_2, 'UMAP1', 'UMAP2', 'age')
define_function.adjust_and_save_plot(adata, folder_path, f'umap_{folder_name}_tissue,', 'umap', min_value_umap_1,
                                     max_value_umap_1, min_value_umap_2, max_value_umap_2, 'UMAP1', 'UMAP2', 'tissue')

adata = define_function.sort_adata_by_attribute(adata, 'tissue')  # sort by tissue
age_dict = define_function.split_adata_by_attribute(adata, 'age')  # split by age
time_points = list(age_dict.keys())  # get the time points
cnsecutive_time_points = []  # Consecutive time points
for i in range(len(time_points) - 1):
    cnsecutive_time_points.append([time_points[i], time_points[i + 1]])

# create the folder to save the images
folder_path_dict = {}  # store the folder path
for time_point in time_points:
    folder_path_dict[time_point] = f'{folder_path}/{time_point}'
    if not os.path.exists(folder_path_dict[time_point]):
        os.makedirs(folder_path_dict[time_point])

for cnsecutive_time_point in cnsecutive_time_points:
    folder_path_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[
        1]] = f'{folder_path}/{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}'
    if not os.path.exists(folder_path_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]]):
        os.makedirs(folder_path_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]])

adata_analyzed_dict = {}  # store the analyzed data
for cnsecutive_time_point in cnsecutive_time_points:
    integrate_adata = ad.concat(
        [age_dict[cnsecutive_time_point[0]], age_dict[cnsecutive_time_point[1]]])  # select consecutive time points
    adata_analyzed_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]] = integrate_adata

    # pca
    define_function.adjust_and_save_plot(integrate_adata,
                             folder_path_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]],
                             f'pca_{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}_age', 'pca', min_value_pca_1,
                             max_value_pca_1, min_value_pca_2, max_value_pca_2, 'PC1', 'PC2', 'age')
    define_function.adjust_and_save_plot(integrate_adata,
                             folder_path_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]],
                             f'pca_{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}_age', 'pca', min_value_pca_1,
                             max_value_pca_1, min_value_pca_2, max_value_pca_2, 'PC1', 'PC2', 'age')

    # umap
    define_function.adjust_and_save_plot(integrate_adata,
                              folder_path_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]],
                              f'umap_{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}_age', 'umap',
                              min_value_umap_1, max_value_umap_1, min_value_umap_2, max_value_umap_2, 'UMAP1', 'UMAP2',
                              'age')
    define_function.adjust_and_save_plot(integrate_adata,
                              folder_path_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]],
                              f'umap_{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}_age', 'umap',
                              min_value_umap_1, max_value_umap_1, min_value_umap_2, max_value_umap_2, 'UMAP1', 'UMAP2',
                              'tissue')