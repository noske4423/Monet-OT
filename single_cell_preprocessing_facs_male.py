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

define_function.plot_mean_expression(adata, folder_path, folder_name)  # plot the mean expression of each gene
define_function.plot_highest_expr_genes(adata, folder_path, folder_name)  # plot the highest expression genes

# pca
first_component_pca = adata.obsm['X_pca'][:, 0]
second_component_pca = adata.obsm['X_pca'][:, 1]

max_value_pca_1 = first_component_pca.max()
min_value_pca_1 = first_component_pca.min()
max_value_pca_2 = second_component_pca.max()
min_value_pca_2 = second_component_pca.min()

define_function.adjust_and_save_plot(adata, folder_path, folder_name, 'pca', min_value_pca_1,
                                     max_value_pca_1,
                                     min_value_pca_2, max_value_pca_2, 'PC1', 'PC2', 'age')
define_function.adjust_and_save_plot(adata, folder_path, folder_name, 'pca', min_value_pca_1,
                                     max_value_pca_1,
                                     min_value_pca_2, max_value_pca_2, 'PC1', 'PC2', 'tissue')

# umap
first_component_umap = adata.obsm['X_umap'][:, 0]
second_component_umap = adata.obsm['X_umap'][:, 1]

max_value_umap_1 = first_component_umap.max()
min_value_umap_1 = first_component_umap.min()
max_value_umap_2 = second_component_umap.max()
min_value_umap_2 = second_component_umap.min()

define_function.adjust_and_save_plot(adata, folder_path, folder_name, 'umap', min_value_umap_1,
                                     max_value_umap_1, min_value_umap_2, max_value_umap_2, 'UMAP1', 'UMAP2', 'age')
define_function.adjust_and_save_plot(adata, folder_path, folder_name, 'umap', min_value_umap_1,
                                     max_value_umap_1, min_value_umap_2, max_value_umap_2, 'UMAP1', 'UMAP2', 'tissue')

adata = define_function.sort_adata_by_attribute(adata, 'tissue')  # sort by tissue
adata_age_dict = define_function.split_adata_by_attribute(adata, 'age')  # split by age
time_points = list(adata_age_dict.keys())  # get the time points

color_list = ['dodgerblue', 'orange', 'limegreen', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan', 'black']
color_age_dict = {}  # store the color for each age
for i in range(len(time_points)):
    color_age_dict[time_points[i]] = color_list[i]

# create the folder to save the images
folder_path_age_dict = {}  # store the folder path
for time_point in time_points:
    folder_path_age_dict[time_point] = f'{folder_path}/{time_point}'
    if not os.path.exists(folder_path_age_dict[time_point]):
        os.makedirs(folder_path_age_dict[time_point])

# plot the pca and umap  for each time point
for time_point in time_points:
    adata_time_point = adata_age_dict[time_point]

    # plot
    define_function.plot_mean_expression(adata_time_point, folder_path_age_dict[time_point],
                                         f'{folder_name}_{time_point}')
    define_function.plot_highest_expr_genes(adata_time_point, folder_path_age_dict[time_point],
                                            f'{folder_name}_{time_point}')

    define_function.adjust_and_save_plot(adata_time_point, folder_path_age_dict[time_point],
                                         f'{folder_name}_{time_point}', 'pca', min_value_pca_1,
                                         max_value_pca_1, min_value_pca_2, max_value_pca_2, 'PC1', 'PC2', 'tissue')
    define_function.adjust_and_save_plot(adata_time_point, folder_path_age_dict[time_point],
                                         f'{folder_name}_{time_point}', 'umap', min_value_umap_1,
                                         max_value_umap_1, min_value_umap_2, max_value_umap_2, 'UMAP1', 'UMAP',
                                         'tissue')

cnsecutive_time_points = []  # Consecutive time points
for i in range(len(time_points) - 1):
    cnsecutive_time_points.append([time_points[i], time_points[i + 1]])

for cnsecutive_time_point in cnsecutive_time_points:
    folder_path_age_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[
        1]] = f'{folder_path}/{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}'
    if not os.path.exists(folder_path_age_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]]):
        os.makedirs(folder_path_age_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]])

adata_analyzed_dict = {}  # store the analyzed data
for cnsecutive_time_point in cnsecutive_time_points:
    integrate_adata = ad.concat(
        [adata_age_dict[cnsecutive_time_point[0]],
         adata_age_dict[cnsecutive_time_point[1]]])  # select consecutive time points
    adata_analyzed_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]] = integrate_adata

    # plot
    define_function.plot_mean_expression(integrate_adata, folder_path_age_dict[
        cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]],
                                         f'{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}')
    define_function.plot_highest_expr_genes(integrate_adata, folder_path_age_dict[
        cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]],
                                            f'{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}')
    # pca
    define_function.adjust_and_save_plot(integrate_adata,
                                         folder_path_age_dict[
                                             cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]],
                                         f'{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}', 'pca',
                                         min_value_pca_1,
                                         max_value_pca_1, min_value_pca_2, max_value_pca_2, 'PC1', 'PC2', 'age')
    define_function.adjust_and_save_plot(integrate_adata,
                                         folder_path_age_dict[
                                             cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]],
                                         f'{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}', 'pca',
                                         min_value_pca_1,
                                         max_value_pca_1, min_value_pca_2, max_value_pca_2, 'PC1', 'PC2', 'tissue')

    # umap
    define_function.adjust_and_save_plot(integrate_adata,
                                         folder_path_age_dict[
                                             cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]],
                                         f'{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}', 'umap',
                                         min_value_umap_1, max_value_umap_1, min_value_umap_2, max_value_umap_2,
                                         'UMAP1', 'UMAP2',
                                         'age')
    define_function.adjust_and_save_plot(integrate_adata,
                                         folder_path_age_dict[
                                             cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]],
                                         f'{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}', 'umap',
                                         min_value_umap_1, max_value_umap_1, min_value_umap_2, max_value_umap_2,
                                         'UMAP1', 'UMAP2',
                                         'tissue')

# split the data by tissue
adata_tissue_age_dict = {}
for time_point in time_points:
    adata_tissue_age_dict[time_point] = define_function.split_adata_by_attribute(adata_age_dict[time_point], 'tissue')

adata_tissue_age_analyzed_dict = {}  # store the analyzed data
for cnsecutive_time_point in cnsecutive_time_points:
    adata_tissue_age_analyzed_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]] = {}
    for tissue in list(adata_tissue_age_dict[cnsecutive_time_point[0]].keys()):
        if tissue in list(adata_tissue_age_dict[cnsecutive_time_point[1]].keys()):
            adata_tissue_age_analyzed_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][
                tissue] = ad.concat(
                [adata_tissue_age_dict[cnsecutive_time_point[0]][tissue],
                 adata_tissue_age_dict[cnsecutive_time_point[1]][tissue]])

folder_path_tissue_age_dict = {}
# plot the pca and umap for each tissue
for time_point in time_points:
    folder_path_tissue_age_dict[time_point] = {}
    adata_tissue_age = adata_tissue_age_dict[time_point]
    for tissue in list(adata_tissue_age.keys()):
        folder_path_tissue_age_dict[time_point][tissue] = f'{folder_path_age_dict[time_point]}/{tissue}'
        if not os.path.exists(folder_path_tissue_age_dict[time_point][tissue]):
            os.makedirs(folder_path_tissue_age_dict[time_point][tissue])

        # plot
        define_function.plot_mean_expression(adata_tissue_age[tissue], folder_path_tissue_age_dict[time_point][tissue],
                                             f'{folder_name}_{time_point}_{tissue}')
        define_function.plot_highest_expr_genes(adata_tissue_age[tissue],
                                                folder_path_tissue_age_dict[time_point][tissue],
                                                f'{folder_name}_{time_point}_{tissue}')

        # pca
        define_function.adjust_and_save_plot(adata_tissue_age[tissue], folder_path_tissue_age_dict[time_point][tissue],
                                             f'{folder_name}_{time_point}_{tissue}', 'pca',
                                             min_value_pca_1,
                                             max_value_pca_1, min_value_pca_2, max_value_pca_2, 'PC1', 'PC2',
                                             'cell_ontology_class')
        define_function.adjust_and_save_plot(adata_tissue_age[tissue], folder_path_tissue_age_dict[time_point][tissue],
                                             f'{folder_name}_{time_point}_{tissue}', 'pca', min_value_pca_1,
                                             max_value_pca_1, min_value_pca_2, max_value_pca_2, 'PC1', 'PC2', None,
                                             color_age_dict[time_point])

        # umap
        define_function.adjust_and_save_plot(adata_tissue_age[tissue], folder_path_tissue_age_dict[time_point][tissue],
                                             f'{folder_name}_{time_point}_{tissue}', 'umap',
                                             min_value_umap_1, max_value_umap_1, min_value_umap_2, max_value_umap_2,
                                             'UMAP1', 'UMAP2', 'cell_ontology_class')
        define_function.adjust_and_save_plot(adata_tissue_age[tissue], folder_path_tissue_age_dict[time_point][tissue],
                                             f'{folder_name}_{time_point}_{tissue}', 'umap', min_value_umap_1,
                                             max_value_umap_1, min_value_umap_2, max_value_umap_2, 'UMAP1', 'UMAP2',
                                             None, color_age_dict[time_point])

for cnsecutive_time_point in cnsecutive_time_points:
    adata_tissue_age_analyzed = adata_tissue_age_analyzed_dict[
        cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]]
    folder_path_tissue_age_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]] = {}

    for tissue in list(adata_tissue_age_analyzed.keys()):
        folder_path_tissue_age_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][
            tissue] = f'{folder_path_age_dict[cnsecutive_time_point[0] + "_" + cnsecutive_time_point[1]]}/{tissue}'
        if not os.path.exists(
                folder_path_tissue_age_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue]):
            os.makedirs(folder_path_tissue_age_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue])

        # plot
        define_function.plot_mean_expression(adata_tissue_age_analyzed[tissue], folder_path_tissue_age_dict[
            cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue],
                                             f'{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}_{tissue}')
        define_function.plot_highest_expr_genes(adata_tissue_age_analyzed[tissue], folder_path_tissue_age_dict[
            cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue],
                                                f'{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}_{tissue}')

        # pca
        define_function.adjust_and_save_plot(adata_tissue_age_analyzed[tissue], folder_path_tissue_age_dict[
            cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue],
                                             f'{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}_{tissue}', 'pca',
                                             min_value_pca_1, max_value_pca_1, min_value_pca_2, max_value_pca_2, 'PC1',
                                             'PC2', 'cell_ontology_class')
        define_function.adjust_and_save_plot(adata_tissue_age_analyzed[tissue], folder_path_tissue_age_dict[
            cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue],
                                             f'{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}_{tissue}', 'pca',
                                             min_value_pca_1, max_value_pca_1, min_value_pca_2, max_value_pca_2, 'PC1',
                                             'PC2', 'age')

        # umap
        define_function.adjust_and_save_plot(adata_tissue_age_analyzed[tissue], folder_path_tissue_age_dict[
            cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue],
                                             f'{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}_{tissue}', 'umap',
                                             min_value_umap_1, max_value_umap_1, min_value_umap_2, max_value_umap_2,
                                             'UMAP1', 'UMAP2', 'cell_ontology_class')
        define_function.adjust_and_save_plot(adata_tissue_age_analyzed[tissue], folder_path_tissue_age_dict[
            cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue],
                                             f'{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}_{tissue}', 'umap',
                                             min_value_umap_1, max_value_umap_1, min_value_umap_2, max_value_umap_2,
                                             'UMAP1', 'UMAP2', 'age')

# split the data by cell ontology class
adata_cell_ontology_class_tissue_age_dict = {}  # adata for each cell ontology class

for time_point in time_points:
    adata_cell_ontology_class_tissue_age_dict[time_point] = {}
    for tissue in list(adata_tissue_age_dict[time_point].keys()):
        adata_cell_ontology_class_tissue_age_dict[time_point][tissue] = define_function.split_adata_by_attribute(
            adata_tissue_age_dict[time_point][tissue], 'cell_ontology_class')

adata_cell_ontology_class_tissue_age_analyzed_dict = {}  # adata for each cell ontology class
for cnsecutive_time_point in cnsecutive_time_points:
    adata_cell_ontology_class_tissue_age_analyzed_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]] = {}
    for tissue in list(adata_cell_ontology_class_tissue_age_dict[cnsecutive_time_point[0]].keys()):
        adata_cell_ontology_class_tissue_age_analyzed_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][
            tissue] = {}
        for cell_ontology_class in list(
                adata_cell_ontology_class_tissue_age_dict[cnsecutive_time_point[0]][tissue].keys()):
            if cell_ontology_class in list(
                    adata_cell_ontology_class_tissue_age_dict[cnsecutive_time_point[1]][tissue].keys()):
                if \
                        adata_cell_ontology_class_tissue_age_dict[cnsecutive_time_point[0]][tissue][
                            cell_ontology_class].shape[
                            0] > 20 and adata_cell_ontology_class_tissue_age_dict[cnsecutive_time_point[1]][tissue][
                            cell_ontology_class].shape[0] > 20:
                    adata_cell_ontology_class_tissue_age_analyzed_dict[
                        cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue][
                        cell_ontology_class] = ad.concat(
                        [adata_cell_ontology_class_tissue_age_dict[cnsecutive_time_point[0]][tissue][
                             cell_ontology_class],
                         adata_cell_ontology_class_tissue_age_dict[cnsecutive_time_point[1]][tissue][
                             cell_ontology_class]])

folder_path_cell_ontology_class_tissue_age_dict = {}  # folder path for each cell ontology class
for time_point in time_points:
    folder_path_cell_ontology_class_tissue_age_dict[time_point] = {}
    for tissue in list(adata_tissue_age_dict[time_point].keys()):
        folder_path_cell_ontology_class_tissue_age_dict[time_point][tissue] = {}
        for cell_ontology_class in list(
                adata_cell_ontology_class_tissue_age_dict[time_point][tissue].keys()):
            folder_path_cell_ontology_class_tissue_age_dict[time_point][tissue][
                cell_ontology_class] = f'{folder_path_tissue_age_dict[time_point][tissue]}/{cell_ontology_class}'
            if not os.path.exists(folder_path_cell_ontology_class_tissue_age_dict[time_point][tissue][
                                      cell_ontology_class]):
                os.makedirs(folder_path_cell_ontology_class_tissue_age_dict[time_point][tissue][
                                cell_ontology_class])
            # plot
            define_function.plot_mean_expression(adata_cell_ontology_class_tissue_age_dict[time_point][tissue][
                                                     cell_ontology_class],
                                                 folder_path_cell_ontology_class_tissue_age_dict[time_point][tissue][
                                                     cell_ontology_class],
                                                 f'{time_point}_{tissue}_{cell_ontology_class}')
            define_function.plot_highest_expr_genes(adata_cell_ontology_class_tissue_age_dict[time_point][tissue][
                                                        cell_ontology_class],
                                                    folder_path_cell_ontology_class_tissue_age_dict[
                                                        time_point][tissue][cell_ontology_class],
                                                    f'{time_point}_{tissue}_{cell_ontology_class}')

            # pca
            define_function.adjust_and_save_plot(adata_cell_ontology_class_tissue_age_dict[time_point][tissue][
                                                     cell_ontology_class],
                                                 folder_path_cell_ontology_class_tissue_age_dict[time_point][tissue][
                                                     cell_ontology_class],
                                                 f'{time_point}_{tissue}_{cell_ontology_class}', 'pca', min_value_pca_1,
                                                 max_value_pca_1, min_value_pca_2, max_value_pca_2, 'PC1', 'PC2', None,
                                                 color_age_dict[time_point])

            # umap
            define_function.adjust_and_save_plot(adata_cell_ontology_class_tissue_age_dict[time_point][tissue][
                                                     cell_ontology_class],
                                                 folder_path_cell_ontology_class_tissue_age_dict[time_point][tissue][
                                                     cell_ontology_class],
                                                 f'{time_point}_{tissue}_{cell_ontology_class}', 'umap',
                                                 min_value_umap_1,
                                                 max_value_umap_1, min_value_umap_2, max_value_umap_2, 'UMAP1',
                                                 'UMAP2', None, color_age_dict[time_point])

for cnsecutive_time_point in cnsecutive_time_points:
    folder_path_cell_ontology_class_tissue_age_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]] = {}
    for tissue in list(adata_cell_ontology_class_tissue_age_analyzed_dict[
                           cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]].keys()):
        folder_path_cell_ontology_class_tissue_age_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][
            tissue] = {}
        for cell_ontology_class in list(
                adata_cell_ontology_class_tissue_age_analyzed_dict[
                    cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue].keys()):
            folder_path_cell_ontology_class_tissue_age_dict[cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][
                tissue][
                cell_ontology_class] = f'{folder_path_tissue_age_dict[cnsecutive_time_point[0] + "_" + cnsecutive_time_point[1]][tissue]}/{cell_ontology_class}'
            if not os.path.exists(folder_path_cell_ontology_class_tissue_age_dict[
                                      cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue][
                                      cell_ontology_class]):
                os.makedirs(folder_path_cell_ontology_class_tissue_age_dict[
                                cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue][
                                cell_ontology_class])

            # plot
            define_function.plot_mean_expression(adata_cell_ontology_class_tissue_age_analyzed_dict[
                                                     cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue][
                                                     cell_ontology_class],
                                                 folder_path_cell_ontology_class_tissue_age_dict[
                                                     cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue][
                                                     cell_ontology_class],
                                                 f'{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}_{tissue}_{cell_ontology_class}')
            define_function.plot_highest_expr_genes(adata_cell_ontology_class_tissue_age_analyzed_dict[
                                                        cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][
                                                        tissue][
                                                        cell_ontology_class],
                                                    folder_path_cell_ontology_class_tissue_age_dict[
                                                        cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][
                                                        tissue][cell_ontology_class],
                                                    f'{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}_{tissue}_{cell_ontology_class}')

            # pca
            define_function.adjust_and_save_plot(adata_cell_ontology_class_tissue_age_analyzed_dict[
                                                     cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue][
                                                     cell_ontology_class],
                                                 folder_path_cell_ontology_class_tissue_age_dict[
                                                     cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][
                                                     tissue][cell_ontology_class],
                                                 f'{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}_{tissue}_{cell_ontology_class}',
                                                 'pca', min_value_pca_1,
                                                 max_value_pca_1, min_value_pca_2, max_value_pca_2, 'PC1', 'PC2', 'age')

            # umap
            define_function.adjust_and_save_plot(adata_cell_ontology_class_tissue_age_analyzed_dict[
                                                     cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue][
                                                     cell_ontology_class],
                                                 folder_path_cell_ontology_class_tissue_age_dict[
                                                     cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][
                                                     tissue][cell_ontology_class],
                                                 f'{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}_{tissue}_{cell_ontology_class}',
                                                 'umap', min_value_umap_1,
                                                 max_value_umap_1, min_value_umap_2, max_value_umap_2, 'UMAP1', 'UMAP2',
                                                 'age')

# adata_cell_ontology_class_tissue_age_analyzed_dict: {cnsecutive_time_point[0]_cnsecutive_time_point[1]: {tissue: {cell_ontology_class: adata}}}

# save adata_cell_ontology_class_tissue_age_analyzed_dict
data_folder_path_cell_ontology_class_tissue_age_analyzed_dict = {}  # {cnsecutive_time_point[0]_cnsecutive_time_point[1]: {tissue: {cell_ontology_class: str}}}
for cnsecutive_time_point in cnsecutive_time_points:
    data_folder_path_cell_ontology_class_tissue_age_analyzed_dict[
        cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]] = {}
    for tissue in list(adata_cell_ontology_class_tissue_age_analyzed_dict[
                           cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]].keys()):
        data_folder_path_cell_ontology_class_tissue_age_analyzed_dict[
            cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue] = {}
        for cell_ontology_class in list(
                adata_cell_ontology_class_tissue_age_analyzed_dict[
                    cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue].keys()):
            data_folder_path_cell_ontology_class_tissue_age_analyzed_dict[
                cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue][
                cell_ontology_class] = f'data/{folder_name}/{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}/{tissue}/{cell_ontology_class}'
            if not os.path.exists(data_folder_path_cell_ontology_class_tissue_age_analyzed_dict[
                                      cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue][
                                      cell_ontology_class]):
                os.makedirs(data_folder_path_cell_ontology_class_tissue_age_analyzed_dict[
                                cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue][
                                cell_ontology_class])
            adata_cell_ontology_class_tissue_age_analyzed_dict[
                cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][
                tissue][cell_ontology_class].write(data_folder_path_cell_ontology_class_tissue_age_analyzed_dict[
                                                       cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][
                                                       tissue][
                                                       cell_ontology_class] + f'/{folder_name}_{cnsecutive_time_point[0]}_{cnsecutive_time_point[1]}_{tissue}_{cell_ontology_class}.h5ad')

cell_ontology_class_tissue_age_analyzed_dict = {}
for cnsecutive_time_point in cnsecutive_time_points:
    cell_ontology_class_tissue_age_analyzed_dict[
        cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]] = {}
    for tissue in list(adata_cell_ontology_class_tissue_age_analyzed_dict[
                           cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]].keys()):
        cell_ontology_class_tissue_age_analyzed_dict[
            cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue] = list(
            adata_cell_ontology_class_tissue_age_analyzed_dict[
                cnsecutive_time_point[0] + '_' + cnsecutive_time_point[1]][tissue].keys())

# save cell_ontology_class_tissue_age_analyzed_dict as csv
cell_ontology_class_tissue_age_analyzed_dict_df = pd.DataFrame.from_dict(cell_ontology_class_tissue_age_analyzed_dict, orient='index')
cell_ontology_class_tissue_age_analyzed_dict_df.to_csv(f'data/{folder_name}/{}_cell_ontology_class_tissue_age_analyzed.csv')

