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
import gc

# Set up the folder
now = datetime.now()
date_str = now.strftime("%Y-%m-%d")
folder_name = 'facs_male'
folder_name_new = 'artificial_data_' + folder_name
data_folder_path = f'./data/{folder_name}'
data_folder_path_new = f'./data/{date_str}/{folder_name_new}'
image_folder_path_new = f'./image/{date_str}/{folder_name_new}'
if not os.path.exists(data_folder_path_new):
    os.makedirs(data_folder_path_new)
if not os.path.exists(image_folder_path_new):
    os.makedirs(image_folder_path_new)

all_files_and_folders = os.listdir(data_folder_path)
print(all_files_and_folders)
cnsecutive_time_points = [name for name in all_files_and_folders if os.path.isdir(os.path.join(data_folder_path, name))]
# split by '_'
cnsecutive_time_points = [name.split('_') for name in cnsecutive_time_points]

print(cnsecutive_time_points)

tissue_age_dict = {}
cell_ontology_class_tissue_age_dict = {}
adata_cell_ontology_class_tissue_age_dict = {}

for cnsecutive_time_point_list in cnsecutive_time_points:
    cnsecutive_time_point = '_'.join(cnsecutive_time_point_list)
    cell_ontology_class_tissue_age_dict[cnsecutive_time_point] = {}
    adata_cell_ontology_class_tissue_age_dict[cnsecutive_time_point] = {}
    tissue_folder_path = f'{data_folder_path}/{cnsecutive_time_point}'
    tissue_files_and_folders = os.listdir(tissue_folder_path)
    tissue_age_dict[cnsecutive_time_point] = [name for name in tissue_files_and_folders if
                                              os.path.isdir(os.path.join(tissue_folder_path, name))]
    for tissue in tissue_age_dict[cnsecutive_time_point]:
        adata_cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue] = {}
        cell_ontology_class_folder_path = f'{tissue_folder_path}/{tissue}'
        cell_ontology_class_files_and_folders = os.listdir(cell_ontology_class_folder_path)
        cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue] = [name for name in
                                                                              cell_ontology_class_files_and_folders if
                                                                              os.path.isdir(os.path.join(
                                                                                  cell_ontology_class_folder_path,
                                                                                  name))]
        for cell_ontology_class in cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue]:
            adata = sc.read_h5ad(
                f'{cell_ontology_class_folder_path}/{cell_ontology_class}/{folder_name}_{cnsecutive_time_point}_{tissue}_{cell_ontology_class}.h5ad')
            adata_cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue][
                cell_ontology_class] = define_function.split_adata_by_attribute(adata, 'age')

# tissue_age_dict :{consecutive_time_point: [tissue1, tissue2, ...], ...}
# cell_ontology_class_tissue_age_dict :{consecutive_time_point: {tissue1: [cell_ontology_class1, cell_ontology_class2, ...], ...}, ...}
# adata_cell_ontology_class_tissue_age_dict :{consecutive_time_point: {tissue1: {cell_ontology_class1: {time_point_yong: adata, time_point_old: adata}, ...}, ...}, ...}

# sort the tissue_age_dict by alphabetical order
for cnsecutive_time_point_list in cnsecutive_time_points:
    cnsecutive_time_point = '_'.join(cnsecutive_time_point_list)
    tissue_age_dict[cnsecutive_time_point].sort()

# sort the cell_ontology_class_tissue_age_dict by alphabetical order
for cnsecutive_time_point_list in cnsecutive_time_points:
    cnsecutive_time_point = '_'.join(cnsecutive_time_point_list)
    # list(cell_ontology_class_tissue_age_dict[cnsecutive_time_point].keys()).sort()
    for tissue in tissue_age_dict[cnsecutive_time_point]:
        cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue].sort()

tissue1_analyzed_list = ['Large_Intestine', 'Limb_Muscle']
tissue2_analyzed_list = ['Brain_Non-Myeloid', 'Heart']

parm_list1_old_dict_list = []
parm_list2_young_dict_list = []
num = 10
mu_list = [i / 5 for i in range(num + 1)]
sigma_list = [i / 5 for i in range(num + 1)]
ratio_list = [i / 10 for i in range(num + 1)]
row_ratio = 0.5
col_ratio = 0.1
mu = 1.0
sigma = 0.5
lambda_dict_list = []
p1_dict_list = []
p2_dict_list = []

for i in range(num):
    print(i)
    parm_list1_old_dict = {}
    parm_list2_young_dict = {}
    parm_list1_old_dict[tissue1_analyzed_list[0]] = [row_ratio, col_ratio, mu_list[i], sigma]
    parm_list1_old_dict[tissue1_analyzed_list[1]] = [row_ratio, col_ratio, 0, sigma]
    parm_list2_young_dict[tissue2_analyzed_list[0]] = [row_ratio, col_ratio, mu_list[i], sigma]
    parm_list2_young_dict[tissue2_analyzed_list[1]] = [row_ratio, col_ratio, 0, sigma]
    parm_list1_old_dict_list.append(parm_list1_old_dict)
    parm_list2_young_dict_list.append(parm_list2_young_dict)

for i in range(num):
    print(i)
    lambda_dict = {}
    p1_dict = {}
    p2_dict = {}
    parm_list1_old_dict = parm_list1_old_dict_list[i]
    parm_list2_young_dict = parm_list2_young_dict_list[i]
    parm_list1_old_dict[tissue1_analyzed_list[0]] = [row_ratio, col_ratio, mu_list[i], sigma]
    parm_list1_old_dict[tissue1_analyzed_list[1]] = [row_ratio, col_ratio, 0, sigma]
    parm_list2_young_dict[tissue2_analyzed_list[0]] = [row_ratio, col_ratio, mu_list[i], sigma]
    parm_list2_young_dict[tissue2_analyzed_list[1]] = [row_ratio, col_ratio, 0, sigma]
    parm_list1_old_dict_list.append(parm_list1_old_dict)
    parm_list2_young_dict_list.append(parm_list2_young_dict)

    for cnsecutive_time_point_list in cnsecutive_time_points:
        cnsecutive_time_point = '_'.join(cnsecutive_time_point_list)
        time_point_young = cnsecutive_time_point_list[0]
        time_point_old = cnsecutive_time_point_list[1]
        lambda_dict[cnsecutive_time_point] = []
        p1_dict[cnsecutive_time_point] = []
        p2_dict[cnsecutive_time_point] = []

        for tissue1 in tissue1_analyzed_list:

            for cell_ontology_class1 in cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue1]:
                adata1_young = adata_cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue1][
                    cell_ontology_class1][time_point_young]
                adata1_old = adata_cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue1][
                    cell_ontology_class1][time_point_old]

                for tissue2 in tissue2_analyzed_list:
                    for cell_ontology_class2 in cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue2]:
                        adata1_young_copy = adata1_young.copy()
                        adata1_old_copy = adata1_old.copy()
                        adata2_young = adata_cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue2][
                            cell_ontology_class2][time_point_young]
                        adata2_old = adata_cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue2][
                            cell_ontology_class2][time_point_old]

                        adata1_young_copy.obs['group'] = f'{tissue1}_{cell_ontology_class1}_yong'
                        adata1_old_copy.obs['group'] = f'{tissue1}_{cell_ontology_class1}_old'
                        adata2_young.obs['group'] = f'{tissue2}_{cell_ontology_class2}_yong'
                        adata2_old.obs['group'] = f'{tissue2}_{cell_ontology_class2}_old'

                        integrated_adata = ad.concat([adata1_young_copy, adata2_young, adata1_old_copy, adata2_old])
                        sc.pp.highly_variable_genes(integrated_adata)
                        integrated_adata = integrated_adata[:, integrated_adata.var['highly_variable']]

                        adata1_young_copy = integrated_adata[integrated_adata.obs['group'] == f'{tissue1}_{cell_ontology_class1}_yong']
                        adata1_old_copy = integrated_adata[integrated_adata.obs['group'] == f'{tissue1}_{cell_ontology_class1}_old']
                        adata2_young = integrated_adata[integrated_adata.obs['group'] == f'{tissue2}_{cell_ontology_class2}_yong']
                        adata2_old = integrated_adata[integrated_adata.obs['group'] == f'{tissue2}_{cell_ontology_class2}_old']

                        artificial_adata1_old, noisy_col_indices = define_function.generate_artificial_adata(adata1_old_copy,
                                                                                                             parm_list1_old_dict[
                                                                                                                 tissue1])
                        artificial_adata2_young, noisy_col_indices = define_function.generate_artificial_adata(
                            adata2_young,
                            parm_list2_young_dict[
                                tissue2],
                            noisy_col_indices)
                        image_folder_tissue_path = f'{image_folder_path_new}/{cnsecutive_time_point}/{mu_list[i]}/{tissue1}/{cell_ontology_class1}/{tissue2}/{cell_ontology_class2}'

                        if not os.path.exists(image_folder_tissue_path):
                            os.makedirs(image_folder_tissue_path)

                        title = f'{folder_name_new}_{cnsecutive_time_point}_{tissue1}_{tissue2}_{cell_ontology_class1}_{cell_ontology_class2}'
                        adata1_young_pca, adata1_old_pca, adata2_young_pca, cumulative_explained_variance_ratio = define_function.process_pca(
                            adata1_young_copy, artificial_adata1_old, artificial_adata2_young, adata2_old,
                            image_folder_tissue_path, title, highly_variable_genes=False)
                        lambda_, p1, p2 = define_function.wproj_adata(adata1_young_pca, adata1_old_pca,
                                                                      adata2_young_pca,
                                                                      image_folder_tissue_path, title)
                        print(f'lambda: {lambda_}, p1: {p1}, p2: {p2}')
