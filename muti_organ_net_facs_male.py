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
folder_name = 'facs_male'
data_folder_path = f'./data/{folder_name}'
image_folder_path = f'./image/{date_str}/{folder_name}'
if not os.path.exists(image_folder_path):
    os.makedirs(image_folder_path)

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
'''
print(tissue_age_dict['3m_18m'])
print(cell_ontology_class_tissue_age_dict['3m_18m']['Brain_Myeloid'])
adata1_young = adata_cell_ontology_class_tissue_age_dict['3m_18m']['Brain_Myeloid']['macrophage']['3m']
adata1_old = adata_cell_ontology_class_tissue_age_dict['3m_18m']['Brain_Myeloid']['macrophage']['18m']
adata2_young = adata_cell_ontology_class_tissue_age_dict['3m_18m']['Brain_Non-Myeloid']['endothelial cell']['3m']

adata1_yong_pca, adata1_old_pca, adata2_yong_pca, cumulative_explained_variance_ratio = define_function.process_pca(
    adata1_young, adata1_old, adata2_young, image_folder_path, 'Brain_Myeloid_macrophage_3m_18m')
print(cumulative_explained_variance_ratio)
lambda_, p1, p2 = define_function.wproj_adata(adata1_yong_pca, adata1_old_pca, adata2_yong_pca,
                                                                      image_folder_path,'Brain_Myeloid_macrophage_3m_18m')
print(lambda_, p1, p2)
'''
lambda_dict = {}
p1_dict = {}
p2_dict = {}
image_folder_path_pca_dict = {}

for cnsecutive_time_point_list in cnsecutive_time_points:
    cnsecutive_time_point = '_'.join(cnsecutive_time_point_list)
    time_point_yong = cnsecutive_time_point_list[0]
    time_point_old = cnsecutive_time_point_list[1]
    lambda_dict[cnsecutive_time_point] = []
    p1_dict[cnsecutive_time_point] = []
    p2_dict[cnsecutive_time_point] = []
    image_folder_path_pca_dict[cnsecutive_time_point] = {}
    for tissue1 in tissue_age_dict[cnsecutive_time_point]:
        image_folder_path_pca_dict[cnsecutive_time_point][tissue1] = {}
        for cell_ontology_class1 in cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue1]:
            adata1_yong = \
                adata_cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue1][cell_ontology_class1][
                    time_point_yong]
            adata1_old = \
                adata_cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue1][cell_ontology_class1][
                    time_point_old]
            lambda_list = []
            p1_list = []
            p2_list = []
            image_folder_path_pca_dict[cnsecutive_time_point][tissue1][cell_ontology_class1] = {}

            for tissue2 in tissue_age_dict[cnsecutive_time_point]:
                image_folder_path_pca_dict[cnsecutive_time_point][tissue1][cell_ontology_class1][tissue2] = {}
                for cell_ontology_class2 in cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue2]:
                    if tissue1 != tissue2 or cell_ontology_class1 != cell_ontology_class2:
                        image_folder_path_pca_dict[cnsecutive_time_point][tissue1][cell_ontology_class1][tissue2][
                            cell_ontology_class2] = f'{image_folder_path}/{cnsecutive_time_point}/{tissue1}/{cell_ontology_class1}/{tissue2}/{cell_ontology_class2}'
                        if not os.path.exists(
                                image_folder_path_pca_dict[cnsecutive_time_point][tissue1][cell_ontology_class1][
                                    tissue2][cell_ontology_class2]):
                            os.makedirs(
                                image_folder_path_pca_dict[cnsecutive_time_point][tissue1][cell_ontology_class1][
                                    tissue2][cell_ontology_class2])
                        adata2_yong = \
                            adata_cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue2][
                                cell_ontology_class2][
                                time_point_yong]
                        adata2_old = \
                            adata_cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue2][
                                cell_ontology_class2][
                                time_point_old]
                        adata1_yong_pca, adata1_old_pca, adata2_yong_pca, cumulative_explained_variance_ratio = define_function.process_pca(
                            adata1_yong, adata1_old, adata2_yong,
                            image_folder_path_pca_dict[cnsecutive_time_point][tissue1][cell_ontology_class1][tissue2][
                                cell_ontology_class2],
                            f'{tissue1}_{cell_ontology_class1}_{tissue2}_{cell_ontology_class2}')
                        print(cumulative_explained_variance_ratio)
                        lambda_, p1, p2 = define_function.wproj_adata(adata1_yong_pca, adata1_old_pca, adata2_yong_pca,
                                                                      image_folder_path_pca_dict[cnsecutive_time_point][
                                                                          tissue1][cell_ontology_class1][tissue2][
                                                                          cell_ontology_class2],
                                                                      f'{tissue1}_{cell_ontology_class1}_{tissue2}_{cell_ontology_class2}')
                        lambda_list.append(lambda_)
                        p1_list.append(p1)
                        p2_list.append(p2)
                        print(lambda_, p1, p2)

            lambda_dict[cnsecutive_time_point].append(lambda_list)
            p1_dict[cnsecutive_time_point].append(p1_list)
            p2_dict[cnsecutive_time_point].append(p2_list)

# save the results
cell_ontology_class_ticks_dict = {}
for cnsecutive_time_point_list in cnsecutive_time_points:
    cnsecutive_time_point = '_'.join(cnsecutive_time_point_list)
    cell_ontology_class_ticks_dict[cnsecutive_time_point] = []
    for tissue in tissue_age_dict[cnsecutive_time_point]:
        for cell_ontology_class in cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue]:
            cell_ontology_class_ticks_dict[cnsecutive_time_point].append('_'.join([tissue, cell_ontology_class]))

for cnsecutive_time_point_list in cnsecutive_time_points:
    cnsecutive_time_point = '_'.join(cnsecutive_time_point_list)
    time_point_yong = cnsecutive_time_point_list[0]
    time_point_old = cnsecutive_time_point_list[1]
    lambda_list = lambda_dict[cnsecutive_time_point]
    p1_list = p1_dict[cnsecutive_time_point]
    p2_list = p2_dict[cnsecutive_time_point]
    lambda_df = pd.DataFrame(lambda_list, index=cell_ontology_class_ticks_dict[cnsecutive_time_point], columns=cell_ontology_class_ticks_dict[cnsecutive_time_point])
    p1_df = pd.DataFrame(p1_list, index=cell_ontology_class_ticks_dict[cnsecutive_time_point], columns=cell_ontology_class_ticks_dict[cnsecutive_time_point])
    p2_df = pd.DataFrame(p2_list, index=cell_ontology_class_ticks_dict[cnsecutive_time_point], columns=cell_ontology_class_ticks_dict[cnsecutive_time_point])
    lambda_df.to_csv(f'{image_folder_path}/{cnsecutive_time_point}/lambda_{time_point_yong}_{time_point_old}.csv')
    p1_df.to_csv(f'{image_folder_path}/{cnsecutive_time_point}/p1_{time_point_yong}_{time_point_old}.csv')
    p2_df.to_csv(f'{image_folder_path}/{cnsecutive_time_point}/p2_{time_point_yong}_{time_point_old}.csv')




