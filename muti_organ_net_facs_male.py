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

for cnsecutive_time_point in cnsecutive_time_points:
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
            adata_cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue][cell_ontology_class] = sc.read_h5ad(
                f'{cell_ontology_class_folder_path}/{cell_ontology_class}/{folder_name}_{cnsecutive_time_point}_{tissue}_{cell_ontology_class}.h5ad')

# tissue_age_dict :{consecutive_time_point: [tissue1, tissue2, ...], ...}
# cell_ontology_class_tissue_age_dict :{consecutive_time_point: {tissue1: [cell_ontology_class1, cell_ontology_class2, ...], ...}, ...}
# adata_cell_ontology_class_tissue_age_dict :{consecutive_time_point: {tissue1: {cell_ontology_class1: adata, cell_ontology_class2: adata, ...}, ...}, ...}

lamnda_dict = {}
p1_dict = {}
p2_dict = {}

for cnsecutive_time_point in cnsecutive_time_points:
    lamnda_dict[cnsecutive_time_point] = {}
    p1_dict[cnsecutive_time_point] = {}
    p2_dict[cnsecutive_time_point] = {}
    for tissue in tissue_age_dict[cnsecutive_time_point]:
        lamnda_dict[cnsecutive_time_point][tissue] = {}
        p1_dict[cnsecutive_time_point][tissue] = {}
        p2_dict[cnsecutive_time_point][tissue] = {}
        for cell_ontology_class in cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue]:
            adata_cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue][cell_ontology_class]