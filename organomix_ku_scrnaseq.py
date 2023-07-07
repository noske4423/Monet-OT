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
import define_function as d_f
import csv
from datetime import datetime
import os
import gc
import multiprocessing as mp
import time
import psutil

organ_list = ['Blood', 'BM', 'Brain', 'Liver', 'LN', 'SmallInt', 'Spleen']
disease_list = ['WT', 'AD']
color_list = ['dodgerblue', 'orange', 'limegreen', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan', 'black']

celltype_ref = {'Blood': 'MCA_broad_cell_type',
                'BM': 'Tabula.Muris_cell_type',
                'Brain': 'Azimuth_predicted.subclass',
                'Liver': 'Tabula.Muris_cell_type',
                'LN': 'MCA_Blood_broad_cell_type',
                'SmallInt': 'MCA_broad_cell_type',
                'Spleen': 'Tabula.Muris_cell_type'
                }
ms_path = '/Users/kataokayuunosuke/MS'
image_path = f'{ms_path}/image'

# load anndata
integrate_adata_path = f'{ms_path}/anndata/integrate_adata_sample.h5ad'
integrate_adata = sc.read_h5ad(integrate_adata_path)

time_point_list = integrate_adata.obs['time.point'].unique().tolist()
ctime_point_list = []
for i in range(len(time_point_list) - 1):
    ctime_point_list.append(f'{time_point_list[i]}_{time_point_list[i + 1]}')

organ_celltype_list = []
k = 0
i = 0
for organ1 in organ_list:
    # load celltype list
    celltype_list_path = f'{ms_path}/{organ1}_cell_type_sample_list.txt'
    with open(celltype_list_path, mode='r') as f:
        celltype_list1 = f.read().splitlines()

    for celltype1 in celltype_list1:
        j = 0

        for organ2 in organ_list:
            # load celltype list
            celltype_list_path = f'{ms_path}/{organ2}_cell_type_sample_list.txt'
            with open(celltype_list_path, mode='r') as f:
                celltype_list2 = f.read().splitlines()
            for celltype2 in celltype_list2:
                organ_celltype_list.append([k, i, j, organ1, celltype1, organ2, celltype2])

                k += 1
                j += 1
        i += 1


if __name__ == '__main__':
    process_list = []
    process_num = 50
    l = 600
    for organ_celltype in organ_celltype_list[l:]:
        k = organ_celltype[0]
        i = organ_celltype[1]
        j = organ_celltype[2]
        organ1 = organ_celltype[3]
        celltype1 = organ_celltype[4]
        organ2 = organ_celltype[5]
        celltype2 = organ_celltype[6]

        p = mp.Process(target=d_f.worker_ku_scrnaseq, args=[[k, i, j, organ1, celltype1, organ2, celltype2]])
        process_list.append(p)
        p.start()
        print(k, i, j)

        if len(process_list) == process_num or k == len(organ_celltype_list) - 1:
            for p in process_list:
                print(p)
                p.join()  # wait for all process to finish

            process_list = []  # reset process list
