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
organomix_path = f'{ms_path}/organomix'
celltype_number = 52
k_list = np.arange(celltype_number*celltype_number).tolist()

# load position_list
with open(f'{ms_path}/organomix/position_list.txt', mode='r') as f:
    position_list = f.read().splitlines()
    position_list = list(map(float, position_list))

# load label_name_list
with open(f'{ms_path}/organomix/label_name_list.txt', mode='r') as f:
    label_name_list = f.read().splitlines()

# load line_position_list
with open(f'{ms_path}/organomix/line_position_list.txt', mode='r') as f:
    line_position_list = f.read().splitlines()
    line_position_list = list(map(float, line_position_list))

# load anndata
integrate_adata_path = f'{ms_path}/anndata/integrate_adata_sample.h5ad'
integrate_adata = sc.read_h5ad(integrate_adata_path)

time_point_list = integrate_adata.obs['time.point'].unique().tolist()
ctime_point_list = []
for i in range(len(time_point_list) - 1):
    ctime_point_list.append(f'{time_point_list[i]}_{time_point_list[i + 1]}')

for ctime_point in ctime_point_list:
    time_point_y = ctime_point.split('_')[0]
    time_point_o = ctime_point.split('_')[1]

    for disease in disease_list:
        result_path = f'{ms_path}/organomix/{disease}_{time_point_y}_{time_point_o}_result.txt'
        with open(result_path, mode='r') as f:
            result = f.read().splitlines()
        result = result[1:]
        results = [i.split('\t') for i in result]

        lambda_matrix = [[0 for i in range(celltype_number)] for j in range(celltype_number)]
        p1_matrix = [[0 for i in range(celltype_number)] for j in range(celltype_number)]
        p2_matrix = [[0 for i in range(celltype_number)] for j in range(celltype_number)]

        for result in results:
            k = int(result[0])
            i = int(result[1])
            j = int(result[2])
            lambda_ = float(result[3])
            p1 = float(result[4])
            p2 = float(result[5])

            lambda_matrix[i][j] = lambda_
            p1_matrix[i][j] = p1
            p2_matrix[i][j] = p2

            if disease == 'AD' and ctime_point == '4.5_6.0' and k in k_list:
                print(k)
                k_list.remove(k)

        # plot lambda matrix
        lambda_matrix = np.array(lambda_matrix)

        fig, ax1 = pl.subplots(figsize=(10, 10))
        img = ax1.imshow(lambda_matrix, cmap='jet', vmin=0, vmax=1)

        for i in range(6):
            ax1.plot([-0.5, 51.5], [line_position_list[i + 1], line_position_list[i + 1]], color='black', linewidth=1)
            ax1.plot([line_position_list[i + 1], line_position_list[i + 1]], [-0.5, 51.5], color='black', linewidth=1)

        xticks = position_list
        yticks = position_list

        ax1.set_xticks(xticks)
        ax1.set_yticks(yticks)

        xlabels = label_name_list
        ylabels = label_name_list

        ax1.set_xticklabels(xlabels, rotation=-90)
        ax1.set_yticklabels(ylabels)

        ax1.xaxis.tick_top()
        ax1.tick_params(axis='both', length=0)

        ax1.set_title(f'Lambda_matrix_{disease}_{time_point_y}_{time_point_o}')
        cbar = pl.colorbar(img, ax=ax1)


        pl.show()
        fig.savefig(f'{ms_path}/organomix/{disease}_{time_point_y}_{time_point_o}_lambda_matrix.png')

        # plot p1 matrix
        p1_matrix = np.array(p1_matrix)

        fig, ax1 = pl.subplots(figsize=(10, 10))
        img = ax1.imshow(p1_matrix, cmap='jet', vmin=0, vmax=1)

        for i in range(6):
            ax1.plot([-0.5, 51.5], [line_position_list[i + 1], line_position_list[i + 1]], color='black', linewidth=1)
            ax1.plot([line_position_list[i + 1], line_position_list[i + 1]], [-0.5, 51.5], color='black', linewidth=1)

        xticks = position_list
        yticks = position_list

        ax1.set_xticks(xticks)
        ax1.set_yticks(yticks)

        xlabels = label_name_list
        ylabels = label_name_list

        ax1.set_xticklabels(xlabels, rotation=-90)
        ax1.set_yticklabels(ylabels)

        ax1.xaxis.tick_top()
        ax1.tick_params(axis='both', length=0)

        ax1.set_title(f'P1_matrix_{disease}_{time_point_y}_{time_point_o}')
        cbar = pl.colorbar(img, ax=ax1)

        pl.show()
        fig.savefig(f'{ms_path}/organomix/{disease}_{time_point_y}_{time_point_o}_P1_matrix.png')



        # plot p2 matrix
        p2_matrix = np.array(p2_matrix)

        fig, ax1 = pl.subplots(figsize=(10, 10))
        img = ax1.imshow(p2_matrix, cmap='jet', vmin=0, vmax=1)

        for i in range(6):
            ax1.plot([-0.5, 51.5], [line_position_list[i + 1], line_position_list[i + 1]], color='black', linewidth=1)
            ax1.plot([line_position_list[i + 1], line_position_list[i + 1]], [-0.5, 51.5], color='black', linewidth=1)

        xticks = position_list
        yticks = position_list

        ax1.set_xticks(xticks)
        ax1.set_yticks(yticks)

        xlabels = label_name_list
        ylabels = label_name_list

        ax1.set_xticklabels(xlabels, rotation=-90)
        ax1.set_yticklabels(ylabels)

        ax1.xaxis.tick_top()
        ax1.tick_params(axis='both', length=0)

        ax1.set_title(f'P2_matrix_{disease}_{time_point_y}_{time_point_o}')
        cbar = pl.colorbar(img, ax=ax1)

        pl.show()
        fig.savefig(f'{ms_path}/organomix/{disease}_{time_point_y}_{time_point_o}_P2_matrix.png')

        # save k_list
        with open(f'{ms_path}/organomix/k_list.txt', mode='w') as f:
            for k in k_list:
                f.write(f'{k}\n')

        # load k_list
        with open(f'{ms_path}/organomix/k_list.txt', mode='r') as f:
            k_list = f.read().splitlines()
        k_list = [int(i) for i in k_list]




print(len(k_list))
