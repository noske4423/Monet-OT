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
import multiprocessing as mp
import time
import psutil

# Set up the folder to save the images
now = datetime.now()
date_str = now.strftime("%Y-%m-%d")
folder_name = 'facs_male'
data_folder_path = f'./data/{folder_name}'
image_folder_path = f'./image/{date_str}/{folder_name}'
if not os.path.exists(image_folder_path):
    os.makedirs(image_folder_path)

all_files_and_folders = os.listdir(data_folder_path)
# print(all_files_and_folders)
cnsecutive_time_points = [name for name in all_files_and_folders if os.path.isdir(os.path.join(data_folder_path, name))]
# split by '_'
cnsecutive_time_points = [name.split('_') for name in cnsecutive_time_points]

# print(cnsecutive_time_points)

tissue_age_dict = {}
cell_ontology_class_tissue_age_dict = {}
adata_cell_ontology_class_tissue_age_dict = {}
mem = psutil.virtual_memory().free / 1e9
print(f'free memory: {mem} GB')

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
        cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue] \
            = [name for name in cell_ontology_class_files_and_folders
               if os.path.isdir(os.path.join(cell_ontology_class_folder_path, name))]

mem = psutil.virtual_memory().free / 1e9
print(f'free memory: {mem} GB')

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

cell_ontology_class_ticks_dict = {}
for cnsecutive_time_point_list in cnsecutive_time_points:
    cnsecutive_time_point = '_'.join(cnsecutive_time_point_list)
    cell_ontology_class_ticks_dict[cnsecutive_time_point] = []
    for tissue in tissue_age_dict[cnsecutive_time_point]:
        for cell_ontology_class in cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue]:
            cell_ontology_class_ticks_dict[cnsecutive_time_point].append('_'.join([tissue, cell_ontology_class]))
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

mem = psutil.virtual_memory().free / 1e9
print(f'free memory: {mem} GB')

task_list = []
for cnsecutive_time_point_list in cnsecutive_time_points:
    cnsecutive_time_point = '_'.join(cnsecutive_time_point_list)
    lambda_dict[cnsecutive_time_point] = [
        [0 for i in range(len(cell_ontology_class_ticks_dict[cnsecutive_time_point]))] for j in
        range(len(cell_ontology_class_ticks_dict[cnsecutive_time_point]))]
    p1_dict[cnsecutive_time_point] = [[0 for i in range(len(cell_ontology_class_ticks_dict[cnsecutive_time_point]))]
                                      for j in range(len(cell_ontology_class_ticks_dict[cnsecutive_time_point]))]
    p2_dict[cnsecutive_time_point] = [[0 for i in range(len(cell_ontology_class_ticks_dict[cnsecutive_time_point]))]
                                      for j in range(len(cell_ontology_class_ticks_dict[cnsecutive_time_point]))]
    i = 0
    for tissue1 in tissue_age_dict[cnsecutive_time_point]:
        for cell_ontology_class1 in cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue1]:
            j = 0
            for tissue2 in tissue_age_dict[cnsecutive_time_point]:
                for cell_ontology_class2 in cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue2]:
                    task_list.append(
                        [i, j, cnsecutive_time_point, tissue1, cell_ontology_class1, tissue2, cell_ontology_class2])
                    j += 1

            i += 1

k = 4150
if __name__ == '__main__':
    process_list = []
    process_num = 100
    for l, task in enumerate(task_list[k:]):
        index = l + k
        i = int(task[0])
        j = int(task[1])
        cnsecutive_time_point = task[2]
        tissue1 = task[3]
        cell_ontology_class1 = task[4]
        tissue2 = task[5]
        cell_ontology_class2 = task[6]
        time_point_yong = cnsecutive_time_point_list[0]
        time_point_old = cnsecutive_time_point_list[1]

        if i != j:
            p = mp.Process(target=define_function.worker,
                           args=[
                               [index, cnsecutive_time_point, tissue1, cell_ontology_class1, tissue2,
                                cell_ontology_class2,
                                folder_name, image_folder_path, data_folder_path]])
            process_list.append(p)
            p.start()
            print(i, j)

            if len(process_list) == process_num or index == len(task_list) - 1:
                for p in process_list:
                    print(p)
                    p.join()  # wait for all process to finish

                process_list = []
                # print(lambda_dict[cnsecutive_time_point])
                mem = psutil.virtual_memory().free / 1e9
                print(i, j, f'free memory: {mem} GB')

'''
if __name__ == '__main__':
    print('start')
    for cnsecutive_time_point_list in cnsecutive_time_points:
        cnsecutive_time_point = '_'.join(cnsecutive_time_point_list)
        time_point_yong = cnsecutive_time_point_list[0]
        time_point_old = cnsecutive_time_point_list[1]
        lambda_dict[cnsecutive_time_point] = []
        p1_dict[cnsecutive_time_point] = []
        p2_dict[cnsecutive_time_point] = []
        for tissue1 in tissue_age_dict[cnsecutive_time_point]:
            for cell_ontology_class1 in cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue1]:
                adata = sc.read_h5ad(
                    f'{data_folder_path}/{cnsecutive_time_point}/{tissue1}/{cell_ontology_class1}/{folder_name}_{cnsecutive_time_point}_{tissue1}_{cell_ontology_class1}.h5ad')
                adata_dict = define_function.split_adata_by_attribute(adata, 'age')
                adata1_young = adata_dict[time_point_yong]
                adata1_old = adata_dict[time_point_old]
                lambda_list = []
                p1_list = []
                p2_list = []

                for tissue2 in tissue_age_dict[cnsecutive_time_point]:
                    for cell_ontology_class2 in cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue2]:
                        if tissue1 != tissue2 or cell_ontology_class1 != cell_ontology_class2:
                            image_folder_path_pca = f'{image_folder_path}/{cnsecutive_time_point}/{tissue1}/{cell_ontology_class1}/{tissue2}/{cell_ontology_class2}'
                            title = f'{tissue1}_{cell_ontology_class1}_{tissue2}_{cell_ontology_class2}_{time_point_yong}_{time_point_old}'
                            if not os.path.exists(image_folder_path_pca):
                                os.makedirs(image_folder_path_pca)
                            mem = psutil.virtual_memory().free / 1e9
                            print(f'free memory: {mem} GB')
                            adata = sc.read_h5ad(
                                f'{data_folder_path}/{cnsecutive_time_point}/{tissue2}/{cell_ontology_class2}/{folder_name}_{cnsecutive_time_point}_{tissue2}_{cell_ontology_class2}.h5ad')
                            adata_dict = define_function.split_adata_by_attribute(adata, 'age')
                            adata2_young = adata_dict[time_point_yong]
                            adata2_old = adata_dict[time_point_old]
                            adata_integrated = define_function.integrate_adata(adata1_young, adata1_old, adata2_young,
                                                                               adata2_old)

                            del adata2_young, adata2_old

                            print(tissue1, cell_ontology_class1, tissue2, cell_ontology_class2)
                            mem = psutil.virtual_memory().free / 1e9
                            print(f'free memory: {mem} GB')

                            p = Pool(5)
                            print('process start')
                            adata_integrated_pca = p.apply(define_function.plot_pca,
                                                           [[adata_integrated, image_folder_path_pca, title]])
                            cumulative_explained_variance_ratio = adata_integrated_pca.uns['pca'][
                                'variance_ratio'].sum()
                            print(3)
                            adata1_young_pca = adata_integrated_pca[
                                adata_integrated_pca.obs['group'] == f'{tissue1}_{cell_ontology_class1}_yong']
                            adata1_old_pca = adata_integrated_pca[
                                adata_integrated_pca.obs['group'] == f'{tissue1}_{cell_ontology_class1}_old']
                            adata2_young_pca = adata_integrated_pca[
                                adata_integrated_pca.obs['group'] == f'{tissue2}_{cell_ontology_class2}_yong']
                            print(cumulative_explained_variance_ratio)
                            print(4)
                            lambda_, p1, p2 = p.apply(define_function.wproj_adata, [
                                [adata1_young_pca, adata1_old_pca, adata2_young_pca, image_folder_path_pca, title]])
                            del adata_integrated_pca

                            p.close()  # no more tasks
                            p.join()  # block until all tasks are done

                            mem = psutil.virtual_memory().free / 1e9
                            print(f'free memory: {mem} GB')

                            lambda_list.append(lambda_)
                            p1_list.append(p1)
                            p2_list.append(p2)
                            print(lambda_, p1, p2)

                            print('process end')
                            print(len(lambda_list))

                            del adata_integrated, adata, adata_dict
                            del lambda_, p1, p2
                            del p, cumulative_explained_variance_ratio, image_folder_path_pca, title
                            gc.collect()

                            mem = psutil.virtual_memory().free / 1e9
                            print(f'free memory: {mem} GB')

                        else:
                            lambda_list.append(0)
                            p1_list.append(0)
                            p2_list.append(0)

                lambda_dict[cnsecutive_time_point].append(lambda_list)
                p1_dict[cnsecutive_time_point].append(p1_list)
                p2_dict[cnsecutive_time_point].append(p2_list)

    for cnsecutive_time_point_list in cnsecutive_time_points:
        cnsecutive_time_point = '_'.join(cnsecutive_time_point_list)
        time_point_yong = cnsecutive_time_point_list[0]
        time_point_old = cnsecutive_time_point_list[1]
        lambda_list = lambda_dict[cnsecutive_time_point]
        p1_list = p1_dict[cnsecutive_time_point]
        p2_list = p2_dict[cnsecutive_time_point]
        lambda_df = pd.DataFrame(lambda_list, index=cell_ontology_class_ticks_dict[cnsecutive_time_point],
                                 columns=cell_ontology_class_ticks_dict[cnsecutive_time_point])
        p1_df = pd.DataFrame(p1_list, index=cell_ontology_class_ticks_dict[cnsecutive_time_point],
                             columns=cell_ontology_class_ticks_dict[cnsecutive_time_point])
        p2_df = pd.DataFrame(p2_list, index=cell_ontology_class_ticks_dict[cnsecutive_time_point],
                             columns=cell_ontology_class_ticks_dict[cnsecutive_time_point])
        lambda_df.to_csv(f'{image_folder_path}/{cnsecutive_time_point}/lambda_{time_point_yong}_{time_point_old}.csv')
        p1_df.to_csv(f'{image_folder_path}/{cnsecutive_time_point}/p1_{time_point_yong}_{time_point_old}.csv')
        p2_df.to_csv(f'{image_folder_path}/{cnsecutive_time_point}/p2_{time_point_yong}_{time_point_old}.csv')
    print('finish')

else:
    print('not start')


if __name__ == '__main__':
    tasks_dict = {}
    for cnsecutive_time_point_list in cnsecutive_time_points:
        cnsecutive_time_point = '_'.join(cnsecutive_time_point_list)
        time_point_yong = cnsecutive_time_point_list[0]
        time_point_old = cnsecutive_time_point_list[1]
        lambda_dict[cnsecutive_time_point] = [
            [0 for i in range(len(cell_ontology_class_ticks_dict[cnsecutive_time_point]))] for j in
            range(len(cell_ontology_class_ticks_dict[cnsecutive_time_point]))]
        p1_dict[cnsecutive_time_point] = [[0 for i in range(len(cell_ontology_class_ticks_dict[cnsecutive_time_point]))]
                                          for j in range(len(cell_ontology_class_ticks_dict[cnsecutive_time_point]))]
        p2_dict[cnsecutive_time_point] = [[0 for i in range(len(cell_ontology_class_ticks_dict[cnsecutive_time_point]))]
                                          for j in range(len(cell_ontology_class_ticks_dict[cnsecutive_time_point]))]
        tasks_dict[cnsecutive_time_point] = {}
        image_folder_path_pca_dict[cnsecutive_time_point] = {}
        i = 0
        for tissue1 in tissue_age_dict[cnsecutive_time_point]:
            image_folder_path_pca_dict[cnsecutive_time_point][tissue1] = {}
            tasks_dict[cnsecutive_time_point][tissue1] = {}
            for cell_ontology_class1 in cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue1]:
                image_folder_path_pca_dict[cnsecutive_time_point][tissue1][cell_ontology_class1] = {}
                tasks_dict[cnsecutive_time_point][tissue1][cell_ontology_class1] = []
                j = 0
                for tissue2 in tissue_age_dict[cnsecutive_time_point]:
                    image_folder_path_pca_dict[cnsecutive_time_point][tissue1][cell_ontology_class1][tissue2] = {}
                    for cell_ontology_class2 in cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue2]:

                        image_folder_path_pca_dict[cnsecutive_time_point][tissue1][cell_ontology_class1][tissue2][
                            cell_ontology_class2] = f'{image_folder_path}/{cnsecutive_time_point}/{tissue1}/{cell_ontology_class1}/{tissue2}/{cell_ontology_class2}'
                        if not os.path.exists(
                                image_folder_path_pca_dict[cnsecutive_time_point][tissue1][cell_ontology_class1][
                                    tissue2][cell_ontology_class2]):
                            os.makedirs(
                                image_folder_path_pca_dict[cnsecutive_time_point][tissue1][cell_ontology_class1][
                                    tissue2][cell_ontology_class2])
                        tasks_dict[cnsecutive_time_point][tissue1][cell_ontology_class1].append(
                            (adata_cell_ontology_class_tissue_age_dict, cnsecutive_time_point,
                             cnsecutive_time_point_list, time_point_yong,
                             time_point_old, tissue1, cell_ontology_class1, tissue2, cell_ontology_class2,
                             image_folder_path_pca_dict, i, j))
                        j += 1
                i += 1

    # multiprocessing pool map function to run the tasks in parallel and save the results
    results_dict = {}
    for cnsecutive_time_point_list in cnsecutive_time_points:
        cnsecutive_time_point = '_'.join(cnsecutive_time_point_list)
        results_dict[cnsecutive_time_point] = {}
        for tissue1 in tissue_age_dict[cnsecutive_time_point]:
            results_dict[cnsecutive_time_point][tissue1] = {}
            for cell_ontology_class1 in cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue1]:
                with Pool(5) as pool:
                    results_dict[cnsecutive_time_point][tissue1][cell_ontology_class1] = pool.map(
                        define_function.worker, tasks_dict[cnsecutive_time_point][tissue1][cell_ontology_class1])
                    print(f'results={results_dict[cnsecutive_time_point][tissue1][cell_ontology_class1]}')
                    pool.close()

    for cnsecutive_time_point_list in cnsecutive_time_points:
        cnsecutive_time_point = '_'.join(cnsecutive_time_point_list)
        for tissue1 in tissue_age_dict[cnsecutive_time_point]:
            for cell_ontology_class1 in cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue1]:
                for result in results_dict[cnsecutive_time_point][tissue1][cell_ontology_class1]:
                    cnsecutive_time_point, lambda_, p1, p2, i, j = result
                    lambda_dict[cnsecutive_time_point][i][j] = lambda_
                    p1_dict[cnsecutive_time_point][i][j] = p1
                    p2_dict[cnsecutive_time_point][i][j] = p2

    print(lambda_dict['3m_18m'])

time_point_list = []
tissue_list = []
cell_ontology_class_list = []

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
            image_folder_path_pca_dict[cnsecutive_time_point][tissue1][cell_ontology_class1] = {}
            time_point_list.append(cnsecutive_time_point)
            tissue_list.append(tissue1)
            cell_ontology_class_list.append(cell_ontology_class1)

i = 0
while i < len(tissue_list):
    tissue1 = tissue_list[i]
    cnsecutive_time_point = time_point_list[i]
    time_point_yong = cnsecutive_time_point.split('_')[0]
    time_point_old = cnsecutive_time_point.split('_')[1]
    cell_ontology_class1 = cell_ontology_class_list[i]
    adata1_yong = \
        adata_cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue1][cell_ontology_class1][
            time_point_yong]
    adata1_old = \
        adata_cell_ontology_class_tissue_age_dict[cnsecutive_time_point][tissue1][cell_ontology_class1][
            time_point_old]
    lambda_list = []
    p1_list = []
    p2_list = []
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
                    adata1_yong, adata1_old, adata2_yong, adata2_old,
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
                del adata1_yong_pca, adata1_old_pca, adata2_yong_pca, cumulative_explained_variance_ratio, lambda_, p1, p2, adata2_yong, adata2_old
                gc.collect() # collect garbage to save memory

            else:
                lambda_list.append(0)
                p1_list.append(0)
                p2_list.append(0)

    del adata1_yong, adata1_old
    gc.collect() # collect garbage to save memory

    lambda_dict[cnsecutive_time_point].append(lambda_list)
    p1_dict[cnsecutive_time_point].append(p1_list)
    p2_dict[cnsecutive_time_point].append(p2_list)

    lambda_list_cosecutive_time_point = lambda_dict[cnsecutive_time_point]
    p1_list_cosecutive_time_point = p1_dict[cnsecutive_time_point]
    p2_list_cosecutive_time_point = p2_dict[cnsecutive_time_point]

    # save lambda, p1, p2
    lambda_list_cosecutive_time_point = np.array(lambda_list_cosecutive_time_point)
    p1_list_cosecutive_time_point = np.array(p1_list_cosecutive_time_point)
    p2_list_cosecutive_time_point = np.array(p2_list_cosecutive_time_point)
    np.save(f'{data_folder_path}/{cnsecutive_time_point}/{tissue1}/lambda_list_{tissue1}_{cnsecutive_time_point}.npy',
            lambda_list_cosecutive_time_point)
    np.save(f'{data_folder_path}/{cnsecutive_time_point}/{tissue1}/p1_list_{tissue1}_{cnsecutive_time_point}.npy',
            p1_list_cosecutive_time_point)
    np.save(f'{data_folder_path}/{cnsecutive_time_point}/{tissue1}/p2_list_{tissue1}_{cnsecutive_time_point}.npy',
            p2_list_cosecutive_time_point)
    print(i)
    i += 1

# load lambda, p1, p2
lambda_list_cosecutive_time_point = np.load(
    f'{data_folder_path}/{cnsecutive_time_point}/{tissue1}/lambda_list_{tissue1}_{cnsecutive_time_point}.npy',
    allow_pickle=True)
p1_list_cosecutive_time_point = np.load(
    f'{data_folder_path}/{cnsecutive_time_point}/{tissue1}/p1_list_{tissue1}_{cnsecutive_time_point}.npy',
    allow_pickle=True)
p2_list_cosecutive_time_point = np.load(
    f'{data_folder_path}/{cnsecutive_time_point}/{tissue1}/p2_list_{tissue1}_{cnsecutive_time_point}.npy',
    allow_pickle=True)

print(lambda_list_cosecutive_time_point)
print(p1_list_cosecutive_time_point)
print(p2_list_cosecutive_time_point)
'''
