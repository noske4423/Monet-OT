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
# adata.var: stores variable annotation (usually one column per gene) as data frame.

adata = adata[adata.obs.sex == 'male']  # select only male cells
print(3)
mean_expression = adata.X.mean(axis=0)  # mean expression of each gene
fig = pl.figure(figsize=(3, 3))
pl.plot(mean_expression.tolist()[0])
pl.xlabel('Genes')
pl.ylabel('Mean expression')
pl.title('Mean expression of each gene')
pl.show()
fig.savefig(f'{folder_path}/facs_male_mean_expression.png')

adata = define_function.sort_adata_by_attribute(adata, 'tissue')  # sort by tissue
age_dict = define_function.split_adata_by_attribute(adata, 'age')  # split by age

os.chdir(folder_path)  # Change the current directory to the desired folder

adata_3_18 = ad.concat([age_dict['3m'], age_dict['18m']])  # select 3m and 18m
sc.tl.pca(adata_3_18, n_comps=20)  # PCA
sc.pl.pca(adata_3_18, components=['1,2', '2,3', '1,3'], color='age', save=f'_facs_male_3m_18m.png',
          title=['facs_male_3m_18m', 'facs_male_3m_18m', 'facs_male_3m_18m'])

adata_18_24 = ad.concat([age_dict['18m'], age_dict['24m']])  # select 18m and 24m
sc.tl.pca(adata_18_24, n_comps=20)  # PCA
sc.pl.pca(adata_18_24, components=['1,2', '2,3', '1,3'], color='age', save=f'_facs_male_18m_24m.png',
          title=['facs_male_3m_18m', 'facs_male_3m_18m', 'facs_male_3m_18m'])  # plot PCA

# umap
sc.pp.neighbors(adata_3_18, n_neighbors=10, n_pcs=20)
sc.tl.umap(adata_3_18)
sc.pl.umap(adata_3_18, color=['age', 'tissue'], save='_facs_male_3m_18m.png',
           title=['facs_male_3m_18m', 'facs_male_3m_18m'])

sc.pp.neighbors(adata_18_24, n_neighbors=10, n_pcs=20)
sc.tl.umap(adata_18_24)
sc.pl.umap(adata_18_24, color=['age', 'tissue'], save=f'_facs_male_18m_24m.png',
           title=['facs_male_18m_24m', 'facs_male_18m_24m'])

tissue_dict = {}
for age in ['3m', '18m', '24m']:
    tissue_dict[age] = define_function.split_adata_by_attribute(age_dict[age], 'tissue')
'''
# save the PCA and umap plot for each tissue
for tissue in list(tissue_dict['3m'].keys()):
    # PCA
    sc.pl.pca(tissue_dict['3m'][tissue], components='1,2', color=['age', 'cell_ontology_class'],
              save=f'_facs_male_3m_{tissue}.png',
              title=['facs_male_3m_' + tissue, 'facs_male_3m_' + tissue])
    sc.pl.pca(tissue_dict['18m'][tissue], components='1,2', color=['age', 'cell_ontology_class'],
              save=f'_facs_male_18m_{tissue}.png',
              title=['facs_male_18m_' + tissue, 'facs_male_18m_' + tissue])
    sc.pl.pca(tissue_dict['24m'][tissue], components='1,2', color=['age', 'cell_ontology_class'],
              save=f'_facs_male_24m_{tissue}.png',
              title=['facs_male_24m_' + tissue, 'facs_male_24m_' + tissue])

    # umap
    sc.pl.umap(tissue_dict['3m'][tissue], color=['age', 'cell_ontology_class'], save=f'_facs_male_3m_{tissue}.png',
               title=['facs_male_3m_' + tissue, 'facs_male_3m_' + tissue])
    sc.pl.umap(tissue_dict['18m'][tissue], color=['age', 'cell_ontology_class'], save=f'_facs_male_18m_{tissue}.png',
               title=['facs_male_18m_' + tissue, 'facs_male_18m_' + tissue])
    sc.pl.umap(tissue_dict['24m'][tissue], color=['age', 'cell_ontology_class'], save=f'_facs_male_24m_{tissue}.png',
               title=['facs_male_24m_' + tissue, 'facs_male_24m_' + tissue])
'''
for tissue in list(tissue_dict['3m'].keys()):
    integrate_adata = ad.concat([tissue_dict['3m'][tissue], tissue_dict['18m'][tissue], tissue_dict['24m'][tissue]])
    integrate_adata_3m_18m = ad.concat([tissue_dict['3m'][tissue], tissue_dict['18m'][tissue]])
    integrate_adata_18m_24m = ad.concat([tissue_dict['18m'][tissue], tissue_dict['24m'][tissue]])

    # PCA
    # integrate_adata
    first_component_pca = integrate_adata.obsm['X_pca'][:, 0]
    second_component_pca = integrate_adata.obsm['X_pca'][:, 1]

    max_value_pca_1 = first_component_pca.max()
    min_value_pca_1 = first_component_pca.min()
    max_value_pca_2 = second_component_pca.max()
    min_value_pca_2 = second_component_pca.min()

    fig = sc.pl.pca(integrate_adata, color='age', return_fig=True)

    ax = fig.gca()
    ax.set_xlim(min_value_pca_1 - (max_value_pca_1 - min_value_pca_1) * 0.05,
                max_value_pca_1 + (max_value_pca_1 - min_value_pca_1) * 0.05)  # adjust the x axis
    ax.set_ylim(min_value_pca_2 - (max_value_pca_2 - min_value_pca_2) * 0.05,
                max_value_pca_2 + (max_value_pca_2 - min_value_pca_2) * 0.05)  # adjust the y axis
    ax.set_title(f'pca_facs_male_{tissue}_age.png')
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')

    # save the plot
    pl.savefig(f'pca_facs_male_{tissue}_age.png')

    fig = sc.pl.pca(integrate_adata, color='cell_ontology_class', return_fig=True)


    ax = fig.gca()
    ax.set_xlim(min_value_pca_1 - (max_value_pca_1 - min_value_pca_1) * 0.05,
                max_value_pca_1 + (max_value_pca_1 - min_value_pca_1) * 0.05)  # adjust the x axis
    ax.set_ylim(min_value_pca_2 - (max_value_pca_2 - min_value_pca_2) * 0.05,
                max_value_pca_2 + (max_value_pca_2 - min_value_pca_2) * 0.05)  # adjust the y axis
    ax.set_title(f'pca_facs_male_{tissue}_cell_ontology_class.png')
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')

    # save the plot
    pl.savefig(f'pca_facs_male_{tissue}_cell_ontology_class.png')

    # integrate_adata_3m_18m
    fig = sc.pl.pca(integrate_adata_3m_18m, color='age', return_fig=True)

    ax = fig.gca()
    ax.set_xlim(min_value_pca_1, max_value_pca_1)
    ax.set_ylim(min_value_pca_2, max_value_pca_2)
    ax.set_title(f'pca_facs_male_{tissue}_3m_18m_age.png')
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')

    # save the plot
    pl.savefig(f'pca_facs_male_{tissue}_3m_18m_age.png')

    fig = sc.pl.pca(integrate_adata_3m_18m, color='cell_ontology_class', return_fig=True)

    ax = fig.gca()
    ax.set_xlim(min_value_pca_1, max_value_pca_1)
    ax.set_ylim(min_value_pca_2, max_value_pca_2)
    ax.set_title(f'pca_facs_male_{tissue}_3m_18m_cell_ontology_class.png')
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')

    # save the plot
    pl.savefig(f'pca_facs_male_{tissue}_3m_18m_cell_ontology_class.png')

    # integrate_adata_18m_24m
    fig = sc.pl.pca(integrate_adata_18m_24m, color='age', return_fig=True)

    ax = fig.gca()
    ax.set_xlim(min_value_pca_1, max_value_pca_1)
    ax.set_ylim(min_value_pca_2, max_value_pca_2)
    ax.set_title(f'pca_facs_male_{tissue}_18m_24m_age.png')
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')

    # save the plot
    pl.savefig(f'pca_facs_male_{tissue}_18m_24m_age.png')

    # tissue_dict['3m'][tissue]
    fig = sc.pl.pca(tissue_dict['3m'][tissue], return_fig=True)

    ax = fig.gca()
    ax.set_xlim(min_value_pca_1, max_value_pca_1)
    ax.set_ylim(min_value_pca_2, max_value_pca_2)
    ax.set_title(f'pca_facs_male_{tissue}_3m.png')
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')

    # save the plot
    pl.savefig(f'pca_facs_male_{tissue}_3m.png')

    fig = sc.pl.pca(tissue_dict['3m'][tissue],color='cell_ontology_class', return_fig=True)

    ax = fig.gca()
    ax.set_xlim(min_value_pca_1, max_value_pca_1)
    ax.set_ylim(min_value_pca_2, max_value_pca_2)
    ax.set_title(f'pca_facs_male_{tissue}_3m_cell_ontology_class.png')
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')

    # save the plot
    pl.savefig(f'pca_facs_male_{tissue}_3m_cell_ontology_class.png')

    # tissue_dict['18m'][tissue]
    fig = sc.pl.pca(tissue_dict['18m'][tissue], return_fig=True)

    ax = fig.gca()
    ax.set_xlim(min_value_pca_1, max_value_pca_1)
    ax.set_ylim(min_value_pca_2, max_value_pca_2)
    ax.set_title(f'pca_facs_male_{tissue}_18m.png')
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')

    # save the plot
    pl.savefig(f'pca_facs_male_{tissue}_18m.png')

    fig = sc.pl.pca(tissue_dict['18m'][tissue],color = 'cell_ontology_class' return_fig=True)

    ax = fig.gca()
    ax.set_xlim(min_value_pca_1, max_value_pca_1)
    ax.set_ylim(min_value_pca_2, max_value_pca_2)
    ax.set_title(f'pca_facs_male_{tissue}_18m_cell_ontology_class.png')
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')

    # save the plot
    pl.savefig(f'pca_facs_male_{tissue}_18m_cell_ontology_class.png')

    # tissue_dict['24m'][tissue]
    fig = sc.pl.pca(tissue_dict['24m'][tissue], return_fig=True)

    ax = fig.gca()
    ax.set_xlim(min_value_pca_1, max_value_pca_1)
    ax.set_ylim(min_value_pca_2, max_value_pca_2)
    ax.set_title(f'pca_facs_male_{tissue}_24m.png')
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')

    # save the plot
    pl.savefig(f'pca_facs_male_{tissue}_24m.png')

    fig = sc.pl.pca(tissue_dict['24m'][tissue],color = 'cell_ontology_class' return_fig=True)

    ax = fig.gca()
    ax.set_xlim(min_value_pca_1, max_value_pca_1)
    ax.set_ylim(min_value_pca_2, max_value_pca_2)
    ax.set_title(f'pca_facs_male_{tissue}_24m_cell_ontology_class.png')
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')

    # save the plot
    pl.savefig(f'pca_facs_male_{tissue}_24m_cell_ontology_class.png')

    # umap
    # integrate_adata
    first_component_umap = integrate_adata.obsm['X_umap'][:, 0]
    second_component_umap = integrate_adata.obsm['X_umap'][:, 1]

    max_value_umap_1 = first_component_umap.max()
    min_value_umap_1 = first_component_umap.min()
    max_value_umap_2 = second_component_umap.max()
    min_value_umap_2 = second_component_umap.min()

    fig = sc.pl.umap(integrate_adata, color=['age', 'cell_ontology_class'], return_fig=True)

    ax = fig.gca()
    ax.set_xlim(min_value_umap_1, max_value_umap_1)
    ax.set_ylim(min_value_umap_2, max_value_umap_2)
    ax.set_title(f'umap_facs_male_{tissue}.png')
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')

    # save the plot
    pl.savefig(f'umap_facs_male_{tissue}.png')

    # integrate_adata_3m_18m
    fig = sc.pl.umap(integrate_adata_3m_18m, color=['age', 'cell_ontology_class'], return_fig=True)

    ax = fig.gca()
    ax.set_xlim(min_value_umap_1, max_value_umap_1)
    ax.set_ylim(min_value_umap_2, max_value_umap_2)
    ax.set_title(f'umap_facs_male_{tissue}_3m_18m.png')
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')

    # save the plot
    pl.savefig(f'umap_facs_male_{tissue}_3m_18m.png')

    # integrate_adata_18m_24m
    fig = sc.pl.umap(integrate_adata_18m_24m, color=['age', 'cell_ontology_class'], return_fig=True)

    ax = fig.gca()
    ax.set_xlim(min_value_umap_1, max_value_umap_1)
    ax.set_ylim(min_value_umap_2, max_value_umap_2)
    ax.set_title(f'umap_facs_male_{tissue}_18m_24m.png')
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')

    # save the plot
    pl.savefig(f'umap_facs_male_{tissue}_18m_24m.png')

tissue_3_18_dict, ontology_3_18_dict, ontology_3_18_list, ontology_3_18_list_dict = define_function.process_adata_by_age(
    age_dict, ['3m', '18m'])
# ontology_3_18_list: the list of ontology that exists in both 3m and 18m
# ontolofy_3_18_list_dict[tissue]: the list of ontology that exists in both 3m and 18m in tissue
# ontology_3_18_dict[age][tissue]: the adata of age and tissue

tissue_18_24_dict, ontology_18_24_dict, ontology_18_24_list, ontology_18_24_list_dict = define_function.process_adata_by_age(
    age_dict, ['18m', '24m'])
# ontology_18_24_list: the list of ontology that exists in both 18m and 24m
# ontolofy_18_24_list_dict[tissue]: the list of ontology that exists in both 18m and 24m in tissue
# ontology_18_24_dict[age][tissue][ontology]: the adata of age and tissue and ontology

# save the PCA and umap plot for each cell ontology
for tissue in list(ontology_3_18_list_dict.keys()):
    for ontology in ontology_3_18_list_dict[tissue]:
        integrate_adata = ad.concat(
            [ontology_3_18_dict['3m'][tissue][ontology], ontology_3_18_dict['18m'][tissue][ontology]])

        # PCA
        # sc.tl.pca(integrate_adata, n_comps=20)
        sc.pl.pca(integrate_adata, components='1,2', color='age', save=f'_facs_male_3m_18m_{tissue}_{ontology}.png',
                  title='facs_male_3m_18m_' + tissue + '_' + ontology)

        # umap
        # sc.pp.neighbors(integrate_adata, n_neighbors=10, n_pcs=20)
        # sc.tl.umap(integrate_adata)
        sc.pl.umap(integrate_adata, color='age', save=f'_facs_male_3m_18m_{tissue}_{ontology}.png',
                   title='facs_male_3m_18m_' + tissue + '_' + ontology)

for tissue in list(ontology_18_24_list_dict.keys()):
    for ontology in ontology_18_24_list_dict[tissue]:
        integrate_adata = ad.concat(
            [ontology_18_24_dict['18m'][tissue][ontology], ontology_18_24_dict['24m'][tissue][ontology]])

        # PCA
        # sc.tl.pca(integrate_adata, n_comps=20)
        sc.pl.pca(integrate_adata, components='1,2', color='age', save=f'_facs_male_18m_24m_{tissue}_{ontology}.png',
                  title='facs_male_18m_24m_' + tissue + '_' + ontology)

        # umap
        # sc.pp.neighbors(integrate_adata, n_neighbors=10, n_pcs=20)
        # sc.tl.umap(integrate_adata)
        sc.pl.umap(integrate_adata, color='age', save=f'_facs_male_18m_24m_{tissue}_{ontology}.png',
                   title='facs_male_18m_24m_' + tissue + '_' + ontology)

os.chdir(original_dir)  # Change the current directory back to the original folder
'''
tissue_3_18_dict = {}
for age in ['3m', '18m']:
    tissue_3_18_dict[age] = define_function.split_adata_by_attribute(age_dict[age], 'tissue')

print(list(tissue_3_18_dict['3m'].keys()))
print(tissue_3_18_dict['3m']['Aorta'].obs)

ontology_3_18_dict = {}
for age in ['3m', '18m']:
    ontology_3_18_dict[age] = {}
    for tissue in list(tissue_3_18_dict[age].keys()):
        print(tissue)
        ontology_3_18_dict[age][tissue] = define_function.split_adata_by_attribute(tissue_3_18_dict[age][tissue], 'cell_ontology_class')
print(ontology_3_18_dict['3m']['Aorta'].keys())

# ontology_3_18_list: the list of ontology that exists in both 3m and 18m
# ontolofy_3_18_list_dict[tissue]: the list of ontology that exists in both 3m and 18m in tissue
ontology_3_18_list = []
ontology_3_18_list_dict = {}
for tissue in list(ontology_3_18_dict['3m'].keys()):
    ontology_3_18_list_dict[tissue] = []
    for ontology in list(ontology_3_18_dict['3m'][tissue].keys()):
        if ontology in list(ontology_3_18_dict['18m'][tissue].keys()):
            if ontology_3_18_dict['3m'][tissue][ontology].n_obs > 20 and ontology_3_18_dict['18m'][tissue][ontology].n_obs > 20:
                ontology_3_18_list_dict[tissue].append(ontology)
                ontology_3_18_list.append(ontology + '(' + tissue + ')')
print(ontology_3_18_list)
print(ontology_3_18_list_dict['Aorta'])

# choose the ontology from ontology_3_18_list
# reconstruct the ontology_3_18_dict
ontology_3_18_dict_new = {}
for age in ['3m', '18m']:
    ontology_3_18_dict_new[age] = {}
    for tissue in list(tissue_3_18_dict[age].keys()):
        print(tissue)
        ontology_3_18_dict_new[age][tissue] = {}
        for ontology in list(ontology_3_18_dict[age][tissue].keys()):
            if ontology in ontology_3_18_list_dict[tissue]:
                ontology_3_18_dict_new[age][tissue][ontology] = ontology_3_18_dict[age][tissue][ontology]
print(ontology_3_18_dict_new['3m']['Aorta'].keys())
ontology_3_18_dict = ontology_3_18_dict_new # update the ontology_3_18_dict

tissue_18_24_dict = {}
for age in ['18m', '24m']:
    tissue_18_24_dict[age] = define_function.split_adata_by_attribute(age_dict[age], 'tissue')

print(list(tissue_18_24_dict['18m'].keys()))

ontology_18_24_dict = {}
for age in ['18m', '24m']:
    ontology_18_24_dict[age] = {}
    for tissue in list(tissue_18_24_dict[age].keys()):
        print(tissue)
        ontology_18_24_dict[age][tissue] = define_function.split_adata_by_attribute(tissue_18_24_dict[age][tissue], 'cell_ontology_class')
print(ontology_18_24_dict['18m']['Aorta'].keys())

# ontology_18_24_list: the list of ontology that exists in both 18m and 24m
# ontolofy_18_24_list_dict[tissue]: the list of ontology that exists in both 18m and 24m in tissue
ontology_18_24_list = []
ontology_18_24_list_dict = {}
for tissue in list(ontology_18_24_dict['18m'].keys()):
    ontology_18_24_list_dict[tissue] = []
    for ontology in list(ontology_18_24_dict['18m'][tissue].keys()):
        if ontology in list(ontology_18_24_dict['24m'][tissue].keys()):
            if ontology_18_24_dict['18m'][tissue][ontology].n_obs > 20 and ontology_18_24_dict['24m'][tissue][ontology].n_obs > 20:
                ontology_18_24_list_dict[tissue].append(ontology)
                ontology_18_24_list.append(ontology + '(' + tissue + ')')
print(ontology_18_24_list)
print(ontology_18_24_list_dict['Aorta'])

# choose the ontology from ontology_18_24_list
# reconstruct the ontology_18_24_dict
ontology_18_24_dict_new = {}
for age in ['18m', '24m']:
    ontology_18_24_dict_new[age] = {}
    for tissue in list(tissue_18_24_dict[age].keys()):
        print(tissue)
        ontology_18_24_dict_new[age][tissue] = {}
        for ontology in list(ontology_18_24_dict[age][tissue].keys()):
            if ontology in ontology_18_24_list_dict[tissue]:
                ontology_18_24_dict_new[age][tissue][ontology] = ontology_18_24_dict[age][tissue][ontology]

print(ontology_18_24_dict_new['18m']['Aorta'].keys())
'''
# save the ontology list
with open('data/facs/ontology_facs_male_3m_18m_list.txt', 'wb') as fp:
    pickle.dump(ontology_3_18_list, fp)

with open('data/facs/ontology_facs_male_18m_24m_list.txt', 'wb') as fp:
    pickle.dump(ontology_18_24_list, fp)

# load the ontology list
with open('data/facs/ontology_facs_male_3m_18m_list.txt', 'rb') as fp:
    ontology_3_18_list = pickle.load(fp)

with open('data/facs/ontology_facs_male_18m_24m_list.txt', 'rb') as fp:
    ontology_18_24_list = pickle.load(fp)

print(ontology_3_18_list)
print(ontology_18_24_list)

# save the ontology list dict as csv
with open('data/facs/ontology_facs_male_3m_18m_list_dict.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    for key, value in ontology_3_18_list_dict.items():
        writer.writerow([key, value])

with open('data/facs/ontology_facs_male_18m_24m_list_dict.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    for key, value in ontology_18_24_list_dict.items():
        writer.writerow([key, value])

# load the ontology list dict as csv
with open('data/facs/ontology_facs_male_3m_18m_list_dict.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    ontology_3_18_list_dict = dict(reader)

with open('data/facs/ontology_facs_male_18m_24m_list_dict.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    ontology_18_24_list_dict = dict(reader)

print(ontology_3_18_list_dict)
print(ontology_18_24_list_dict)

# save the ontology dict
with open('data/facs/ontology_facs_male_3m_18m_dict.txt', 'wb') as fp:
    pickle.dump(ontology_3_18_dict, fp)

with open('data/facs/ontology_facs_male_18m_24m_dict.txt', 'wb') as fp:
    pickle.dump(ontology_18_24_dict, fp)

# load the ontology dict
with open('data/facs/ontology_facs_male_3m_18m_dict.txt', 'rb') as fp:
    ontology_3_18_dict = pickle.load(fp)

with open('data/facs/ontology_facs_male_18m_24m_dict.txt', 'rb') as fp:
    ontology_18_24_dict = pickle.load(fp)

print(ontology_3_18_dict['3m']['Aorta']['aortic endothelial cell'].obs)
print(ontology_18_24_dict['18m']['Aorta']['aortic endothelial cell'].obs)
