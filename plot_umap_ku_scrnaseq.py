import os
import scanpy as sc
import define_function as d_f
import anndata as ad
import numpy as np

organ_list = ['Blood', 'BM', 'Brain', 'Liver', 'LN', 'SmallInt', 'Spleen']
disease_list = ['WT', 'AD']
color_list = ['dodgerblue', 'orange', 'limegreen', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan', 'black']

cell_type_ref = {'Blood': 'MCA_broad_cell_type',
                 'BM': 'Tabula.Muris_cell_type',
                 'Brain': 'Azimuth_predicted.subclass',
                 'Liver': 'Tabula.Muris_cell_type',
                 'LN': 'MCA_Blood_broad_cell_type',
                 'SmallInt': 'MCA_broad_cell_type',
                 'Spleen': 'Tabula.Muris_cell_type'
                 }
ms_path = '/Users/kataokayuunosuke/MS'
image_path = f'{ms_path}/image'

for organ in organ_list:
    organ_path = f'{ms_path}/anndata/{organ}_all.integrated_doublet.removed_annotated.h5ad'
    organ_adata = sc.read_h5ad(organ_path)  # read the data
    time_point_list = organ_adata.obs['time.point'].unique()
    color_age_dict = {}  # store the color for each age
    organ_image_path = f'{image_path}/{organ}'
    if not os.path.exists(organ_image_path):
        os.mkdir(organ_image_path)

    for i in range(len(time_point_list)):
        color_age_dict[time_point_list[i]] = color_list[i]

    # plot umap
    print(organ_adata.obsm['X_umap'])
    first_component_umap = organ_adata.obsm['X_umap'][:, 0]
    second_component_umap = organ_adata.obsm['X_umap'][:, 1]

    max_value_umap_1 = first_component_umap.max()
    min_value_umap_1 = first_component_umap.min()
    max_value_umap_2 = second_component_umap.max()
    min_value_umap_2 = second_component_umap.min()

    d_f.adjust_and_save_plot(organ_adata, organ_image_path, organ, 'umap', min_value_umap_1,
                             max_value_umap_1, min_value_umap_2, max_value_umap_2, 'UMAP1', 'UMAP2',
                             cell_type_ref[organ])

    d_f.adjust_and_save_plot(organ_adata, organ_image_path, organ, 'umap', min_value_umap_1,
                             max_value_umap_1, min_value_umap_2, max_value_umap_2, 'UMAP1', 'UMAP2',
                             'samples4')

    for disease in disease_list:
        for time_point in time_point_list:
            adata = organ_adata[(organ_adata.obs['time.point'] == time_point) & (organ_adata.obs['disease'] == disease), :]

            d_f.adjust_and_save_plot(adata, organ_image_path, f'{disease}_{time_point}_{organ}', 'umap', min_value_umap_1,
                                     max_value_umap_1, min_value_umap_2, max_value_umap_2, 'UMAP1', 'UMAP2',
                                     cell_type_ref[organ])

