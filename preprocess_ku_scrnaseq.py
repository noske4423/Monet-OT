import os
import scanpy as sc
import define_function as d_f
import anndata as ad
import numpy as np

organ_list = ['Blood', 'BM', 'Brain', 'Liver', 'LN', 'SmallInt', 'Spleen']
disease_list = ['WT', 'AD']
color_list = ['dodgerblue', 'orange', 'limegreen', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan', 'black']

cell_type_ref ={'Blood': 'MCA_broad_cell_type',
                'BM': 'Tabula.Muris_cell_type',
                'Brain': 'Azimuth_predicted.subclass',
                'Liver': 'Tabula.Muris_cell_type',
                'LN': 'MCA_Blood_broad_cell_type',
                'SmallInt': 'MCA_broad_cell_type',
                'Spleen': 'Tabula.Muris_cell_type'
                }
ms_path = '/Users/kataokayuunosuke/MS'
image_path = f'{ms_path}/image'

organ_adata_list = []
for organ in organ_list:
    organ_path = f'{ms_path}/anndata/{organ}_all.integrated_doublet.removed_annotated.h5ad'
    organ_adata = sc.read_h5ad(organ_path)  # read the data
    organ_adata_list.append(organ_adata)

    print(organ, organ_adata.n_obs, organ_adata.n_vars)
    # print(organ + '_Azimuth_predicted.celltype.l1=' + str(organ_adata.obs['Azimuth_predicted.celltype.l1'].unique()))
    # print(organ + '_Azimuth_predicted.celltype.l2=' + str(organ_adata.obs['Azimuth_predicted.celltype.l2'].unique()))
    # print(organ + '_MCA_broad_cell_type=' + str(organ_adata.obs['MCA_broad_cell_type'].unique()))
    print(organ + '_samples2=' + str(organ_adata.obs['samples2'].unique()))
    print(organ + '_samples3=' + str(organ_adata.obs['samples3'].unique()))
    print(organ + '_samples4=' + str(organ_adata.obs['samples4'].unique()))
    print(organ + '_time.point=' + str(organ_adata.obs['time.point'].unique()))
    print(organ + '_cell_type' + str(organ_adata.obs[cell_type_ref[organ]].unique()))

    cell_type_list = organ_adata.obs[cell_type_ref[organ]].unique()
    organ_adata.obs['cell_type'] = organ_adata.obs[cell_type_ref[organ]]
    time_point_list = organ_adata.obs['time.point'].unique()
    color_age_dict = {}  # store the color for each age
    for i in range(len(time_point_list)):
        color_age_dict[time_point_list[i]] = color_list[i]

    # save cell type list
    cell_type_list_path = ms_path
    if not os.path.exists(cell_type_list_path):
        os.makedirs(cell_type_list_path)

    cell_type_list_file = f'{cell_type_list_path}/{organ}_cell_type_list.txt'
    with open(cell_type_list_file, mode='w') as f:
        f.write('\n'.join(cell_type_list))

    # load cell type list
    cell_type_list_file = f'{cell_type_list_path}/{organ}_cell_type_list.txt'
    with open(cell_type_list_file, mode='r') as f:
        cell_type_list = f.read().splitlines()



integrate_adata = ad.concat(organ_adata_list, label='organ', keys=organ_list)
sc.pp.normalize_total(integrate_adata, target_sum=1e4, exclude_highly_expressed=True)
sc.pp.log1p(integrate_adata)

# sc.pp.filter_cells(integrate_adata, min_genes=250)
# sc.pp.filter_cells(integrate_adata, min_counts=2500)
# sc.pp.filter_genes(integrate_adata, min_cells=3)
print(integrate_adata.n_obs)
print(integrate_adata.n_vars)

# save integrated data
integrate_adata_path = f'{ms_path}/anndata/integrate_adata.h5ad'
integrate_adata.write(integrate_adata_path)

adata_sample_list = []
for organ in organ_list:
    # load celltype list
    celltype_list_path = f'{ms_path}/{organ}_cell_type_list.txt'
    with open(celltype_list_path) as f:
        celltype_list = f.read().splitlines()
    celltype_sample_list = []
    for celltype in celltype_list:
        cellnumber_list = []

        for disease in disease_list:
            for time_point in time_point_list:
                print(f'{disease}_{time_point}_{organ}_{celltype}')
                split_adata = \
                integrate_adata[(integrate_adata.obs['disease'] == disease)
                                & (integrate_adata.obs['time.point'] == time_point)
                                & (integrate_adata.obs['organ'] == organ)
                                & (integrate_adata.obs['cell_type'] == celltype)]
                print(split_adata.n_obs)
                cellnumber_list.append(split_adata.n_obs)
        if min(cellnumber_list) >= 50:
            celltype_sample_list.append(celltype)
            celltype_adata = integrate_adata[integrate_adata.obs['cell_type'] == celltype]

            celltype_adata_copy = celltype_adata.copy()
            sc.pp.highly_variable_genes(celltype_adata_copy, n_top_genes=2000)
            celltype_adata_copy = celltype_adata_copy[:, celltype_adata_copy.var['highly_variable']]
            high_variable_genes_list = celltype_adata_copy.var_names.tolist()
            # save high variable genes list
            high_variable_genes_list_path = f'{ms_path}/high_variable_genes/{organ}_{celltype}_high_variable_genes_list.txt'
            if not os.path.exists(f'{ms_path}/high_variable_genes'):
                os.makedirs(f'{ms_path}/high_variable_genes')
            with open(high_variable_genes_list_path, mode='w') as f:
                f.write('\n'.join(high_variable_genes_list))

            # sampling cell
            for disease in disease_list:
                for time_point in time_point_list:
                    print(f'{disease}_{time_point}_{organ}_{celltype}')
                    split_adata = \
                    integrate_adata[(integrate_adata.obs['disease'] == disease)
                                    & (integrate_adata.obs['time.point'] == time_point)
                                    & (integrate_adata.obs['organ'] == organ)
                                    & (integrate_adata.obs['cell_type'] == celltype)]
                    split_adata_sample = split_adata[np.random.choice(split_adata.obs.index, 50, replace=False)]

                    adata_sample_list.append(split_adata_sample)

    # save celltype sample list
    celltype_sample_list_path = f'{ms_path}/{organ}_cell_type_sample_list.txt'
    with open(celltype_sample_list_path, mode='w') as f:
        f.write('\n'.join(celltype_sample_list))

integrate_adata_sample = ad.concat(adata_sample_list)

# save anndata
integrate_adata_sample_path = f'{ms_path}/anndata/integrate_adata_sample.h5ad'
integrate_adata_sample.write(integrate_adata_sample_path)

