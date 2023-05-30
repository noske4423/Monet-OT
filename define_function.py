import anndata as ad
import numpy as np
import scanpy as sc
import random
import ot
import matplotlib.pylab as pl


def sort_adata_by_attribute(adata, attribute):
    """Sorts an AnnData object by the specified attribute.

    Args:
        adata (AnnData): The input AnnData object.
        attribute (str): The attribute to sort by.

    Returns:
        AnnData: The sorted AnnData object.
    """
    groups = adata.obs.groupby(attribute).indices
    return ad.concat([adata[inds] for inds in groups.values()], merge="same")


def split_adata_by_attribute(adata, attribute):
    """Splits an AnnData object into a dictionary of AnnData objects based on the specified attribute.

    Args:
        adata (AnnData): The input AnnData object.
        attribute (str): The attribute to split by.

    Returns:
        dict: A dictionary with keys as unique values of the attribute and values as the corresponding AnnData objects.
    """
    groups = adata.obs.groupby(attribute).indices
    return {key: adata[adata.obs[attribute] == key] for key in groups.keys()}


def process_adata_by_age(age_dict, ages):
    """
    This function processes anndata objects of two different ages, separates them by tissue and cell ontology class,
    and then filters the ontology classes based on their presence in both ages and the number of observations.

    Parameters:
    age_dict (dict): A dictionary containing anndata objects. The keys are age identifiers (like '3m', '18m')
                     and the values are the corresponding anndata objects.
    ages (list): A list of two age identifiers (like ['3m', '18m']). These are the ages that will be processed.

    Returns:
    tissue_dict (dict): A dictionary where keys are age identifiers and values are dictionaries.
                        Each inner dictionary has tissue types as keys and sub-anndata objects as values, split by tissue.
    ontology_dict_new (dict): A similar structure as `tissue_dict`, but split further by cell ontology classes.
    ontology_list (list): A list of ontology classes that are present in both ages and have more than 20 observations.
    ontology_list_dict (dict): A dictionary where keys are tissue types and values are lists of ontology classes
                               that are present in both ages and have more than 20 observations.
    """
    tissue_dict = {}
    ontology_dict = {}
    ontology_list_dict = {}
    ontology_list = []

    for age in ages:
        tissue_dict[age] = split_adata_by_attribute(age_dict[age], 'tissue')

    for age in ages:
        ontology_dict[age] = {}
        for tissue in list(tissue_dict[age].keys()):
            ontology_dict[age][tissue] = split_adata_by_attribute(tissue_dict[age][tissue],
                                                                  'cell_ontology_class')

    for tissue in list(ontology_dict[ages[0]].keys()):
        ontology_list_dict[tissue] = []
        for ontology in list(ontology_dict[ages[0]][tissue].keys()):
            if ontology in list(ontology_dict[ages[1]][tissue].keys()):
                if ontology_dict[ages[0]][tissue][ontology].n_obs > 20 and ontology_dict[ages[1]][tissue][
                    ontology].n_obs > 20:
                    ontology_list_dict[tissue].append(ontology)
                    ontology_list.append(ontology + '(' + tissue + ')')

    ontology_dict_new = {}
    for age in ages:
        ontology_dict_new[age] = {}
        for tissue in list(tissue_dict[age].keys()):
            ontology_dict_new[age][tissue] = {}
            for ontology in list(ontology_dict[age][tissue].keys()):
                if ontology in ontology_list_dict[tissue]:
                    ontology_dict_new[age][tissue][ontology] = ontology_dict[age][tissue][ontology]

    return tissue_dict, ontology_dict_new, ontology_list, ontology_list_dict


def plot_highest_expr_genes(adata, folder_path, title, n_genes=10):
    """Plots the highest expressed genes in the AnnData object.

    Args:
        adata (AnnData): The input AnnData object.
        folder_path (str): The path to the folder where the plot will be saved.
        title (str): The title of the plot.
        n_genes (int): The number of genes to plot.
    """
    ax = sc.pl.highest_expr_genes(adata, n_top=n_genes, show=False)
    ax.set_title(f'highest_expr_genes_{title}')
    pl.savefig(f'{folder_path}/highest_expr_genes_{title}.png')


def plot_mean_expression(adata, folder_path, title):
    """Plots the mean expression of each gene in the AnnData object.

    Args:
        adata (AnnData): The input AnnData object.
    """
    mean_expression = adata.X.mean(axis=0)  # mean expression of each gene
    fig = pl.figure(figsize=(5, 5))
    pl.plot(mean_expression.tolist()[0])
    pl.xlabel('genes')
    pl.ylabel('mean_expression')
    pl.title(f'mean_expression_{title}')
    fig.savefig(f'{folder_path}/mean_expression_{title}.png')


def adjust_and_save_plot(adata, folder_path, title, method, min_value_1, max_value_1, min_value_2, max_value_2, xlabel,
                         ylabel, color, color_name='blue'):
    """Adjusts and saves a plot.

    Args:
        adata (AnnData): The input AnnData object.
        folder_path (str): The path to the folder where the plot will be saved.
        title (str): The title of the plot.
        method (str): The method used to generate the plot.
        xlabel (str): The label of the x-axis.
        ylabel (str): The label of the y-axis.
        color (str): The color of the plot.
    """
    if method == 'pca':
        if color == None:
            title = f'{method}_{title}'
            pl.figure()
            pl.scatter(adata.obsm['X_pca'][:, 0], adata.obsm['X_pca'][:, 1], c=color_name)
            pl.xlim([min_value_1 - (max_value_1 - min_value_1) * 0.05,
                     max_value_1 + (max_value_1 - min_value_1) * 0.05])
            pl.ylim([min_value_2 - (max_value_2 - min_value_2) * 0.05,
                     max_value_2 + (max_value_2 - min_value_2) * 0.05])
            pl.title(title)
            pl.xlabel(xlabel)
            pl.ylabel(ylabel)
            pl.xticks([])  # remove xticks
            pl.yticks([])  # remove yticks

            pl.savefig(f'{folder_path}/{title}.png')
            pl.close()

        else:
            title = f'{method}_{title}_{color}'
            # sc.set_figure_params(figsize=[12, 6])
            fig = sc.pl.pca(adata, color=color, return_fig=True)
            ax = fig.gca()
            ax.set_xlim(min_value_1 - (max_value_1 - min_value_1) * 0.05,
                        max_value_1 + (max_value_1 - min_value_1) * 0.05)
            ax.set_ylim(min_value_2 - (max_value_2 - min_value_2) * 0.05,
                        max_value_2 + (max_value_2 - min_value_2) * 0.05)
            ax.set_title(title)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            legend = ax.get_legend()

            pl.savefig(f'{folder_path}/{title}.png', bbox_extra_artists=(legend,), bbox_inches='tight')

    elif method == 'umap':
        if color == None:
            title = f'{method}_{title}'
            pl.figure()
            pl.scatter(adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1], c=color_name)
            pl.xlim([min_value_1 - (max_value_1 - min_value_1) * 0.05,
                     max_value_1 + (max_value_1 - min_value_1) * 0.05])
            pl.ylim([min_value_2 - (max_value_2 - min_value_2) * 0.05,
                     max_value_2 + (max_value_2 - min_value_2) * 0.05])
            pl.title(title)
            pl.xlabel(xlabel)
            pl.ylabel(ylabel)
            pl.xticks([])  # remove xticks
            pl.yticks([])  # remove yticks

            pl.savefig(f'{folder_path}/{title}.png')
            pl.close()

        else:
            title = f'{method}_{title}_{color}'
            fig = sc.pl.umap(adata, color=color, return_fig=True)
            ax = fig.gca()
            ax.set_xlim(min_value_1 - (max_value_1 - min_value_1) * 0.05,
                        max_value_1 + (max_value_1 - min_value_1) * 0.05)
            ax.set_ylim(min_value_2 - (max_value_2 - min_value_2) * 0.05,
                        max_value_2 + (max_value_2 - min_value_2) * 0.05)
            ax.set_title(title)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            legend = ax.get_legend()

            pl.savefig(f'{folder_path}/{title}.png', bbox_extra_artists=(legend,), bbox_inches='tight')


def extract_edges_above_threshold(g, thr):
    """Extracts edges from a matrix with values above a specified threshold.

    Args:
        g (np.array): The input 2D array (matrix) from which to extract edges.
        thr (float): The threshold value.

    Returns:
        list: A list of edges in the form [row, column, value] where value > thr.
    """
    edge_list = [[i, j, g[i][j]] for i in range(g.shape[0]) for j in range(g.shape[1]) if g[i][j] > thr]
    return edge_list


def interpolate_edges(x1, x2, edge_list, t):
    """Interpolates edge values based on input weights and interpolation factor.

    Args:
        x1 (np.array): The first input array of weights.
        x2 (np.array): The second input array of weights.
        edge_list (list): A list of edges in the form [row, column, value].
        t (float): The interpolation factor (0 <= t <= 1).

    Returns:
        tuple: A tuple containing the interpolated edge values and the normalized weights.
    """
    x_list = [(1 - t) * x1[edge[0]] + t * x2[edge[1]] for edge in edge_list]
    w_list = [edge[2] for edge in edge_list]
    normalized_w_list = np.array(w_list) / np.sum(w_list)

    return np.array(x_list), normalized_w_list


def preprocess_data(a3, a18, b3, celltype_a, tissue_a, tissue_b):
    """Preprocesses the data for the interpolation.

    :param a3: List of AnnData objects for 3 months group A
    :param a18: List of AnnData objects for 18 months group A
    :param b3: List of AnnData objects for 3 months group B
    :param celltype_a: List of cell types for group A
    :param tissue_a: List of tissues for group A
    :param tissue_b: List of tissues for group B
    :return: adata_list
    """
    adata_list = []

    for i in range(len(a3)):
        adata1 = a3[i]
        adata3 = a18[i]

        for j in range(len(b3)):
            if celltype_a[i] == celltype_a[j] and tissue_a[i] == tissue_b[j]:
                adata_list.append([i, j, adata1, adata1, adata3])
            else:
                adata2 = b3[j]
                adata_new = ad.concat([adata1, adata2, adata3])
                sc.pp.highly_variable_genes(adata_new)
                adata_new = adata_new[:, adata_new.var.highly_variable]  # filter highly variable genes
                sc.tl.pca(adata_new, n_comps=50)  # run PCA

                w = adata_new.varm['PCs'].T
                adata1_new = adata_new[adata1.obs.index]  # filter cells
                adata2_new = adata_new[adata2.obs.index]
                adata3_new = adata_new[adata3.obs.index]

                x1 = np.array(adata1_new.obsm['X_pca'])  # get PCA coordinates
                x2 = np.array(adata2_new.obsm['X_pca'])
                x3 = np.array(adata3_new.obsm['X_pca'])

                adata_list.append([i, j, x1, x2, x3, w, adata1.obs.index, adata2.obs.index, adata3.obs.index])
                print(i, j)

    return adata_list


def apply_noise(adata_list, p1, p2, mu1, mu2, sigma1, sigma2, noise_idx):
    """Applies noise to the data.

    :param adata_list: List of arrays containing processed data
    :param p1: Proportion of noise for X2
    :param p2: Proportion of noise for X3
    :param mu1: Mean for noise in X2
    :param mu2: Mean for noise in X3
    :param sigma1: Standard deviation for noise in X2
    :param sigma2: Standard deviation for noise in X3
    :param noise_idx: Index for the noise dimension
    :return: adata_list_new
    """
    adata_list_new = []

    for adata_array in adata_list:
        x1, x2, x3, w = adata_array[2], adata_array[3], adata_array[4], adata_array[5]
        num2, num3 = int(p1 * x2.shape[0]), int(p2 * x3.shape[0])

        y2 = np.zeros_like(x2)
        y3 = np.zeros_like(x3)

        y2[random.sample(range(x2.shape[0]), num2), noise_idx] = np.random.normal(loc=mu1, scale=sigma1, size=num2)
        y3[random.sample(range(x3.shape[0]), num3), noise_idx] = np.random.normal(loc=mu2, scale=sigma2, size=num3)
        x2 += y2
        x3 += y3

        adata_list_new.append(
            [adata_array[0], adata_array[1], x1, x2, x3, w, adata_array[6], adata_array[7], adata_array[8]])

    return adata_list_new


def reconstruct_adata(adata_list):
    """Reconstructs the AnnData objects from the processed data.

    :param adata_list: List of arrays containing processed data
    :return: adata_list_new
    """
    adata_list_new = []
    for adata_array in adata_list:
        x1, x2, x3, w = adata_array[2], adata_array[3], adata_array[4], adata_array[5]

        x1_recon, x2_recon, x3_recon = np.dot(x1, w), np.dot(x2, w), np.dot(x3, w)

        adata1, adata2, adata3 = ad.AnnData(x1_recon), ad.AnnData(x2_recon), ad.AnnData(x3_recon)
        adata1.obs.index, adata2.obs.index, adata3.obs.index = adata_array[6], adata_array[7], adata_array[8]

        adata_new = ad.concat([adata1, adata2, adata3])
        sc.tl.pca(adata_new, n_comps=20)

        adata1_new = adata_new[adata1.obs.index]
        adata2_new = adata_new[adata2.obs.index]
        adata3_new = adata_new[adata3.obs.index]

        x1_new, x2_new, x3_new = np.array(adata1_new.obsm['X_pca']), np.array(
            adata2_new.obsm['X_pca']), np.array(adata3_new.obsm['X_pca'])

        adata_list_new.append([adata_array[0], adata_array[1], x1_new, x2_new, x3_new, np.array(adata2_new.X)])

    return adata_list_new


def process_pca(adata1, adata1_next, adata2, folder_name, title):
    tissue1 = adata1.obs['tissue'].unique()[0]
    tissue2 = adata2.obs['tissue'].unique()[0]
    cell_ontology_class1 = adata1.obs['cell_ontology_class'].unique()[0]
    cell_ontology_class2 = adata2.obs['cell_ontology_class'].unique()[0]
    adata1.obs['group'] = f'{tissue1}_{cell_ontology_class1}_yong'
    adata1_next.obs['group'] = f'{tissue1}_{cell_ontology_class1}_old'
    adata2.obs['group'] = f'{tissue2}_{cell_ontology_class2}_yong'

    integrate_adata = ad.concat([adata1, adata1_next, adata2])
    sc.pp.highly_variable_genes(integrate_adata)
    integrate_adata = integrate_adata[:, integrate_adata.var.highly_variable]  # filter highly variable genes

    adata1 = integrate_adata[integrate_adata.obs['group'] == f'{tissue1}_{cell_ontology_class1}_yong']
    adata1_next = integrate_adata[integrate_adata.obs['group'] == f'{tissue1}_{cell_ontology_class1}_old']
    adata2 = integrate_adata[integrate_adata.obs['group'] == f'{tissue2}_{cell_ontology_class2}_yong']

    # Set column average to 0
    adata1.X = adata1.X - adata1.X.mean(axis=0)
    adata1_next.X = adata1_next.X - adata1_next.X.mean(axis=0)
    adata2.X = adata2.X - adata2.X.mean(axis=0)

    integrate_adata = ad.concat([adata1, adata1_next, adata2])
    sc.tl.pca(integrate_adata, n_comps=50)  # run PCA
    fig = sc.pl.pca(integrate_adata, color='group', return_fig=True)
    legend = fig.get_axes()[0].get_legend()
    pl.title(title)
    fig.savefig(f'{folder_name}/pca_{title}.png', bbox_extra_artists=(legend,), bbox_inches='tight')

    sc.settings.figdir = f'{folder_name}'
    sc.pl.pca_variance_ratio(integrate_adata, log=True, n_pcs=50, save=f'{title}.png')
    sc.pl.pca_loadings(integrate_adata, components=[1, 2, 3], save=f'{title}.png')
    cumulative_explained_variance_ratio = integrate_adata.uns['pca']['variance_ratio'].sum()
    adata1.obsm['X_pca'] = integrate_adata[adata1.obs.index].obsm['X_pca']
    adata1_next.obsm['X_pca'] = integrate_adata[adata1_next.obs.index].obsm['X_pca']
    adata2.obsm['X_pca'] = integrate_adata[adata2.obs.index].obsm['X_pca']

    return adata1, adata1_next, adata2, cumulative_explained_variance_ratio


def wproj_adata(adata1_young, adata1_old, adata2_young, folder_name, title):
    lambda_list = []
    W_list = []
    x1_young = np.array(adata1_young.obsm['X_pca'])
    x1_old = np.array(adata1_old.obsm['X_pca'])
    x2_young = np.array(adata2_young.obsm['X_pca'])

    # x1_young = x1_young - np.mean(x1_young, axis=0)
    # x1_old = x1_old - np.mean(x1_old, axis=0)
    # x2_young = x2_young - np.mean(x2_young, axis=0)

    w1_young = np.ones((x1_young.shape[0],)) / x1_young.shape[0]
    w1_old = np.ones((x1_old.shape[0],)) / x1_old.shape[0]
    w2_young = np.ones((x2_young.shape[0],)) / x2_young.shape[0]

    M_1y_1o = ot.dist(x1_young, x1_old, p=2)  # Euclidean distance matrix
    M = M_1y_1o
    M_1y_1o = M_1y_1o / M.sum()
    # OT_1y_1o = ot.emd(w1_young, w1_old, np.array(M_1y_1o))  # OT matrix
    W_1y_1o = ot.emd2(w1_young, w1_old, np.array(M_1y_1o))  # Wasserstein distance between young and old

    M_1y_2y = ot.dist(x1_young, x2_young, p=2)  # Euclidean distance matrix
    print(M_1y_2y.sum())
    M_1y_2y = M_1y_2y / M.sum()
    OT_1y_2y = ot.emd(w1_young, w2_young, np.array(M_1y_2y))  # OT matrix
    W_1y_2y = ot.emd2(w1_young, w2_young,
                      np.array(M_1y_2y))  # Wasserstein distance between adata1_young and adata2_young

    M_2y_1o = ot.dist(x2_young, x1_old, p=2)  # Euclidean distance matrix
    M_2y_1o = M_2y_1o / M.sum()
    # OT_2y_1o = ot.emd(w2_young, w1_old, np.array(M_2y_1o))  # OT matrix
    W_2y_1o = ot.emd2(w2_young, w1_old, np.array(M_2y_1o))  # Wasserstein distance between adata2_young and adata1_old

    edge_list_1y_2y = extract_edges_above_threshold(OT_1y_2y, 0)  # extract edges above threshold
    print(x1_young.shape[0], x2_young.shape[0], x1_old.shape[0], len(edge_list_1y_2y))

    lambda_list.append(0.0)
    W_list.append(W_1y_1o)

    lambda_list.append(1.0)
    W_list.append(W_2y_1o)

    W_bc_1o = W_1y_1o

    m = 0  # number of iterations
    while m < 19:
        lambda_ = (m + 1) / 20
        x_bc, w_bc = interpolate_edges(x1_young, x2_young, edge_list_1y_2y, lambda_)

        M_bc_1o = ot.dist(x_bc, x1_old, p=2)
        print(M_bc_1o.sum())
        M_bc_1o = M_bc_1o / M.sum()
        W_bc_1o_new = ot.emd2(w_bc, w1_old, M_bc_1o, numItermax=1000000)
        lambda_list.append(lambda_)
        W_list.append(W_bc_1o_new)

        if W_bc_1o_new > W_bc_1o:
            break

        W_bc_1o = W_bc_1o_new  # update W_bc_1o
        m += 1  # increase m by 1

    n = 1  # number of iterations
    while lambda_ > -2:
        lambda_ = (m + 1) / 20 - n / 100
        if lambda_ <= 0:
            lambda_ = 0
            W_1y_bc = W_1y_2y * lambda_ * abs(lambda_)
            p1 = W_1y_bc / W_1y_1o
            p2 = W_1y_bc / (abs(W_1y_bc) + W_bc_1o)
            break

        x_bc, w_bc = interpolate_edges(x1_young, x2_young, edge_list_1y_2y, lambda_)

        M_bc_1o = ot.dist(x_bc, x1_old, p=2)
        M_bc_1o = M_bc_1o / M.sum()
        W_bc_1o_new = ot.emd2(w_bc, w1_old, M_bc_1o, numItermax=1000000)
        lambda_list.append(lambda_)
        W_list.append(W_bc_1o_new)

        if W_bc_1o_new > W_bc_1o:
            lambda_ = round((m + 1) / 20 - (n - 1) / 100, 2)
            W_1y_bc = W_1y_2y * lambda_ * abs(lambda_)
            p1 = W_1y_bc / W_1y_1o
            p2 = W_1y_bc / (abs(W_1y_bc) + W_bc_1o)
            break

        W_bc_1o = W_bc_1o_new
        n += 1  # increase n by 1

    pl.plot(lambda_list, W_list, 'x')
    pl.xlabel('lambda')
    pl.ylabel('W_bc_A_old')
    pl.savefig(f'{folder_name}/lambda_D_{title}.png')

    print(w_bc.sum())

    return lambda_, p1, p2


def wproj_adata_list(adata_list, m, n):
    """Computes projection weights and distances for each pair of samples.

    Args:
        adata_list (list): A list of AnnData objects containing the data.
        m (int): Number of rows in the output matrices.
        n (int): Number of columns in the output matrices.

    Returns:
        tuple: A tuple containing three 2D lists (matrices) T, P1, and P2.
    """
    T = [[0 for _ in range(m)] for _ in range(n)]
    P1 = [[0 for _ in range(m)] for _ in range(n)]
    P2 = [[0 for _ in range(m)] for _ in range(n)]

    for adata_array in adata_list:
        x1 = adata_array[2]
        x2 = adata_array[3]
        x3 = adata_array[4]

        x1_list = [x1.mean(0).tolist() for _ in range(x1.shape[0])]
        x2_list = [x2.mean(0).tolist() for _ in range(x2.shape[0])]
        x3_list = [x3.mean(0).tolist() for _ in range(x3.shape[0])]

        x1 = np.array(x1) - np.array(x1_list)
        x2 = np.array(x2) - np.array(x2_list)
        x3 = np.array(x3) - np.array(x3_list)

        w1 = np.ones((x1.shape[0],)) / x1.shape[0]
        w2 = np.ones((x2.shape[0],)) / x2.shape[0]
        w3 = np.ones((x3.shape[0],)) / x3.shape[0]

        m = 0
        M = ot.dist(x1, x2, p=2)
        M = M / M.sum()
        G = ot.emd(w1, w2, np.array(M))
        L = ot.emd2(w1, w2, np.array(M))

        n = 1
        edge_list = extract_edges_above_threshold(G, 1 / x1.shape[0] / x2.shape[0] / 2)

        M = ot.dist(x1, x3, p=2)
        M = M / M.sum()
        D_init = ot.emd2(w1, w3, np.array(M), numItermax=1000000) * M.sum() * M.sum()
        D = D_init

        while m < 19:
            t = (m + 1) / 20
            x_new, w_new = interpolate_edges(x1, x2, edge_list, t)

            M = ot.dist(x_new, x3, p=2)
            M = M / M.sum()
            D_new = ot.emd2(w_new, w3, M, numItermax=1000000) * M.sum() * M.sum()

            if D_new > D:
                break

            x = x_new
            w = w_new
            D = D_new
            m += 1

        x = x_new
        w = w_new
        D = D_new
        while t > -2:
            t = (m + 1) / 20 - n / 100
            if t <= 0:
                t = 0
                L_new = L * t * abs(t)
                p1 = L_new / D_init
                p2 = L_new / (abs(L_new) + D)
                break

            x_new, w_new = interpolate_edges(x1, x2, edge_list, t)

            M = ot.dist(x_new, x3, p=2)
            M = M / M.sum()
            D_new = ot.emd2(w_new, w3, M, numItermax=1000000) * M.sum() * M.sum()

            if D_new > D:
                t = round((m + 1) / 20 - (n - 1) / 100, 2)
                L_new = L * t * abs(t)
                p1 = L_new / D_init
                p2 = L_new / (abs(L_new) + D)
                break

            x = x_new
            w = w_new
            D = D_new
            n += 1

        T[adata_array[0]][adata_array[1]] = t
        P1[adata_array[0]][adata_array[1]] = p1
        P2[adata_array[0]][adata_array[1]] = p2

    return T, P1, P2
