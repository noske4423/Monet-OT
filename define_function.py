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


def adjust_and_save_plot(adata, title, method, min_value_1, max_value_1, min_value_2, max_value_2, xlabel, ylabel, color):
    """Adjusts and saves a plot.

    Args:
        adata (AnnData): The input AnnData object.
        title (str): The title of the plot.
        method (str): The method used to generate the plot.
        xlabel (str): The label of the x-axis.
        ylabel (str): The label of the y-axis.
        color (str): The color of the plot.
    """
    if method == 'pca':
        fig = sc.pl.pca(adata, color=color, return_fig=True)
        ax = fig.gca()
        ax.set_xlim(min_value_1 - (max_value_1 - min_value_1) * 0.05,
                    max_value_1 + (max_value_1 - min_value_1) * 0.05)
        ax.set_ylim(min_value_2 - (max_value_2 - min_value_2) * 0.05,
                    max_value_2 + (max_value_2 - min_value_2) * 0.05)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        pl.savefig(title + '.png')

    elif method == 'umap':
        fig = sc.pl.umap(adata, color=color, return_fig=True)
        ax = fig.gca()
        ax.set_xlim(min_value_1 - (max_value_1 - min_value_1) * 0.05,
                    max_value_1 + (max_value_1 - min_value_1) * 0.05)
        ax.set_ylim(min_value_2 - (max_value_2 - min_value_2) * 0.05,
                    max_value_2 + (max_value_2 - min_value_2) * 0.05)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        pl.savefig(title + '.png')


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


def wproj(adata_list, m, n):
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
