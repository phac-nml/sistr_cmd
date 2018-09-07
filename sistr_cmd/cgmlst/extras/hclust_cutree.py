from fastcluster import linkage
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import fcluster
import pandas as pd
import numpy as np


def profiles_to_np_array(profiles_csv_path):
    """

    """
    df = pd.read_csv(profiles_csv_path, index_col=0)
    arr = np.array(df, dtype=np.float64)
    genomes = df.index
    markers = df.columns
    return arr, genomes, markers


def nr_profiles(arr, genomes):
    """
    Get a condensed cgMLST pairwise distance matrix for specified Genomes_
    where condensed means redundant cgMLST profiles are only represented once in the distance matrix.

    Args:
        user_name (list): List of Genome_ names to retrieve condensed distance matrix for

    Returns:
        (numpy.array, list): tuple of condensed cgMLST distance matrix and list of grouped Genomes_
    """
    gs_collapse = []
    genome_idx_dict = {}
    indices = []
    patt_dict = {}
    for i, g in enumerate(genomes):
        p = arr[i, :].tostring()
        if p in patt_dict:
            parent = patt_dict[p]
            idx = genome_idx_dict[parent]
            gs_collapse[idx].append(g)
        else:
            indices.append(i)
            patt_dict[p] = g
            genome_idx_dict[g] = len(gs_collapse)
            gs_collapse.append([g])
    return arr[indices, :], gs_collapse


def dist_matrix_hamming(arr):
    return pdist(arr, metric='hamming')


def complete_linkage(dm):
    """
    Perform complete linkage hierarchical clustering on a distance matrix.

    Args:
        dm (numpy.array): Distance matrix

    Returns:
        (object): fastcluster complete linkage hierarchical clustering object
    """
    return linkage(dm, 'complete')


def cutree(Z, thresholds):
    out = {}
    for t in thresholds:
        out[t] = fcluster(Z, t, criterion='distance')
    return pd.DataFrame(out)

def expand_clusters_dataframe(df_clusters, genome_groups):
    lens_genome_groups = [len(xs) for xs in genome_groups]
    idxs = np.repeat(df_clusters.index.values, lens_genome_groups)
    df_cl_exp = df_clusters.reindex(idxs, method='ffill')
    df_cl_exp.index = [g for gs in genome_groups for g in gs]
    return df_cl_exp