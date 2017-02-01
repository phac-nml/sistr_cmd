import zlib
from collections import defaultdict
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import fcluster, linkage


NT_TO_INT = {'A':1,'C':2,'G':3,'T':4}
INT_TO_NT = {1:'A',2:'C',3:'G',4:'T'}


def group_alleles_by_size(alleles):
    allele_size_seqs = defaultdict(list)
    for allele in alleles:
        allele_size_seqs[len(allele)].append(allele)
    return allele_size_seqs


def seq_int_arr(seqs):
    """Convert list of ACGT strings to matix of 1-4 ints

    Args:
        seqs (list of str): nucleotide sequences with only 'ACGT' characters

    Returns:
        numpy.array of int: matrix of integers from 1 to 4 inclusive representing A, C, G, and T
        str: nucleotide sequence string
    """
    return np.array([[NT_TO_INT[c] for c in x.upper()] for x in seqs])


def group_alleles_by_start_end_Xbp(arr, bp=28):
    """Group alleles by matching ends

    Args:
        arr (numpy.array): 2D int matrix of alleles
        bp (int): length of ends to group by

    Returns:
        dict of lists: key of start + end strings to list of indices of alleles with matching ends
    """
    starts = arr[:,0:bp]
    ends = arr[:,-bp:]
    starts_ends_idxs = defaultdict(list)
    l, seq_len = arr.shape
    for i in range(l):
        start_i = starts[i]
        end_i = ends[i]
        start_i_str = ''.join([str(x) for x in start_i])
        end_i_str = ''.join([str(x) for x in end_i])
        starts_ends_idxs[start_i_str + end_i_str].append(i)
    return starts_ends_idxs


def allele_clusters(dists, t=0.025):
    """Flat clusters from distance matrix

    Args:
        dists (numpy.array): pdist distance matrix
        t (float): fcluster (tree cutting) distance threshold

    Returns:
        dict of lists: cluster number to list of indices of distances in cluster
    """
    clusters = fcluster(linkage(dists), 0.025, criterion='distance')
    cluster_idx = defaultdict(list)
    for idx, cl in enumerate(clusters):
        cluster_idx[cl].append(idx)
    return cluster_idx


def dm_subset(dm_sq, idxs):
    """Get subset of distance matrix given list of indices

    Args:
        dm_sq (numpy.array): squareform distance matrix from pdist
        idxs (list of int): list of indices
    Returns:
        numpy.array: subset of `dm_sq` with `shape == (len(idxs), len(idxs))`
    """
    return dm_sq[idxs][:,idxs]


def min_row_dist_sum_idx(dists):
    """Find the index of the row with the minimum row distance sum 

    This should return the index of the row index with the least distance overall 
    to all other rows. 

    Args:
        dists (np.array): must be square distance matrix

    Returns:
        int: index of row with min dist row sum
    """
    row_sums = np.apply_along_axis(arr=dists, axis=0, func1d=np.sum)
    return row_sums.argmin()


def seq_int_arr_to_nt(arr):
    """Convert 1-4 int array to ACGT string

    Args:
        arr (numpy.array of int): array of integers from 1 to 4 inclusive representing A, C, G, and T

    Returns:
        str: nucleotide sequence string
    """
    return ''.join([INT_TO_NT[x] for x in arr])


def find_centroid_alleles(alleles, bp=28, t=0.025):
    """Reduce list of alleles to set of centroid alleles based on size grouping, ends matching and hierarchical clustering

    Workflow for finding centroid alleles:

    - grouping by size (e.g. 100bp, 101bp, 103bp, etc)
    - then grouped by `bp` nucleotides at ends matching
    - size and ends grouped alleles hierarchically clustered (Hamming distance, complete linkage)
    - tree cutting at threshold `t`
    - select allele with minimum distance to other alleles in cluster as centroid

    Args:
        alleles (iterable): collection of allele nucleotide sequences
        bp (int): number of bp matching at allele ends for size grouping (default=28 due to default blastn megablast word size)
        t (float): cluster generation (tree cutting) distance threshold for size grouped alleles

    Returns:
        set of str: centroid alleles
    """
    centroid_alleles = set()
    len_allele = group_alleles_by_size(alleles)
    for length, seqs in len_allele.items():
        # if only one alelle of a particular size, add as centroid, move onto next size group
        if len(seqs) == 1:
            centroid_alleles.add(seqs[0])
            continue
        # convert allele nucleotide sequences to integer matrix
        seq_arr = seq_int_arr(seqs)
        # group alleles by matching ends
        starts_ends_idxs = group_alleles_by_start_end_Xbp(seq_arr, bp=bp)
        for k, idxs in starts_ends_idxs.items():
            # if only one allele for a particular matching ends group, then add as centroid and move onto next ends group
            if len(idxs) == 1:
                centroid_alleles.add(seqs[idxs[0]])
                continue
            # fetch subset of int allele sequences for a matching ends group
            seq_arr_subset = seq_arr[idxs]
            # Hamming distances between alleles
            dists = pdist(seq_arr_subset, 'hamming')
            # create flat clusters (tree cut) at t threshold
            cl = allele_clusters(dists, t=t)
            # for each allele cluster
            dm_sq = squareform(dists)
            for cl_key, cl_idxs in cl.items():
                # if only 1 or 2 alleles in cluster then return first
                if len(cl_idxs) == 1 or len(cl_idxs) == 2:
                    # get first cluster index and get nt seq for that index
                    centroid_alleles.add(seq_int_arr_to_nt(seq_arr_subset[cl_idxs[0]]))
                    continue
                # else find allele with min distances to all other alleles in cluster
                dm_sub = dm_subset(dm_sq, cl_idxs)
                min_idx = min_row_dist_sum_idx(dm_sub)
                # add nucleotide seq for cluster centroid allele to centroids set
                centroid_alleles.add(seq_int_arr_to_nt(seq_arr_subset[min_idx]))
            #end for cl_key, cl_idxs in cl.iteritems():
        #end for k, idxs in starts_ends_idxs.iteritems():
    #end for length, seqs in alleles.iteritems():
    return centroid_alleles
