from sistr.src.cgmlst.extras.hclust_cutree import profiles_to_np_array, dist_matrix_hamming, complete_linkage, cutree, nr_profiles, expand_clusters_dataframe

from sistr.src.cgmlst import CGMLST_PROFILES_PATH


def test_read_profiles():
    M, genomes, markers = profiles_to_np_array(CGMLST_PROFILES_PATH)
    print(M.shape)
    print(genomes)
    print(markers)
    assert M.shape == (len(genomes), len(markers))


def test_cutree():
    print(f"Reading cgMLST profiles file {CGMLST_PROFILES_PATH} ...")
    M, genomes, markers = profiles_to_np_array(CGMLST_PROFILES_PATH)
    print(M.shape)
    M_nr, genome_groups = nr_profiles(M, genomes)
    print(M_nr.shape)
    dm = dist_matrix_hamming(M_nr)
    print(dm.shape)
    Z = complete_linkage(dm)
    print(Z.shape)
    df_clusters = cutree(Z, [0, 0.1, 0.2, 0.3, 0.4, 0.5])
    print(df_clusters)
    df_clusters = expand_clusters_dataframe(df_clusters, genome_groups)
    print(df_clusters)
    # check that each genome group of size >= 2 is part of the same 0.0 dist cluster
    for gs in genome_groups:
        gs=list(set(gs))
        if len(gs) == 1:
            continue
        df = df_clusters.loc[gs]
        assert df[0.5].unique().size == 1
