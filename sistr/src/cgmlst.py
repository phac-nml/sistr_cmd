import os
from pkg_resources import resource_filename, Requirement
import numpy as np
import pandas as pd

# base_dir = os.path.join(os.path.dirname(__file__), '../')
# here = lambda x: os.path.abspath(os.path.join(base_dir, x))

CGMLST_FASTA_PATH = resource_filename('sistr', 'data/cgmlst/cgmlst330.fasta')
CGMLST_PROFILES_PATH = resource_filename('sistr', 'data/cgmlst/cgmlst330-profiles-nr.csv')

GENOMES_TO_SEROVAR_PATH = resource_filename('sistr', 'data/cgmlst/genomes-to-serovar.txt')

# CGMLST_FASTA_PATH = 'data/cgmlst/cgmlst330.fasta'
# CGMLST_PROFILES_PATH = 'data/cgmlst/cgmlst330-profiles-nr.csv'
#
# GENOMES_TO_SEROVAR_PATH = 'data/cgmlst/genomes-to-serovar.txt'

def cgmlst_profiles():
    return pd.read_csv(CGMLST_PROFILES_PATH, index_col=0)

def genomes_to_serovar():
    rtn = {}
    with open(GENOMES_TO_SEROVAR_PATH) as f:
        for l in f:
            l = l.strip()
            genome, serovar = l.split('\t')
            rtn[genome] = serovar
    return rtn


def perfect_matches_to_marker_results(df):
    """Perfect BLAST matches to marker results dict

    Parse perfect BLAST matches to marker results dict.
    The `qseqid` column must contain cgMLST330 query IDs with `{marker name}|{allele number}` format.

    Args:
        df (pandas.DataFrame): DataFrame of perfect BLAST matches

    Returns:
        dict: cgMLST330 marker names to matching allele numbers
    """
    marker_results = {}
    for x in df['qseqid']:
        marker, result = x.split('|')
        marker_results[marker] = int(result)
    return marker_results


def find_closest_related_genome(marker_results, df_genome_profiles):
    """

    Args:
        df_genome_profiles (pandas.DataFrame):

    Returns:
        (dict, list): Most closely related Genome and list of other related Genomes_ in order of relatedness
    """
    marker_names = df_genome_profiles.columns
    n_markers = len(marker_names)

    profile = [marker_results[marker_name] if marker_name in marker_results else None for marker_name in marker_names]
    genome_profile = np.array(profile, dtype=np.float64)
    profiles_matrix = np.array(df_genome_profiles, dtype=np.float64)
    genome_profile_similarity_counts = np.apply_along_axis(lambda x: (x == genome_profile).sum(), 1, profiles_matrix)

    df_relatives = pd.DataFrame()
    df_relatives['matching'] = genome_profile_similarity_counts
    df_relatives['distance'] = 1.0 - (df_relatives['matching'] / float(n_markers))
    df_relatives.index = df_genome_profiles.index
    df_relatives.sort_values(by='distance', inplace=True)
    return df_relatives

