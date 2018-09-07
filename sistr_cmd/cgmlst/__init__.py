# -*- coding: utf-8 -*-

from typing import Dict, Optional, Tuple
import logging
from collections import Counter

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from sistr_cmd.blast_wrapper.reader import BlastReader
from sistr_cmd.blast_wrapper.runner import BlastRunner
from sistr_cmd.cgmlst import hash_allele_seq, allele_result_dict
from sistr_cmd.cgmlst.allele_processing import alleles_to_retrieve, allele_result_dict, get_allele_sequences
from sistr_cmd.cgmlst.const import CGMLST_CENTROID_FASTA_PATH, CGMLST_FULL_FASTA_PATH
from sistr_cmd.cgmlst.msa import msa_ref_vs_novel, number_gapped_ungapped, MSA_GAP_PROP_THRESHOLD
from sistr_cmd.cgmlst.results import CgMLSTPrediction

from sistr_cmd.cgmlst.util import hash_allele_seq, ref_cgmlst_profiles
from sistr_cmd.serovar_prediction.constants import CGMLST_SUBSPECIATION_DISTANCE_THRESHOLD, genomes_to_subspecies, \
    genomes_to_serovar
from sistr_cmd.util import first_row


def matches_to_marker_results(df: pd.DataFrame) -> Dict[str, Dict]:
    """Perfect BLAST matches to marker results dict

    Parse perfect BLAST matches to marker results dict.


    Args:
        df (pandas.DataFrame): DataFrame of perfect BLAST matches

    Returns:
        dict: cgMLST330 marker names to matching allele numbers
    """
    assert isinstance(df, pd.DataFrame)
    from collections import defaultdict
    marker_blast_results = defaultdict(list)
    for idx, row in df.iterrows():
        marker = row['marker']
        marker_blast_results[marker].append(row)

    marker_results = {}
    for marker, blast_results in marker_blast_results.items():
        if len(blast_results) > 1:
            logging.debug(f'Multiple potential cgMLST allele matches (n={len(blast_results)}) found for marker '
                          f'{marker}. Selecting match on longest contig.')
            _process_multiple_blast_results_for_marker(blast_results)
        elif len(blast_results) == 1:
            marker_blast_results[marker] = _process_single_blast_result_for_marker(blast_results)
        else:
            err_msg = 'Empty list of matches for marker {}'.format(marker)
            logging.error(err_msg)
            raise Exception(err_msg)
    return marker_results


def _process_single_blast_result_for_marker(blast_results):
    row = blast_results[0]
    seq = row['sseq']
    if '-' in seq:
        logging.warning('Gaps found in allele. Removing gaps. %s', row)
        seq = seq.replace('-', '').upper()
    allele_match_name = hash_allele_seq(seq)
    return allele_result_dict(allele_match_name, seq, row.to_dict())


def _process_multiple_blast_results_for_marker(blast_results):
    df_marker_blast_results = pd.DataFrame(blast_results)
    df_marker_blast_results.sort_values('slen', ascending=False, inplace=True)
    blast_result_row = first_row(df_marker_blast_results)
    allele_match_name = blast_result_row['allele_name']
    subject_length = blast_result_row['slen']
    logging.debug('Selecting allele %s from contig with length %s', allele_match_name, subject_length)
    seq = blast_result_row['sseq']
    if '-' in seq:
        logging.info(
            f'Gaps found in allele {blast_result_row["qseqid"]}. Removing gaps. BLAST+ result: {blast_result_row.to_dict()}')
        seq = seq.replace('-', '').upper()
    allele_hash = hash_allele_seq(seq)
    if allele_hash != allele_match_name:
        logging.info(
            f'New allele name assigned for {blast_result_row["qseqid"]}. NEW={allele_hash}, OLD={allele_match_name}')
    return allele_result_dict(allele_hash, seq, blast_result_row.to_dict())


def find_relatives(marker_results: Dict[str, Optional[float]], df_genome_profiles: pd.DataFrame) -> pd.DataFrame:
    """Find the closest related reference genomes given some cgMLST typing results

    Args:
        marker_results: Dict of marker name to result
        df_genome_profiles:

    Returns:
        DataFrame of distances and number of alleles matching reference genomes
    """
    marker_names = df_genome_profiles.columns
    n_markers = len(marker_names)
    ref_profiles = df_genome_profiles.values  # type: np.ndarray
    profile = [marker_results[marker_name] if marker_name in marker_results else None for marker_name in marker_names]
    genome_profile = np.array([profile], dtype=np.float64)
    distances = cdist(genome_profile, ref_profiles, metric='hamming')

    df_relatives = pd.DataFrame(dict(matching=(np.round(1.0 - distances) * n_markers),
                                     distance=distances),
                                index=df_genome_profiles.index)
    df_relatives.sort_values(by='distance', inplace=True)
    return df_relatives


def cgmlst_subspecies_call(df_relatives) -> Optional[Tuple[str, float, dict]]:
    """Call Salmonella subspecies based on cgMLST results

    This method attempts to find the majority subspecies type within curated
    public genomes above a cgMLST allelic profile distance threshold.

    Note:
        ``CGMLST_SUBSPECIATION_DISTANCE_THRESHOLD`` is the cgMLST distance
        threshold used to determine the subspecies by cgMLST. It is set at a
        distance of 0.9 which translates to a cgMLST allelic similarity of 10%.
        A threshold of 0.9 is generous and reasonable given the congruence
        between subspecies designations and 10% cgMLST clusters by Adjusted
        Rand (~0.850) and Adjusted Wallace metrics (~0.850 both ways).

    Args:
        df_relatives (pandas.DataFrame): Table of genomes related by cgMLST to input genome

    Returns:
        None: if no curated public genomes found to have a cgMLST profile similarity of 10% or greater
        (string, float, dict): most common subspecies, closest related public genome distance, subspecies frequencies
    """

    closest_distance = df_relatives['distance'].min()

    if closest_distance > CGMLST_SUBSPECIATION_DISTANCE_THRESHOLD:
        logging.warning('Min cgMLST distance (%s) above subspeciation distance threshold (%s)',
                        closest_distance,
                        CGMLST_SUBSPECIATION_DISTANCE_THRESHOLD)
        return None
    else:
        df_relatives = df_relatives.loc[df_relatives.distance <= CGMLST_SUBSPECIATION_DISTANCE_THRESHOLD, :]
        df_relatives = df_relatives.sort_values('distance', ascending=True)
        logging.debug('df_relatives by cgmlst %s', df_relatives.head())
        genome_spp = genomes_to_subspecies()
        subspecies_below_threshold = [genome_spp[member_genome] if member_genome in genome_spp else None for
                                      member_genome in df_relatives.index]
        subspecies_below_threshold = filter(None, subspecies_below_threshold)
        subspecies_counter = Counter(subspecies_below_threshold)
        logging.debug('Subspecies counter: %s', subspecies_counter)
        return subspecies_counter.most_common(1)[0][0], closest_distance, dict(subspecies_counter)





def marker_result(all_marker_results: Dict[str, Dict]) -> Dict[str, Optional[int]]:
    out = {}
    for marker, res in all_marker_results.items():
        if res and 'name' in res:
            out[marker] = int(res['name'])
        else:
            out[marker] = None
    return out


def get_cgmlst_marker_results(df_cgmlst_profiles, marker_match_results, retrieved_marker_alleles):
    out = marker_match_results.copy()
    for marker, res in retrieved_marker_alleles.items():
        out[marker] = res
    for marker in df_cgmlst_profiles.columns:
        if marker not in out:
            out[marker] = {'blast_result': None,
                           'name': None,
                           'seq': None, }
    return out


def cgmlst_sequence_type(all_marker_results):
    cgmlst_st = None
    cgmlst_markers_sorted = sorted(all_marker_results.keys())
    cgmlst_allele_names = []
    marker = None
    for marker in cgmlst_markers_sorted:
        try:
            aname = all_marker_results[marker]['name']
            if aname:
                cgmlst_allele_names.append(str(aname))
            else:
                break
        except:
            break
    if len(cgmlst_allele_names) == len(cgmlst_markers_sorted):
        cgmlst_st = hash_allele_seq('-'.join(cgmlst_allele_names))
        logging.info('cgMLST330 Sequence Type=%s', cgmlst_st)
    else:
        logging.warning('Could not compute cgMLST330 Sequence Type due to missing data (marker %s)', marker)
    return cgmlst_st


def _add_ref_genome_serovar_info(df_relatives: pd.DataFrame) -> pd.DataFrame:
    """Add reference genome serovar information to cgMLST results DataFrame under `serovar` column.

    Args:
        df_relatives: cgMLST results

    Returns:
        cgMLST results with `serovar` column added containing reference genome serovar designations
    """
    genome_serovar = genomes_to_serovar()
    df_relatives['serovar'] = [genome_serovar[genome] for genome in df_relatives.index]
    return df_relatives


def top_serovar(df_relatives: pd.DataFrame, threshold: float = 1.0) -> Tuple[float, int, Optional[str], Optional[str]]:
    """Get the serovar for the most related reference genome by cgMLST distance.

    Args:
        df_relatives: Related genomes by cgMLST; sorted by distance in ascending order
        threshold: Distance threshold for relatedness

    Returns:
        Tuple of cgMLST distance, number of matching cgMLST alleles, closest matching reference genome name, serovar of closest matching reference genome
    """
    serovar = None
    matching_genome = None
    matching_alleles = 0
    distance = threshold
    for idx, row in df_relatives.iterrows():
        distance = row['distance']
        matching_alleles = row['matching']
        serovar = row['serovar'] if distance <= threshold else None
        matching_genome = idx if distance <= threshold else None
        logging.info(f'Top serovar by cgMLST profile matching: "{serovar}" with {matching_alleles} matching alleles, '
                     f'distance={distance:.1%}')
        break
    return distance, matching_alleles, matching_genome, serovar


def subspeciation(df_relatives: pd.DataFrame) -> Optional[str]:
    """Subspecies prediction by cgMLST relatedness to reference genomes.

    Args:
        df_relatives: Related genomes by cgMLST

    Returns:
        Most likely Salmonella subspecies by cgMLST relatedness to reference genomes.
    """
    spp = None
    subspeciation_tuple = cgmlst_subspecies_call(df_relatives)
    if subspeciation_tuple:
        spp, distance, spp_counter = subspeciation_tuple
        logging.info(f'Top subspecies by cgMLST is "{spp}" (min dist={distance}, Counter={spp_counter})')
    else:
        logging.warning('Subspeciation by cgMLST was not possible!')
    return spp
