#!/usr/bin/env python
import logging
import os
import argparse
import numpy as np

from sistr.src.logger import init_console_logger
from sistr.src.cgmlst.extras.hclust_cutree import profiles_to_np_array, nr_profiles, dist_matrix_hamming, complete_linkage, cutree, expand_clusters_dataframe


def init_arg_parser():
    prog_desc = 'Hierarchical clustering flat clusters from cgMLST profiles'
    parser = argparse.ArgumentParser(prog='hc_clusters_from_cgmlst_profiles',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=prog_desc)

    parser.add_argument('-i', '--input',
                        required=True,
                        help='cgMLST allelic profiles. Column names == marker names; row names == genome names. First row is column header. First column is genome names.')
    parser.add_argument('-o', '--output',
                        required=True,
                        help='Output path for CSV of hierarchical clustering flat clusters table.')
    parser.add_argument('-b', '--begin-threshold',
                        type=float,
                        default=0.0,
                        help='Flat cluster distance threshold range start (default: 0.0)')
    parser.add_argument('-e', '--end-threshold',
                        type=float,
                        default=1.0,
                        help='Flat cluster distance threshold range end inclusive (default: 1.0)')
    parser.add_argument('-s', '--step-threshold',
                        type=float,
                        default=0.01,
                        help='Flat cluster distance threshold range step value (default: 0.01).')
    parser.add_argument('-v', '--verbose',
                        action='count',
                        default=0,
                        help='Logging verbosity (-v to log warnings; -vvv to log debug info)')
    return parser


def main():
    parser = init_arg_parser()
    args = parser.parse_args()
    init_console_logger(args.verbose)
    input_profiles = args.input
    if not os.path.exists(input_profiles):
        err_msg = 'Input cgMLST profiles file does not exist at {}. cgMLST profiles file required!'.format(input_profiles)
        logging.error(err_msg)
        raise Exception(err_msg)
    output_path = args.output
    if os.path.exists(output_path):
        err_msg = 'File already exists at the output path "{}"! Specify a different output path.'.format(output_path)
        logging.error(err_msg)
        raise Exception(err_msg)

    begin_threshold = args.begin_threshold
    end_threshold = args.end_threshold
    step_threshold = args.step_threshold
    assert begin_threshold >= 0.0, "-b/--begin-threshold must be positive number!"
    assert begin_threshold <= 1.0, "-b/--begin-threshold must be less than or equal to 1.0"
    assert begin_threshold < end_threshold, "-b/--begin-threshold must be less than -e/--end-threshold"
    assert end_threshold >= 0.0, "-e/--end-threshold must be positive number!"
    assert end_threshold <= 1.0, "-e/--end-threshold must be less than or equal to 1.0"
    assert step_threshold <= 1.0, "-s/--step-threshold must be less than or equal to 1.0"
    assert step_threshold >= 0.0, "-s/--step-threshold must be positive number!"

    logging.info('Reading profiles from %s', input_profiles)
    profiles_matrix, genomes, markers = profiles_to_np_array(input_profiles)
    logging.info('Profiles matrix shape: %s', profiles_matrix.shape)
    logging.debug('Genomes: %s', genomes)
    logging.debug('Markers: %s', markers)
    logging.info('Finding non-redundant profiles. Grouping genomes by distinct profiles.')
    nr_profiles_matrix, genome_groups = nr_profiles(profiles_matrix, genomes)
    logging.info('Non-redundant profiles matrix shape: %s', nr_profiles_matrix.shape)
    logging.debug('Genome groups: %s', genome_groups)
    logging.info('Computing Hamming distance matrix from profiles')
    dm = dist_matrix_hamming(nr_profiles_matrix)
    logging.info('Complete linkage of Hamming distance matrix')
    hc = complete_linkage(dm)
    thresholds = np.arange(begin_threshold, end_threshold + step_threshold, step_threshold)
    logging.info('Generating flat clusters from %s to %s with step %s',
                 begin_threshold,
                 end_threshold,
                 step_threshold)
    logging.debug('Thresholds: %s', thresholds)
    df_clusters = cutree(hc, thresholds)
    logging.info('Flat clusters generated for non-redundant profiles. Expanding to all profiles.')
    logging.debug('df_clusters: %s', df_clusters)
    df_clusters = expand_clusters_dataframe(df_clusters, genome_groups)
    df_clusters.to_csv(output_path)
    logging.info('HC flat clusters written to "%s"', output_path)
    logging.info('Done!')


if __name__ == '__main__':
    main()
