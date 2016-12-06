#!/usr/bin/env python
from collections import defaultdict

import logging
import os
import argparse
import re
from sistr.src.cgmlst import allele_name

from sistr.src.logger import init_console_logger
from sistr.src.parsers import parse_fasta
from sistr.src.cgmlst.extras.centroid_cgmlst_alleles import find_centroid_alleles


def init_arg_parser():
    prog_desc = 'Output centroid alleles for cgMLST markers'
    parser = argparse.ArgumentParser(prog='reduce_to_centroid_alleles',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=prog_desc)

    parser.add_argument('-i', '--input',
                        required=True,
                        help='cgMLST FASTA file with header format ">{marker name}|{allele name}"')
    parser.add_argument('-o', '--output',
                        help='Reduced centroid cgMLST alleles FASTA output path (default: "<input_filename>.centroid.fasta"')
    parser.add_argument('-w', '--word_size',
                        type=int,
                        default=28,
                        help='Number of bp at ends of each allele to group by. Related to blastn word size.')
    parser.add_argument('-t', '--cluster-threshold',
                        type=float,
                        default=0.025,
                        help='Flat cluster generation threshold from hierarchically clustered alleles.')
    parser.add_argument('--threads',
                        type=int,
                        default=1)
    parser.add_argument('-v', '--verbose',
                        action='count',
                        default=0,
                        help='Logging verbosity (-v to log warnings; -vvv to log debug info)')
    return parser


def parse_cgmlst_alleles(cgmlst_fasta):
    """Parse cgMLST alleles from fasta file
    cgMLST FASTA file must have a header format of ">{marker name}|{allele name}"

    Args:
        cgmlst_fasta (str): cgMLST fasta file path

    Returns:
        dict of list: Marker name to list of allele sequences
    """
    out = defaultdict(list)
    for header, seq in parse_fasta(cgmlst_fasta):
        if not '|' in header:
            raise Exception('Unexpected format for cgMLST fasta file header. No "|" (pipe) delimiter present! Header="{}"'.format(header))
        marker_name, allele_name = header.split('|')
        out[marker_name].append(seq)
    return out


def find_all_centroid_alleles(marker_alleles, threads=1, word_size=28, cluster_threshold=0.025):
    out = {}
    if threads == 1:
        logging.info('Running centroid allele finding in serial single-threaded mode')
        for marker_name, alleles in marker_alleles.items():
            centroids = find_centroid_alleles(alleles, bp=word_size, t=cluster_threshold)
            logging.info('Found {} centroid alleles for marker {} ({}% reduction from {})'.format(
                         len(centroids),
                         marker_name,
                         (1.0 - len(centroids) / float(len(alleles))) * 100),
                         len(alleles))
            out[marker_name] = centroids
    else:
        logging.info('Running centroid allele finding in parallel mode with %s threads', threads)
        from multiprocessing import Pool
        pool = Pool(processes=threads)
        jobs = []
        markers = sorted(marker_alleles.keys())
        for marker_name in markers:
            alleles = marker_alleles[marker_name]
            job = pool.apply_async(find_centroid_alleles,
                             (alleles,),
                             {'bp': word_size, 't': cluster_threshold})
            jobs.append(job)
        logging.info('Queued all markers for centroid allele finding; awaiting results')
        for marker_name, job in zip(markers, jobs):
            out[marker_name] = job.get()
        for marker_name, centroids in out.items():
            alleles = marker_alleles[marker_name]
            logging.info('Found {} centroid alleles for marker {} ({}% reduction from {})'.format(
                     len(centroids),
                     marker_name,
                     (1.0 - len(centroids) / float(len(alleles))) * 100,
                     len(alleles)))
    return out


def write_alleles(marker_alleles, output_path):
    with open(output_path, 'w') as fout:
        for marker_name, seqs in marker_alleles.items():
            for seq in seqs:
                aname = allele_name(seq)
                fout.write('>{}|{}\n{}\n'.format(marker_name, aname, seq))


def run_allele_reduction(input_fasta, output_path, threads=1, word_size=28, cluster_threshold=0.025):
    logging.info('Parsing alleles from input fasta %s', input_fasta)
    marker_alleles = parse_cgmlst_alleles(input_fasta)
    n_total_alleles = sum([len(v) for k,v in marker_alleles.items()])
    logging.info('Parsed %s alleles for %s markers from %s',
                 n_total_alleles,
                 len(marker_alleles),
                 input_fasta)

    logging.info('Finding centroid alleles for each marker')
    marker_centroids = find_all_centroid_alleles(marker_alleles,
                                                 threads=threads,
                                                 word_size=word_size,
                                                 cluster_threshold=cluster_threshold)
    n_centroids = sum([len(v) for k,v in marker_centroids.items()])
    logging.info('Found {} centroid alleles ({}% of input)'.format(
                 n_centroids,
                 (n_centroids / float(n_total_alleles)) * 100))

    logging.info('Outputting centroid alleles to output path "%s"', output_path)
    write_alleles(marker_centroids, output_path)
    logging.info('Centroid alleles written to "%s"', output_path)
    logging.info('Done!')


def main():
    parser = init_arg_parser()
    args = parser.parse_args()
    init_console_logger(args.verbose)
    input_fasta = args.input
    if not os.path.exists(input_fasta):
        err_msg = 'Input fasta file does not exist at {}. cgMLST fasta file required!'.format(input_fasta)
        logging.error(err_msg)
        raise Exception(err_msg)
    output_path = args.output
    if output_path is None:
        input_fasta_dir = os.path.dirname(input_fasta)
        input_filename = os.path.basename(input_fasta)
        output_filename = re.sub(r'(\.fasta$)|(\.\w{1,}$)', '', input_filename) + '-centroid.fasta'
        output_path = os.path.join(input_fasta_dir, output_filename)
        logging.info('No output path specified. Using "%s"', output_path)
    if os.path.exists(output_path):
        err_msg = 'File already exists at the output path "{}"! Specify a different output path.'.format(output_path)
        logging.error(err_msg)
        raise Exception(err_msg)

    run_allele_reduction(input_fasta, output_path,
                         threads=args.threads,
                         word_size=args.word_size,
                         cluster_threshold=args.cluster_threshold)


if __name__ == '__main__':
    main()