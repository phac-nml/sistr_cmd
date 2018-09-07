# -*- coding: utf-8 -*-

import argparse
import logging

from sistr_cmd.util.logger import init_console_logger

from sistr_cmd.sistr_cmd import genome_name_from_fasta_path, sistr_predict, write_cgmlst_profiles, \
    write_cgmlst_results_json, write_novel_alleles
from sistr_cmd.version import __version__


def init_parser():

    parser = argparse.ArgumentParser(prog='sistr_cmd',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=prog_desc)

    parser.add_argument('fastas',
                        metavar='F',
                        nargs='*',
                        help='Input genome FASTA file')
    parser.add_argument('-i',
                        '--input-fasta-genome-name',
                        nargs=2,
                        metavar=('fasta_path', 'genome_name'),
                        action='append',
                        help='fasta file path to genome name pair')
    parser.add_argument('-f',
                        '--output-format',
                        default='json',
                        choices=['json', 'csv', 'pickle'],
                        help='Output format (json, csv, pickle)')
    parser.add_argument('-o',
                        '--output-prediction',
                        help='SISTR serovar prediction output path')
    parser.add_argument('-M',
                        '--more-results',
                        action='count',
                        default=0,
                        help='Output more detailed results (-M) and all antigen search blastn results (-MM)')
    parser.add_argument('-p',
                        '--cgmlst-profiles',
                        help='Output CSV file destination for cgMLST allelic profiles')
    parser.add_argument('-n',
                        '--novel-alleles',
                        help='Output FASTA file destination of novel cgMLST alleles from input genomes')
    parser.add_argument('-a',
                        '--alleles-output',
                        help='Output path of allele sequences and info to JSON')
    parser.add_argument('-T',
                        '--tmp-dir',
                        default='/tmp',
                        help='Base temporary working directory for intermediate analysis files.')
    parser.add_argument('-K',
                        '--keep-tmp',
                        action='store_true',
                        help='Keep temporary analysis files.')
    parser.add_argument('--use-full-cgmlst-db',
                        action='store_true',
                        help='Use the full set of cgMLST alleles which can include highly similar alleles. By default the smaller "centroid" alleles or representative alleles are used for each marker. ')
    parser.add_argument('--no-cgmlst',
                        action='store_true',
                        help='Do not run cgMLST serovar prediction')
    parser.add_argument('-m', '--run-mash',
                        action='store_true',
                        help='Determine Mash MinHash genomic distances to Salmonella genomes with trusted serovar designations. Mash binary must be in accessible via $PATH (e.g. /usr/bin).')
    parser.add_argument('--qc',
                        action='store_true',
                        help='Perform basic QC to provide level of confidence in serovar prediction results.')
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=1,
                        help='Number of parallel threads to run sistr_cmd analysis.')
    parser.add_argument('-v',
                        '--verbose',
                        action='count',
                        default=0,
                        help='Logging verbosity level (-v == show warnings; -vvv == show debug info)')
    parser.add_argument('-V', '--version', action='version', version='%(prog)s {}'.format(__version__))
    return parser


def main():
    parser = init_parser()
    args = parser.parse_args()
    init_console_logger(args.verbose)
    logging.info('Running sistr_cmd {}'.format(__version__))
    input_fastas = args.fastas
    paths_names = args.input_fasta_genome_name
    if len(input_fastas) == 0 and (paths_names is None or len(paths_names) == 0):
        raise Exception('No FASTA files specified!')
    if paths_names is None:
        genome_names = [genome_name_from_fasta_path(x) for x in input_fastas]
    else:
        if len(input_fastas) == 0 and len(paths_names) > 0:
            input_fastas = [x for x,y in paths_names]
            genome_names = [y for x,y in paths_names]
        elif len(input_fastas) > 0 and len(paths_names) > 0:
            tmp = input_fastas
            input_fastas = [x for x,y in paths_names] + tmp
            genome_names = [y for x,y in paths_names] + [genome_name_from_fasta_path(x) for x in tmp]
        else:
            raise Exception('Unhandled fasta input args: input_fastas="{}" | input_fasta_genome_name="{}"'.format(
                input_fastas,
                paths_names))

    tmp_dir = args.tmp_dir
    keep_tmp = args.keep_tmp
    output_format = args.output_format
    output_path = args.output_prediction

    n_threads = args.threads
    if n_threads == 1:
        logging.info('Serial single threaded run mode on %s genomes', len(input_fastas))
        outputs = [sistr_predict(input_fasta, genome_name, tmp_dir, keep_tmp, args) for input_fasta, genome_name in zip(input_fastas, genome_names)]
    else:
        from multiprocessing import Pool
        logging.info('Initializing thread pool with %s threads', n_threads)
        pool = Pool(processes=n_threads)
        logging.info('Running SISTR analysis asynchronously on %s genomes', len(input_fastas))
        res = [pool.apply_async(sistr_predict, (input_fasta, genome_name, tmp_dir, keep_tmp, args)) for input_fasta, genome_name in zip(input_fastas, genome_names)]

        logging.info('Getting SISTR analysis results')
        outputs = [x.get() for x in res]

    prediction_outputs = [x for x,y in outputs]
    cgmlst_results = [y for x,y in outputs]

    if output_path:
        from sistr_cmd.writers import write
        logging.info('Writing results with %s verbosity',
                     args.more_results)
        write(output_path, output_format, prediction_outputs, more_results=args.more_results)
    else:
        import json
        from sistr_cmd.writers import to_dict
        logging.warning('No prediction results output file written! Writing results summary to stdout as JSON')
        exclude_keys_in_output = {'blast_results', 'sseq'}
        if args.more_results >= 2:
            exclude_keys_in_output.remove('blast_results')
            exclude_keys_in_output.remove('sseq')
        elif args.more_results == 1:
            exclude_keys_in_output.remove('sseq')
        outs = [to_dict(x, 0, exclude_keys=exclude_keys_in_output) for x in prediction_outputs]
        print(json.dumps(outs))
    if args.cgmlst_profiles:
        write_cgmlst_profiles(genome_names, cgmlst_results, args.cgmlst_profiles)
        logging.info('cgMLST allelic profiles written to %s', args.cgmlst_profiles)
    if args.alleles_output:
        write_cgmlst_results_json(genome_names, cgmlst_results, args.alleles_output)
        logging.info('JSON of allele data written to %s for %s cgMLST allele results', args.alleles_output, len(cgmlst_results))
    if args.novel_alleles:
        count = write_novel_alleles(cgmlst_results, args.novel_alleles)
        logging.info('Wrote %s alleles to %s', count, args.novel_alleles)