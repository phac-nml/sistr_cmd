#!/usr/bin/env python
import argparse
import logging
import os

from sistr.src.blast_wrapper import BlastRunner, BlastReader
from sistr.src.logger import init_console_logger
from sistr.src.serovar_prediction import SerovarPredictor, serovar_from_cgmlst_and_antigen_blast


def init_parser():
    prog_desc = '''
SISTR (Salmonella In Silico Typing Resource) Command-line Tool
==============================================================
Serovar predictions from whole-genome sequence assemblies by determination of antigen gene and cgMLST gene alleles using BLAST.

If you find this program useful in your research, please cite as:

The Salmonella In Silico Typing Resource (SISTR): an open web-accessible tool for rapidly typing and subtyping draft Salmonella genome assemblies.
Catherine Yoshida, Peter Kruczkiewicz, Chad R. Laing, Erika J. Lingohr, Victor P.J. Gannon, John H.E. Nash, Eduardo N. Taboada.
PLoS ONE 11(1): e0147101. doi: 10.1371/journal.pone.0147101
'''

    parser = argparse.ArgumentParser(prog='sistr',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=prog_desc)

    parser.add_argument('fastas',
                        metavar='F',
                        nargs='+',
                        help='Input genome FASTA file')
    parser.add_argument('-f',
                        '--output-format',
                        default='json',
                        help='Output format (json, csv, pickle)')
    parser.add_argument('-o',
                        '--output-dest',
                        help='Output')
    parser.add_argument('-T',
                        '--tmp-dir',
                        default='tmp',
                        help='Base temporary working directory for intermediate analysis files.')
    parser.add_argument('-K',
                        '--keep-tmp',
                        action='store_true',
                        help='Keep temporary analysis files.')
    parser.add_argument('--no-cgmlst',
                        action='store_true',
                        help='Do not run cgMLST serovar prediction')
    parser.add_argument('-m', '--run-mash',
                        action='store_true',
                        help='''Determine Mash MinHash genomic distances to Salmonella genomes with trusted serovar designations.
                        Mash binary must be in accessible via $PATH (e.g. /usr/bin).
                        ''')
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=1,
                        help='Number of parallel threads to run sistr_cmd analysis.')
    parser.add_argument('-v',
                        '--verbose',
                        action='count',
                        default=2,
                        help='Logging verbosity level (-v == show warnings; -vvv == show debug info)')
    return parser


def run_cgmlst(prediction, blast_runner):
    from sistr.src.cgmlst import\
        cgmlst_profiles, \
        CGMLST_FASTA_PATH, \
        genomes_to_serovar, \
        find_closest_related_genome, \
        perfect_matches_to_marker_results

    df_cgmlst_profiles = cgmlst_profiles()

    logging.debug('{} distinct cgMLST330 profiles'.format(df_cgmlst_profiles.shape[0]))

    logging.info('Running BLAST on serovar predictive cgMLST330 alleles')
    blast_outfile = blast_runner.blast_against_query(CGMLST_FASTA_PATH)
    logging.info('Reading BLAST output file "{}"'.format(blast_outfile))
    blast_reader = BlastReader(blast_outfile)
    if blast_reader.df is None:
        logging.error('No cgMLST330 alleles found!')
        return
    logging.info('Found {} cgMLST330 allele BLAST results'.format(blast_reader.df.shape[0]))
    df_perfect_cgmlst_results = blast_reader.perfect_matches()
    if df_perfect_cgmlst_results is None:
        logging.warning('No cgMLST330 allele matches found for 330 markers!')
        return
    logging.info('Found {} perfect matches to cgMLST330 alleles'.format(df_perfect_cgmlst_results.shape[0]))
    cgmlst_results = perfect_matches_to_marker_results(df_perfect_cgmlst_results)
    logging.info('Calculating number of matching alleles to serovar predictive cgMLST330 profiles')
    df_relatives = find_closest_related_genome(cgmlst_results, df_cgmlst_profiles)
    genome_serovar_dict = genomes_to_serovar()
    df_relatives['serovar'] = [genome_serovar_dict[genome] for genome in df_relatives.index]
    logging.debug('Top 5 serovar predictive cgMLST profiles:\n{}'.format(df_relatives.head()))

    cgmlst_serovar = None
    cgmlst_matching_genome = None
    cgmlst_matching_alleles = 0
    cgmlst_distance = 1.0
    for idx, row in df_relatives.iterrows():
        cgmlst_distance = row['distance']
        cgmlst_matching_alleles = row['matching']
        cgmlst_serovar = row['serovar'] if cgmlst_distance <= 1.0 else None
        cgmlst_matching_genome = idx if cgmlst_distance <= 1.0 else None
        logging.info('Top serovar by cgMLST profile matching: "{}" with {} matching alleles, distance={:.1%}'.format(
        cgmlst_serovar,
        cgmlst_matching_alleles,
        cgmlst_distance
    ))
        break
    prediction.cgmlst_distance = cgmlst_distance
    prediction.cgmlst_genome_match = cgmlst_matching_genome
    prediction.serovar_cgmlst = cgmlst_serovar
    prediction.cgmlst_matching_alleles = cgmlst_matching_alleles


def run_mash(input_fasta, prediction):
    from sistr.src.mash import mash_dist_trusted, mash_output_to_pandas_df

    mash_out = mash_dist_trusted(input_fasta)
    df_mash = mash_output_to_pandas_df(mash_out)
    logging.debug('Mash top 5 results:\n{}'.format(df_mash[['ref', 'dist', 'n_match', 'serovar']].head()))
    df_mash_top_5 = df_mash[['ref', 'dist', 'n_match', 'serovar']].head(n=5)
    prediction.__dict__['mash_top_5'] = df_mash_top_5.to_dict()
    for idx, row in df_mash_top_5.iterrows():
        mash_genome = row['ref']
        mash_serovar = row['serovar']
        mash_distance = row['dist']
        mash_match = row['n_match']

        log_msg = 'Top serovar by Mash: "{}" with dist={}, # matching sketches={}, matching genome={}'
        logging.info(log_msg.format(mash_serovar, mash_distance, mash_match, mash_genome))

        mash_result_dict = {
            'mash_genome': mash_genome,
            'mash_serovar': mash_serovar,
            'mash_distance': mash_distance,
            'mash_match': mash_match
        }
        for k in mash_result_dict:
            prediction.__dict__[k] = mash_result_dict[k]
        break


def sistr_predict(input_fasta, tmp_dir, keep_tmp, should_run_mash, no_cgmlst):
    blast_runner = None
    try:
        assert os.path.exists(input_fasta), "Input fasta file '%s' must exist!" % input_fasta
        fasta_filename = os.path.basename(input_fasta)
        genome_tmp_dir = tmp_dir + '-' + fasta_filename
        blast_runner = BlastRunner(input_fasta, genome_tmp_dir)
        logging.info('Initializing temporary analysis directory "%s" and preparing for BLAST searching.', genome_tmp_dir)
        blast_runner.prep_blast()
        logging.info('Temporary FASTA file copied to %s', blast_runner.tmp_fasta_path)
        serovar_predictor = SerovarPredictor(blast_runner)
        serovar_predictor.predict_serovar_from_antigen_blast()
        prediction = serovar_predictor.get_serovar_prediction()
        prediction.genome = fasta_filename
        if should_run_mash:
            run_mash(input_fasta, prediction)
        if not no_cgmlst:
            run_cgmlst(prediction, blast_runner)
        serovar_from_cgmlst_and_antigen_blast(prediction, serovar_predictor)
        logging.info('%s | Antigen gene BLAST serovar prediction: "%s" serogroup=%s:H1=%s:H2=%s',
                     fasta_filename,
                     prediction.serovar_antigen,
                     prediction.serogroup,
                     prediction.h1,
                     prediction.h2)
        logging.info('%s | Overall serovar prediction: %s',
                     fasta_filename,
                     prediction.serovar)
    finally:
        if not keep_tmp:
            logging.info('Deleting temporary working directory at %s', blast_runner.tmp_work_dir)
            blast_runner.cleanup()
        else:
            logging.info('Keeping temp dir at %s', blast_runner.tmp_work_dir)
    return prediction


def main():
    parser = init_parser()
    args = parser.parse_args()
    init_console_logger(args.verbose)
    input_fastas = args.fastas
    if len(input_fastas) == 0:
        logging.error('No FASTA files specified!')
        exit(1)

    tmp_dir = args.tmp_dir
    keep_tmp = args.keep_tmp
    output_format = args.output_format
    output_path = args.output_dest

    from multiprocessing import Pool
    n_threads = args.threads
    logging.info('Initializing thread pool with %s threads', n_threads)
    pool = Pool(processes=n_threads)
    logging.info('Running SISTR analysis asynchronously on %s genomes', len(input_fastas))
    res = [pool.apply_async(sistr_predict, (input_fasta, tmp_dir, keep_tmp, args.run_mash, args.no_cgmlst)) for input_fasta in input_fastas]
    logging.info('Getting SISTR analysis results')
    outputs = [x.get() for x in res]

    if output_path:
        from sistr.src.writers import write
        write(output_path, output_format, outputs)
    else:
        logging.warning('No output file written!')


if __name__ == '__main__':
    main()
