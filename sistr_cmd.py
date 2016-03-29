import argparse
import logging

from src.blast_wrapper import BlastRunner, BlastReader
from src.logger import init_console_logger
from src.serovar_prediction import SerovarPredictor, serovar_from_cgmlst_and_antigen_blast

prog_desc = '''
SISTR (Salmonella In Silico Typing Resource) Command-line Tool
==============================================================
Serovar predictions from whole-genome sequence assemblies by determination of antigen gene and cgMLST gene alleles using BLAST.
'''

parser = argparse.ArgumentParser(prog='predict_serovar',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=prog_desc)

parser.add_argument('-i',
                    '--input',
                    help='Input genome FASTA file')
parser.add_argument('-f',
                    '--output_format',
                    default='json',
                    help='Output format (json, csv, pickle)')
parser.add_argument('-o',
                    '--output_dest',
                    help='Output')
parser.add_argument('-T',
                    '--tmp_dir',
                    default='tmp',
                    help='Base temporary working directory for intermediate analysis files.')
parser.add_argument('-K',
                    '--keep_tmp',
                    action='store_true',
                    help='Keep temporary analysis files.')
# parser.add_argument('--antigen_only', action='store_true', help='Antigen gene serovar prediction only')
parser.add_argument('--no_cgmlst',
                    action='store_true',
                    help='Do not run cgMLST serovar prediction')
parser.add_argument('-m', '--run_mash',
                    action='store_true',
                    help='''Determine Mash MinHash genomic distances to Salmonella genomes with trusted serovar designations.
                    Mash binary must be in accessible via $PATH (e.g. /usr/bin).
                    ''')
parser.add_argument('-v',
                    '--verbose',
                    action='count',
                    default=2,
                    help='Logging verbosity level (-v == show warnings; -vvv == show debug info)')


def run_mash():
    from src.mash import mash_dist_trusted, mash_output_to_pandas_df

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
        for k, v in mash_result_dict.iteritems():
            prediction.__dict__[k] = v
        break


if __name__ == '__main__':
    args = parser.parse_args()
    init_console_logger(args.verbose)

    input_fasta = args.input
    tmp_dir = args.tmp_dir
    keep_tmp = args.keep_tmp
    output_format = args.output_format
    output_path = args.output_dest

    if input_fasta is None:
        logging.error('You need to specify a genome fasta path!')
        parser.print_help()
        exit(1)

    import os
    assert os.path.exists(input_fasta), "Input fasta file '{}' must exist!".format(input_fasta)

    fasta_filename = os.path.basename(input_fasta)


    blast_runner = BlastRunner(input_fasta, tmp_dir)
    logging.info('Initializing temporary analysis directory and preparing for BLAST searching.')
    blast_runner.prep_blast()
    logging.info('Temporary FASTA file copied to {}'.format(blast_runner.tmp_fasta_path))

    serovar_predictor = SerovarPredictor(blast_runner)
    serovar_predictor.predict_serovar_from_antigen_blast()

    prediction = serovar_predictor.get_serovar_prediction()

    if not args.no_cgmlst:
        # TODO: refactor cgmlst into cgmlst.py
        from src.cgmlst import *

        df_cgmlst_profiles = cgmlst_profiles()

        logging.debug('{} distinct cgMLST330 profiles'.format(df_cgmlst_profiles.shape[0]))

        logging.info('Running BLAST on serovar predictive cgMLST330 alleles')
        blast_outfile = blast_runner.blast_against_query(CGMLST_FASTA_PATH)
        logging.info('Reading BLAST output file "{}"'.format(blast_outfile))
        blast_reader = BlastReader(blast_outfile)
        logging.info('Found {} cgMLST330 allele BLAST results'.format(blast_reader.df.shape[0]))
        df_perfect_cgmlst_results = blast_reader.perfect_matches()
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

    serovar_from_cgmlst_and_antigen_blast(prediction, serovar_predictor)

    if not keep_tmp:
        logging.info('Deleting temporary working directory at {}'.format(blast_runner.tmp_work_dir))
        blast_runner.cleanup()
    else:
        logging.info('Keeping temp dir at {}'.format(blast_runner.tmp_work_dir))

    if args.run_mash:
        run_mash()

    logging.info('Antigen gene BLAST serovar prediction: "{}" serogroup={}:H1={}:H2={}'.format(prediction.serovar_antigen, prediction.serogroup, prediction.h1, prediction.h2))
    logging.info('Overall serovar prediction: {}'.format(prediction.serovar))

    if output_path:
        from src.writers import write
        write(output_path, output_format, prediction)
