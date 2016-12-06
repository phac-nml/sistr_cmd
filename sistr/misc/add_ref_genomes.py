#!/usr/bin/env python
import argparse
from collections import defaultdict
import logging
import os
from subprocess import Popen
import re
from datetime import datetime

import pandas as pd
import numpy as np

from sistr.misc.reduce_to_centroid_alleles import run_allele_reduction
from sistr.sistr_cmd import genome_name_from_fasta_path
from sistr.src.blast_wrapper import BlastRunner
from sistr.src.logger import init_console_logger
from sistr.src.parsers import parse_fasta
from sistr.src.serovar_prediction import SerovarPredictor, overall_serovar_call
from sistr.src.cgmlst import CGMLST_PROFILES_PATH, run_cgmlst, allele_name, CGMLST_FULL_FASTA_PATH
from sistr.src.serovar_prediction.constants import GENOMES_TO_SEROVAR_PATH, GENOMES_TO_SPP_PATH, SEROVAR_TABLE_PATH
from sistr.src.mash import MASH_SKETCH_FILE


def init_parser():
    prog_desc = '''Add reference genomes to sistr_cmd
Supply genome FASTA files and a table with genome name to serovar (and subspecies). If genome not present in table or table not supplied then the serovar and subspecies predictions will be used instead.

sistr_cmd ref genome info files will be written to an output directory
'''

    parser = argparse.ArgumentParser(prog='predict_serovar',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=prog_desc)

    parser.add_argument('fastas',
                        metavar='F',
                        nargs='+',
                        help='Input genome FASTA file(s). Genome names in filenames before file extension (e.g. for "g1.fasta" genome name is "g1")')
    parser.add_argument('-o',
                        '--outdir',
                        required=True,
                        help='Output destination')
    parser.add_argument('-s',
                        '--serovar-table',
                        help='Table with serovar (and subspecies). CSV expected if extension is .csv; otherwise tab delimited expected. Columns=[genome,serovar, subspecies(optional)]')
    parser.add_argument('--force',
                        action='store_true',
                        help='Force overwrite of output directory if it exists!')
    parser.add_argument('-T',
                        '--tmp-dir',
                        default='/tmp',
                        help='Base temporary working directory for intermediate analysis files.')
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=1,
                        help='Number of parallel threads to run sistr_cmd analysis.')
    parser.add_argument('-v',
                        '--verbose',
                        action='count',
                        default=2,
                        help='Logging verbosity level (-v to show warnings; -vvv to show debug info)')
    return parser


def sketch_fasta(fasta_path, outdir):
    """Create a Mash sketch from an input fasta file

    Args:
        fasta_path (str): input fasta file path. Genome name in fasta filename
        outdir (str): output directory path to write Mash sketch file to

    Returns:
        str: output Mash sketch file path
    """
    genome_name = genome_name_from_fasta_path(fasta_path)
    outpath = os.path.join(outdir, genome_name)
    args = ['mash', 'sketch', '-o', outpath, fasta_path]
    logging.info('Running Mash sketch with command: %s', ' '.join(args))
    p = Popen(args)
    p.wait()
    sketch_path = outpath + '.msh'
    assert os.path.exists(sketch_path), 'Mash sketch for genome {} was not created at {}'.format(
        genome_name,
        sketch_path)
    return sketch_path


def merge_sketches(outdir, sketch_paths):
    """Merge new Mash sketches with current Mash sketches

    Args:
        outdir (str): output directory to write merged Mash sketch file
        sketch_paths (list of str): Mash sketch file paths for input fasta files

    Returns:
        str: output path for Mash sketch file with new and old sketches
    """
    merge_sketch_path = os.path.join(outdir, 'sistr.msh')
    args = ['mash', 'paste', merge_sketch_path]
    for x in sketch_paths:
        args.append(x)
    args.append(MASH_SKETCH_FILE)
    logging.info('Running Mash paste with command: %s', ' '.join(args))
    p = Popen(args)
    p.wait()
    assert os.path.exists(merge_sketch_path), 'Merged sketch was not created at {}'.format(merge_sketch_path)
    return merge_sketch_path


def create_subdirs(outdir, *args):
    subdir = os.path.join(outdir, *args)
    try:
        os.makedirs(subdir)
        return subdir
    except Exception as ex:
        if os.path.exists(subdir):
            return subdir
        logging.error(ex)


def merge_cgmlst_prediction(serovar_prediction, cgmlst_prediction):
    serovar_prediction.cgmlst_distance = cgmlst_prediction['distance']
    serovar_prediction.cgmlst_genome_match = cgmlst_prediction['genome_match']
    serovar_prediction.serovar_cgmlst = cgmlst_prediction['serovar']
    serovar_prediction.cgmlst_matching_alleles = cgmlst_prediction['matching_alleles']
    serovar_prediction.cgmlst_subspecies = cgmlst_prediction['subspecies']
    return serovar_prediction


def run_sistr(input_fasta, tmp_dir):
    blast_runner = None
    try:
        assert os.path.exists(input_fasta), "Input fasta file '%s' must exist!" % input_fasta
        fasta_filename = os.path.basename(input_fasta)
        genome_name = genome_name_from_fasta_path(input_fasta)
        dtnow = datetime.now()
        genome_tmp_dir = os.path.join(tmp_dir, dtnow.strftime("%Y%m%d%H%M%S") + '-' + 'SISTR' + '-' + genome_name)
        blast_runner = BlastRunner(input_fasta, genome_tmp_dir)
        logging.info('Initializing temporary analysis directory "%s" and preparing for BLAST searching.',
                     genome_tmp_dir)
        blast_runner.prep_blast()
        logging.info('Temporary FASTA file copied to %s', blast_runner.tmp_fasta_path)

        cgmlst_prediction, cgmlst_results = run_cgmlst(blast_runner)

        spp = cgmlst_prediction['subspecies']

        serovar_predictor = SerovarPredictor(blast_runner, spp)
        serovar_predictor.predict_serovar_from_antigen_blast()
        prediction = serovar_predictor.get_serovar_prediction()
        merge_cgmlst_prediction(prediction, cgmlst_prediction)
        overall_serovar_call(prediction, serovar_predictor)
        logging.info('%s | Antigen gene BLAST serovar prediction: "%s" serogroup=%s:H1=%s:H2=%s',
                     fasta_filename,
                     prediction.serovar_antigen,
                     prediction.serogroup,
                     prediction.h1,
                     prediction.h2)
        logging.info('%s | Subspecies prediction: %s',
                     fasta_filename,
                     spp)
        logging.info('%s | Overall serovar prediction: %s',
                     fasta_filename,
                     prediction.serovar)
    finally:
        logging.info('Deleting temporary working directory at %s', blast_runner.tmp_work_dir)
        blast_runner.cleanup()
    return prediction, cgmlst_results


def cgmlst_profiles_df(fastas, cgmlst_results):
    genome_marker_cgmlst_result = {}
    for fasta, res in zip(fastas, cgmlst_results):
        genome = genome_name_from_fasta_path(fasta)
        tmp = {}
        for marker, res_dict in res.items():
            aname = res_dict['name']
            tmp[marker] = int(aname) if aname is not None else None
        genome_marker_cgmlst_result[genome] = tmp
    return pd.DataFrame(genome_marker_cgmlst_result).transpose()


def write_cgmlst_fasta(outdir, cgmlst_results):
    marker_allele_seqs = defaultdict(set)

    allowed_nts = set('ATGCatgc')
    for h, s in parse_fasta(CGMLST_FULL_FASTA_PATH):
        marker, allele = h.split('|')
        s = s.replace('-', '')
        forbidden_char = set(s) - allowed_nts
        if len(forbidden_char) > 0:
            logging.warning('Forbidden nucleotide characters %s in allele "%s". Skipping this allele!',
                            forbidden_char,
                            h)
            continue
        marker_allele_seqs[marker].add(s)

    # init default dict with int where values start as int 0
    new_allele_count = defaultdict(int)

    for x in cgmlst_results:
        for marker, res in x.items():
            seq = res['seq']
            if seq is not None:
                if seq not in marker_allele_seqs[marker]:
                    new_allele_count[marker] += 1
                if '-' in seq:
                    logging.error('marker %s | result %s', marker, res)
                marker_allele_seqs[marker].add(seq)

    for marker in sorted(new_allele_count.keys()):
        logging.info('Added %s new alleles for marker %s', new_allele_count[marker], marker)

    new_cgmlst_fasta_path = os.path.join(outdir, 'cgmlst-full.fasta')
    with open(new_cgmlst_fasta_path, 'w') as fout:
        for marker in sorted(marker_allele_seqs.keys()):
            seqs = marker_allele_seqs[marker]
            for seq in seqs:
                fout.write('>{}|{}\n{}\n'.format(marker, allele_name(seq), seq))
    logging.info('cgMLST FASTA written to "%s" with %s novel alleles',
                 new_cgmlst_fasta_path,
                 sum([v for k, v in new_allele_count.items()]))
    return new_cgmlst_fasta_path


def write_cgmlst_profiles_csv(outdir, cgmlst_results, genome_names):
    df_profiles_old = pd.read_csv(CGMLST_PROFILES_PATH, index_col=0)
    markers = df_profiles_old.columns

    genome_marker_allele_results = defaultdict(dict)
    for genome, cgmlst_result in zip(genome_names, cgmlst_results):
        for marker in markers:
            allele = None
            if marker in cgmlst_result:
                r = cgmlst_result[marker]
                if 'name' in r:
                    allele = int(r['name']) if r['name'] is not None else None
                else:
                    allele = None
            genome_marker_allele_results[genome][marker] = allele
    df_profiles_new = pd.DataFrame(genome_marker_allele_results).transpose()
    df_all_profiles = pd.concat([df_profiles_new, df_profiles_old])
    profiles_output_path = os.path.join(outdir, 'cgmlst-profiles.csv')
    df_all_profiles.to_csv(profiles_output_path, float_format='%.0f')
    assert os.path.exists(profiles_output_path), 'cgMLST profiles CSV file was not written to "{}"'.format(
        profiles_output_path)
    logging.info('cgMLST profiles (dim=%s) CSV written to "%s"',
                 df_all_profiles.shape,
                 profiles_output_path)


def read_genomes_to_x(path):
    out = {}
    with open(path) as f:
        for l in f:
            l = l.strip()
            g, s = l.split('\t')
            out[g] = s
    return out


def write_genomes_to_x_table(path, genome_to_x):
    with open(path, 'w') as fout:
        for k, v in genome_to_x.items():
            fout.write('{}\t{}\n'.format(k, v))


def write_serovar_and_spp_tables(outdir, df_serovar, predictions, genome_names):
    genome_serovar = read_genomes_to_x(GENOMES_TO_SEROVAR_PATH)
    genome_spp = read_genomes_to_x(GENOMES_TO_SPP_PATH)

    # prediction serovars and subspecies
    pred_genome_serovar = {}
    pred_genome_spp = {}
    for genome, prediction in zip(genome_names, predictions):
        pred_dict = prediction.__dict__
        pred_genome_serovar[genome] = pred_dict['serovar']
        if 'cgmlst_subspecies' in pred_dict:
            pred_genome_spp[genome] = pred_dict['cgmlst_subspecies']
        else:
            pred_genome_spp[genome] = None

    if df_serovar is not None:
        for i, row in df_serovar.iterrows():
            genome = row['genome']
            serovar = row['serovar']
            if not serovar in pred_genome_serovar[genome]:
                logging.warning('Genome "%s" user specified serovar "%s" not in serovar prediction "%s"',
                                genome,
                                serovar,
                                pred_genome_serovar[genome])
            if 'subspecies' in df_serovar:
                spp = row['subspecies']
                if spp != pred_genome_spp[genome]:
                    logging.warning('Genome "%s" provided subspecies of "%s" does not match prediction of "%s"',
                                    genome,
                                    spp,
                                    pred_genome_spp[genome])
            else:
                spp = pred_genome_spp[genome]
                logging.warning('Genome "%s" subspecies info not provided. Using subspecies prediction of "%s"',
                                genome,
                                spp)
            genome_serovar[genome] = serovar
            genome_spp[genome] = spp
    else:
        logging.warning(
            'User did not specify serovar/subspecies table! Using SISTR serovar and subspecies predictions for all genomes.')
        for genome in genome_names:
            genome_serovar[genome] = pred_genome_serovar[genome]
            genome_spp[genome] = pred_genome_spp[genome]
    genomes_to_serovar_path = os.path.join(outdir, 'genomes-to-serovar.txt')
    genomes_to_spp_path = os.path.join(outdir, 'genomes-to-subspecies.txt')
    write_genomes_to_x_table(genomes_to_serovar_path, genome_serovar)
    assert os.path.exists(genomes_to_serovar_path), '{} file could not be written!'.format(
        genomes_to_serovar_path)
    logging.info('Wrote genomes to serovars table at %s', genomes_to_serovar_path)
    write_genomes_to_x_table(genomes_to_spp_path, genome_spp)
    assert os.path.exists(genomes_to_spp_path), '{} file could not be written!'.format(
        genomes_to_spp_path)
    logging.info('Wrote genomes to subspecies table at %s', genomes_to_spp_path)


def create_merge_mash_sketches(input_fastas, data_outdir, sketch_outdir):
    sketch_paths = [sketch_fasta(fasta, sketch_outdir) for fasta in input_fastas]
    merge_sketches(data_outdir, sketch_paths)


def main():
    parser = init_parser()
    args = parser.parse_args()
    init_console_logger(args.verbose)
    logging.debug(args)
    input_fastas = args.fastas
    outdir = args.outdir
    tmp_dir = args.tmp_dir
    serovar_table_path = args.serovar_table
    threads = args.threads
    force = args.force

    assert len(input_fastas) > 0, 'No FASTA files specified!'
    for input_fasta in input_fastas:
        assert os.path.exists(input_fasta), 'Genome FASTA file does not exist at "{}"'.format(input_fasta)
    genome_names = [genome_name_from_fasta_path(x) for x in input_fastas]
    logging.info('You have specified %s genomes to add to current sistr_cmd data files! %s',
                 len(genome_names),
                 genome_names)

    if os.path.exists(outdir):
        if not force:
            raise Exception('Output directory already exists at {}!'.format(outdir))
        else:
            logging.warning('Using existing output directory at %s', outdir)
    try:
        os.makedirs(outdir)
    except:
        pass
    assert os.path.exists(outdir), 'Output directory could not be created!'

    if serovar_table_path:
        assert os.path.exists(serovar_table_path), 'Provided serovar table path does not exist! {}'.format(
            serovar_table_path)
        logging.info('Parsing serovar table from "%s"', serovar_table_path)
        if re.match(r'.*.csv$', serovar_table_path):
            logging.info('Trying to read serovar table "%s" as CSV', serovar_table_path)
            df_serovar = pd.read_csv(serovar_table_path)
        else:
            logging.info('Trying to read serovar table "%s" as tab-delimited', serovar_table_path)
            df_serovar = pd.read_table(serovar_table_path)
        expected_columns = ['genome', 'serovar']
        assert np.all(
            df_serovar.columns.isin(expected_columns)), 'User serovar table did not contain expected columns {}'.format(
            expected_columns)
        if 'subspecies' not in df_serovar.columns:
            logging.warning(
                'User serovar table did not contain "subspecies" column so the sistr_cmd subspecies prediction will be used!')
        genome_names_series = pd.Series(genome_names)
        genomes_in_serovar_table = genome_names_series.isin(df_serovar.genome)
        if not np.all(genomes_in_serovar_table):
            missing_genomes = ', '.join([x for x in genome_names_series[~genomes_in_serovar_table]])
            logging.error('The following genomes were not found in the serovar table: %s', missing_genomes)
            raise Exception('Not all user provided genome FASTA files in the provided serovar table!')

        df_wklm = pd.read_csv(SEROVAR_TABLE_PATH)
        logging.info('Checking for non-standard serovar designations')
        serovars_not_in_wklm = df_serovar.serovar[~df_serovar.serovar.isin(df_wklm.Serovar)]
        for row_idx, serovar in serovars_not_in_wklm.iteritems():
            logging.warning('Non-standard serovar %s at row %s for genome %s!', serovar, row_idx,
                            df_serovar.ix[row_idx]['genome'])
    else:
        logging.warning('No genome to serovar table specified! Using SISTR serovar predictions')
        df_serovar = None

    if threads == 1:
        logging.info('Serial single threaded run mode on %s genomes', len(input_fastas))
        outputs = [run_sistr(input_fasta, tmp_dir) for input_fasta in input_fastas]
    else:
        from multiprocessing import Pool

        logging.info('Initializing thread pool with %s threads', threads)
        pool = Pool(processes=threads)
        logging.info('Running SISTR analysis asynchronously on %s genomes', len(input_fastas))
        res = [pool.apply_async(run_sistr, (input_fasta, tmp_dir)) for input_fasta in input_fastas]

        logging.info('Getting SISTR analysis results')
        outputs = [x.get() for x in res]
    # collect results from sistr analysis
    prediction_outputs = [x for x, y in outputs]
    cgmlst_results = [y for x, y in outputs]
    # create some output dirs
    data_outdir = create_subdirs(outdir, 'data')
    cgmlst_outdir = create_subdirs(outdir, 'data', 'cgmlst')
    sketch_outdir = create_subdirs(outdir, 'mash-sketches')
    # write files with new and old data
    cgmlst_fasta = write_cgmlst_fasta(cgmlst_outdir, cgmlst_results)
    write_cgmlst_profiles_csv(cgmlst_outdir, cgmlst_results, genome_names)
    write_serovar_and_spp_tables(data_outdir, df_serovar, prediction_outputs, genome_names)
    create_merge_mash_sketches(input_fastas, data_outdir, sketch_outdir)

    centroid_alleles_path = os.path.join(cgmlst_outdir, 'cgmlst-centroid.fasta')
    run_allele_reduction(cgmlst_fasta, centroid_alleles_path, threads=threads)

    logging.info('Done!')


if __name__ == '__main__':
    main()
