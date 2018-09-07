# -*- coding: utf-8 -*-

from __future__ import print_function

from collections import Counter
from datetime import datetime
import logging
import os
import re

from sistr_cmd.blast_wrapper.runner import BlastRunner
from sistr_cmd.cgmlst import cgmlst_serovar_prediction
from sistr_cmd.mash import MashPrediction
from sistr_cmd.qc import quality_check
from sistr_cmd.serovar_prediction import SerovarPrediction, serovar_table


def merge_mash_prediction(prediction, mash_prediction):
    for k in mash_prediction:
        prediction.__dict__[k] = mash_prediction[k]
    return prediction


def merge_cgmlst_prediction(serovar_prediction, cgmlst_prediction: ):
    serovar_prediction.cgmlst_distance = cgmlst_prediction.distance
    serovar_prediction.cgmlst_genome_match = cgmlst_prediction['genome_match']
    serovar_prediction.serovar_cgmlst = cgmlst_prediction['serovar']
    serovar_prediction.cgmlst_matching_alleles = cgmlst_prediction['matching_alleles']
    serovar_prediction.cgmlst_subspecies = cgmlst_prediction['subspecies']
    serovar_prediction.cgmlst_ST = cgmlst_prediction['cgmlst330_ST']
    return serovar_prediction


def infer_o_antigen(prediction):
    df_serovar = serovar_table()
    predicted_serovars = prediction.serovar.split('|') if '|' in prediction.serovar else [prediction.serovar]
    series_o_antigens = df_serovar.O_antigen[df_serovar.Serovar.isin(predicted_serovars)]
    if series_o_antigens.size == 0:
        sg = prediction.serogroup
        if sg is None or sg == '' or sg == '-':
            prediction.o_antigen = '-'
        else:
            series_o_antigens = df_serovar.O_antigen[df_serovar.Serogroup == sg]
            counter_o_antigens = Counter(series_o_antigens)
            prediction.o_antigen = counter_o_antigens.most_common(1)[0][0]
    else:
        counter_o_antigens = Counter(series_o_antigens)
        prediction.o_antigen = counter_o_antigens.most_common(1)[0][0]


def sistr_predict(input_fasta, genome_name=None, tmp_dir='/tmp', keep_tmp=False, run_mash=False, qc=True, no_cgmlst=False, use_full_cgmlst_db=False):
    blast_runner = None
    try:
        assert os.path.exists(input_fasta), "Input fasta file '%s' must exist!" % input_fasta
        if genome_name is None or genome_name == '':
            genome_name = genome_name_from_fasta_path(input_fasta)
        dtnow = datetime.now()
        genome_name_no_spaces = re.sub(r'\W', '_', genome_name)
        genome_tmp_dir = os.path.join(tmp_dir, dtnow.strftime("%Y%m%d%H%M%S") + '-' + 'SISTR' + '-' + genome_name_no_spaces)
        blast_runner = BlastRunner(input_fasta, genome_tmp_dir)
        logging.info('Initializing temporary analysis directory "%s" and preparing for BLAST searching.', genome_tmp_dir)
        blast_runner.prep_blast()
        logging.info('Temporary FASTA file copied to %s', blast_runner.tmp_fasta_path)
        spp = None
        mash_prediction = None
        if run_mash:
            mash_prediction = MashPrediction.predict(input_fasta)
            spp = mash_prediction['mash_subspecies']

        cgmlst_prediction = None
        if not no_cgmlst:
            cgmlst_prediction = cgmlst_serovar_prediction(blast_runner, use_full_cgmlst_db=use_full_cgmlst_db)
            spp = cgmlst_prediction.subspecies

        serovar_predictor = SerovarPrediction.predict(blast_runner)
        serovar_predictor.predict_serovar_from_antigen_blast()
        prediction = serovar_predictor.get_serovar_prediction()
        prediction.genome = genome_name
        prediction.fasta_filepath = os.path.abspath(input_fasta)
        if cgmlst_prediction:
            merge_cgmlst_prediction(prediction, cgmlst_prediction)
        if mash_prediction:
            merge_mash_prediction(prediction, mash_prediction)
        overall_serovar_call(prediction, serovar_predictor)
        infer_o_antigen(prediction)
        logging.info('%s | Antigen gene BLAST serovar prediction: "%s" serogroup=%s %s:%s:%s',
                     genome_name,
                     prediction.serovar_antigen,
                     prediction.serogroup,
                     prediction.o_antigen,
                     prediction.h1,
                     prediction.h2)
        logging.info('%s | Subspecies prediction: %s',
                     genome_name,
                     spp)
        logging.info('%s | Overall serovar prediction: %s',
                     genome_name,
                     prediction.serovar)
        if qc:
            qc_status, qc_msgs = quality_check(blast_runner.tmp_fasta_path, cgmlst_results, prediction)
            prediction.qc_status = qc_status
            prediction.qc_messages = ' | '.join(qc_msgs)
    finally:
        if not keep_tmp:
            logging.info('Deleting temporary working directory at %s', blast_runner.tmp_work_dir)
            blast_runner.cleanup()
        else:
            logging.info('Keeping temp dir at %s', blast_runner.tmp_work_dir)
    return prediction, cgmlst_results


def genome_name_from_fasta_path(fasta_path):
    """Extract genome name from fasta filename

    Get the filename without directory and remove the file extension.

    Example:
        With fasta file path ``/path/to/genome_1.fasta``::

            fasta_path = '/path/to/genome_1.fasta'
            genome_name = genome_name_from_fasta_path(fasta_path)
            print(genome_name)
            # => "genome_1"

    Args:
        fasta_path (str): fasta file path

    Returns:
        str: genome name
    """
    filename = os.path.basename(fasta_path)
    return re.sub(r'(\.fa$)|(\.fas$)|(\.fasta$)|(\.fna$)|(\.\w{1,}$)', '', filename)


def write_cgmlst_profiles(fastas, cgmlst_results, output_path):
    genome_marker_cgmlst_result = {}
    for genome, res in zip(fastas, cgmlst_results):
        tmp = {}
        for marker, res_dict in res.items():
            aname = res_dict['name']
            tmp[marker] = int(aname) if aname is not None else None
        genome_marker_cgmlst_result[genome] = tmp
    import pandas as pd
    df = pd.DataFrame(genome_marker_cgmlst_result).transpose()
    df.to_csv(output_path, float_format='%.0f')


def write_cgmlst_results_json(input_fastas, cgmlst_results, output_path):
    import json
    with open(output_path, 'w') as fout:
        json.dump({x:y for x,y in zip(input_fastas, cgmlst_results)}, fout)


def write_novel_alleles(cgmlst_results, output_path):
    count = 0
    with open(output_path, 'w') as fout:
        for x in cgmlst_results:
            for marker, res in x.items():
                name = res['name']
                seq = res['seq']
                br = res['blast_result']
                if br is not None and isinstance(br, dict):
                    trunc = br['trunc']
                    if not trunc:
                        fout.write('>{}|{}\n{}\n'.format(marker, name, seq))
                        count += 1
    return count


