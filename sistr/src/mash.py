import logging
import re, os
from pkg_resources import resource_filename
from subprocess import Popen, PIPE
import pandas as pd
from sistr.src.serovar_prediction.constants import genomes_to_serovar, genomes_to_subspecies
from sistr.src.serovar_prediction.constants import MASH_SUBSPECIATION_DISTANCE_THRESHOLD


MASH_BIN = 'mash'
MASH_SKETCH_FILE = resource_filename('sistr', 'data/sistr.msh')


def mash_dist_trusted(fasta_path):
    """
    Compute Mash distances of sketch file of genome fasta to RefSeq sketch DB.

    Args:
        mash_bin (str): Mash binary path

    Returns:
        (str): Mash STDOUT string
    """
    args = [MASH_BIN,
            'dist',
            MASH_SKETCH_FILE,
            fasta_path]
    p = Popen(args, stderr=PIPE, stdout=PIPE)
    (stdout, stderr) = p.communicate()
    retcode = p.returncode
    if retcode != 0:
        raise Exception('Could not run Mash dist {}'.format(stderr))

    return stdout


def mash_output_to_pandas_df(mash_out):
    from io import BytesIO
    df = pd.read_table(BytesIO(mash_out), header=None)
    df.columns = ['ref', 'query', 'dist', 'pval', 'matching']
    refs = [re.sub(r'(\.fa$)|(\.fas$)|(\.fasta$)|(\.fna$)', '', os.path.basename(r)) for r in df['ref']]
    df['ref'] = refs
    genome_serovar_dict = genomes_to_serovar()
    df['serovar'] = [genome_serovar_dict[genome] for genome in refs]
    df['n_match'] = [int(x.split('/')[0]) for x in df['matching']]
    df.sort_values(by='dist', inplace=True)
    df = df[df['dist'] < 1.0]
    return df


def mash_subspeciation(df_mash):
    closest_distance = df_mash['dist'].min()
    if closest_distance > MASH_SUBSPECIATION_DISTANCE_THRESHOLD:
        logging.warning('Min Mash distance (%s) above subspeciation distance threshold (%s)',
                        closest_distance,
                        MASH_SUBSPECIATION_DISTANCE_THRESHOLD)
        return None
    else:
        df_mash_spp = df_mash[df_mash['dist'] <= MASH_SUBSPECIATION_DISTANCE_THRESHOLD]
        genomes = df_mash_spp['ref']
        from collections import Counter
        genome_spp = genomes_to_subspecies()
        subspecies_below_threshold = [genome_spp[g] if g in genome_spp else None for g in genomes]
        subspecies_below_threshold = filter(None, subspecies_below_threshold)
        subspecies_counter = Counter(subspecies_below_threshold)
        logging.info('Mash subspecies counter: %s', subspecies_counter)
        return (subspecies_counter.most_common(1)[0][0], closest_distance, dict(subspecies_counter))
