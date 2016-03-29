import os
from subprocess import Popen, PIPE
import pandas as pd
from src.cgmlst import genomes_to_serovar

MASH_BIN = 'mash'

base_dir = os.path.join(os.path.dirname(__file__), '..')
here = lambda x: os.path.abspath(os.path.join(base_dir, x))

MASH_SKETCH_FILE = here('data/mash-sistr-3967-trusted.msh')


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
    from StringIO import StringIO
    mash_header_list = ['ref', 'query', 'dist', 'pval', 'matching']
    mash_header_str = '\t'.join(mash_header_list)
    sio_o = StringIO('\n'.join([mash_header_str, mash_out]))
    df = pd.read_table(sio_o)
    refs = [r.replace('.fasta', '') for r in df['ref']]
    df['ref'] = refs
    genome_serovar_dict = genomes_to_serovar()
    df['serovar'] = [genome_serovar_dict[genome] for genome in refs]
    df['n_match'] = [int(x.split('/')[0]) for x in df['matching']]
    df.sort_values(by='dist', inplace=True)
    return df
