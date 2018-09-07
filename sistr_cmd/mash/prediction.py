import logging
import os
import re
from io import BytesIO
from subprocess import Popen, PIPE

import attr
import pandas as pd
from pkg_resources import resource_filename
from typing import Optional

from sistr_cmd.serovar_prediction.constants import genomes_to_serovar, genomes_to_subspecies
from sistr_cmd.util import first_row


@attr.s
class MashPrediction:
    serovar = attr.ib(default=None)
    subspecies = attr.ib(default=None)
    distance = attr.ib(default=1.0)
    n_matching = attr.ib(default=0)
    matching_genome = attr.ib(default=None)
    matching_genome_serovar = attr.ib(default=None)
    matching_genome_subspecies = attr.ib(default=None)
    top_matches = attr.ib(default=None)

    MASH_DIST_COLUMNS = ['matching_genome', 'query', 'distance', 'p_value', 'matching_sketches']
    MASH_BIN = 'mash'
    REF_SKETCH_FILE = resource_filename('sistr_cmd', 'data/sistr_cmd.msh')
    SUBSPECIATION_DISTANCE_THRESHOLD = 0.01
    DISTANCE_THRESHOLD = 0.005
    N_TOP_RESULTS = 5

    @staticmethod
    def run_mash_dist(fasta_path: str) -> bytes:
        """Run `mash dist` with genome FASTA against Salmonella Mash Sketch DB.

        Args:
            fasta_path: Genome FASTA file path

        Returns:
            Mash dist standard output
        """
        args = [MashPrediction.MASH_BIN,
                'dist',
                MashPrediction.REF_SKETCH_FILE,
                fasta_path]
        p = Popen(args, stderr=PIPE, stdout=PIPE)
        (stdout, stderr) = p.communicate()
        retcode = p.returncode
        if retcode != 0:
            raise Exception('Could not run Mash dist {}'.format(stderr))
        return stdout

    @staticmethod
    def mash_dataframe(mash_out: bytes) -> pd.DataFrame:
        df = pd.read_table(BytesIO(mash_out), header=None)
        df.columns = MashPrediction.MASH_DIST_COLUMNS
        # get base filenames with file extension removed from reference genome column
        matching_genomes = [re.sub(r'(\.fa$)|(\.fas$)|(\.fasta$)|(\.fna$)', '', os.path.basename(r)) for r in df['matching_genome']]
        df['matching_genome'] = matching_genomes
        # add serovar and subspecies info to Mash results dataframe
        genome_serovar = genomes_to_serovar()
        df['serovar'] = [genome_serovar[genome] for genome in matching_genomes]
        genome_spp = genomes_to_subspecies()
        df['subspecies'] = [genome_spp[genome] for genome in matching_genomes]
        df['n_matching'] = [int(x.split('/')[0]) for x in df['matching_sketches']]
        df.sort_values(by='distance', inplace=True)
        # filter for rows with distance below 1.0 since many rows may have a distance of 1.0 which is uninformative
        df = df[df['distance'] < 1.0]
        return df



    @staticmethod
    def mash_subspeciation(df_mash: pd.DataFrame) -> Optional[str]:
        min_dist = df_mash['distance'].min()
        if min_dist > MashPrediction.SUBSPECIATION_DISTANCE_THRESHOLD:
            logging.warning(f'Minimum Mash distance of {min_dist} was over subspeciation distance threshold of'
                            f' {MashPrediction.SUBSPECIATION_DISTANCE_THRESHOLD}. Cannot assign subspecies with '
                            f'confidence.')
            return None
        else:
            df_mash_spp = df_mash[df_mash['distance'] <= MashPrediction.SUBSPECIATION_DISTANCE_THRESHOLD]
            spp_counts = df_mash_spp['subspecies'].value_counts()
            most_freq_spp = spp_counts.index[0]
            most_freq_spp_count = spp_counts[0]
            logging.info(f'Mash subspecies="{most_freq_spp}" with {most_freq_spp_count} reference genomes matching '
                         f'within a distance of {MashPrediction.SUBSPECIATION_DISTANCE_THRESHOLD}. All counts: '
                         f'{spp_counts}')
            return most_freq_spp

    @classmethod
    def predict(cls, input_fasta: str) -> 'MashPrediction':
        mash_out = MashPrediction.run_mash_dist(input_fasta)
        df_mash = MashPrediction.mash_dataframe(mash_out)
        df_top_mash_matches = df_mash.head(n=MashPrediction.N_TOP_RESULTS) # type: pd.DataFrame
        subspecies = MashPrediction.mash_subspeciation(df_mash)
        row = first_row(df_top_mash_matches)
        top_serovar = row['serovar']
        top_distance = row['distance']

        obj = cls()
        obj.serovar = top_serovar if top_distance <= MashPrediction.DISTANCE_THRESHOLD else None
        obj.subspecies = subspecies
        obj.distance = top_distance
        obj.n_matching = row['n_matching']
        obj.matching_genome = row['matching_genome']
        obj.matching_genome_serovar = top_serovar
        obj.matching_genome_subspecies = row['subspecies']
        obj.top_matches = df_top_mash_matches.to_dict('records')
        logging.info(f'{obj} for "{input_fasta}"')
        return obj
