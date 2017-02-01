from datetime import datetime
import logging
import shutil

from subprocess import Popen, PIPE
import os
import pandas as pd
import numpy as np
from pandas.io.common import EmptyDataError


BLAST_TABLE_COLS = '''
qseqid
stitle
pident
length
mismatch
gapopen
qstart
qend
sstart
send
evalue
bitscore
qlen
slen
sseq
'''.strip().split('\n')


class BlastRunner:
    blast_db_created = False

    def __init__(self, fasta_path, tmp_work_dir):
        self.tmp_work_dir = tmp_work_dir
        self.fasta_path = fasta_path


    def _create_tmp_folder(self):
        try:
            os.makedirs(self.tmp_work_dir)
        except:
            pass
        return self.tmp_work_dir

    def _copy_fasta_to_work_dir(self):
        filename = os.path.basename(self.fasta_path)
        dest_path = os.path.join(self.tmp_work_dir, filename)
        if self.fasta_path == dest_path:
            self.tmp_fasta_path = dest_path
            return dest_path
        shutil.copyfile(self.fasta_path, dest_path)
        self.tmp_fasta_path = dest_path
        return dest_path

    def _run_makeblastdb(self):
        work_dir = os.path.dirname(self.tmp_fasta_path)
        filename = os.path.basename(self.tmp_fasta_path)
        nin_filepath = os.path.join(work_dir, filename + '.nin')
        if os.path.exists(nin_filepath):
            self.blast_db_created = True
            return self.tmp_fasta_path

        p = Popen(['makeblastdb',
                   '-in', '{}'.format(self.tmp_fasta_path),
                   '-dbtype', 'nucl'],
                  stdout=PIPE,
                  stderr=PIPE)
        p.wait()
        stdout = p.stdout.read()
        stderr = p.stderr.read()
        if stdout is not None and stdout != '':
            logging.debug('makeblastdb on {0} STDOUT: {1}'.format(self.tmp_fasta_path, stdout))
        if stderr is not None and stderr != '':
            logging.debug('makeblastdb on {0} STDERR: {1}'.format(self.tmp_fasta_path, stderr))

        if os.path.exists(nin_filepath):
            self.blast_db_created = True
            return self.tmp_fasta_path
        else:
            ex_msg = 'makeblastdb was not able to create a BLAST DB for {0}. STDERR: {1}'.format(filename, stderr)
            logging.error(ex_msg)
            raise Exception(ex_msg)

    def blast_against_query(self, query_fasta_path, blast_task='megablast', evalue=1e-20, min_pid=85):

        if not self.blast_db_created:
            self.prep_blast()

        gene_filename = os.path.basename(query_fasta_path)
        genome_filename = os.path.basename(self.tmp_fasta_path)
        timestamp = '{:%Y%b%d_%H_%M_%S}'.format(datetime.now())
        outfile = os.path.join(self.tmp_work_dir, '{}-{}-{}.blast'.format(gene_filename,
                                                                          genome_filename,
                                                                          timestamp))
        p = Popen(['blastn',
                   '-task', blast_task,
                   '-query', query_fasta_path,
                   '-db', '{}'.format(self.tmp_fasta_path),
                   '-evalue', '{}'.format(evalue),
                   '-dust', 'no',
                   '-perc_identity', '{}'.format(min_pid),
                   '-out', outfile,
                   '-outfmt', '6 {}'.format(' '.join(BLAST_TABLE_COLS))],
                  stdout=PIPE,
                  stderr=PIPE)

        p.wait()

        stdout = p.stdout.read()
        stderr = p.stderr.read()
        if stdout is not None and stdout != '':
            logging.debug('blastn on db {} and query {} STDOUT: {}'.format(genome_filename, gene_filename, stdout))
        if stderr is not None and stderr != '':
            logging.debug('blastn on db {} and query {} STDERR: {}'.format(genome_filename, gene_filename, stderr))

        if os.path.exists(outfile):
            return outfile
        else:
            ex_msg = 'blastn on db {} and query {} did not produce expected output file at {}'.format(genome_filename,
                                                                                                      gene_filename,
                                                                                                      outfile)
            logging.error(ex_msg)
            raise Exception(ex_msg)

    def cleanup(self):
        self.blast_db_created = False
        shutil.rmtree(self.tmp_work_dir)

    def prep_blast(self):
        self._create_tmp_folder()
        self._copy_fasta_to_work_dir()
        self._run_makeblastdb()

    def run_blast(self, query_fasta_path):
        self.prep_blast()
        blast_outfile = self.blast_against_query(query_fasta_path)
        return blast_outfile


class BlastReader:
    is_missing = True
    is_perfect_match = False
    is_trunc = False
    df = None


    def __init__(self, blast_outfile):
        """Read BLASTN output file into a pandas DataFrame
        Sort the DataFrame by BLAST bitscore.
        If there are no BLASTN results, then no results can be returned.

        Args:
            blast_outfile (str): `blastn` output file path

        Raises:
            EmptyDataError: No data could be parsed from the `blastn` output file
        """
        self.blast_outfile = blast_outfile
        try:
            self.df = pd.read_table(self.blast_outfile, header=None)
            self.df.columns = BLAST_TABLE_COLS
            # calculate the coverage for when results need to be validated
            self.df.loc[:, 'coverage'] = self.df.length / self.df.qlen
            self.df.sort_values(by='bitscore', ascending=False, inplace=True)
            self.df.loc[:, 'is_trunc'] = BlastReader.trunc(qstart=self.df.qstart,
                                                           qend=self.df.qend,
                                                           qlen=self.df.qlen,
                                                           sstart=self.df.sstart,
                                                           send=self.df.send,
                                                           slen=self.df.slen)

            logging.debug(self.df.head())
            self.is_missing = False
        except EmptyDataError as exc:
            logging.warning('No BLASTN results to parse from file %s', blast_outfile)
            self.is_missing = True


    def df_dict(self):
        if not self.is_missing:
            return self.df.to_dict()

    @staticmethod
    def df_first_row_to_dict(df):
        """First DataFrame row to list of dict

        Args:
            df (pandas.DataFrame): A DataFrame with at least one row

        Returns:
            A list of dict that looks like:

                [{'C1': 'x'}, {'C2': 'y'}, {'C3': 'z'}]

            from a DataFrame that looks like:

                    C1  C2  C3
                1   x   y   z

            Else if `df` is `None`, returns `None`
        """
        if df is not None:
            return [dict(r) for i, r in df.head(1).iterrows()][0]

    @staticmethod
    def is_blast_result_trunc(qstart, qend, sstart, send, qlen, slen):
        """Check if a query sequence is truncated by the end of a subject sequence

        Args:
            qstart (int): Query sequence start index
            qend (int): Query sequence end index
            sstart (int): Subject sequence start index
            send (int): Subject sequence end index
            qlen (int): Query sequence length
            slen (int): Subject sequence length

        Returns:
            bool: Result truncated by subject sequence end?
        """
        q_match_len = abs(qstart - qend) + 1
        s_max = max(sstart, send)
        s_min = min(sstart, send)
        return (q_match_len < qlen) and (s_max >= slen or s_min <= 1)

    @staticmethod
    def trunc(qstart, qend, sstart, send, qlen, slen):
        """Check if a query sequence is truncated by the end of a subject sequence

        Args:
            qstart (int pandas.Series): Query sequence start index
            qend (int pandas.Series): Query sequence end index
            sstart (int pandas.Series): Subject sequence start index
            send (int pandas.Series): Subject sequence end index
            qlen (int pandas.Series): Query sequence length
            slen (int pandas.Series): Subject sequence length

        Returns:
            Boolean pandas.Series: Result truncated by subject sequence end?
        """
        ssum2 = (send + sstart) / 2.0
        sabs2 = np.abs(send - sstart) / 2.0
        smax = ssum2 + sabs2
        smin = ssum2 - sabs2
        q_match_len = np.abs(qstart - qend) + 1
        return (q_match_len < qlen) & ((smax >= slen) |  (smin <= 1))


    def perfect_matches(self):
        """
        Return pandas DataFrame with perfect BLAST matches (100% identity and coverage)

        Returns:
            pandas.DataFrame or None: DataFrame of perfect BLAST matches or None if no perfect matches exist
        """
        if self.is_missing:
            return None

        df_perfect_matches = self.df[(self.df['coverage'] == 1.0) & (self.df['pident'] == 100.0)]
        if df_perfect_matches.shape[0] == 0:
            return None
        return df_perfect_matches



    def top_result(self):
        """Return top `blastn` result
        Try to find a 100% identity and coverage result (perfect match).
        If one does not exist, then retrieve the result with the highest bitscore.

        Returns:
            Ordered dict of BLASTN results or None if no BLASTN results generated
        """

        if self.is_missing:
            return None

        df_perfect_matches = self.df[(self.df['coverage'] == 1.0) & (self.df['pident'] == 100.0)]
        if df_perfect_matches.shape[0]:
            self.is_perfect_match = True
            return BlastReader.df_first_row_to_dict(df_perfect_matches)

        # Return the result with the highest bitscore.
        # This is the first result in dataframe since the df is ordered by
        # bitscore in descending order.
        result_dict = BlastReader.df_first_row_to_dict(self.df)
        result_trunc = BlastReader.is_blast_result_trunc(qstart=result_dict['qstart'],
                                                         qend=result_dict['qend'],
                                                         sstart=result_dict['sstart'],
                                                         send=result_dict['send'],
                                                         qlen=result_dict['qlen'],
                                                         slen=result_dict['slen'])
        self.is_trunc = result_trunc
        return result_dict