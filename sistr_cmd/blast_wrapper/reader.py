import logging

import attr
import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError
from typing import Dict, Optional, Mapping

from sistr_cmd.blast_wrapper.const import BLAST_TABLE_COLS
from sistr_cmd.util import first_row


@attr.s
class BlastReader:
    blast_outfile = attr.ib()  # type: str
    df = attr.ib()  # type: pd.DataFrame
    is_missing = attr.ib(default=True)  # type: bool
    has_perfect_match = attr.ib(default=False)  # type: bool
    is_trunc = attr.ib(default=False)  # type: bool

    @df.default
    def read_blast_outfile(self):
        """Read BLASTN output file into a pandas DataFrame
        Sort the DataFrame by BLAST bitscore.
        If there are no BLASTN results, then no results can be returned.

        Args:
            blast_outfile (str): `blastn` output file path

        Raises:
            EmptyDataError: No data could be parsed from the `blastn` output file
        """
        df = None
        try:
            df = pd.read_table(self.blast_outfile, header=None, names=BLAST_TABLE_COLS)
            # calculate the coverage for when results need to be validated
            df.loc[:, 'coverage'] = df.length / df.qlen
            df.sort_values(by='bitscore', ascending=False, inplace=True)
            truncations = BlastReader.trunc(qstart=df.qstart,
                                            qend=df.qend,
                                            qlen=df.qlen,
                                            sstart=df.sstart,
                                            send=df.send,
                                            slen=df.slen)
            df.loc[:, 'is_trunc'] = truncations
            self.is_missing = False
        except EmptyDataError:
            self.is_missing = True
            logging.warning('No BLASTN results to parse from file %s', self.blast_outfile)
        return df

    def df_dict(self) -> Optional[Dict]:
        if not self.is_missing:
            return self.df.to_dict('records')

    @staticmethod
    def df_first_row_to_dict(df: pd.DataFrame):
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
        return (q_match_len < qlen) & ((smax >= slen) | (smin <= 1))

    def perfect_matches(self) -> Optional[pd.DataFrame]:
        """Get perfect BLAST matches (100% identity and coverage)

        Returns:
            DataFrame of perfect BLAST matches if any
        """
        if self.is_missing:
            return None
        df_perfect_matches = self.df[(self.df['coverage'] == 1.0) & (self.df['pident'] == 100.0)]
        if df_perfect_matches.shape[0] == 0:
            return None
        return df_perfect_matches

    def top_result(self) -> Optional[Mapping]:
        """Get top BLASTN result.

        Try to find a 100% identity and coverage, perfect match result. If one does not exist, then retrieve the result
        with the highest bitscore.

        Returns:
            Dict of top BLASTN result match info or `None` if no BLASTN results generated
        """

        if self.is_missing:
            return None

        df_perfect_matches = self.df[(self.df['coverage'] == 1.0) & (self.df['pident'] == 100.0)]
        if df_perfect_matches.shape[0]:
            self.has_perfect_match = True
            return first_row(df_perfect_matches).to_dict()

        # No perfect match? Then return the result with the highest bitscore. This is the first row since the DataFrame
        # is ordered by BLAST bitscore in descending order.
        top_result_row = first_row(self.df)
        self.is_trunc = top_result_row['is_trunc']
        return top_result_row.to_dict()
