# -*- coding: utf-8 -*-
from typing import Tuple, Union

import numpy as np
import pandas as pd

#: dict for nucleotide complement substitution in `revcomp` function
NUCLEOTIDE_SUBSTITION = {x:y for x, y in zip('acgtrymkswhbvdnxACGTRYMKSWHBVDNX', 'tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX')}


def revcomp(nt_seq):
    """Reverse complement a nucleotide sequence

    Args:
        nt_seq: Nucleotide sequence

    Returns:
        Reverse complemented nucleotide sequence
    """
    return ''.join([NUCLEOTIDE_SUBSTITION[c] for c in nt_seq[::-1]])


def extend_matches(df: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray, pd.Series, pd.Series, pd.Series]:
    """Get the extended clipped (clamped) start and end subject sequence indices

    Also get whether each match needs to be reverse complemented and whether each extended match would be truncated by
    the end of the subject sequence.

    Args:
        df: BLAST+ results

    Returns:
        Tuple of:
            - integers for extended and clipped start indices
            - integers for extended and clipped end indices
            - booleans for "does extracted seq need to be reverse complemented?"
            - booleans for "would the extended seq be truncated by the ends of the subject sequence?"
            - booleans for "was the subject seq extended?"
    """
    needs_revcomp = df.sstart > df.send # type: pd.Series
    add_to_end = df.qlen - df.qend # type: pd.Series
    add_to_start = df.qstart - 1 # type: pd.Series
    ssum2 = (df.send + df.sstart) / 2.0
    sabs2 = np.abs(df.send - df.sstart) / 2.0
    end_idx = ssum2 + sabs2 - 1
    start_idx = ssum2 - sabs2 - 1
    start_idx[needs_revcomp] -= add_to_end
    start_idx[~needs_revcomp] -= add_to_start
    end_idx[needs_revcomp] += add_to_start
    end_idx[~needs_revcomp] += add_to_end
    clipped_start_idx = np.clip(start_idx, 0, (df.slen - 1))
    clipped_end_idx = np.clip(end_idx, 0, (df.slen - 1))
    trunc = (clipped_start_idx != start_idx) | (clipped_end_idx != end_idx)
    is_extended = (add_to_start > 0) | (add_to_end > 0) # type: pd.Series
    return clipped_start_idx, clipped_end_idx, needs_revcomp, trunc, is_extended


def extract(seq: str, start: int, end: int, needs_revcomp: bool) -> str:
    """

    :param seq:
    :param start:
    :param end:
    :param needs_revcomp:
    :return:
    """
    start = int(start)
    end = int(end)
    out_seq = seq[start:(end + 1)]
    if needs_revcomp:
        out_seq = revcomp(out_seq)
    return out_seq
