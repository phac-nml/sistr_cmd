import numpy as np


NT_SUB = {x:y for x,y in zip('acgtrymkswhbvdnxACGTRYMKSWHBVDNX', 'tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX')}


revcomp = lambda s: ''.join([NT_SUB[c] for c in s[::-1]])


def extend_subj_match_vec(df):
    """
    Get the extended clipped (clamped) start and end subject sequence indices

    Also get whether each match needs to be reverse complemented and whether each extended match would be truncated by
    the end of the subject sequence.

    Args:
        df (pandas.DataFrame): blastn results dataframe

    Returns:
        int pandas.Series: extended and clipped start indices
        int pandas.Series: extended and clipped end indices
        bool pandas.Series: does extracted seq need to be reverse complemented?
        bool pandas.Series: would the extended seq be truncated by the ends of the subject sequence?
        bool pandas.Series: was the subject seq extended?
    """
    needs_revcomp = df.sstart > df.send
    add_to_end = df.qlen - df.qend
    add_to_start = df.qstart - 1
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
    is_extended = (add_to_start > 0) | (add_to_end > 0)
    return clipped_start_idx, clipped_end_idx, needs_revcomp, trunc, is_extended


def retrieve_seq(seq, start, end, needs_revcomp):
    start = int(start)
    end = int(end)
    out_seq = seq[start:(end + 1)]
    if needs_revcomp:
        out_seq = revcomp(out_seq)
    return out_seq
