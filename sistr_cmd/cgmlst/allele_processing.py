import logging
import zlib
from collections import defaultdict
from typing import Dict, List, Optional, Mapping, Union, Tuple

import pandas as pd

from sistr_cmd.blast_wrapper.helpers import extract
from sistr_cmd.cgmlst import MSA_GAP_PROP_THRESHOLD, msa_ref_vs_novel, number_gapped_ungapped
from sistr_cmd.io.parsers import parse_fasta


def hash_allele_seq(seq: str) -> int:
    """CRC32 unsigned integer from allele nucleotide sequence.
    The "& 0xffffffff" is used to generate an unsigned integer and
    will generate the same number for both Python 2 and 3
    (https://docs.python.org/2/library/zlib.html#zlib.crc32).

    Args:
        seq: nucleotide string

    Returns:
        CRC32 checksum as unsigned 32bit integer of `seq`
    """
    seq = str(seq).encode()
    return zlib.crc32(seq) & 0xffffffff


def alleles_to_retrieve(df: pd.DataFrame) -> Dict[str, List[pd.Series]]:
    """Alleles to retrieve from the input genome sequence.

    Get a dict of the genome fasta contig title to a list of blastn results of the allele sequences that must be
    retrieved from the genome contig.

    Args:
        df (pandas.DataFrame): BLASTN results in a Pandas DataFrame

    Returns:
        Top BLASTN result records for each marker grouped by the input contig they belong to
    """
    contig_blastn_records = defaultdict(list)
    markers = df.marker.unique()
    for m in markers:
        dfsub = df[df.marker == m]
        for i, r in dfsub.iterrows():
            if r.coverage < 1.0:
                contig_blastn_records[r.stitle].append(r)
            break
    return contig_blastn_records


def allele_result_dict(name: Optional[int], seq: Optional[str], blast_result: Mapping) -> Dict[
    str, Optional[Union[str, int, Mapping]]]:
    """Get cgMLST allele result dict

    Args:
        name: Allele name as number derived from hashing allele sequence
        seq: Allele sequence
        blast_result: Top BLAST result dict

    Returns:
        Allele information dict
    """
    return {
        'name': name,
        'seq': seq,
        'blast_result': blast_result
    }


def get_allele_sequences(genome_fasta_path: str,
                         cgmlst_fasta_path: str,
                         contig_blastn_records: Dict[str, List[pd.Series]]) -> Dict[str, Dict]:
    """Get allele sequences from input genome sequence file and assign allele name from hashing of allele sequence.

    For each allele, check:
    - that this allele does not introduce more than MSA_GAP_PROP_THRESHOLD % gaps when aligned against its closest matching reference allele


    Args:
        genome_fasta_path: Input genome FASTA path
        contig_blastn_records:
        using_full_cgmlst: Is the full cgMLST DB being used? Read from the full cgMLST FASTA file, otherwise, read from the centroid cgMLST allele DB file.

    Returns:
        Dict of marker name to dict of allele result
    """
    out = {}
    for contig_name, contig_seq in parse_fasta(genome_fasta_path):
        if contig_name not in contig_blastn_records:
            continue
        for blast_result in contig_blastn_records[contig_name]:  # type: pd.Series
            start_idx = blast_result['start_idx']
            end_idx = blast_result['end_idx']
            needs_revcomp = blast_result['needs_revcomp']
            marker = blast_result['marker']
            subject_title = blast_result['stitle']
            query_length = blast_result['qlen']
            query_id = blast_result['qseqid']
            logging.debug(f'seq len {len(contig_seq)}| start {start_idx}| end {end_idx}| revcomp? {needs_revcomp}')
            allele_seq = extract(contig_seq, start_idx, end_idx, needs_revcomp)
            ref_seqid = blast_result['qseqid']
            ref_seq = get_ref_cgmlst_sequence(cgmlst_fasta_path, ref_seqid)
            if ref_seq is None:
                raise IOError(f'Could not retrieve allele "{ref_seqid}" from "{cgmlst_fasta_path}"')
            trimmed_msa_novel, gapped, ungapped, p_gapped = check_full_allele_alignment(blast_result, ref_seq,
                                                                                        allele_seq)
            # if there are too many gaps within the trimmed extracted allele seq then result is equivalent
            # to missing or contig trunc
            if p_gapped > MSA_GAP_PROP_THRESHOLD:
                logging.warning(f'Too many gapped sites in extracted allele seq for marker {marker} '
                                f'contained {gapped} gaps out of {gapped + ungapped} bp ({p_gapped} > '
                                f'{MSA_GAP_PROP_THRESHOLD}); contig name/stitle: {subject_title} in file '
                                f'"{genome_fasta_path}". This allele will not be used in the final results.')
                blast_result['too_many_gaps'] = True
                out[blast_result.marker] = allele_result_dict(None, None, blast_result.to_dict())
                continue
            # otherwise if there are an acceptable number of gaps then remove gap characters and uppercase
            # Mafft MSA extracted and trimmed seq
            allele_seq = trimmed_msa_novel.replace('-', '').upper()
            new_allele_name = hash_allele_seq(allele_seq)
            logging.info(f'Marker {marker} | Recovered novel allele with gaps (n={gapped}) of length '
                         f'{len(allele_seq)} vs length {query_length} for ref allele {query_id}. Novel allele '
                         f'name={new_allele_name}')
            out[blast_result.marker] = allele_result_dict(new_allele_name, allele_seq, blast_result.to_dict())
    return out


def get_ref_cgmlst_sequence(cgmlst_fasta_path: str, seq_id: str) -> Optional[str]:
    """Try to find the reference cgMLST sequence for a given sequence id

    Args:
        cgmlst_fasta_path: Reference cgMLST alleles FASTA file path
        seq_id: cgMLST sequence id

    Returns:
        cgMLST sequence if it exists in the cgMLST FASTA file
    """
    for header, fasta_seq in parse_fasta(cgmlst_fasta_path):
        if header == seq_id:
            return fasta_seq


def check_full_allele_alignment(blast_result: pd.Series, ref_seq: str, allele_seq: str) -> Tuple[str, int, int, float]:
    msa_ref, msa_novel = msa_ref_vs_novel(ref_seq, allele_seq)
    # if there are gaps at the start or end of the ref allele MSA then trim those from both MSAs
    trim_left = 0
    while (msa_ref[trim_left] == '-'):
        trim_left += 1
    trim_right = len(msa_ref)
    while (msa_ref[trim_right - 1] == '-'):
        trim_right -= 1
    trimmed_msa_ref = msa_ref[trim_left:trim_right]
    trimmed_msa_novel = msa_novel[trim_left:trim_right]
    logging.debug(msa_ref)
    logging.debug(msa_novel)
    logging.debug('%s:%s', trim_left, trim_right)
    logging.debug(trimmed_msa_ref)
    logging.debug(trimmed_msa_novel)
    gapped, ungapped = number_gapped_ungapped(trimmed_msa_ref, trimmed_msa_novel)
    p_gapped = gapped / float((gapped + ungapped))
    blast_result['qseq_msa'] = msa_ref
    blast_result['qseq_msa_trimmed'] = trimmed_msa_ref
    blast_result['sseq_msa'] = msa_novel
    blast_result['sseq_msa_trimmed'] = trimmed_msa_novel
    blast_result['sseq_msa_gaps'] = gapped
    blast_result['sseq_msa_p_gaps'] = p_gapped
    return trimmed_msa_novel, gapped, ungapped, p_gapped