# -*- coding: utf-8 -*-
"""
File parsing functions and helpers
"""

#: set: valid IUPAC nucleotide characters for checking FASTA format
from typing import Tuple, Iterable

VALID_NUCLEOTIDES = {'A', 'a',
                     'C', 'c',
                     'G', 'g',
                     'T', 't',
                     'R', 'r',
                     'Y', 'y',
                     'S', 's',
                     'W', 'w',
                     'K', 'k',
                     'M', 'm',
                     'B', 'b',
                     'D', 'd',
                     'H', 'h',
                     'V', 'v',
                     'N', 'n',
                     'X', 'x', } # X for masked nucleotides


def parse_fasta(filepath: str) -> Iterable[Tuple[str, str]]:
    '''
    Parse a fasta file returning a generator yielding tuples of fasta headers to sequences.

    Note:
        This function should give equivalent results to SeqIO from BioPython

        .. code-block:: python

            from Bio import SeqIO
            # biopython to dict of header-seq
            hseqs_bio = {r.description:str(r.seq) for r in SeqIO.parse(fasta_path, 'fasta')}
            # this func to dict of header-seq
            hseqs = {header:seq for header, seq in parse_fasta(fasta_path)}
            # both methods should return the same dict
            assert hseqs == hseqs_bio

    Args:
        filepath (str): Fasta file path

    Returns:
        generator: yields tuples of (<fasta header>, <fasta sequence>)
    '''
    with open(filepath, 'r') as f:
        seqs = ''
        header = ''
        for line in f:
            line = line.strip()
            if line == '':
                continue
            if line.startswith('>'):
                if header == '':
                    header = line.replace('>','')
                else:
                    yield header, seqs
                    seqs = ''
                    header = line.replace('>','')
            else:
                seqs += line
        yield header, seqs
