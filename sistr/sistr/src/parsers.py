# TODO: check format of file? guess format maybe; use BioPython to parse variety of formats?

#: set: valid IUPAC nucleotide characters for checking FASTA format
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


def parse_fasta(filepath):
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
        seqs = []
        header = ''
        for line in f:
            line = line.strip()
            if line == '':
                continue
            if line[0] == '>':
                if header == '':
                    header = line.replace('>','')
                else:
                    yield header, ''.join(seqs)
                    seqs = []
                    header = line.replace('>','')
            else:
                seqs.append(line)
        yield header, ''.join(seqs)


def fasta_format_check(fasta_path, logger):
    """
    Check that a file is valid FASTA format.

     - First non-blank line needs to begin with a '>' header character.
     - Sequence can only contain valid IUPAC nucleotide characters

    Args:
        fasta_str (str): FASTA file contents string

    Raises:
        Exception: If invalid FASTA format
    """

    header_count = 0
    line_count = 1
    nt_count = 0
    with open(fasta_path) as f:
        for l in f:
            l = l.strip()
            if l == '':
                continue
            if l[0] == '>':
                header_count += 1
                continue
            if header_count == 0 and l[0] != '>':
                error_msg = 'First non-blank line (L:{line_count}) does not contain FASTA header. Line beginning with ">" expected.' \
                    .format(line_count=line_count)
                logger.error(error_msg)
                raise Exception(error_msg)
            non_nucleotide_chars_in_line = set(l) - VALID_NUCLEOTIDES

            if len(non_nucleotide_chars_in_line) > 0:
                error_msg = 'Line {line} contains the following non-nucleotide characters: {non_nt_chars}' \
                    .format(line=line_count,
                            non_nt_chars=', '.join([x for x in non_nucleotide_chars_in_line]))
                logger.error(error_msg)
                raise Exception(error_msg)
            nt_count += len(l)
            line_count += 1

        if nt_count == 0:
            error_msg = 'File "{}" does not contain any nucleotide sequence.'.format(fasta_path)
            logger.error(error_msg)
            raise Exception(error_msg)

        logger.info('Valid FASTA format "{}" ({} bp)'.format(fasta_path, nt_count))
