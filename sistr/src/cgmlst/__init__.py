from pkg_resources import resource_filename
import logging
import zlib
from collections import Counter, defaultdict
import numpy as np
import pandas as pd
from sistr.src.blast_wrapper import BlastReader
from sistr.src.blast_wrapper.helpers import extend_subj_match_vec, retrieve_seq
from sistr.src.cgmlst.msa import msa_ref_vs_novel, number_gapped_ungapped, MSA_GAP_PROP_THRESHOLD
from sistr.src.parsers import parse_fasta
from sistr.src.serovar_prediction.constants import CGMLST_SUBSPECIATION_DISTANCE_THRESHOLD, genomes_to_subspecies


CGMLST_CENTROID_FASTA_PATH = resource_filename('sistr', 'data/cgmlst/cgmlst-centroid.fasta')
CGMLST_FULL_FASTA_PATH = resource_filename('sistr', 'data/cgmlst/cgmlst-full.fasta')
CGMLST_PROFILES_PATH = resource_filename('sistr', 'data/cgmlst/cgmlst-profiles.csv')
BLASTN_PIDENT_THRESHOLD = 90.0


def allele_name(seq):
    """CRC32 unsigned integer from allele nucleotide sequence.
    The "& 0xffffffff" is used to generate an unsigned integer and
    will generate the same number for both Python 2 and 3
    (https://docs.python.org/2/library/zlib.html#zlib.crc32).

    Args:
        seq (str): nucleotide string

    Returns:
        int: CRC32 checksum as unsigned 32bit integer
    """
    seq = str(seq).encode()
    return zlib.crc32(seq) & 0xffffffff


def ref_cgmlst_profiles():
    return pd.read_csv(CGMLST_PROFILES_PATH, index_col=0)


def process_cgmlst_results(df):
    """Append informative fields to cgMLST330 BLAST results DataFrame

    The `qseqid` column must contain cgMLST330 query IDs with `{marker name}|{allele number}` format.
    The `qseqid` parsed allele numbers and marker names are appended as new fields.

    `is_perfect` column contains boolean values for whether an allele result is 100% identity and coverage.
    `has_perfect_match` denotes if a cgMLST330 marker has a perfect allele match.
    The top result with the largest bitscore for a marker with no perfect match is used to retrieve the allele present
    at that marker locus.

    Args:
        df (pandas.DataFrame): DataFrame of cgMLST330 BLAST results

    Returns:
        pandas.DataFrame: cgMLST330 BLAST results DataFrame with extra fields (`marker`, `allele`, `is_perfect`, `has_perfect_match`)
    """
    assert isinstance(df, pd.DataFrame)
    markers = []
    alleles = []
    for x in  df['qseqid']:
        marker, allele = x.split('|')
        markers.append(marker)
        alleles.append(int(allele))
    df.loc[:, 'marker'] = markers
    df.loc[:, 'allele'] = alleles
    df.loc[:, 'is_match'] = (df['coverage'] >= 1.0) & (df['pident'] >= 90.0) & ~(df['is_trunc'])
    df.loc[:, 'allele_name'] = df.apply(lambda x: allele_name(x.sseq.replace('-', '')), axis=1)
    df.loc[:, 'is_perfect'] = (df['coverage'] == 1.0) & (df['pident'] == 100.0)
    df_perf = df[df['is_perfect']]
    perf_markers = df_perf['marker'].unique()
    df.loc[:, 'has_perfect_match'] = df['marker'].isin(perf_markers)
    start_idxs, end_idxs, needs_revcomps, trunc, is_extended = extend_subj_match_vec(df)
    df.loc[:, 'start_idx'] = start_idxs
    df.loc[:, 'end_idx'] = end_idxs
    df.loc[:, 'needs_revcomp'] = needs_revcomps
    df.loc[:, 'trunc'] = trunc
    df.loc[:, 'is_extended'] = is_extended
    df.loc[:, 'sseq_msa_gaps'] = np.zeros(df.shape[0], dtype=np.int64)
    df.loc[:, 'sseq_msa_p_gaps'] = np.zeros(df.shape[0], dtype=np.float64)
    df.loc[:, 'too_many_gaps'] = trunc

    return df


def alleles_to_retrieve(df):
    """Alleles to retrieve from genome fasta

    Get a dict of the genome fasta contig title to a list of blastn results of the allele sequences that must be
    retrieved from the genome contig.

    Args:
        df (pandas.DataFrame): blastn results dataframe

    Returns:
        {str:[pandas.Series]}: dict of contig title (header name) to list of top blastn result records for each marker
            for which the allele sequence must be retrieved from the original sequence.
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


def allele_result_dict(name, seq, blast_result):
    return {'name': name,
            'seq': seq,
            'blast_result': blast_result, }


def get_allele_sequences(genome_fasta_path, contig_blastn_records, full=False):
    cgmlst_fasta_path = CGMLST_CENTROID_FASTA_PATH if not full else CGMLST_FULL_FASTA_PATH
    out = {}
    for header, seq in parse_fasta(genome_fasta_path):
        if header in contig_blastn_records:
            for r in contig_blastn_records[header]:
                start_idx = r['start_idx']
                end_idx = r['end_idx']
                needs_revcomp = r['needs_revcomp']
                logging.debug('seq len {}| start {}| end {}| revcomp? {}'.format(len(seq), start_idx, end_idx, needs_revcomp))
                allele_seq = retrieve_seq(seq, start_idx, end_idx, needs_revcomp)
                ref_seqid = r['qseqid']
                ref_seq = None
                for h, s in parse_fasta(cgmlst_fasta_path):
                    if h == ref_seqid:
                        ref_seq = s
                        break
                if ref_seq is None:
                    raise Exception('Could not retrieve allele %s from %s', ref_seqid, cgmlst_fasta_path)
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
                r['qseq_msa'] = msa_ref
                r['qseq_msa_trimmed'] = trimmed_msa_ref
                r['sseq_msa'] = msa_novel
                r['sseq_msa_trimmed'] = trimmed_msa_novel
                r['sseq_msa_gaps'] = gapped
                r['sseq_msa_p_gaps'] = p_gapped
                # if there are too many gaps within the trimmed extracted allele seq then result is equivalent
                # to missing or contig trunc
                if p_gapped > MSA_GAP_PROP_THRESHOLD:
                    logging.error('Too many gapped sites in extracted allele seq for marker %s contained %s gaps out of %s bp (%s > %s); stitle: %s',
                                  r['marker'],
                                  gapped,
                                  (gapped + ungapped),
                                  p_gapped,
                                  MSA_GAP_PROP_THRESHOLD,
                                  r['stitle'])
                    r['too_many_gaps'] = True
                    out[r.marker] = allele_result_dict(None, None, r.to_dict())
                    continue
                # otherwise if there are an acceptable number of gaps then remove gap characters and uppercase
                # Mafft MSA extracted and trimmed seq
                allele_seq = trimmed_msa_novel.replace('-', '').upper()
                new_allele_name = allele_name(allele_seq)
                logging.info('Marker %s | Recovered novel allele with gaps (n=%s) of length %s vs length %s for ref allele %s. Novel allele name=%s',
                             r['marker'],
                             gapped,
                             len(allele_seq),
                             r['qlen'],
                             r['qseqid'],
                             new_allele_name)
                out[r.marker] = allele_result_dict(new_allele_name, allele_seq, r.to_dict())
    return out


def matches_to_marker_results(df):
    """Perfect BLAST matches to marker results dict

    Parse perfect BLAST matches to marker results dict.


    Args:
        df (pandas.DataFrame): DataFrame of perfect BLAST matches

    Returns:
        dict: cgMLST330 marker names to matching allele numbers
    """
    assert isinstance(df, pd.DataFrame)
    from collections import defaultdict
    d = defaultdict(list)
    for idx, row in df.iterrows():
        marker = row['marker']
        d[marker].append(row)

    marker_results = {}
    for k,v in d.items():
        if len(v) > 1:
            logging.debug('Multiple potential cgMLST allele matches (n=%s) found for marker %s. Selecting match on longest contig.', len(v), k)
            df_marker = pd.DataFrame(v)
            df_marker.sort_values('slen', ascending=False, inplace=True)
            for i,r in df_marker.iterrows():
                allele = r['allele_name']
                slen = r['slen']
                logging.debug('Selecting allele %s from contig with length %s', allele, slen)
                seq = r['sseq']
                if '-' in seq:
                    logging.warning('Gaps found in allele. Removing gaps. %s', r)
                    seq = seq.replace('-', '').upper()
                    allele = allele_name(seq)
                marker_results[k] = allele_result_dict(allele, seq, r.to_dict())
                break
        elif len(v) == 1:
            row = v[0]
            seq = row['sseq']
            if '-' in seq:
                logging.warning('Gaps found in allele. Removing gaps. %s', r)
                seq = seq.replace('-', '').upper()
                allele = allele_name(seq)
            marker_results[k] = allele_result_dict(allele, seq, row.to_dict())
        else:
            err_msg = 'Empty list of matches for marker {}'.format(k)
            logging.error(err_msg)
            raise Exception(err_msg)
    return marker_results


def find_closest_related_genome(marker_results, df_genome_profiles):
    """

    Args:
        df_genome_profiles (pandas.DataFrame):

    Returns:
        (dict, list): Most closely related Genome and list of other related Genomes_ in order of relatedness
    """
    marker_names = df_genome_profiles.columns
    n_markers = len(marker_names)

    profile = [marker_results[marker_name] if marker_name in marker_results else None for marker_name in marker_names]
    genome_profile = np.array(profile, dtype=np.float64)
    profiles_matrix = np.array(df_genome_profiles, dtype=np.float64)
    genome_profile_similarity_counts = np.apply_along_axis(lambda x: (x == genome_profile).sum(), 1, profiles_matrix)

    df_relatives = pd.DataFrame()
    df_relatives['matching'] = genome_profile_similarity_counts
    df_relatives['distance'] = 1.0 - (df_relatives['matching'] / float(n_markers))
    df_relatives.index = df_genome_profiles.index
    df_relatives.sort_values(by='distance', inplace=True)
    return df_relatives


def cgmlst_subspecies_call(df_relatives):
    """Call Salmonella subspecies based on cgMLST results

    This method attempts to find the majority subspecies type within curated
    public genomes above a cgMLST allelic profile distance threshold.

    Note:
        ``CGMLST_SUBSPECIATION_DISTANCE_THRESHOLD`` is the cgMLST distance
        threshold used to determine the subspecies by cgMLST. It is set at a
        distance of 0.9 which translates to a cgMLST allelic similarity of 10%.
        A threshold of 0.9 is generous and reasonable given the congruence
        between subspecies designations and 10% cgMLST clusters by Adjusted
        Rand (~0.850) and Adjusted Wallace metrics (~0.850 both ways).

    Args:
        df_relatives (pandas.DataFrame): Table of genomes related by cgMLST to input genome

    Returns:
        None: if no curated public genomes found to have a cgMLST profile similarity of 10% or greater
        (string, float, dict): most common subspecies, closest related public genome distance, subspecies frequencies
    """

    closest_distance = df_relatives['distance'].min()

    if closest_distance > CGMLST_SUBSPECIATION_DISTANCE_THRESHOLD:
        logging.warning('Min cgMLST distance (%s) above subspeciation distance threshold (%s)',
            closest_distance,
            CGMLST_SUBSPECIATION_DISTANCE_THRESHOLD)
        return None
    else:
        df_relatives = df_relatives[df_relatives.distance <= CGMLST_SUBSPECIATION_DISTANCE_THRESHOLD]
        df_relatives.sort_values('distance', ascending=True, inplace=True)
        logging.debug('df_relatives by cgmlst %s', df_relatives.head())
        genome_spp = genomes_to_subspecies()
        subspecies_below_threshold = [genome_spp[member_genome] if member_genome in genome_spp else None for member_genome in df_relatives.index]
        subspecies_below_threshold = filter(None, subspecies_below_threshold)
        subspecies_counter = Counter(subspecies_below_threshold)
        logging.debug('Subspecies counter: %s', subspecies_counter)
        return (subspecies_counter.most_common(1)[0][0], closest_distance, dict(subspecies_counter))


def run_cgmlst(blast_runner, full=False):
    """Perform in silico cgMLST on an input genome

    Args:
        blast_runner (sistr.src.blast_wrapper.BlastRunner): blastn runner object with genome fasta initialized

    Returns:
        dict: cgMLST ref genome match, distance to closest ref genome, subspecies and serovar predictions
        dict: marker allele match results (seq, allele name, blastn results)
    """
    from sistr.src.serovar_prediction.constants import genomes_to_serovar

    df_cgmlst_profiles = ref_cgmlst_profiles()

    logging.debug('{} distinct cgMLST330 profiles'.format(df_cgmlst_profiles.shape[0]))

    logging.info('Running BLAST on serovar predictive cgMLST330 alleles')
    cgmlst_fasta_path = CGMLST_CENTROID_FASTA_PATH if not full else CGMLST_FULL_FASTA_PATH
    blast_outfile = blast_runner.blast_against_query(cgmlst_fasta_path)
    logging.info('Reading BLAST output file "{}"'.format(blast_outfile))
    blast_reader = BlastReader(blast_outfile)
    if blast_reader.df is None:
        logging.error('No cgMLST330 alleles found!')
        return ({'distance': 1.0,
            'genome_match': None,
            'serovar': None,
            'matching_alleles': 0,
            'subspecies': None,
            'cgmlst330_ST': None,},
                {}, )
    logging.info('Found {} cgMLST330 allele BLAST results'.format(blast_reader.df.shape[0]))


    df_cgmlst_blastn = process_cgmlst_results(blast_reader.df)

    marker_match_results = matches_to_marker_results(df_cgmlst_blastn[df_cgmlst_blastn.is_match])
    contig_blastn_records = alleles_to_retrieve(df_cgmlst_blastn)
    retrieved_marker_alleles = get_allele_sequences(blast_runner.fasta_path,
                                                    contig_blastn_records,
                                                    full=full)
    logging.info('Type retrieved_marker_alleles %s', type(retrieved_marker_alleles))
    all_marker_results = marker_match_results.copy()
    for marker, res in retrieved_marker_alleles.items():
        all_marker_results[marker] = res
    for marker in df_cgmlst_profiles.columns:
        if marker not in all_marker_results:
            all_marker_results[marker] = {'blast_result': None,
                                          'name': None,
                                          'seq': None,}
    cgmlst_results = {}
    for marker, res in all_marker_results.items():
        try:
            cgmlst_results[marker] = int(res['name'])
        except:
            logging.error('Missing cgmlst_results for %s', marker)
            logging.debug(res)

    logging.info('Calculating number of matching alleles to serovar predictive cgMLST330 profiles')
    df_relatives = find_closest_related_genome(cgmlst_results, df_cgmlst_profiles)
    genome_serovar_dict = genomes_to_serovar()
    df_relatives['serovar'] = [genome_serovar_dict[genome] for genome in df_relatives.index]
    logging.debug('Top 5 serovar predictive cgMLST profiles:\n{}'.format(df_relatives.head()))
    spp = None
    subspeciation_tuple = cgmlst_subspecies_call(df_relatives)
    if subspeciation_tuple is not None:
        spp, distance, spp_counter = subspeciation_tuple
        logging.info('Top subspecies by cgMLST is "{}" (min dist={}, Counter={})'.format(spp, distance, spp_counter))
    else:
        logging.warning('Subspeciation by cgMLST was not possible!')

    cgmlst_serovar = None
    cgmlst_matching_genome = None
    cgmlst_matching_alleles = 0
    cgmlst_distance = 1.0
    for idx, row in df_relatives.iterrows():
        cgmlst_distance = row['distance']
        cgmlst_matching_alleles = row['matching']
        cgmlst_serovar = row['serovar'] if cgmlst_distance <= 1.0 else None
        cgmlst_matching_genome = idx if cgmlst_distance <= 1.0 else None
        logging.info('Top serovar by cgMLST profile matching: "{}" with {} matching alleles, distance={:.1%}'.format(
        cgmlst_serovar,
        cgmlst_matching_alleles,
        cgmlst_distance
    ))
        break

    cgmlst_st = None
    cgmlst_markers_sorted = sorted(all_marker_results.keys())
    cgmlst_allele_names = []
    marker = None
    for marker in cgmlst_markers_sorted:
        try:
            aname = all_marker_results[marker]['name']
            if aname:
                cgmlst_allele_names.append(str(aname))
            else:
                break
        except:
            break
    if len(cgmlst_allele_names) == len(cgmlst_markers_sorted):
        cgmlst_st = allele_name('-'.join(cgmlst_allele_names))
        logging.info('cgMLST330 Sequence Type=%s', cgmlst_st)
    else:
        logging.warning('Could not compute cgMLST330 Sequence Type due to missing data (marker %s)', marker)
    return ({'distance': cgmlst_distance,
            'genome_match': cgmlst_matching_genome,
            'serovar': cgmlst_serovar,
            'matching_alleles': cgmlst_matching_alleles,
            'subspecies': spp,
            'cgmlst330_ST': cgmlst_st,},
           all_marker_results, )
