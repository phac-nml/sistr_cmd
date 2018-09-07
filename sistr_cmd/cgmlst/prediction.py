import logging
from typing import Optional, Dict
from pkg_resources import resource_filename

import attr
import pandas as pd
import numpy as np

from sistr_cmd.blast_wrapper import BlastRunner, BlastReader

N_ALLELES = 330

@attr.s
class CgMLSTPrediction:
    serovar = attr.ib(default=None)
    subspecies = attr.ib(default=None)
    distance = attr.ib(default=1.0)
    matching_genome = attr.ib(default=None)
    matching_genome_serovar = attr.ib(default=None)
    matching_genome_subspecies = attr.ib(default=None)
    top_matches = attr.ib(default=None)
    matching_alleles = attr.ib(default=0) # type: int
    n_missing_alleles = attr.ib(default=N_ALLELES) # type: int
    n_present_alleles = attr.ib(default=0) # type: int
    cgmlst_ST = attr.ib(default=None) # type: Optional[int]
    marker_results = attr.ib(default=None) # type: Optional[Dict]
    ref_genomes_with_serovar = attr.ib(default=None) # type: int
    input_genome = attr.ib(default=None) # type: str

    CGMLST_CENTROID_FASTA_PATH = resource_filename('sistr_cmd', 'data/cgmlst/cgmlst-centroid.fasta')
    CGMLST_FULL_FASTA_PATH = resource_filename('sistr_cmd', 'data/cgmlst/cgmlst-full.fasta')
    CGMLST_PROFILES_PATH = resource_filename('sistr_cmd', 'data/cgmlst/cgmlst-profiles.hdf')
    BLASTN_PIDENT_THRESHOLD = 90.0

    # cgMLST330 distance threshold for refining overall serovar prediction
    CGMLST_DISTANCE_THRESHOLD = 0.1
    CGMLST_SUBSPECIATION_DISTANCE_THRESHOLD = 0.9

    @classmethod
    def cgmlst_serovar_prediction(cls, blast_runner: BlastRunner, use_full_cgmlst_db: bool = False) -> 'CgMLSTPrediction':
        """Perform in silico cgMLST on an input genome

        Args:
            blast_runner (sistr.src.blast_wrapper.BlastRunner): blastn runner object with genome fasta initialized

        Returns:
            dict: cgMLST ref genome match, distance to closest ref genome, subspecies and serovar predictions
            dict: marker allele match results (seq, allele name, blastn results)
        """
        df_cgmlst_profiles = cls.ref_cgmlst_profiles()
        logging.info(f'Running BLAST. Querying "{blast_runner.fasta_path}" serovar predictive cgMLST330 alleles')
        cgmlst_fasta_path = cls.CGMLST_CENTROID_FASTA_PATH if not use_full_cgmlst_db else cls.CGMLST_FULL_FASTA_PATH
        blast_outfile = blast_runner.blast_against_query(cgmlst_fasta_path)
        logging.info('Reading BLAST output file "{}"'.format(blast_outfile))
        blast_reader = BlastReader(blast_outfile)
        if blast_reader.df is None:
            logging.error('No cgMLST330 alleles found!')
            obj = cls()
            obj.input_genome = blast_runner.filename()
            return  obj
        logging.info('Found {} cgMLST330 allele BLAST results'.format(blast_reader.df.shape[0]))

        df_cgmlst_blastn = CgMLSTPrediction.process_cgmlst_results(blast_reader.df)
        marker_match_results = matches_to_marker_results(df_cgmlst_blastn[df_cgmlst_blastn.is_match])
        contig_blastn_records = alleles_to_retrieve(df_cgmlst_blastn)
        retrieved_marker_alleles = get_allele_sequences(blast_runner.fasta_path,
                                                        contig_blastn_records,
                                                        using_full_cgmlst=use_full_cgmlst_db)
        logging.info('Type retrieved_marker_alleles %s', type(retrieved_marker_alleles))
        all_marker_results = get_cgmlst_marker_results(df_cgmlst_profiles, marker_match_results,
                                                       retrieved_marker_alleles)

        logging.info('Calculating number of matching alleles to serovar predictive cgMLST330 profiles')
        cgmlst_marker_results = marker_result(all_marker_results)
        df_relatives = find_relatives(cgmlst_marker_results, df_cgmlst_profiles)
        _add_ref_genome_serovar_info(df_relatives)
        logging.debug('Top 5 serovar predictive cgMLST profiles:\n{}'.format(df_relatives.head()))
        cgmlst_subspecies = subspeciation(df_relatives)

        cgmlst_distance, cgmlst_matching_alleles, cgmlst_matching_genome, cgmlst_serovar = top_serovar(df_relatives)

        cgmlst_st = cgmlst_sequence_type(all_marker_results)

        n_missing_alleles = sum([1 for marker, result in cgmlst_marker_results.items() if result is None])
        n_present_alleles = len(cgmlst_marker_results) - n_missing_alleles

        return cls(
            serovar=cgmlst_serovar,
            subspecies=cgmlst_subspecies,
            genome_match=cgmlst_matching_genome,
            distance=cgmlst_distance,
            matching_alleles=cgmlst_matching_alleles,
            n_missing_alleles=n_missing_alleles,
            n_present_alleles=n_present_alleles,
            cgmlst_ST=cgmlst_st,
            marker_results=all_marker_results,
            ref_genomes_with_serovar=(df_relatives.serovar == cgmlst_serovar).sum())

    @staticmethod
    def ref_cgmlst_profiles() -> pd.DataFrame:
        """Read in reference genome cgMLST profiles.

        Returns:
            DataFrame of cgMLST profiles for SISTR reference genomes
        """
        return pd.read_hdf(CgMLSTPrediction.CGMLST_PROFILES_PATH, key='cgmlst')

    @staticmethod
    def process_cgmlst_results(df: pd.DataFrame) -> pd.DataFrame:
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
        for x in df['qseqid']:
            marker, allele = x.split('|')
            markers.append(marker)
            alleles.append(int(allele))
        df.loc[:, 'marker'] = markers
        df.loc[:, 'allele'] = alleles
        df.loc[:, 'is_match'] = (df['coverage'] >= 1.0) & (df['pident'] >= 90.0) & ~(df['is_trunc'])
        df.loc[:, 'allele_name'] = df.apply(lambda x: hash_allele_seq(x.sseq.replace('-', '')), axis=1)
        df.loc[:, 'is_perfect'] = (df['coverage'] == 1.0) & (df['pident'] == 100.0)
        df_perf = df[df['is_perfect']]
        perf_markers = df_perf['marker'].unique()
        df.loc[:, 'has_perfect_match'] = df['marker'].isin(perf_markers)
        start_idxs, end_idxs, needs_revcomps, trunc, is_extended = extend_matches(df)
        df.loc[:, 'start_idx'] = start_idxs
        df.loc[:, 'end_idx'] = end_idxs
        df.loc[:, 'needs_revcomp'] = needs_revcomps
        df.loc[:, 'trunc'] = trunc
        df.loc[:, 'is_extended'] = is_extended
        df.loc[:, 'sseq_msa_gaps'] = np.zeros(df.shape[0], dtype=np.int64)
        df.loc[:, 'sseq_msa_p_gaps'] = np.zeros(df.shape[0], dtype=np.float64)
        df.loc[:, 'too_many_gaps'] = trunc

        return df