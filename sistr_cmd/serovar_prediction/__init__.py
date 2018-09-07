# -*- coding: utf-8 -*-

import logging
from typing import Optional, Dict, Type, ClassVar

import attr
import pandas as pd

from sistr_cmd.blast_wrapper.reader import BlastReader
from sistr_cmd.blast_wrapper.runner import BlastRunner
from sistr_cmd.cgmlst.const import CGMLST_DISTANCE_THRESHOLD
from sistr_cmd.mash.const import MASH_DISTANCE_THRESHOLD
from sistr_cmd.serovar_prediction.constants import \
    FLJB_FASTA_PATH, \
    FLIC_FASTA_PATH, \
    H2_FLJB_SIMILARITY_GROUPS, \
    H1_FLIC_SIMILARITY_GROUPS, \
    WZY_FASTA_PATH, \
    WZX_FASTA_PATH, \
    SEROGROUP_SIMILARITY_GROUPS, \
    SEROVAR_TABLE_PATH

SPP_NAME_TO_ROMAN = {'enterica': 'I',
                     'salamae': 'II',
                     'arizonae': 'IIIa',
                     'diarizonae': 'IIIb',
                     'houtenae': 'IV',
                     'bongori': 'V',
                     'indica': 'VI'}


def get_antigen_name(qseqid: Optional[str]) -> Optional[str]:
    """Get the antigen name from the BLASTN result query ID.

    The last item delimited by `|` characters is the antigen name for all antigens (H1, H2, serogroup).

    Args:
        qseqid: BLASTN result query ID

    Returns:
        The antigen name
    """
    if qseqid:
        return qseqid.split('|')[-1]


def serovar_table() -> pd.DataFrame:
    """Get the WHO 2007 Salmonella enterica serovar table with serogroup, H1 and H2 antigen info as a Pandas DataFrame.

    Returns:
        The serovar info table
    """
    return pd.read_csv(SEROVAR_TABLE_PATH)

@attr.s
class BlastResultMixin(object):
    blast_results: pd.DataFrame = attr.ib(default=None)
    top_result: Optional[Dict] = attr.ib(default=None)
    result: Optional[str] = attr.ib(default=None)
    is_trunc: bool = attr.ib(default=False)
    is_missing: bool = attr.ib(default=False)
    is_perfect_match: bool = attr.ib(default=False)

    @classmethod
    def predict(cls: Type['BlastResultMixin'], blast_runner: BlastRunner) -> 'BlastResultMixin':
        obj = get_antigen_gene_blast_results(blast_runner, cls, cls.FASTA_PATH)
        if not obj.is_missing:
            top_result = obj.top_result
            top_result_pident = top_result['pident']
            top_result_length = top_result['length']

            if top_result_pident < cls.PIDENT_THRESHOLD:
                obj.is_missing = True
                obj.result = None
                return obj

            if top_result_length < cls.LENGTH_THRESHOLD:
                obj.is_missing = True
                obj.result = None
                return obj

            obj.result = get_antigen_name(top_result['qseqid'])
        return obj


def get_antigen_gene_blast_results(blast_runner: BlastRunner, cls: Type[BlastResultMixin], antigen_gene_fasta: str) -> BlastResultMixin:
    blast_outfile = blast_runner.blast(antigen_gene_fasta)
    blast_reader = BlastReader(blast_outfile)
    is_missing = blast_reader.is_missing
    obj = cls()
    if not is_missing:
        obj.is_missing = is_missing
        obj.blast_results = blast_reader.df
        obj.top_result = blast_reader.top_result()
        obj.is_perfect_match = blast_reader.has_perfect_match
        obj.is_trunc = blast_reader.is_trunc
    return obj

@attr.s
class WzxPrediction(BlastResultMixin):
    PIDENT_THRESHOLD: ClassVar[float] = 88.0
    LENGTH_THRESHOLD: ClassVar[int] = 400
    FASTA_PATH: ClassVar[str] = WZX_FASTA_PATH


@attr.s
class WzyPrediction(BlastResultMixin):
    PIDENT_THRESHOLD = 88.0
    LENGTH_THRESHOLD = 600
    FASTA_PATH = WZY_FASTA_PATH


@attr.s
class SerogroupPrediction():
    serogroup: Optional[str] = attr.ib(default=None)
    wzx_prediction: WzxPrediction = attr.ib(default=None)
    wzy_prediction: WzyPrediction = attr.ib(default=None)

    @classmethod
    def predict(cls: Type["SerogroupPrediction"], blast_runner: BlastRunner) -> "SerogroupPrediction":
        obj = cls(wzx_prediction=WzxPrediction.predict(blast_runner),
                   wzy_prediction=WzyPrediction.predict(blast_runner))

        if obj.wzy_prediction.is_perfect_match or obj.wzx_prediction.is_perfect_match:
            if obj.wzy_prediction.is_perfect_match:
                obj.serogroup = obj.wzy_prediction.result
            if obj.wzx_prediction.is_perfect_match:
                obj.serogroup = obj.wzx_prediction.result
        elif obj.wzx_prediction.is_missing and obj.wzy_prediction.is_missing:
            obj.serogroup = None
        elif obj.wzx_prediction.is_missing and not obj.wzy_prediction.is_missing:
            obj.serogroup = obj.wzy_prediction.result
        elif obj.wzy_prediction.is_missing and not obj.wzx_prediction.is_missing:
            obj.serogroup = obj.wzx_prediction.result
        elif obj.wzy_prediction.result == obj.wzx_prediction.result:
            obj.result = obj.wzx_prediction.result
        else:
            top_wzy_result = obj.wzy_prediction.top_result
            top_wzx_result = obj.wzx_prediction.top_result

            wzx_cov = top_wzx_result['coverage']
            wzx_pident = top_wzx_result['pident']

            wzy_cov = top_wzy_result['coverage']
            wzy_pident = top_wzy_result['pident']

            if wzx_cov >= wzy_cov and wzx_pident >= wzy_pident:
                obj.serogroup = obj.wzx_prediction.result
            elif wzx_cov < wzy_cov and wzx_pident < wzy_pident:
                obj.serogroup = obj.wzy_prediction.result
            else:
                obj.serogroup = obj.wzx_prediction.result
        return obj


@attr.s
class H1FliCPrediction(BlastResultMixin):
    """

    """
    MAX_MISMATCHES: ClassVar[int] = 25
    MAX_MISMATCHES_HIGH_QUALITY: ClassVar[int] = 5
    MIN_LENGTH: ClassVar[int] = 700
    MIN_LENGTH_HIGH_QUALITY: ClassVar[int] = 1000
    FASTA_PATH: ClassVar[str] = FLIC_FASTA_PATH

    @classmethod
    def predict(cls: Type['H1FliCPrediction'], blast_runner: BlastRunner) ->  'H1FliCPrediction':
        obj = get_antigen_gene_blast_results(blast_runner, cls, cls.FASTA_PATH) # type: H1FliCPrediction
        if not obj.is_missing:
            if not obj.h1_prediction.is_perfect_match:
                df_blast_results = obj.blast_results
                df_blast_results = df_blast_results[(df_blast_results['mismatch'] <= H1FliCPrediction.MAX_MISMATCHES)
                                                    & (df_blast_results['length'] >= H1FliCPrediction.MIN_LENGTH)]

                if df_blast_results.shape[0] == 0:
                    obj.is_missing = True
                    obj.top_result = None
                    obj.result = None
                    return obj

                df_larger_blast_results = df_blast_results[
                    (df_blast_results['mismatch'] <= 5) & (df_blast_results['length'] >= 1000)]

                if df_larger_blast_results.shape[0] > 0:
                    df_blast_results = df_larger_blast_results.sort_values(by='mismatch')
                else:
                    df_blast_results = df_blast_results.sort_values(by='bitscore', ascending=False)

                result_dict = BlastReader.df_first_row_to_dict(df_blast_results)
                result_trunc = BlastReader.is_blast_result_trunc(qstart=result_dict['qstart'],
                                                                 qend=result_dict['qend'],
                                                                 sstart=result_dict['sstart'],
                                                                 send=result_dict['send'],
                                                                 qlen=result_dict['qlen'],
                                                                 slen=result_dict['slen'])
                obj.top_result = result_dict
                obj.is_trunc = result_trunc
            else:
                obj.result = get_antigen_name(obj.top_result['qseqid'])
        return obj


@attr.s
class H2FljBPrediction(BlastResultMixin):
    """

    """
    MAX_MISMATCHES: ClassVar[int] = 50
    MAX_MISMATCHES_HIGH_QUALITY: ClassVar[int] = 5
    MIN_LENGTH: ClassVar[int] = 600
    MIN_LENGTH_HIGH_QUALITY: ClassVar[int] = 1000
    MIN_PIDENT: ClassVar[float] = 88.0
    FASTA_PATH: ClassVar[str] = FLJB_FASTA_PATH

    @classmethod
    def predict(cls: Type['H2FljBPrediction'], blast_runner: BlastRunner) -> 'H2FljBPrediction':
        obj = get_antigen_gene_blast_results(blast_runner, cls, cls.FASTA_PATH) # type: H2FljBPrediction

        if obj.is_missing:
            obj.result = '-'
        else:
            if not obj.is_perfect_match:
                top_result = obj.top_result
                match_len = top_result['length']
                pident = top_result['pident']

                # short lower %ID matches are treated as missing or '-' for H2
                if match_len <= H2FljBPrediction.MIN_LENGTH and pident < H2FljBPrediction.MIN_PIDENT:
                    obj.result = '-'
                    obj.is_missing = True
                    return obj
                # short matches are treated as missing or '-' for H2
                if match_len <= H2FljBPrediction.MIN_LENGTH and not obj.h2_prediction.is_trunc:
                    obj.result = '-'
                    obj.is_missing = True
                    return obj

                df_blast_results = obj.blast_results
                df_blast_results = df_blast_results[
                    (df_blast_results['mismatch'] <= H2FljBPrediction.MAX_MISMATCHES)
                    & (df_blast_results['length'] >= H2FljBPrediction.MIN_LENGTH)]

                if df_blast_results.shape[0] == 0:
                    obj.is_missing = True
                    obj.top_result = None
                    obj.result = '-'
                    return obj

                df_larger_blast_results = df_blast_results[
                    (df_blast_results['mismatch'] <= H2FljBPrediction.MAX_MISMATCHES_HIGH_QUALITY)
                    & (df_blast_results['length'] >= H2FljBPrediction.MIN_LENGTH_HIGH_QUALITY)]

                result_dict = BlastReader.df_first_row_to_dict(df_larger_blast_results.sort_values(by='mismatch') if df_larger_blast_results.shape[0] > 0 else df_blast_results.sort_values(by='bitscore', ascending=False))
                obj.top_result = result_dict
                obj.is_trunc = BlastReader.is_blast_result_trunc(**result_dict)
            obj.result = get_antigen_name(obj.top_result['qseqid'])
        return obj



@attr.s
class SerovarPrediction():
    genome = attr.ib(default=None)  # type: str
    serovar = attr.ib(default=None)  # type: str
    serovar_antigen = attr.ib(default=None)  # type: str
    serogroup = attr.ib(default=None)  # type: str
    serogroup_prediction = attr.ib(default=None)  # type: SerogroupPrediction
    h1 = attr.ib(default=None)  # type: str
    h1_flic_prediction = attr.ib(default=None)  # type: H1FliCPrediction
    h2 = attr.ib(default=None)  # type: str
    h2_fljb_prediction = attr.ib(default=None)  # type: H2FljBPrediction
    serovar_cgmlst = attr.ib(default=None)  # type: str
    cgmlst_distance = attr.ib(default=1.0)  # type: float
    cgmlst_matching_alleles = attr.ib(default=0)  # type: int
    cgmlst_genome_match = attr.ib(default=None)  # type: str
    cgmlst_subspecies = attr.ib(default=None)  # type: str



    @classmethod
    def predict(cls: Type['SerovarPrediction'], blast_runner: BlastRunner) -> 'SerovarPrediction':
        obj = cls()
        obj.h1_flic_prediction = H1FliCPrediction.predict(blast_runner)
        obj.h1 = obj.h1_flic_prediction.result
        obj.h2_fljb_prediction = H2FljBPrediction.predict(blast_runner)
        obj.h2 = obj.h2_fljb_prediction.result
        obj.serogroup_prediction = SerogroupPrediction.predict(blast_runner)
        obj.serogroup = obj.serogroup_prediction.serogroup

        return obj

    @staticmethod
    def get_serovar(df, sg, h1, h2, spp):
        h2_is_missing = '-' in h2
        b_sg = df['Serogroup'].isin(sg)
        b_h1 = df['H1'].isin(h1)
        if h2_is_missing:
            b_h2 = df['can_h2_be_missing']
        else:
            b_h2 = df['H2'].isin(h2)

        if spp is not None:
            b_spp = df['subspecies'] == spp
        else:
            b_spp = b_sg
        df_prediction = df[(b_spp & b_sg & b_h1 & b_h2)]
        logging.debug('Serovar prediction for %s %s:%s:%s is %s', spp, sg, h1, h2, list(df_prediction['Serovar']))
        if df_prediction.shape[0] > 0:
            return '|'.join(list(df_prediction['Serovar']))

    def predict_serovar_from_antigen_blast(self):
        if not self.serogroup or not self.h2 or not self.h1:
            self.predict_antigens()

        df = serovar_table()
        sg = self.serogroup
        h1 = self.h1
        h2 = self.h2

        # no antigen results then serovar == '-:-:-'
        if sg is None \
                and h1 is None \
                and h2 == '-':
            self.serovar = '-:-:-'
            return self.serovar

        for sg_groups in SEROGROUP_SIMILARITY_GROUPS:
            if sg in sg_groups:
                sg = sg_groups
                break
        if sg is None:
            sg = list(df['Serogroup'].unique())
        if not isinstance(sg, list):
            sg = [sg]

        for h1_groups in H1_FLIC_SIMILARITY_GROUPS:
            if h1 in h1_groups:
                h1 = h1_groups
                break
        if h1 is None:
            h1 = list(df['H1'].unique())
        if not isinstance(h1, list):
            h1 = [h1]

        for h2_groups in H2_FLJB_SIMILARITY_GROUPS:
            if h2 in h2_groups:
                h2 = h2_groups
                break

        if not isinstance(h2, list):
            h2 = [h2]

        self.serovar = SerovarPrediction.get_serovar(df, sg, h1, h2, self.subspecies)
        if self.serovar is None:
            try:
                spp_roman = SPP_NAME_TO_ROMAN[self.subspecies]
            except:
                spp_roman = None
            from collections import Counter
            c = Counter(df.O_antigen[df.Serogroup.isin(sg)])
            o_antigen = c.most_common()[0][0]
            h1_first = h1[0]
            h2_first = h2[0]
            if spp_roman:
                self.serovar = '{} {}:{}:{}'.format(spp_roman, o_antigen, h1_first, h2_first)
            else:
                self.serovar = '{}:{}:{}'.format(o_antigen, h1_first, h2_first)
        return self.serovar

    def overall_serovar_call(self):
        """
        Predict serovar from cgMLST cluster membership analysis and antigen BLAST results.
        SerovarPrediction object is assigned H1, H2 and Serogroup from the antigen BLAST results.
        Antigen BLAST results will predict a particular serovar or list of serovars, however,
        the cgMLST membership may be able to help narrow down the list of potential serovars.

        Notes:
            If the cgMLST predicted serovar is within the list of antigen BLAST predicted serovars,
            then the serovar is assigned the cgMLST predicted serovar.


            If all antigens are found, but an antigen serovar is not found then the serovar is assigned
            a pseudo-antigenic formula (Serogroup:H1:H2), otherwise the serovar is assigned the cgMLST prediction.


            If the antigen predicted serovar does not match the cgMLST predicted serovar,

            - the serovar is the cgMLST serovar if the cgMLST cluster level is <= 0.1 (10% or less)
            - otherwise, the serovar is antigen predicted serovar(s)

        Args:
            prediction (SerovarPrediction): Serovar prediction results (antigen+cgMLST[+Mash])
            antigen_results (SerovarPredictor): Antigen search results

        Returns:
            SerovarPrediction: Serovar prediction results with overall prediction from antigen + cgMLST
        """

        h1 = self.h1 or self.h1_flic_prediction.result
        h2 = self.h2 or self.h2_fljb_prediction.result
        sg = self.serogroup or self.serogroup_prediction.serogroup
        spp = self.cgmlst_subspecies or self.__dict__['mash_subspecies'] if 'mash_match' in self.__dict__ else None

        self.serovar_antigen = self.serovar

        cgmlst_serovar = self.serovar_cgmlst
        cgmlst_distance = float(self.cgmlst_distance)

        null_result = '-:-:-'

        try:
            spp_roman = SPP_NAME_TO_ROMAN[spp]
        except:
            spp_roman = None

        is_antigen_null = lambda x: (x is None or x == '' or x == '-')

        if self.serovar is None:
            if is_antigen_null(sg) and is_antigen_null(h1) and is_antigen_null(h2):
                if spp_roman is not None:
                    self.serovar = '{} {}:{}:{}'.format(spp_roman, sg, h1, h2)
                else:
                    self.serovar = '{}:{}:{}'.format(spp_roman, sg, h1, h2)
            elif cgmlst_serovar is not None and cgmlst_distance <= CGMLST_DISTANCE_THRESHOLD:
                self.serovar = cgmlst_serovar
            else:
                self.serovar = null_result
                if 'mash_match' in self.__dict__:
                    spd = self.__dict__
                    mash_dist = float(spd['mash_distance'])
                    if mash_dist <= MASH_DISTANCE_THRESHOLD:
                        self.serovar = spd['mash_serovar']
        else:
            serovars_from_antigen = self.serovar_antigen.split('|')
            if not isinstance(serovars_from_antigen, list):
                serovars_from_antigen = [serovars_from_antigen]
            if cgmlst_serovar is not None:
                if cgmlst_serovar in serovars_from_antigen:
                    self.serovar = cgmlst_serovar
                else:
                    if float(cgmlst_distance) <= CGMLST_DISTANCE_THRESHOLD:
                        self.serovar = cgmlst_serovar
            elif 'mash_match' in self.__dict__:
                spd = self.__dict__
                mash_serovar = spd['mash_serovar']
                mash_dist = float(spd['mash_distance'])
                if mash_serovar in serovars_from_antigen:
                    self.serovar = mash_serovar
                else:
                    if mash_dist <= MASH_DISTANCE_THRESHOLD:
                        self.serovar = mash_serovar

            if self.serovar is None:
                self.serovar = self.serovar_antigen

        if self.h1 is None:
            self.h1 = '-'
        if self.h2 is None:
            self.h2 = '-'
        if self.serogroup is None:
            self.serogroup = '-'
        if self.serovar_antigen is None:
            if spp_roman is not None:
                self.serovar_antigen = '{} -:-:-'.format(spp_roman)
            else:
                self.serovar_antigen = '-:-:-'
        if self.serovar is None:
            self.serovar = self.serovar_antigen
        return self

