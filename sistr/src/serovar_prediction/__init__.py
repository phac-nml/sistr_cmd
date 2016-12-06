import logging
import pandas as pd

from sistr.src.blast_wrapper import BlastReader
from sistr.src.serovar_prediction.constants import \
    FLJB_FASTA_PATH, \
    FLIC_FASTA_PATH, \
    H2_FLJB_SIMILARITY_GROUPS, \
    H1_FLIC_SIMILARITY_GROUPS, \
    WZY_FASTA_PATH, \
    WZX_FASTA_PATH, \
    SEROGROUP_SIMILARITY_GROUPS, \
    SEROVAR_TABLE_PATH, \
    CGMLST_DISTANCE_THRESHOLD, MASH_DISTANCE_THRESHOLD

spp_name_to_roman = {'enterica': 'I',
                     'salamae': 'II',
                     'arizonae': 'IIIa',
                     'diarizonae': 'IIIb',
                     'houtenae': 'IV',
                     'bongori': 'V',
                     'indica': 'VI'}


class BlastResultMixin(object):
    blast_results = None
    top_result = None
    is_trunc = False
    is_missing = False
    is_perfect_match = False


class WzxPrediction(BlastResultMixin):
    serogroup = None


class WzyPrediction(BlastResultMixin):
    serogroup = None


class SerogroupPrediction():
    serogroup = None

    wzx_prediction = None
    wzy_prediction = None


class H1FliCPrediction(BlastResultMixin):
    h1 = None


class H2FljBPrediction(BlastResultMixin):
    h2 = None


class SerovarPrediction():
    genome = None
    serovar = None
    serovar_cgmlst = None
    cgmlst_distance = 1.0
    cgmlst_matching_alleles = 0
    cgmlst_genome_match = None
    cgmlst_subspecies = None

    serovar_antigen = None
    serogroup = None
    serogroup_prediction = None
    h1 = None
    h1_flic_prediction = None
    h2 = None
    h2_fljb_prediction = None


def get_antigen_name(qseqid):
    """
    Get the antigen name from the BLASTN result query ID.
    The last item delimited by | characters is the antigen name for all
    antigens (H1, H2, serogroup)

    @type qseqid: str
    @param qseqid: BLASTN result query ID
    @return: antigen name
    """
    if qseqid:
        return qseqid.split('|')[-1]


def serovar_table():
    """
    Get the WHO 2007 Salmonella enterica serovar table with serogroup, H1 and
    H2 antigen info as a Pandas DataFrame.
    @return: Pandas DataFrame of serovar table
    """
    return pd.read_csv(SEROVAR_TABLE_PATH)


class BlastAntigenGeneMixin:
    def get_antigen_gene_blast_results(self, model_obj, antigen_gene_fasta):
        blast_outfile = self.blast_runner.blast_against_query(antigen_gene_fasta)
        blast_reader = BlastReader(blast_outfile)
        is_missing = blast_reader.is_missing
        model_obj.is_missing = is_missing
        if not is_missing:
            model_obj.blast_results = blast_reader.df_dict()

            model_obj.top_result = blast_reader.top_result()
            model_obj.is_perfect_match = blast_reader.is_perfect_match
            model_obj.is_trunc = blast_reader.is_trunc

        return model_obj


class SerogroupPredictor(BlastAntigenGeneMixin):
    def __init__(self, blast_runner):
        """
        SerogroupPredictor takes a initialized BlastRunner object where the
        temp work folder has been created and the genome fasta has been copied
        over and a BLASTN DB has been made for it.
        This class then queries wzx and wzy against the genome using the
        BlastRunner to get the wzx and wzy serogroup predictions.

        @type blast_runner: app.blast_wrapper.BlastRunner
        @param blast_runner: Initialized BlastRunner object
        """
        self.blast_runner = blast_runner

        self.wzx_prediction = WzxPrediction()
        self.wzy_prediction = WzyPrediction()
        self.serogroup_prediction = SerogroupPrediction()


    def search_for_wzx(self):
        self.wzx_prediction = self.get_antigen_gene_blast_results(self.wzx_prediction, WZX_FASTA_PATH)
        if not self.wzx_prediction.is_missing:
            top_result = self.wzx_prediction.top_result
            top_result_pident = top_result['pident']
            top_result_length = top_result['length']

            if top_result_pident < 88.0:
                self.wzx_prediction.is_missing = True
                self.wzx_prediction.serogroup = None
                return

            if top_result_length < 400:
                self.wzx_prediction.is_missing = True
                self.wzx_prediction.serogroup = None
                return

            self.wzx_prediction.serogroup = get_antigen_name(top_result['qseqid'])


    def search_for_wzy(self):
        self.wzy_prediction = self.get_antigen_gene_blast_results(self.wzy_prediction, WZY_FASTA_PATH)
        if not self.wzy_prediction.is_missing:
            top_result = self.wzy_prediction.top_result
            top_result_pident = top_result['pident']
            top_result_length = top_result['length']

            if top_result_pident < 88.0:
                self.wzy_prediction.is_missing = True
                self.wzy_prediction.serogroup = None
                return

            if top_result_length < 600:
                self.wzy_prediction.is_missing = True
                self.wzy_prediction.serogroup = None
                return

            self.wzy_prediction.serogroup = get_antigen_name(top_result['qseqid'])


    def predict(self):
        self.search_for_wzx()
        self.search_for_wzy()
        self.serogroup_prediction.wzx_prediction = self.wzx_prediction
        self.serogroup_prediction.wzy_prediction = self.wzy_prediction
        if self.wzx_prediction.is_perfect_match:
            self.serogroup_prediction.serogroup = self.wzx_prediction.serogroup
        if self.wzy_prediction.is_perfect_match:
            self.serogroup_prediction.serogroup = self.wzy_prediction.serogroup
        if self.wzy_prediction.is_perfect_match or self.wzx_prediction.is_perfect_match:
            return

        if self.wzx_prediction.is_missing and self.wzy_prediction.is_missing:
            return
        if self.wzx_prediction.is_missing:
            self.serogroup_prediction.serogroup = self.wzy_prediction.serogroup
            return
        if self.wzy_prediction.is_missing:
            self.serogroup_prediction.serogroup = self.wzx_prediction.serogroup
            return

        if self.wzy_prediction.serogroup == self.wzx_prediction.serogroup:
            self.serogroup_prediction.serogroup = self.wzx_prediction.serogroup
            return

        top_wzy_result = self.wzy_prediction.top_result
        top_wzx_result = self.wzx_prediction.top_result

        wzx_cov = top_wzx_result['coverage']
        wzx_pident = top_wzx_result['pident']

        wzy_cov = top_wzy_result['coverage']
        wzy_pident = top_wzy_result['pident']

        if wzx_cov >= wzy_cov and wzx_pident >= wzy_pident:
            self.serogroup_prediction.serogroup = self.wzx_prediction.serogroup
        elif wzx_cov < wzy_cov and wzx_pident < wzy_pident:
            self.serogroup_prediction.serogroup = self.wzy_prediction.serogroup
        else:
            self.serogroup_prediction.serogroup = self.wzx_prediction.serogroup


class H1Predictor(BlastAntigenGeneMixin):
    def __init__(self, blast_runner):
        self.blast_runner = blast_runner
        self.h1_prediction = H1FliCPrediction()

    def predict(self):
        self.h1_prediction = self.get_antigen_gene_blast_results(self.h1_prediction, FLIC_FASTA_PATH)
        if not self.h1_prediction.is_missing:
            if not self.h1_prediction.is_perfect_match:
                df_blast_results = pd.DataFrame(self.h1_prediction.blast_results)
                df_blast_results = df_blast_results[
                    (df_blast_results['mismatch'] <= 25) & (df_blast_results['length'] >= 700)]

                if df_blast_results.shape[0] == 0:
                    self.h1_prediction.is_missing = True
                    self.h1_prediction.top_result = None
                    self.h1_prediction.h1 = None
                    return

                df_blast_results_over1000 = df_blast_results[
                    (df_blast_results['mismatch'] <= 5) & (df_blast_results['length'] >= 1000)]

                if df_blast_results_over1000.shape[0] > 0:
                    df_blast_results = df_blast_results_over1000.sort_values(by='mismatch')
                else:
                    df_blast_results = df_blast_results.sort_values(by='bitscore', ascending=False)

                result_dict = BlastReader.df_first_row_to_dict(df_blast_results)
                result_trunc = BlastReader.is_blast_result_trunc(qstart=result_dict['qstart'],
                                                                 qend=result_dict['qend'],
                                                                 sstart=result_dict['sstart'],
                                                                 send=result_dict['send'],
                                                                 qlen=result_dict['qlen'],
                                                                 slen=result_dict['slen'])
                self.h1_prediction.top_result = result_dict
                self.h1_prediction.is_trunc = result_trunc
            self.h1_prediction.h1 = get_antigen_name(self.h1_prediction.top_result['qseqid'])


class H2Predictor(BlastAntigenGeneMixin):
    def __init__(self, blast_runner):
        self.blast_runner = blast_runner
        self.h2_prediction = H2FljBPrediction()

    def predict(self):
        self.h2_prediction = self.get_antigen_gene_blast_results(self.h2_prediction, FLJB_FASTA_PATH)
        if not self.h2_prediction.is_missing:
            if not self.h2_prediction.is_perfect_match:
                top_result = self.h2_prediction.top_result
                match_len = top_result['length']
                pident = top_result['pident']

                # short lower %ID matches are treated as missing or '-' for H2
                if match_len <= 600 and pident < 88.0:
                    self.h2_prediction.h2 = '-'
                    self.h2_prediction.is_missing = True
                    return

                if match_len <= 600 and not self.h2_prediction.is_trunc:
                    self.h2_prediction.h2 = '-'
                    self.h2_prediction.is_missing = True
                    return

                df_blast_results = pd.DataFrame(self.h2_prediction.blast_results)
                df_blast_results = df_blast_results[
                    (df_blast_results['mismatch'] <= 50) & (df_blast_results['length'] >= 700)]

                if df_blast_results.shape[0] == 0:
                    self.h2_prediction.is_missing = True
                    self.h2_prediction.top_result = None
                    self.h2_prediction.h2 = '-'
                    return

                df_blast_results_over1000 = df_blast_results[
                    (df_blast_results['mismatch'] <= 5) & (df_blast_results['length'] >= 1000)]

                if df_blast_results_over1000.shape[0] > 0:
                    df_blast_results = df_blast_results_over1000.sort_values(by='mismatch')
                else:
                    df_blast_results = df_blast_results.sort_values(by='bitscore', ascending=False)

                result_dict = BlastReader.df_first_row_to_dict(df_blast_results)
                result_trunc = BlastReader.is_blast_result_trunc(qstart=result_dict['qstart'],
                                                                 qend=result_dict['qend'],
                                                                 sstart=result_dict['sstart'],
                                                                 send=result_dict['send'],
                                                                 qlen=result_dict['qlen'],
                                                                 slen=result_dict['slen'])
                self.h2_prediction.top_result = result_dict
                self.h2_prediction.is_trunc = result_trunc
            self.h2_prediction.h2 = get_antigen_name(self.h2_prediction.top_result['qseqid'])

        if self.h2_prediction.is_missing:
            self.h2_prediction.h2 = '-'


class SerovarPredictor:
    serogroup = None
    h1 = None
    h2 = None
    serovar = None
    subspecies = None

    def __init__(self, blast_runner, subspecies):
        """

        """
        self.blast_runner = blast_runner
        self.subspecies = subspecies
        self.serogroup_predictor = SerogroupPredictor(self.blast_runner)
        self.h1_predictor = H1Predictor(self.blast_runner)
        self.h2_predictor = H2Predictor(self.blast_runner)

    def predict_antigens(self):
        self.h1_predictor.predict()
        self.h2_predictor.predict()
        self.serogroup_predictor.predict()
        self.h1 = self.h1_predictor.h1_prediction.h1
        self.h2 = self.h2_predictor.h2_prediction.h2
        self.serogroup = self.serogroup_predictor.serogroup_prediction.serogroup
        return self.serogroup, self.h1, self.h2

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

        self.serovar = SerovarPredictor.get_serovar(df, sg, h1, h2, self.subspecies)
        if self.serovar is None:
            try:
                spp_roman = spp_name_to_roman[self.subspecies]
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

    def get_serovar_prediction(self):
        serovar_pred = SerovarPrediction()
        sg_pred = self.serogroup_predictor.serogroup_prediction

        h1_pred = self.h1_predictor.h1_prediction
        h2_pred = self.h2_predictor.h2_prediction

        serovar_pred.serogroup_prediction = sg_pred
        serovar_pred.serogroup = self.serogroup

        serovar_pred.h1_flic_prediction = h1_pred
        serovar_pred.h1 = self.h1

        serovar_pred.h2_fljb_prediction = h2_pred
        serovar_pred.h2 = self.h2

        return serovar_pred


def overall_serovar_call(serovar_prediction, antigen_predictor):
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
        serovar_prediction (src.serovar_prediction.SerovarPrediction): Serovar prediction results (antigen+cgMLST[+Mash])
        antigen_predictor (src.serovar_prediction.SerovarPredictor): Antigen search results

    Returns:
        src.serovar_prediction.SerovarPrediction: Serovar prediction results with overall prediction from antigen + cgMLST
    """
    assert isinstance(serovar_prediction, SerovarPrediction)
    assert isinstance(antigen_predictor, SerovarPredictor)

    h1 = antigen_predictor.h1
    h2 = antigen_predictor.h2
    sg = antigen_predictor.serogroup
    spp = serovar_prediction.cgmlst_subspecies
    if spp is None:
        if 'mash_match' in serovar_prediction.__dict__:
            spp = serovar_prediction.__dict__['mash_subspecies']

    serovar_prediction.serovar_antigen = antigen_predictor.serovar

    cgmlst_serovar = serovar_prediction.serovar_cgmlst
    cgmlst_distance = float(serovar_prediction.cgmlst_distance)

    null_result = '-:-:-'


    try:
        spp_roman = spp_name_to_roman[spp]
    except:
        spp_roman = None

    is_antigen_null = lambda x: (x is None or x == '' or x == '-')


    if antigen_predictor.serovar is None:
        if is_antigen_null(sg) and is_antigen_null(h1) and is_antigen_null(h2):
            if spp_roman is not None:
                serovar_prediction.serovar = '{} {}:{}:{}'.format(spp_roman, sg, h1, h2)
            else:
                serovar_prediction.serovar = '{}:{}:{}'.format(spp_roman, sg, h1, h2)
        elif cgmlst_serovar is not None and cgmlst_distance <= CGMLST_DISTANCE_THRESHOLD:
            serovar_prediction.serovar = cgmlst_serovar
        else:
            serovar_prediction.serovar = null_result
            if 'mash_match' in serovar_prediction.__dict__:
                spd = serovar_prediction.__dict__
                mash_dist = float(spd['mash_distance'])
                if mash_dist <= MASH_DISTANCE_THRESHOLD:
                    serovar_prediction.serovar = spd['mash_serovar']
    else:
        serovars_from_antigen = antigen_predictor.serovar.split('|')
        if not isinstance(serovars_from_antigen, list):
            serovars_from_antigen = [serovars_from_antigen]
        if cgmlst_serovar is not None:
            if cgmlst_serovar in serovars_from_antigen:
                serovar_prediction.serovar = cgmlst_serovar
            else:
                if float(cgmlst_distance) <= CGMLST_DISTANCE_THRESHOLD:
                    serovar_prediction.serovar = cgmlst_serovar
        elif 'mash_match' in serovar_prediction.__dict__:
            spd = serovar_prediction.__dict__
            mash_serovar = spd['mash_serovar']
            mash_dist = float(spd['mash_distance'])
            if mash_serovar in serovars_from_antigen:
                serovar_prediction.serovar = mash_serovar
            else:
                if mash_dist <= MASH_DISTANCE_THRESHOLD:
                    serovar_prediction.serovar = mash_serovar

        if serovar_prediction.serovar is None:
            serovar_prediction.serovar = serovar_prediction.serovar_antigen

    if serovar_prediction.h1 is None:
        serovar_prediction.h1 = '-'
    if serovar_prediction.h2 is None:
        serovar_prediction.h2 = '-'
    if serovar_prediction.serogroup is None:
        serovar_prediction.serogroup = '-'
    if serovar_prediction.serovar_antigen is None:
        if spp_roman is not None:
            serovar_prediction.serovar_antigen = '{} -:-:-'.format(spp_roman)
        else:
            serovar_prediction.serovar_antigen = '-:-:-'
    if serovar_prediction.serovar is None:
        serovar_prediction.serovar = serovar_prediction.serovar_antigen
    return serovar_prediction