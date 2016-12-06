import logging
from sistr.src.parsers import parse_fasta
from sistr.src.qc.constants import ERR_MISSING_CGMLST_MARKERS, WARN_MISSING_CGMLST_MARKERS, SALMONELLA_GENOME_SIZE_MBP
from sistr.src.serovar_prediction import CGMLST_DISTANCE_THRESHOLD


def qc(fasta_path, cgmlst_results, prediction):
    qc_status = 'PASS'
    qc_msgs = []
    genome_size = sum([len(s) for h, s in parse_fasta(fasta_path)])
    lb_salm_gsize, ub_salm_gsize = SALMONELLA_GENOME_SIZE_MBP
    is_gsize_acceptable = (genome_size >= lb_salm_gsize and genome_size <= ub_salm_gsize)
    logging.info('Genome size=%s (within gsize thresholds? %s)', genome_size, is_gsize_acceptable)
    if not is_gsize_acceptable:
        qc_status = 'WARNING'
        qc_msgs.append('WARNING: Input genome size ({} bp) not within expected range of {}-{} (bp) for Salmonella'.format(genome_size, lb_salm_gsize, ub_salm_gsize))
    if len(cgmlst_results) == 0:
        missing_cgmlst_count = 330
    else:
        missing_cgmlst_count = 0
        for marker, results in cgmlst_results.items():
            if results['name'] is None:
                missing_cgmlst_count += 1

    if missing_cgmlst_count >= ERR_MISSING_CGMLST_MARKERS:
        qc_status = 'FAIL'
        qc_msgs.append('FAIL: Large number of cgMLST330 loci missing (n={} > {})'.format(
            missing_cgmlst_count,
            ERR_MISSING_CGMLST_MARKERS))
    elif missing_cgmlst_count >= WARN_MISSING_CGMLST_MARKERS:
        qc_status = 'WARNING'
        qc_msgs.append('FAIL: Moderate number of cgMLST330 loci missing (n={} > {})'.format(
            missing_cgmlst_count,
            WARN_MISSING_CGMLST_MARKERS))

    matching_cgmlst_alleles_threshold = (1.0 - CGMLST_DISTANCE_THRESHOLD) * 330
    if prediction.cgmlst_matching_alleles < matching_cgmlst_alleles_threshold:
        if qc_status != 'FAIL':
            qc_status = 'WARNING'
        qc_msgs.append('WARNING: Only matched {} cgMLST330 loci. Min threshold for confident serovar prediction from cgMLST is {}'.format(prediction.cgmlst_matching_alleles, matching_cgmlst_alleles_threshold))

    if prediction.h1 is None or prediction.h1 == '-':
        if qc_status != 'FAIL':
            qc_status = 'WARNING'
        qc_msgs.append('WARNING: H1 antigen gene (fliC) missing. Cannot determine H1 antigen. Cannot accurately predict serovar from antigen genes.')

    if prediction.serogroup is None \
        or prediction.serogroup == '-' \
        or prediction.serogroup == '':
        qc_status = 'FAIL'
        qc_msgs.append('FAIL: Wzx/Wzy genes missing. Cannot determine O-antigen group/serogroup. Cannot accurately predict serovar from antigen genes.')

    qc_msgs.sort()
    return qc_status, qc_msgs
