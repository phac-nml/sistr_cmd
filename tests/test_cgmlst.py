'''
Test running cgMLST

This is a longer test (>20 sec) due to running blastn of all alleles of all 330 cgMLST loci against a subject genome.
'''
from __future__ import print_function

from sistr.sistr_cmd import run_cgmlst
from sistr.src.cgmlst import allele_name


def test_hash_allele():
    seq = 'ATGC'
    chksum = allele_name(seq)
    assert chksum == 968142784

    seq_unicode = u'ATGC'
    chksum = allele_name(seq_unicode)
    assert chksum == 968142784


def test_run_cgmlst(blast_runner):
    """
    Expected perfect match against test genome 00_0163 since it is one of the ref genomes.
    """
    cgmlst_serovar_pred, marker_res = run_cgmlst(blast_runner)
    assert isinstance(cgmlst_serovar_pred, dict)
    spp = cgmlst_serovar_pred['subspecies']
    serovar = cgmlst_serovar_pred['serovar']
    genome_match = cgmlst_serovar_pred['genome_match']
    matching_alleles = cgmlst_serovar_pred['matching_alleles']
    distance = cgmlst_serovar_pred['distance']
    assert matching_alleles == 330
    assert distance == 0.0
    assert spp == 'salamae'
    assert serovar == 'II 58:l,z13,z28:z6'
    assert isinstance(marker_res, dict)
    # 330 marker results for 330 cgMLST loci
    assert len(marker_res) == 330
