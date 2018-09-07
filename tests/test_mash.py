# -*- coding: utf-8 -*-

from sistr_cmd.mash import mash_serovar_prediction


def test_run_mash(fasta_path):
    out = mash_serovar_prediction(fasta_path)
    assert 'mash_match' in out
    assert 'mash_genome' in out
    assert 'mash_distance' in out
    assert 'mash_subspecies' in out
    spp = out['mash_subspecies']
    serovar = out['mash_serovar']
    genome = out['mash_genome']
    assert spp == 'salamae'
    assert serovar == 'II 58:l,z13,z28:z6'
    assert genome == '00_0163'
