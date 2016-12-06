from sistr.sistr_cmd import run_mash


def test_run_mash(fasta_path):
    out = run_mash(fasta_path)
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
