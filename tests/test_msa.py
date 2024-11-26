from __future__ import print_function
from sistr.src.cgmlst.msa import msa_mafft, msa_ref_vs_novel, parse_aln_out


def test_parse_aln_out_string():
    test_aln_str = """
>1
atgc-
>2
-tgca
>3
atgca
"""
    aln_dict = {h:s for h,s in parse_aln_out(test_aln_str)}
    assert "1" in aln_dict
    assert "2" in aln_dict
    assert "3" in aln_dict
    assert aln_dict['1'] == 'atgc-'
    assert aln_dict['2'] == '-tgca'
    assert aln_dict['3'] == 'atgca'
    print(aln_dict)


def test_msa_mafft_str():
    test_input_str = """
>1
atgc
>2
tgca
>3
atgca
"""
    aln_dict = msa_mafft(test_input_str)
    assert isinstance(aln_dict, dict)
    assert "1" in aln_dict
    assert "2" in aln_dict
    assert "3" in aln_dict
    assert aln_dict['1'] == 'atgc-'
    assert aln_dict['2'] == '-tgca'
    assert aln_dict['3'] == 'atgca'
    print(aln_dict)


def test_msa_ref_vs_novel():
    ref_nt = 'atgtgc'
    novel_nt = 'atgcatgc'
    ref_msa, novel_msa = msa_ref_vs_novel(ref_nt, novel_nt)
    assert novel_msa == novel_nt
    assert ref_msa == 'atgtgc--' 
    print(ref_nt, novel_nt)
    print(ref_msa, novel_msa)
