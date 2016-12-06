import os
import pandas as pd

from sistr.src.blast_wrapper import BlastRunner, BLAST_TABLE_COLS, BlastReader
from sistr.src.serovar_prediction import WZX_FASTA_PATH, get_antigen_name, SerovarPredictor


def test_BlastRunner(fasta_path, tmp_dir):
    br = BlastRunner(fasta_path, tmp_dir)
    br.prep_blast()
    assert os.path.exists(tmp_dir)
    fasta_copy_path = os.path.join(tmp_dir, os.path.basename(fasta_path))
    nin_path = fasta_copy_path + '.nin'
    assert os.path.exists(fasta_copy_path)
    assert os.path.exists(nin_path)

    blast_outfile = br.run_blast(WZX_FASTA_PATH)
    assert os.path.exists(blast_outfile)
    # Pandas should be able to parse in the table
    df = pd.read_table(blast_outfile)
    # the number of columns needs to be the expected number based on blastn outfmt specs
    assert len(df.columns) == len(BLAST_TABLE_COLS)

    br.cleanup()
    assert not os.path.exists(tmp_dir)


def test_BlastReader(blast_runner):
    blast_outfile = blast_runner.run_blast(WZX_FASTA_PATH)
    blast_reader = BlastReader(blast_outfile)
    top_result = blast_reader.top_result()
    assert get_antigen_name(top_result['qseqid']) == 'O58'


def test_BlastReader_is_blast_result_trunc():

    # not truncated; match found in the middle of the subject sequence
    assert not BlastReader.is_blast_result_trunc(qstart=1,
                                                 qend=100,
                                                 sstart=101,
                                                 send=200,
                                                 qlen=100,
                                                 slen=1000)

    # not truncated; shorter match (-10bp) found in the middle of the subject
    # sequence
    assert not BlastReader.is_blast_result_trunc(qstart=1,
                                                 qend=90,
                                                 sstart=101,
                                                 send=190,
                                                 qlen=100,
                                                 slen=1000)

    # not truncated; shorter match (-20bp) found in the middle of the subject
    # sequence
    assert not BlastReader.is_blast_result_trunc(qstart=1,
                                                 qend=80,
                                                 sstart=101,
                                                 send=180,
                                                 qlen=100,
                                                 slen=1000)


    # truncated at the start of the subject
    assert BlastReader.is_blast_result_trunc(qstart=51,
                                             qend=100,
                                             sstart=1,
                                             send=50,
                                             qlen=100,
                                             slen=1000)

    # truncated at the end of the subject
    assert BlastReader.is_blast_result_trunc(qstart=51,
                                             qend=100,
                                             sstart=951,
                                             send=1000,
                                             qlen=100,
                                             slen=1000)


def test_SerovarPredictor(blast_runner):
    sp = SerovarPredictor(blast_runner, "salamae")
    sp.predict_serovar_from_antigen_blast()
    assert sp.serogroup == 'O58'
    assert sp.h1 == 'l,z13,z28'
    assert sp.h2 == 'z6'
    assert sp.serovar == 'II 58:l,z13,z28:z6'
