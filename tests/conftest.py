import pytest

from sistr.src.blast_wrapper import BlastRunner

@pytest.fixture(scope='module')
def fasta_path():
    return 'tests/00_0163.fasta'

@pytest.fixture(scope='module')
def tmp_dir():
    return '/tmp/test-blast-runner/'

@pytest.yield_fixture(scope='module')
def blast_runner(fasta_path, tmp_dir):
    _blast_runner = BlastRunner(fasta_path, tmp_dir)
    _blast_runner.prep_blast()

    yield _blast_runner

    _blast_runner.cleanup()
