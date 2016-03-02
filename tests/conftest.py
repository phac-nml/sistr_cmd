
from src.blast_wrapper import BlastRunner

__author__ = 'peter'

import pytest



@pytest.yield_fixture(scope='module')
def blast_runner():
    _blast_runner = BlastRunner()

    yield _blast_runner

    _blast_runner.cleanup()
