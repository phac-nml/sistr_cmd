# -*- coding: utf-8 -*-

import sys
import io
from distutils.core import setup
from glob import glob
from os.path import basename, splitext, dirname, join
from setuptools import find_packages
from setuptools.command.test import test as TestCommand

from sistr_cmd.version import __version__

classifiers = """
Development Status :: 4 - Beta
Environment :: Console
License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)
Intended Audience :: Science/Research
Topic :: Scientific/Engineering :: Bio-Informatics
Programming Language :: Python :: 3.6
Programming Language :: Python :: Implementation :: CPython
Operating System :: POSIX :: Linux
""".strip().split('\n')


def read(*names, **kwargs):
    return io.open(
        join(dirname(__file__), *names),
        encoding=kwargs.get('encoding', 'utf8')
    ).read()


# The following taken from https://docs.pytest.org/en/latest/goodpractices.html#manual-integration
class PyTest(TestCommand):
    """Required to be able to run tests with `python setup.py test`"""
    user_options = [("pytest-args=", "a", "Arguments to pass to pytest")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = ""

    def run_tests(self):
        import shlex
        # import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(shlex.split(self.pytest_args))
        sys.exit(errno)


setup(
    name='sistr_cmd',
    version=__version__,
    packages=find_packages(exclude=['tests']),
    url='https://github.com/peterk87/sistr_cmd',
    license='GPLv3',
    author='Peter Kruczkiewicz',
    author_email='peter.kruczkiewicz@gmail.com',
    description=('Serovar predictions from Salmonella whole-genome sequence assemblies by determination of antigen gene'
                 'and cgMLST gene alleles using BLAST. Mash MinHash can also be used for serovar prediction.'),
    long_description=read('README.rst'),
    keywords=[
        'Salmonella',
        'serotype',
        'serovar',
        'serotyping',
        'genotyping',
        'cgMLST',
        'BLAST',
        'Mash',
        'MinHash'
    ],
    classifiers=classifiers,
    package_dir={'sistr_cmd': 'sistr_cmd'},
    package_data={'sistr_cmd': ['data/*.msh',
                                'data/*.csv',
                                'data/*.txt',
                                'data/antigens/*.fasta',
                                'data/cgmlst/*.fasta',
                                'data/cgmlst/*.txt',
                                'data/cgmlst/*.csv',
                                'data/cgmlst/*.hdf'
                                ]},
    py_modules=[splitext(basename(path))[0] for path in glob('sistr_cmd/*.py')],
    install_requires=[
        'attrs>=18.1.0',
        'numpy>=1.15.0',
        'scipy>=1.1.0',
        'pandas>=0.23.4',
        'tables>=3.4.4',
    ],
    tests_require=['pytest', ],
    cmdclass={'test': PyTest},  # run tests with `python setup.py test`
    entry_points={
        'console_scripts': [
            'sistr=sistr_cmd.main:main',
        ],
    },
)
