import os
from distutils.core import setup
from setuptools import find_packages

classifiers = """
Development Status :: 3 - Alpha
Environment :: Console
License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)
Intended Audience :: Science/Research
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bio-Informatics
Programming Language :: Python
Programming Language :: Python :: 2.7
Programming Language :: Python :: 3.4
Programming Language :: Python :: 3.5
Programming Language :: Python :: 3.6
Programming Language :: Python :: Implementation :: CPython
Operating System :: POSIX :: Linux
""".strip().split('\n')

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='sistr_cmd',
    version='0.1.0',
    packages=find_packages(exclude=['tests']),
    url='https://github.com/peterk87/sistr_cmd',
    license='GPLv3',
    author='Peter Kruczkiewicz',
    author_email='peter.kruczkiewicz@gmail.com',
    description=('Serovar predictions from Salmonella whole-genome sequence assemblies by determination of antigen gene'
                 'and cgMLST gene alleles using BLAST. Mash MinHash can also be used for serovar prediction.'),
    keywords='Salmonella serotyping genotyping cgMLST BLAST Mash MinHash',
    classifiers=classifiers,
    package_dir={'sistr':'sistr'},
    package_data={'sistr': ['data/*.msh',
                            'data/*.csv',
                            'data/antigens/*.fasta',
                            'data/cgmlst/*.fasta',
                            'data/cgmlst/*.txt',
                            'data/cgmlst/*.csv',]},
    install_requires=['numpy', 'pandas'],
    extras_require={
        'test': ['pytest',],
    },
    entry_points={
        'console_scripts': [
            'sistr=sistr.sistr_cmd:main',
        ],
    },
)
