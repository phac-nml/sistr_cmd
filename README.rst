*************************************************************
Salmonella In Silico Typing Resource (SISTR) commandline tool
*************************************************************


Serovar predictions from whole-genome sequence assemblies by determination of antigen gene and cgMLST gene alleles using BLAST.
`Mash MinHash <https://mash.readthedocs.io/en/latest/>`_ can also be used for serovar prediction.

*Don't want to use a command-line app?* Try the `SISTR web app <https://lfz.corefacility.ca/sistr-app/>`_

Citation
========

If you find this tool useful, please cite as:

.. epigraph::

   The *Salmonella In Silico* Typing Resource (SISTR): an open web-accessible tool for rapidly typing and subtyping draft *Salmonella* genome assemblies. Catherine Yoshida, Peter Kruczkiewicz, Chad R. Laing, Erika J. Lingohr, Victor P.J. Gannon, John H.E. Nash, Eduardo N. Taboada. PLoS ONE 11(1): e0147101. doi: 10.1371/journal.pone.0147101

   -- Paper Link: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0147101

BibTeX
------

.. code-block:: none

    @article{Yoshida2016,
      doi = {10.1371/journal.pone.0147101},
      url = {http://dx.doi.org/10.1371/journal.pone.0147101},
      year  = {2016},
      month = {jan},
      publisher = {Public Library of Science ({PLoS})},
      volume = {11},
      number = {1},
      pages = {e0147101},
      author = {Catherine E. Yoshida and Peter Kruczkiewicz and Chad R. Laing and Erika J. Lingohr and Victor P. J. Gannon and John H. E. Nash and Eduardo N. Taboada},
      editor = {Michael Hensel},
      title = {The Salmonella In Silico Typing Resource ({SISTR}): An Open Web-Accessible Tool for Rapidly Typing and Subtyping Draft Salmonella Genome Assemblies},
      journal = {{PLOS} {ONE}}
    }


Installation
============

You can install ``sistr_cmd`` using ``pip``:

.. code-block:: bash

    pip install sistr_cmd


``sistr_cmd`` is available from PYPI at https://pypi.python.org/pypi/sistr-cmd

Dependencies
============

- BLAST+ (>= v2.2.30)
- Python (>= v2.7 OR >= v3.4)
- `Mash v1.0+ <https://github.com/marbl/Mash/releases>`_ [optional]

Python Dependencies
-------------------

``sistr_cmd`` requires the following Python libraries:

- numpy (>=1.10.4)
- pandas (>=0.17.1)


You can run the following commands to get up-to-date versions of ``numpy`` and ``pandas``

.. code-block:: bash

    pip install --upgrade pip
    pip install wheel
    pip install numpy pandas

Usage
=====


.. code-block:: none

    usage: predict_serovar [-h] [-f OUTPUT_FORMAT] [-o OUTPUT_DEST] [-T TMP_DIR]
                           [-K] [--no-cgmlst] [-m] [-t THREADS] [-v]
                           F [F ...]

    SISTR (Salmonella In Silico Typing Resource) Command-line Tool
    ==============================================================
    Serovar predictions from whole-genome sequence assemblies by determination of antigen gene and cgMLST gene alleles using BLAST.

    If you find this program useful in your research, please cite as:

    The Salmonella In Silico Typing Resource (SISTR): an open web-accessible tool for rapidly typing and subtyping draft Salmonella genome assemblies.
    Catherine Yoshida, Peter Kruczkiewicz, Chad R. Laing, Erika J. Lingohr, Victor P.J. Gannon, John H.E. Nash, Eduardo N. Taboada.
    PLoS ONE 11(1): e0147101. doi: 10.1371/journal.pone.0147101

    positional arguments:
      F                     Input genome FASTA file

    optional arguments:
      -h, --help            show this help message and exit
      -f OUTPUT_FORMAT, --output-format OUTPUT_FORMAT
                            Output format (json, csv, pickle)
      -o OUTPUT_DEST, --output-dest OUTPUT_DEST
                            Output
      -T TMP_DIR, --tmp-dir TMP_DIR
                            Base temporary working directory for intermediate
                            analysis files.
      -K, --keep-tmp        Keep temporary analysis files.
      --no-cgmlst           Do not run cgMLST serovar prediction
      -m, --run-mash        Determine Mash MinHash genomic distances to Salmonella
                            genomes with trusted serovar designations. Mash binary
                            must be in accessible via $PATH (e.g. /usr/bin).
      -t THREADS, --threads THREADS
                            Number of parallel threads to run sistr_cmd analysis.
      -v, --verbose         Logging verbosity level (-v == show warnings; -vvv ==
                            show debug info)

Example output
==============

By running the following command on a FASTA file of *Salmonella enterica* strain LT2 (https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP014051.1):

.. code-block:: bash

    sistr -f csv -o out.csv -m LT2.fasta


You should see some log messages like so:

.. code-block:: none

    <TIME> INFO: Initializing temporary analysis directory and preparing for BLAST searching. [in sistr_cmd.py:152]
    <TIME> INFO: Temporary FASTA file copied to tmp/LT2.fasta [in sistr_cmd.py:154]
    <TIME> INFO: Running BLAST on serovar predictive cgMLST330 alleles [in sistr_cmd.py:64]
    <TIME> INFO: Reading BLAST output file "tmp/cgmlst330.fasta-LT2.fasta.blast" [in sistr_cmd.py:66]
    <TIME> INFO: Found 39203 cgMLST330 allele BLAST results [in sistr_cmd.py:71]
    <TIME> INFO: Found 330 perfect matches to cgMLST330 alleles [in sistr_cmd.py:76]
    <TIME> INFO: Calculating number of matching alleles to serovar predictive cgMLST330 profiles [in sistr_cmd.py:78]
    <TIME> INFO: Top serovar by cgMLST profile matching: "Typhimurium" with 330 matching alleles, distance=0.0% [in sistr_cmd.py:96]
    <TIME> INFO: Top serovar by Mash: "Typhimurium" with dist=0.0, # matching sketches=1000, matching genome=LT2 [in sistr_cmd.py:120]
    <TIME> INFO: Antigen gene BLAST serovar prediction: "Typhimurium" serogroup=B:H1=i:H2=1,2 [in sistr_cmd.py:169]
    <TIME> INFO: Overall serovar prediction: Typhimurium [in sistr_cmd.py:170]
    <TIME> INFO: Writing output "csv" file to "out.csv" [in src/writers.py:33]
    <TIME> INFO: Deleting temporary working directory at tmp [in sistr_cmd.py:178]


CSV Output
----------

.. code-block:: csv

    cgmlst_distance,cgmlst_genome_match,cgmlst_matching_alleles,genome,h1,h2,mash_distance,mash_genome,mash_match,mash_serovar,serogroup,serovar,serovar_antigen,serovar_cgmlst
    0.0,LT2,330,LT2.fasta,i,"1,2",0.0,LT2,1000,Typhimurium,B,Typhimurium,Typhimurium,Typhimurium

.. csv-table:: 

   cgmlst_distance,cgmlst_genome_match,cgmlst_matching_alleles,genome,h1,h2,mash_distance,mash_genome,mash_match,mash_serovar,serogroup,serovar,serovar_antigen,serovar_cgmlst
    0.0,LT2,330,LT2.fasta,i,"1,2",0.0,LT2,1000,Typhimurium,B,Typhimurium,Typhimurium,Typhimurium



JSON output
-----------

.. code-block:: json

    [{
      "mash_serovar": "Typhimurium",
      "serovar_cgmlst": "Typhimurium",
      "cgmlst_matching_alleles": 330,
      "mash_match": 1000,
      "serovar_antigen": "Typhimurium",
      "cgmlst_distance": 0.0,
      "mash_distance": 0.0,
      "h2": "1,2",
      "cgmlst_genome_match": "LT2",
      "h1": "i",
      "mash_genome": "LT2",
      "serovar": "Typhimurium",
      "genome": "LT2.fasta",
      "serogroup": "B"
    }]


TODO
====

- add more genomes to improve cgMLST and Mash serovar calling (7511 -> infinity)
- output of cgMLST allele results
    + output novel alleles?


Issues
======

If you encounter any problems or have any questions feel free to create an issue anonymously or not to let us know so we can address it!

License
=======

Copyright 2016 Public Health Agency of Canada

Distributed under the GNU Public License version 3.0
