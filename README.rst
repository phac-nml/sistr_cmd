****************************************************************************
``sistr_cmd``: Salmonella In Silico Typing Resource (SISTR) commandline tool
****************************************************************************


Serovar predictions from whole-genome sequence assemblies by determination of antigen gene and cgMLST gene alleles using BLAST.
`Mash MinHash <https://mash.readthedocs.io/en/latest/>`_ can also be used for serovar prediction.

**Latest stable version: `v0.3.4 <https://github.com/peterk87/sistr_cmd/releases/tag/v0.3.4>`_ **


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

Using Conda [Recommended]
-----------

You can install ``sistr_cmd`` using `Conda <https://conda.io/miniconda.html>`_ from the `BioConda channel <https://bioconda.github.io/>`_:

.. code-block:: bash

	# Install conda. Miniconda is recommended if you don't have Conda installed already
	# see https://conda.io/miniconda.html
	# Add Bioconda channel and other channels https://bioconda.github.io/#set-up-channels
	conda config --add channels conda-forge
	conda config --add channels defaults
	conda config --add channels r
	conda config --add channels bioconda
	# Install sistr_cmd and its dependencies
	conda install sistr_cmd
	# sistr_cmd should be installed in your $PATH
	sistr --help

Installing ``sistr_cmd`` is recommended for the least amount of headache since Conda will ensure that all necessary external dependencies are installed along with ``sistr_cmd`` (i.e. ``blast+``, ``mafft``, ``mash``). This will also help you get ``sistr_cmd`` running on older systems (e.g. CentOS 5) or where you may not have many user privileges. 


Using ``pip``
-------------

You can install ``sistr_cmd`` using ``pip``:

.. code-block:: bash

	pip install sistr_cmd


``sistr_cmd`` is available from PYPI at https://pypi.python.org/pypi/sistr-cmd

**NOTE:** You will need to ensure that external dependencies are installed (i.e. ``blast+``, ``mafft``, ``mash`` [optionally])


Dependencies
============

These are the external dependencies required for ``sistr_cmd``:

- Python (>= v2.7 OR >= v3.4)
- BLAST+ (>= v2.2.30)
- MAFFT (>=v7.271 (2016/1/6))
- `Mash v1.0+ <https://github.com/marbl/Mash/releases>`_ [optional]

Python Dependencies
-------------------

``sistr_cmd`` requires the following Python libraries:

- numpy (>=1.11.1)
- pandas (>=0.18.1)


You can run the following commands to get up-to-date versions of ``numpy`` and ``pandas``

.. code-block:: bash

	pip install --upgrade pip
	pip install wheel
	pip install numpy pandas

Usage
=====

If you run ``sistr -h``, you should see the following usage info:

.. code-block:: none

	usage: sistr_cmd [-h] [-i fasta_path genome_name] [-f OUTPUT_FORMAT]
	                 [-o OUTPUT_PREDICTION] [-p CGMLST_PROFILES]
	                 [-n NOVEL_ALLELES] [-a ALLELES_OUTPUT] [-T TMP_DIR] [-K]
	                 [--use-full-cgmlst-db] [--no-cgmlst] [-m] [--qc] [-t THREADS]
	                 [-v] [-V]
	                 [F [F ...]]

	SISTR (Salmonella In Silico Typing Resource) Command-line Tool
	==============================================================
	Serovar predictions from whole-genome sequence assemblies by determination of antigen gene and cgMLST gene alleles using BLAST.

	Note about using the "--use-full-cgmlst-db" flag:
	    The "centroid" allele database is ~10% the size of the full set so analysis is much quicker with the "centroid" vs "full" set of alleles. Results between 2 cgMLST allele sets should not differ.

	If you find this program useful in your research, please cite as:

	The Salmonella In Silico Typing Resource (SISTR): an open web-accessible tool for rapidly typing and subtyping draft Salmonella genome assemblies.
	Catherine Yoshida, Peter Kruczkiewicz, Chad R. Laing, Erika J. Lingohr, Victor P.J. Gannon, John H.E. Nash, Eduardo N. Taboada.
	PLoS ONE 11(1): e0147101. doi: 10.1371/journal.pone.0147101

	positional arguments:
	  F                     Input genome FASTA file

	optional arguments:
	  -h, --help            show this help message and exit
	  -i fasta_path genome_name, --input-fasta-genome-name fasta_path genome_name
	                        fasta file path to genome name pair
	  -f OUTPUT_FORMAT, --output-format OUTPUT_FORMAT
	                        Output format (json, csv, pickle)
	  -o OUTPUT_PREDICTION, --output-prediction OUTPUT_PREDICTION
	                        SISTR serovar prediction output path
	  -p CGMLST_PROFILES, --cgmlst-profiles CGMLST_PROFILES
	                        Output CSV file destination for cgMLST allelic
	                        profiles
	  -n NOVEL_ALLELES, --novel-alleles NOVEL_ALLELES
	                        Output FASTA file destination of novel cgMLST alleles
	                        from input genomes
	  -a ALLELES_OUTPUT, --alleles-output ALLELES_OUTPUT
	                        Output path of allele sequences and info to JSON
	  -T TMP_DIR, --tmp-dir TMP_DIR
	                        Base temporary working directory for intermediate
	                        analysis files.
	  -K, --keep-tmp        Keep temporary analysis files.
	  --use-full-cgmlst-db  Use the full set of cgMLST alleles which can include
	                        highly similar alleles. By default the smaller
	                        "centroid" alleles or representative alleles are used
	                        for each marker.
	  --no-cgmlst           Do not run cgMLST serovar prediction
	  -m, --run-mash        Determine Mash MinHash genomic distances to Salmonella
	                        genomes with trusted serovar designations. Mash binary
	                        must be in accessible via $PATH (e.g. /usr/bin).
	  --qc                  Perform basic QC to provide level of confidence in
	                        serovar prediction results.
	  -t THREADS, --threads THREADS
	                        Number of parallel threads to run sistr_cmd analysis.
	  -v, --verbose         Logging verbosity level (-v == show warnings; -vvv ==
	                        show debug info)
	  -V, --version         show program's version number and exit



Example Usage
-------------

By running the following command on a FASTA file of *Salmonella enterica* strain LT2 (https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP014051.1):

.. code-block:: bash

	sistr --qc -vv --alleles-output allele-results.json --novel-alleles novel-alleles.fasta --cgmlst-profiles cgmlst-profiles.csv -f tab -o sistr-output.tab LT2.fasta


You should see some log messages like so:

.. code-block:: none

	<time> INFO: Running sistr_cmd 0.3.4 [in /usr/lib/python2.7/site-packages/sistr/sistr_cmd.py:290]
	<time> INFO: Serial single threaded run mode on 1 genomes [in /usr/lib/python2.7/site-packages/sistr/sistr_cmd.py:319]
	<time> INFO: Initializing temporary analysis directory "/tmp/20170309104912-SISTR-LT2" and preparing for BLAST searching. [in /usr/lib/python2.7/site-packages/sistr/sistr_cmd.py:175]
	<time> INFO: Temporary FASTA file copied to /tmp/20170309104912-SISTR-LT2/LT2_fasta [in /usr/lib/python2.7/site-packages/sistr/sistr_cmd.py:177]
	<time> INFO: Running BLAST on serovar predictive cgMLST330 alleles [in /usr/lib/python2.7/site-packages/sistr/src/cgmlst/__init__.py:319]
	<time> INFO: Reading BLAST output file "/tmp/20170309104912-SISTR-LT2/cgmlst-centroid.fasta-LT2_fasta-2017Mar09_10_49_13.blast" [in /usr/lib/python2.7/site-packages/sistr/src/cgmlst/__init__.py:322]
	<time> INFO: Found 6525 cgMLST330 allele BLAST results [in /usr/lib/python2.7/site-packages/sistr/src/cgmlst/__init__.py:333]
	<time> INFO: Marker NZ_AOXE01000081.1_201 | Recovered novel allele with gaps (n=0) of length 477 vs length 477 for ref allele NZ_AOXE01000081.1_201|2823059714. Novel allele name=3250876267 [in /usr/lib/python2.7/site-packages/sistr/src/cgmlst/__init__.py:181]
	<time> INFO: Type retrieved_marker_alleles <type 'dict'> [in /usr/lib/python2.7/site-packages/sistr/src/cgmlst/__init__.py:343]
	<time> INFO: Calculating number of matching alleles to serovar predictive cgMLST330 profiles [in /usr/lib/python2.7/site-packages/sistr/src/cgmlst/__init__.py:360]
	<time> INFO: Top subspecies by cgMLST is "enterica" (min dist=0.00909090909091, Counter={'enterica': 11532}) [in /usr/lib/python2.7/site-packages/sistr/src/cgmlst/__init__.py:369]
	<time> INFO: Top serovar by cgMLST profile matching: "Typhimurium" with 327 matching alleles, distance=0.9% [in /usr/lib/python2.7/site-packages/sistr/src/cgmlst/__init__.py:385]
	<time> INFO: cgMLST330 Sequence Type=660408169 [in /usr/lib/python2.7/site-packages/sistr/src/cgmlst/__init__.py:404]
	<time> INFO: LT2 | Antigen gene BLAST serovar prediction: "Typhimurium" serogroup=B 1,4,[5],12:i:1,2 [in /usr/lib/python2.7/site-packages/sistr/sistr_cmd.py:207]
	<time> INFO: LT2 | Subspecies prediction: enterica [in /usr/lib/python2.7/site-packages/sistr/sistr_cmd.py:210]
	<time> INFO: LT2 | Overall serovar prediction: Typhimurium [in /usr/lib/python2.7/site-packages/sistr/sistr_cmd.py:213]
	<time> INFO: Genome size=4857473 (within gsize thresholds? True) [in /usr/lib/python2.7/site-packages/sistr/src/qc/__init__.py:13]
	<time> INFO: Deleting temporary working directory at /tmp/20170309104912-SISTR-LT2 [in /usr/lib/python2.7/site-packages/sistr/sistr_cmd.py:220]
	<time> INFO: Writing output "tab" file to "sistr-output.tab" [in /usr/lib/python2.7/site-packages/sistr/src/writers.py:38]
	<time> INFO: cgMLST allelic profiles written to cgmlst-profiles.csv [in /usr/lib/python2.7/site-packages/sistr/sistr_cmd.py:340]
	<time> INFO: JSON of allele data written to allele-results.json for 1 cgMLST allele results [in /usr/lib/python2.7/site-packages/sistr/sistr_cmd.py:343]
	<time> INFO: Wrote 330 alleles to novel-alleles.fasta [in /usr/lib/python2.7/site-packages/sistr/sistr_cmd.py:346]


``sistr_cmd`` Output
====================

``sistr_cmd`` has several output options. The primary output is the serovar prediction and in silico typing results output (e.g. ``-o sistr-results.tab``).

Summary of output options:

- primary results output 
	+ serovar prediction, cgMLST results, Mash results
	+ format (``-f <format>``): ``tab``, ``csv``, ``json``, ``pickle``
	+ ``-o sistr-results``
- cgMLST allele results
	+ in-depth allele search results for each input genome for each cgMLST locus (330 loci in total)
	+ includes extracted allele sequences, top ``blastn`` results and summarized ``mafft`` results
	+ format: JSON
	+ ``-a allele-results.json``
- cgMLST allelic profiles
	+ table of allele designations for each genome for each cgMLST locus
	+ row names: genome names
	+ column names: cgMLST marker names
	+ format: CSV
	+ ``--cgmlst-profiles cgmlst-profiles.csv``


Primary results output (``-o sistr-results``)
------------------------------------------

Tab-delimited results output (``-f tab``):

.. code-block:: tab
	
	cgmlst_ST	cgmlst_distance	cgmlst_genome_match	cgmlst_matching_alleles	cgmlst_subspecies	fasta_filepath	genome	h1	h2	o_antigen	qc_messages	qc_status	serogroup	serovar	serovar_antigen	serovar_cgmlst
	660408169	0.00909090909091	LT2	327	enterica	/home/peter/Downloads/sistr-LT2-example/LT2.fasta	LT2	i	1,2	1,4,[5],12		PASS	B	Typhimurium	Typhimurium	Typhimurium

CSV results output (``-f csv``):

.. code-block:: csv

	cgmlst_ST,cgmlst_distance,cgmlst_genome_match,cgmlst_matching_alleles,cgmlst_subspecies,fasta_filepath,genome,h1,h2,o_antigen,qc_messages,qc_status,serogroup,serovar,serovar_antigen,serovar_cgmlst
	660408169,0.00909090909091,LT2,327,enterica,/home/peter/Downloads/sistr-LT2-example/LT2.fasta,LT2,i,"1,2","1,4,[5],12",,PASS,B,Typhimurium,Typhimurium,Typhimurium

How the results should look in a table:

.. csv-table:: 

	cgmlst_ST,cgmlst_distance,cgmlst_genome_match,cgmlst_matching_alleles,cgmlst_subspecies,fasta_filepath,genome,h1,h2,o_antigen,qc_messages,qc_status,serogroup,serovar,serovar_antigen,serovar_cgmlst
	660408169,0.00909090909091,LT2,327,enterica,/home/peter/Downloads/sistr-LT2-example/LT2.fasta,LT2,i,"1,2","1,4,[5],12",,PASS,B,Typhimurium,Typhimurium,Typhimurium


JSON results output:

.. code-block:: json

	[
	  {
	    "serovar_cgmlst": "Typhimurium",
	    "cgmlst_matching_alleles": 327,
	    "h1": "i",
	    "serovar_antigen": "Typhimurium",
	    "cgmlst_distance": 0.009090909090909038,
	    "h2": "1,2",
	    "cgmlst_genome_match": "LT2",
	    "cgmlst_ST": 660408169,
	    "serovar": "Typhimurium",
	    "fasta_filepath": "/full/path/to/LT2.fasta",
	    "genome": "LT2",
	    "serogroup": "B",
	    "qc_messages": "",
	    "qc_status": "PASS",
	    "o_antigen": "1,4,[5],12",
	    "cgmlst_subspecies": "enterica"
	  }
	]

cgMLST allele search results
-------------------------------------

You can produce in-depth allele search results with the ``-a``/``--alleles-output`` commandline argument.
These results may be useful for understanding unexpected or low confidence serovar predictions.

Schema:
~~~~~~~

.. code-block:: json
	
	{
		<genome name>: {
			// for each 
			<cgMLST marker id>: {
				// top blast result on largest contig
				blast_result: {
					// perfect match to a previously identified allele?
					"is_perfect": boolean,
					// blastn subject sequence length
					"slen": integer,
					// blastn percent identity
					"pident": numeric,
					// cgMLST marker name
					"marker": string,
					// blastn query sequence id
					"qseqid": string,
					// blastn query sequence start index
					"qstart": integer,
					// is match truncated by end of sequence? 
					"is_trunc": boolean,
					// number of MSA gaps in subject sequence
					"sseq_msa_gaps": integer,
					// blastn subject sequence
					"sseq": string,
					// blastn bitscore
					"bitscore": numeric,
					// proportion of subject sequence MSA with gaps
					"sseq_msa_p_gaps": numeric,
					// blastn E-value
					"evalue": numeric,
					// blastn gap open
					"gapopen": integer,
					// blastn subject sequence end index
					"send": integer,
					// does this allele have a perfect match?
					"has_perfect_match": boolean,
					// matching allele name
					"allele": integer,
					// subject sequence start index
					"sstart": integer,
					// extracted allele name (CRC32 of subject nucleotide sequence)
					"allele_name": integer,
					// adjusted subject sequence start index
					"start_idx": numeric,
					// blastn query end index
					"qend": integer,
					// did the extracted allele sequence need to be reverse complemented?
					"needs_revcomp": boolean,
					// did the extracted allele sequence need to be extended to match the length of the query sequence?
					"is_extended": boolean,
					// blastn number of mismatches
					"mismatch": integer,
					// extracted allele coverage i.e. (length of extracted allele) / (length of closest matching allele)
					"coverage": numeric,
					// too many gaps within the MSA of extracted allele sequence and closest matching allele?
					"too_many_gaps": boolean,
					// adjusted subject end index
					"end_idx": numeric,
					// is extracted allele truncated by end of sequence? 
					"trunc": boolean,
					// blastn subject sequence title
					"stitle": string,
					// blastn query sequence length
					"qlen": integer,
					// valid allele match found?
					"is_match": true,
					// blastn alignment length
					"length": integer
				},
				// CRC32 unsigned 32-bit integer allele name from allele sequence
				"name": integer,
				// extracted allele sequence
				"seq": string
			}
			
		}}

Example:
~~~~~~~~

Here's some truncated example allele search results output:

.. code-block:: json

	{
	  "LT2": {
	    "NZ_AOXE01000034.1_82": {
	      "blast_result": {
	        "is_perfect": false,
	        "slen": 4857473,
	        "pident": 99.479,
	        "marker": "NZ_AOXE01000034.1_82",
	        "qseqid": "NZ_AOXE01000034.1_82|340989631",
	        "qstart": 1,
	        "is_trunc": false,
	        "sseq_msa_gaps": 0,
	        "sseq": "ATGCCAACCAGACCACCTTATCCGCGGGAAGCTTATATCGTCACCATTGAAAAAGGCACGCCGGGCCAGACGGTGACGTGGTATCAGCTACGGGCTGACCATCCGAAACCTGATTCGCTCATCAGCGAGCATCCGACCGCAGAAGAAGCGATGGATGCGAAAAATCGTTACGAAGATCCGGATAAATCATAG",
	        "bitscore": 350.0,
	        "sseq_msa_p_gaps": 0.0,
	        "evalue": 3.289999999999999E-97,
	        "gapopen": 0,
	        "send": 358277,
	        "has_perfect_match": false,
	        "allele": 340989631,
	        "sstart": 358468,
	        "allele_name": 1204520418,
	        "start_idx": 358276.0,
	        "qend": 192,
	        "needs_revcomp": true,
	        "is_extended": false,
	        "mismatch": 1,
	        "coverage": 1.0,
	        "too_many_gaps": false,
	        "end_idx": 358467.0,
	        "trunc": false,
	        "stitle": "NZ_CP014051.1 Salmonella enterica strain LT2, complete genome",
	        "qlen": 192,
	        "is_match": true,
	        "length": 192
	      },
	      "name": 1204520418,
	      "seq": "ATGCCAACCAGACCACCTTATCCGCGGGAAGCTTATATCGTCACCATTGAAAAAGGCACGCCGGGCCAGACGGTGACGTGGTATCAGCTACGGGCTGACCATCCGAAACCTGATTCGCTCATCAGCGAGCATCCGACCGCAGAAGAAGCGATGGATGCGAAAAATCGTTACGAAGATCCGGATAAATCATAG"
	    },
	    // 329 other cgMLST allele results
	  },
	  "another-genome": { /* allele results */}
	}


cgMLST allelic profiles output (``--cgmlst-profiles cgmlst-profiles.csv``)
--------------------------------------------------------------------------

With the ``-p``/``--cgmlst-profiles`` commandline argument, you can output the 330 loci cgMLST allelic profiles for your input genomes (i.e. the allele designation for each cgMLST locus for each input genome). 
You can use this information to construct phylogenetic trees from this data using a tool such as `Phyloviz Online <https://online.phyloviz.net/index>`_. 
This type of analysis may be useful to explore why unexpected serovar prediction results were generated (e.g. your genomes are genetically very different from each other). 

Example cgMLST profiles output:

.. csv-table::

	,NC_003198.1_3005,NC_006905.1_2841,NC_011149.1_467,NC_017623.1_3300,NZ_ABFH02000002.1_1303,NZ_AOXE01000003.1_37,NZ_AOXE01000003.1_39,NZ_AOXE01000003.1_57,NZ_AOXE01000003.1_7,NZ_AOXE01000003.1_70,NZ_AOXE01000004.1_10,NZ_AOXE01000004.1_101,NZ_AOXE01000004.1_12,NZ_AOXE01000004.1_134,NZ_AOXE01000004.1_135,NZ_AOXE01000004.1_14,NZ_AOXE01000004.1_140,NZ_AOXE01000004.1_154,NZ_AOXE01000004.1_35,NZ_AOXE01000004.1_36,NZ_AOXE01000004.1_39,NZ_AOXE01000004.1_59,NZ_AOXE01000004.1_68,NZ_AOXE01000004.1_74,NZ_AOXE01000004.1_87,NZ_AOXE01000007.1_13,NZ_AOXE01000007.1_18,NZ_AOXE01000007.1_20,NZ_AOXE01000007.1_48,NZ_AOXE01000008.1_59,NZ_AOXE01000008.1_63,NZ_AOXE01000009.1_17,NZ_AOXE01000011.1_101,NZ_AOXE01000011.1_77,NZ_AOXE01000011.1_82,NZ_AOXE01000011.1_83,NZ_AOXE01000011.1_85,NZ_AOXE01000016.1_13,NZ_AOXE01000016.1_8,NZ_AOXE01000017.1_117,NZ_AOXE01000017.1_118,NZ_AOXE01000017.1_130,NZ_AOXE01000017.1_4,NZ_AOXE01000017.1_40,NZ_AOXE01000017.1_43,NZ_AOXE01000017.1_54,NZ_AOXE01000017.1_59,NZ_AOXE01000017.1_80,NZ_AOXE01000017.1_82,NZ_AOXE01000017.1_96,NZ_AOXE01000019.1_13,NZ_AOXE01000019.1_14,NZ_AOXE01000019.1_24,NZ_AOXE01000021.1_10,NZ_AOXE01000021.1_11,NZ_AOXE01000021.1_165,NZ_AOXE01000021.1_29,NZ_AOXE01000021.1_38,NZ_AOXE01000021.1_49,NZ_AOXE01000021.1_6,NZ_AOXE01000021.1_61,NZ_AOXE01000021.1_79,NZ_AOXE01000023.1_11,NZ_AOXE01000023.1_25,NZ_AOXE01000023.1_30,NZ_AOXE01000024.1_3,NZ_AOXE01000024.1_35,NZ_AOXE01000024.1_38,NZ_AOXE01000025.1_13,NZ_AOXE01000025.1_14,NZ_AOXE01000025.1_20,NZ_AOXE01000031.1_102,NZ_AOXE01000031.1_106,NZ_AOXE01000031.1_70,NZ_AOXE01000031.1_80,NZ_AOXE01000033.1_11,NZ_AOXE01000033.1_12,NZ_AOXE01000033.1_14,NZ_AOXE01000033.1_17,NZ_AOXE01000033.1_19,NZ_AOXE01000033.1_2,NZ_AOXE01000033.1_21,NZ_AOXE01000033.1_26,NZ_AOXE01000033.1_3,NZ_AOXE01000033.1_30,NZ_AOXE01000033.1_34,NZ_AOXE01000033.1_38,NZ_AOXE01000033.1_43,NZ_AOXE01000033.1_51,NZ_AOXE01000034.1_103,NZ_AOXE01000034.1_106,NZ_AOXE01000034.1_111,NZ_AOXE01000034.1_112,NZ_AOXE01000034.1_113,NZ_AOXE01000034.1_119,NZ_AOXE01000034.1_126,NZ_AOXE01000034.1_127,NZ_AOXE01000034.1_133,NZ_AOXE01000034.1_134,NZ_AOXE01000034.1_164,NZ_AOXE01000034.1_173,NZ_AOXE01000034.1_53,NZ_AOXE01000034.1_82,NZ_AOXE01000035.1_13,NZ_AOXE01000035.1_21,NZ_AOXE01000036.1_108,NZ_AOXE01000036.1_116,NZ_AOXE01000036.1_15,NZ_AOXE01000036.1_157,NZ_AOXE01000036.1_16,NZ_AOXE01000036.1_2,NZ_AOXE01000036.1_3,NZ_AOXE01000036.1_31,NZ_AOXE01000036.1_39,NZ_AOXE01000036.1_43,NZ_AOXE01000036.1_58,NZ_AOXE01000036.1_66,NZ_AOXE01000036.1_98,NZ_AOXE01000040.1_19,NZ_AOXE01000040.1_28,NZ_AOXE01000040.1_31,NZ_AOXE01000041.1_33,NZ_AOXE01000041.1_73,NZ_AOXE01000041.1_75,NZ_AOXE01000041.1_76,NZ_AOXE01000041.1_84,NZ_AOXE01000041.1_85,NZ_AOXE01000041.1_87,NZ_AOXE01000043.1_4,NZ_AOXE01000047.1_56,NZ_AOXE01000047.1_57,NZ_AOXE01000050.1_18,NZ_AOXE01000050.1_44,NZ_AOXE01000052.1_115,NZ_AOXE01000052.1_128,NZ_AOXE01000052.1_131,NZ_AOXE01000052.1_137,NZ_AOXE01000052.1_141,NZ_AOXE01000052.1_23,NZ_AOXE01000052.1_36,NZ_AOXE01000052.1_38,NZ_AOXE01000052.1_41,NZ_AOXE01000052.1_43,NZ_AOXE01000052.1_78,NZ_AOXE01000052.1_92,NZ_AOXE01000053.1_113,NZ_AOXE01000053.1_128,NZ_AOXE01000053.1_130,NZ_AOXE01000053.1_166,NZ_AOXE01000053.1_173,NZ_AOXE01000053.1_180,NZ_AOXE01000053.1_190,NZ_AOXE01000053.1_217,NZ_AOXE01000053.1_86,NZ_AOXE01000059.1_11,NZ_AOXE01000059.1_129,NZ_AOXE01000059.1_133,NZ_AOXE01000059.1_15,NZ_AOXE01000059.1_174,NZ_AOXE01000059.1_182,NZ_AOXE01000059.1_184,NZ_AOXE01000059.1_189,NZ_AOXE01000059.1_229,NZ_AOXE01000059.1_31,NZ_AOXE01000059.1_32,NZ_AOXE01000059.1_325,NZ_AOXE01000059.1_328,NZ_AOXE01000059.1_333,NZ_AOXE01000059.1_335,NZ_AOXE01000059.1_336,NZ_AOXE01000059.1_338,NZ_AOXE01000059.1_35,NZ_AOXE01000059.1_353,NZ_AOXE01000059.1_363,NZ_AOXE01000059.1_37,NZ_AOXE01000059.1_370,NZ_AOXE01000059.1_372,NZ_AOXE01000059.1_38,NZ_AOXE01000059.1_395,NZ_AOXE01000059.1_396,NZ_AOXE01000059.1_408,NZ_AOXE01000059.1_411,NZ_AOXE01000059.1_418,NZ_AOXE01000059.1_42,NZ_AOXE01000059.1_427,NZ_AOXE01000059.1_430,NZ_AOXE01000059.1_433,NZ_AOXE01000059.1_435,NZ_AOXE01000059.1_437,NZ_AOXE01000059.1_440,NZ_AOXE01000059.1_442,NZ_AOXE01000059.1_49,NZ_AOXE01000059.1_60,NZ_AOXE01000059.1_66,NZ_AOXE01000059.1_67,NZ_AOXE01000059.1_68,NZ_AOXE01000059.1_69,NZ_AOXE01000059.1_72,NZ_AOXE01000059.1_79,NZ_AOXE01000059.1_9,NZ_AOXE01000059.1_94,NZ_AOXE01000061.1_12,NZ_AOXE01000061.1_20,NZ_AOXE01000061.1_22,NZ_AOXE01000061.1_3,NZ_AOXE01000064.1_26,NZ_AOXE01000064.1_27,NZ_AOXE01000064.1_36,NZ_AOXE01000068.1_19,NZ_AOXE01000068.1_20,NZ_AOXE01000068.1_27,NZ_AOXE01000068.1_29,NZ_AOXE01000068.1_37,NZ_AOXE01000068.1_38,NZ_AOXE01000068.1_45,NZ_AOXE01000068.1_46,NZ_AOXE01000068.1_5,NZ_AOXE01000068.1_52,NZ_AOXE01000068.1_58,NZ_AOXE01000068.1_65,NZ_AOXE01000068.1_67,NZ_AOXE01000068.1_70,NZ_AOXE01000068.1_72,NZ_AOXE01000068.1_76,NZ_AOXE01000072.1_100,NZ_AOXE01000072.1_104,NZ_AOXE01000072.1_12,NZ_AOXE01000072.1_13,NZ_AOXE01000072.1_3,NZ_AOXE01000072.1_41,NZ_AOXE01000072.1_42,NZ_AOXE01000072.1_60,NZ_AOXE01000072.1_65,NZ_AOXE01000072.1_73,NZ_AOXE01000072.1_8,NZ_AOXE01000072.1_82,NZ_AOXE01000072.1_83,NZ_AOXE01000072.1_86,NZ_AOXE01000072.1_93,NZ_AOXE01000073.1_11,NZ_AOXE01000073.1_130,NZ_AOXE01000073.1_144,NZ_AOXE01000073.1_15,NZ_AOXE01000073.1_19,NZ_AOXE01000073.1_48,NZ_AOXE01000073.1_79,NZ_AOXE01000073.1_85,NZ_AOXE01000073.1_98,NZ_AOXE01000077.1_25,NZ_AOXE01000077.1_28,NZ_AOXE01000077.1_29,NZ_AOXE01000077.1_33,NZ_AOXE01000077.1_35,NZ_AOXE01000079.1_15,NZ_AOXE01000079.1_4,NZ_AOXE01000080.1_12,NZ_AOXE01000080.1_13,NZ_AOXE01000080.1_20,NZ_AOXE01000081.1_103,NZ_AOXE01000081.1_105,NZ_AOXE01000081.1_124,NZ_AOXE01000081.1_136,NZ_AOXE01000081.1_179,NZ_AOXE01000081.1_186,NZ_AOXE01000081.1_190,NZ_AOXE01000081.1_193,NZ_AOXE01000081.1_195,NZ_AOXE01000081.1_200,NZ_AOXE01000081.1_201,NZ_AOXE01000081.1_209,NZ_AOXE01000081.1_210,NZ_AOXE01000081.1_211,NZ_AOXE01000081.1_212,NZ_AOXE01000081.1_214,NZ_AOXE01000081.1_215,NZ_AOXE01000081.1_220,NZ_AOXE01000081.1_223,NZ_AOXE01000081.1_249,NZ_AOXE01000081.1_251,NZ_AOXE01000081.1_262,NZ_AOXE01000081.1_264,NZ_AOXE01000081.1_267,NZ_AOXE01000081.1_272,NZ_AOXE01000081.1_282,NZ_AOXE01000081.1_283,NZ_AOXE01000081.1_286,NZ_AOXE01000081.1_294,NZ_AOXE01000081.1_40,NZ_AOXE01000081.1_48,NZ_AOXE01000081.1_49,NZ_AOXE01000081.1_52,NZ_AOXE01000081.1_55,NZ_AOXE01000081.1_59,NZ_AOXE01000081.1_62,NZ_AOXE01000081.1_64,NZ_AOXE01000081.1_76,NZ_AOXE01000081.1_79,NZ_AOXE01000081.1_83,NZ_AOXE01000081.1_87,NZ_AOXE01000081.1_92,NZ_AOXE01000081.1_97,NZ_AOXE01000083.1_45,NZ_AOXE01000083.1_47,NZ_AOXE01000083.1_53,NZ_AOXE01000083.1_74,NZ_AOXE01000083.1_86,NZ_AOXE01000085.1_10,NZ_AOXE01000085.1_17,NZ_AOXE01000085.1_20,NZ_AOXE01000085.1_34,NZ_AOXE01000085.1_57,NZ_AOXE01000085.1_58,NZ_AOXE01000085.1_60,NZ_AOXE01000085.1_62,NZ_AOXE01000085.1_63,NZ_AOXE01000085.1_65,NZ_AOXI01000002.1_306,NZ_AOXI01000005.1_72,NZ_AOXI01000016.1_73,NZ_AOYI01000008.1_9,NZ_AOYL01000006.1_89,NZ_AOYO01000084.1_456,NZ_AOYX01000009.1_43,NZ_AOYX01000031.1_11,NZ_AOYX01000060.1_42,NZ_AOYX01000075.1_47,NZ_AOYX01000092.1_135,NZ_APAO01000014.1_55,NZ_AYDA01000043.1_275,NZ_CM001471.1_3941
	LT2,419666160,2853045644,161888011,3634146466,4104237653,1415645483,3477025620,1721939526,3058650972,222074798,1268872669,1327657400,2932852483,35548503,1766268720,1547995867,16813142,802448923,1415808768,2274790611,2777872337,2099395296,2748382047,2612985220,388687368,2980577262,3237112777,2730940254,1884122039,2106993987,1303226364,1885567356,1334856249,2176838220,4150577466,2927589129,2375168468,1280527088,33744532,2586709922,2487576818,1551963583,1499568936,2864140755,254929250,3756506314,1302755412,827295360,1515943483,3460855762,3546715191,1818309659,4174276796,3219814449,4276130360,2623667813,3429770389,2533606710,636455403,2919337619,961839877,1115148926,1917632219,1898052981,246770399,3308702850,1262642563,919676378,1859197400,3656308000,1868074335,2767245212,2545200385,4234954887,2955003506,2353669245,2323569073,2605175336,1116230143,2675920931,3411211333,908026377,4265398385,383184955,122484534,835140794,4287331645,1667877667,3641670366,2040936010,2052683636,3745633893,20195024,3868995612,1305201076,1110352487,1508329141,3448126424,1015550385,2725387546,1944654920,2538869957,1204520418,4178717449,1225424198,3005760350,2832724106,519139030,1267870073,192899265,1781476728,1436385457,1781281347,3925631537,1292522932,3000213477,3045375618,1848379870,371825940,1772513075,3141100755,1558476698,2400355742,1080346469,3606169715,551814723,708436169,3287828105,2463803910,510160738,1204841711,3641973773,1013441200,2514330295,656883055,694036224,4237193625,1196686387,2598101487,761097923,1319443396,82932838,3506531371,1477881595,398274060,683371912,1852780082,3908070206,830225365,4090514893,890190767,2206608985,1925258705,29321829,2223490890,271857015,3092958502,3178611537,3536902160,1817113462,1649013537,663662898,1217930253,3092689021,4044543395,2652912325,3799959984,629769704,574581854,1226961299,876735054,3264162588,2113055247,1796868112,2467894309,3206518183,128195820,99201986,1011589736,302830379,257452301,3352853114,1498827863,943409602,3605841992,2006818648,3689427398,2536338461,400314092,3233075956,2409457196,1409117107,3371366009,1674121234,4239265450,4135632021,3520112373,2229984332,1300322297,3320173704,1175502179,1923155005,1044890724,2920049295,1129751402,3371174263,774642033,1114480804,1854288580,3832041954,3349781307,3517434139,332066877,3565621154,3987728523,3952282231,3078292395,1426386553,181287415,702142682,3713645811,484067521,739645922,594076509,3766280094,1236556171,916567331,2034737222,2838102363,196575977,2615307626,1822361792,1882995852,3419283395,1375971183,428510990,1391862775,125247390,3458312089,340459710,621463099,2734837303,1982008089,4105928256,106320492,2486472784,904262353,2754698655,1017068241,766393622,849365577,4213771231,352954861,1741286159,1417246238,2758229229,1743327588,1917216667,4140525988,4012682490,1812733644,1179850848,2753336226,3221170065,1590132301,2639318025,1024177630,3445065778,3250876267,2308303834,10121991,3254127745,3746351055,3555078814,2116920669,148325698,1739708597,3137569212,808545123,1710702657,2203341188,4254642704,2240911600,2037995616,3625135283,1407835461,4169948067,2414391967,3858965322,488347769,1610885895,2560486208,4157124496,2938422205,411913884,3422434776,233023645,3273153987,976767449,2887404622,2763485497,2649675872,4206264737,1856369146,3972685647,1484137762,1613239859,574411238,3106239735,227014981,2712433588,3346526922,2414482651,1458595125,3308323309,2769765934,728160368,1187450835,340611593,305671293,491087478,3161492497,1671956993,2182926251,486088087,2151397774,1342450405,4071108027,3364789160,2201764151



QC by ``sistr_cmd`` (``--qc``)
-------------------

If you are running ``sistr_cmd`` with the ``--qc`` commandline argument, ``sistr_cmd`` will run some basic QC to determine the level of confidence in the serovar prediction. 

The ``qc_status`` field should contain a value of ``PASS`` if your genome passes all QC checks, otherwise, it will be ``WARNING`` or ``FAIL`` if there are issues with your results and/or input genome sequence.

The ``qc_messages`` field will contain useful information about why you may have a low confidence serovar prediction result. The QC messages will be delimited by `` | ``.

For example, here are the QC messages for an unusually small *Salmonella* assembly where the predicted serovar was "-:-:-":

.. code-block::

	FAIL: Large number of cgMLST330 loci missing (n=272 > 30)
	FAIL: Wzx/Wzy genes missing. Cannot determine O-antigen group/serogroup. Cannot accurately predict serovar from antigen genes.
	WARNING: H1 antigen gene (fliC) missing. Cannot determine H1 antigen. Cannot accurately predict serovar from antigen genes.
	WARNING: Input genome size (699860 bp) not within expected range of 4000000-6000000 (bp) for Salmonella
	WARNING: Only matched 57 cgMLST330 loci. Min threshold for confident serovar prediction from cgMLST is 297.0

The QC messages produced by ``sistr_cmd`` should help you understand your serovar prediction results.


Issues
======

If you encounter any problems or have any questions feel free to create an issue anonymously or not to let us know so we can address it!

Feature requests and pull requests are welcome!


Want to help improve this tool?
===============================

Do you have any *Salmonella* genomes with trustworthy serovar info? Would you like SISTR to provide better serovar predictions? You can help by contributing those genomes along with their serovar info!

SISTR relies on a database of cgMLST allelic profiles from *Salmonella* genomes with validated serovar info to make accurate serovar predictions (since antigenic determinations from a handful of genes like wzx or fliC can only get you so far). So the more genomes there are in the SISTR database, the more accurate the serovar predictions, especially if those genomes belong to uncommon or rare serovars or lineages.

Help us improve SISTR serovar predictions! Contribute *Salmonella* genomes to SISTR!


You can contribute by:

- let us know here: https://github.com/peterk87/sistr_cmd/issues/15
- linking to your genome on NCBI SRA/BioSample/Assembly
- sending us an email at sistr.salmonella@gmail.com
- contacting the authors of SISTR


Development
===========

Getting started::
	
	git clone https://github.com/peterk87/sistr_cmd.git
	cd sistr_cmd/
	export PYTHONPATH=$(pwd)
	# run tests
	py.test tests/

Pull requests for feature additions and bug fixes welcome!


Using ``sistr_cmd`` in your Python application
----------------------------------------------

Want to use ``sistr_cmd`` directly in your Python application?

Install ``sistr_cmd`` using pip or Conda.

You can run SISTR serovar predictions like so:

.. code-block:: python

	from sistr.sistr_cmd import sistr_predict
	# create mock commandline arguments class
	class SistrCmdMockArgs:
	    run_mash = True
	    no_cgmlst = False
	    qc = True
	    use_full_cgmlst_db = False
	# run SISTR serovar prediction
	sistr_results, allele_results = sistr_predict(genome_fasta_path, genome_name, keep_tmp=False, tmp_dir='/tmp/sistr_cmd', args=SistrCmdMockArgs)
	# use sistr_cmd results for something


License
=======

Copyright 2016 Public Health Agency of Canada

Distributed under the GNU Public License version 3.0
