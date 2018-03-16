Analysis code and summary data for:

Mandy S. Wong, Justin B. Kinney\*, Adrian R. Krainer\* (2018) Quantitative activity provide and context dependence of all human 5' splice sites. Submitted. \*Equal contribution.

The computational pipeline used in this study was split into two parts: a "cluster" component run on the Black N Blue High (BNB) Performance Compute Cluster at Cold Spring Harbor Laboratory (CSHL), and a "local" component suitable for execution on a standard laptop computer. 

Code for the cluster component is provided in the cluster/ directory. This code was configured specifically for the BNB cluster at CSHL, and will likely have to be modified before it is run on a different cluster. To run this  pipeline, first download all Illumina sequencing data from SRA BioProject number PRJNA420342 and deposit into the directory cluster/data/illumina_runs/. These file names should match the file names listed in cluster/data/metadata.xlsx. The pipeline can then be run by executing

$ cd cluster

$ python2 run_pipeline.py

Note that this pipeline should be run under Python 2.7.11. 

Code for the local component is in the local/ directory. Unlike the cluster pipeline, this code should be run using Python >=3.6.3. First, copy the results from the cluster/saved directory to the local/from_pipeline directory. Next, download the human genome (hg38.fa) from the UCSC genome browser and decompress it. Then do

$ cd local

$ python3 run_pipeline.py

This generates textual output (deposited in local/output/) as well as graphical output (deposited in local/plots). 
