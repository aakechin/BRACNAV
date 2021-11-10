# BRACNAV
BRACNAV (BRAc Copy Number Alteration Viewer)  is a new tool for calling copy number variations (CNVs) in BRCA1 and BRCA2 genes. 
# Quickstart
BRACNAV can be used with graphical interface or in command-line version. As an input files both versions use:
* Tab-separated file (TSV) with coverage for each target region (see the detailed description below).
* TSV-file with coordinates for each target region. They should correspond to target regions of the coverage file.
Other options and input files are optional.

# Installation

As an input arguments they use:
* Reference version (hg19 or hg38). BRCA1 and BRCA2 genome coordinates depend on this version.
* Not clust - by default, one of normalization steps includes clusterization of sample coverage values, and normalization is additionally carried out inside sample clusters.
* Hard score threshold
