# BRACNAV
BRACNAV (BRAc Copy Number Alteration Viewer)  is a new tool for calling copy number variations (CNVs) in BRCA1 and BRCA2 genes. 
# Installation
## Installation from zip file
* Download [last BRACNAV release](https://github.com/aakechin/BRACNAV.zip):
* Unzip:
 ```
unzip BRACNAV.zip
cd BRACNAV
 ```
## Installation from github
```
git clone  https://github.com/aakechin/BRACNAV.git
cd BRACNAV
```
# Usage
## Windows 7/10
Run executable file BRACNAV.exe
## Ubuntu
To use graphical interface run:
```
python3 main_BRACNAV.py
```
or in command-line version:
```
python bracnav.py -in file_with_coverage.csv -af file_with_coordinates.csv -out output_file
```
BRACNAV can be used with graphical interface or in command-line version. As an input files both versions use:
* Tab-separated file (TSV) with coverage for each target region (see the detailed description below).
* TSV-file with coordinates for each target region. They should correspond to target regions of the coverage file.
Other options and input files are optional.

## Advanced parameters
As an input arguments they use:
* `-ref` Reference version (hg19 or hg38). BRCA1 and BRCA2 genome coordinates depend on this version.
* `-notclust` Not clust - by default, one of normalization steps includes clusterization of sample coverage values, and normalization is additionally carried out inside sample clusters.
* `-hard_score` Hard score threshold - minimal hard score for large rearrangement detection (default: 9.9)
* `-hard_pvalue` Hard p-value threshold - maximal p-value for hard filtering mutations (default: 0.001)
* `-score` Minimal score for large rearrangement detection (default: 2)
* 
* 
