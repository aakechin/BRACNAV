# BRACNAV
BRACNAV (BRAc Copy Number Alteration Viewer)  is a new tool for calling copy number variations (CNVs) in BRCA1 and BRCA2 genes. 
## Installation
### Installation from .zip file
* Download [last BRACNAV release](https://github.com/aakechin/BRACNAV/releases/tag/v1.0):
* Unzip:
 ```
unzip BRACNAV.zip
cd BRACNAV
 ```
### Installation from github
```
git clone  https://github.com/aakechin/BRACNAV.git
cd BRACNAV
```
## Usage
### Windows 7/10
Run executable file **BRACNAV.exe**
### Ubuntu
To use graphical interface run:
```
main_BRACNAV
```
or in command-line version:
```
python bracnav.py -in file_with_coverage.csv -af file_with_coordinates.csv -out output_file
```
BRACNAV can be used with graphical interface or in command-line version. As an input files both versions use:
* Tab-separated file (TSV) with coverage for each target region (see the detailed description below).
* TSV-file with coordinates for each target region. They should correspond to target regions of the coverage file.

*Other options and input files are optional.*
### Input file format
* TSV-file with amplicon coverages should have following columns: 

    Patient#, amplicon#1...
                   
*Patient_ID, Barcodes, Median_Coverage, Number_<30 are optional columns between Patient# and amplicon#1*
### Advanced options
* `-ref` Reference version (hg19 or hg38). BRCA1 and BRCA2 genome coordinates depend on this version (default: hg19).

* `-notclust` Not clust - by default, one of normalization steps includes clusterization of sample coverage values, and normalization is additionally carried out inside sample clusters.
* `-hard_score` Hard score threshold - minimal hard score for large rearrangement detection (default: 9.9)
* `-hard_pvalue` Hard p-value threshold - maximal p-value for hard filtering mutations (default: 0.001)
* `-score` Minimal score for large rearrangement detection (default: 2)
* `-pvalue` Maximal p-value for filtering mutations (default: 0.05)
* `-cov` Minimal median coverage for patients (default: 100)
* `-perm` Number of permutations for calculating p-value (default: 1000)
* `-whole` Value for calculating score when whole exon is covered (default: 1)
* `-part` Value for calculating score when part of exon is covered (default: 0.5)
* `-non` Value for calculating score when exon is not covered (default: 0)
* `-del1` Normalized value of coverage for considering an amplicon as likely deleted (default: 1.3)
* `-del2` Normalized value of coverage for considering an amplicon as probably deleted (default: 1.8)
* `-dupl1` Normalized value of coverage for considering an amplicon as likely duplicated (default: 2.7)
* `-dupl2` Normalized value of coverage for considering an amplicon as probably duplicated (default: 2.4)
* `-delta` Minimal relative difference between two values to be considered as significant (default: 0.05)
