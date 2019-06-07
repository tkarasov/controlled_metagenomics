### Max Planck Institute for Developmental Biology
##### Talia Karasov
This repository contains scripts related to the controlled metagenomic infections performed in 2018 and 2019

## DEPENDENCIES
### Read Mapping of Plant and Metagenome
* python3: https://www.python.org/download/releases/3.0/
* bwa (version 0.7.17-r1188): https://github.com/lh3/bwa
* picard: https://broadinstitute.github.io/picard/
* samtools (version 1.6-19-g1c03df6 (using htslib 1.6-53-ge539b32)): http://www.htslib.org/
* centrifuge (version 1.0.4): https://ccb.jhu.edu/software/centrifuge/manual.shtml

### Analysis of metagenomes post-centrifuge
* R (version 3.5.3: https://www.r-project.org/
* Major R packages
  * Vegan (version 2.5.5): https://cran.r-project.org/web/packages/vegan/index.html

## SCRIPTS
#### The general pipeline for assessing metagenomic load can be found in this repository
* https://github.com/tkarasov/metagenomics_pipeline

The manuscript used the pipeline relevant for centrifuge. 
* [Script to run the whole pipeline](https://github.com/tkarasov/metagenomics_pipeline/blob/master/centrifuge/centrifuge_total_pipeline.sh)

#### Script-by-Script rundown of the Pipeline to build the metagenome table
* [Script to remove plant reads from the metagenomes](https://github.com/tkarasov/metagenomics_pipeline/blob/master/centrifuge/run_plantRemoval_tlk_centrifuge.sh)
* [Script to run centrifuge](https://github.com/tkarasov/metagenomics_pipeline/blob/master/centrifuge/centrifuge_db.sh)
* [Script to separate fungal, oomycete and bacterial reads output from centrfiuge](https://github.com/tkarasov/metagenomics_pipeline/blob/master/centrifuge/classify_eukaryote_prokaryote.py)
* [Script to recalibrate metagenomic output from centrifuge](https://github.com/tkarasov/metagenomics_pipeline/blob/master/centrifuge/recalibrate_metagenome_table_centrifuge.py)

#### Scripts for Figure generation

## DATA 
#### For metagenome classification and normalization
* [The ncbi taxonomy tree used for family-wise calculation of load](https://github.com/tkarasov/metagenomics_pipeline/blob/master/data/megan_genus_tree_10_2_2018.tre)
* [The ncbi genome size file](https://github.com/tkarasov/metagenomics_pipeline/blob/master/data/genomes.csv)

#### LABORATORY INFECTION TRIALS
* Pst DC3000 infections
* Hpa infections
