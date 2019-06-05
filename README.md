# Controlled Metagenomics
### Max Planck Institute for Developmental Biology
##### Talia Karasov
This repository contains scripts related to the controlled metagenomic infections performed in 2018 and 2019

## SCRIPTS
#### The general pipeline for assessing metagenomic load can be found in this repository
* https://github.com/tkarasov/metagenomics_pipeline

The manuscript used the pipeline relevant for centrifuge. 
* [Script to run the whole pipeline](https://github.com/tkarasov/metagenomics_pipeline/blob/master/centrifuge/centrifuge_total_pipeline.sh)

#### Script-by-Script rundown of the Pipeline
* [Script to remove plant reads from the metagenomes](https://github.com/tkarasov/metagenomics_pipeline/blob/master/centrifuge/run_plantRemoval_tlk_centrifuge.sh)
* [Script to run centrifuge](https://github.com/tkarasov/metagenomics_pipeline/blob/master/centrifuge/centrifuge_db.sh)
* [Script to separate fungal, oomycete and bacterial reads output from centrfiuge](https://github.com/tkarasov/metagenomics_pipeline/blob/master/centrifuge/classify_eukaryote_prokaryote.py)
* [Script to recalibrate metagenomic output from centrifuge](https://github.com/tkarasov/metagenomics_pipeline/blob/master/centrifuge/recalibrate_metagenome_table_centrifuge.py)

## DATA
