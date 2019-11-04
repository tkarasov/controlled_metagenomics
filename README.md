### Max Planck Institute for Developmental Biology
#### Talia Karasov
This repository contains scripts related to the controlled metagenomic infections performed in 2018 and 2019

## DEPENDENCIES
### Read Mapping of Plant and Metagenome
* [python3](https://www.python.org/download/releases/3.0/)
* [bwa (version 0.7.17-r1188)](https://github.com/lh3/bwa)
* [picard (version 2.17.3)](https://broadinstitute.github.io/picard/)
* [samtools (version 1.6-19-g1c03df6(using htslib 1.6-53-ge539b32)](http://www.htslib.org/)
* [centrifuge (version 1.0.4)](https://ccb.jhu.edu/software/centrifuge/manual.shtml)

### Analysis of metagenomes post-centrifuge
* [R (version 3.5.3)](https://www.r-project.org/)
* Major R packages
  * [vegan (version 2.5.5)](https://cran.r-project.org/web/packages/vegan/index.html)

## SCRIPTS
#### The general pipeline for assessing metagenomic load can be found in this repository though there is a lot of unrelated code here.
* https://github.com/tkarasov/metagenomics_pipeline

The manuscript used the pipeline for centrifuge. 
* [Script to run the whole pipeline](https://github.com/tkarasov/metagenomics_pipeline/blob/master/centrifuge/final_centrifuge_total_pipeline.sh)

#### Script-by-Script rundown of the Pipeline to build the metagenome table
* [Script to remove plant reads from the metagenomes](https://github.com/tkarasov/metagenomics_pipeline/blob/master/centrifuge/run_plantRemoval_tlk_centrifuge.sh)
* [Script to run centrifuge](https://github.com/tkarasov/metagenomics_pipeline/blob/master/centrifuge/centrifuge_total_step1_v2.sh)
* [Script to make kreport](https://github.com/tkarasov/metagenomics_pipeline/blob/master/centrifuge/centrifuge_total_step2_v2.sh)
* [Script to recalibrate kreport](https://github.com/tkarasov/metagenomics_pipeline/blob/master/centrifuge/centrifuge_total_step3_v2.sh)
* [Script to subsample centrifuge output to look at normalized data](https://github.com/tkarasov/metagenomics_pipeline/blob/master/centrifuge/subsample_centrifuge_table.py)
* [Script to separate fungal, oomycete and bacterial reads output from centrfiuge](https://github.com/tkarasov/metagenomics_pipeline/blob/master/centrifuge/classify_eukaryote_prokaryote.py)
* [The actual analysis script to recalibrate the output from centrifuge to genome coverage](https://github.com/tkarasov/metagenomics_pipeline/blob/master/centrifuge/recalibrate_metagenome_table_centrifuge_v2.py)

#### Scripts for Figure generation and statistics
* [Folder for scripts for figure generation] (https://github.com/tkarasov/controlled_metagenomics/tree/master/scripts/keep_used_in_publication)
* [Figure 1: Stacked barplots](https://github.com/tkarasov/controlled_metagenomics/blob/master/scripts/keep_used_in_publication/figure1_stacked_barplots.Rmd)
* [Figure 2: Diversity measurements](https://github.com/tkarasov/controlled_metagenomics/blob/master/scripts/keep_used_in_publication/figure2_diversity_real.Rmd)
* [Figure 2: PCoA](https://github.com/tkarasov/controlled_metagenomics/blob/master/scripts/keep_used_in_publication/figure1_2_pcoa_statistics.Rmd)
* [Figure 3: Lab comparisons (DC3000)](https://github.com/tkarasov/controlled_metagenomics/blob/master/scripts/keep_used_in_publication/figure3_dc3000_reads_metagenome.Rmd)
* [Figure 3: Lab comparisons (HpA)](https://github.com/tkarasov/controlled_metagenomics/blob/master/scripts/keep_used_in_publication/figure3_hpa_reads_metagenome.Rmd)
* [Deseq: Variance stabilization between Sweden and Germany](https://github.com/tkarasov/controlled_metagenomics/blob/master/scripts/keep_used_in_publication/figureX_deseq_fin.Rmd)


## DATA 
#### For metagenome classification and normalization
* [The ncbi taxonomy tree used for family-wise calculation of load](https://github.com/tkarasov/metagenomics_pipeline/blob/master/data/megan_genus_tree_10_2_2018.tre)
* [The ncbi genome size file](https://github.com/tkarasov/metagenomics_pipeline/blob/master/data/genomes.csv)

### Field trials
* Swedish metagenome table for oomycete:
  * [Bacteria](https://github.com/tkarasov/controlled_metagenomics/blob/master/data/swedish_meta_family_corrected_per_plant_bacteria.csv)
  * [Oomycetes](https://github.com/tkarasov/controlled_metagenomics/blob/master/data/swedish_meta_family_corrected_per_plant_oomycete.csv)
  * [Fungi](https://github.com/tkarasov/controlled_metagenomics/blob/master/data/swedish_meta_family_corrected_per_plant_fungi.csv)

* German metagenome table
  * [Bacteria](https://github.com/tkarasov/controlled_metagenomics/blob/master/data/german_meta_family_corrected_per_plant_bacteria.csv)
  * [Oomycetes](https://github.com/tkarasov/controlled_metagenomics/blob/master/data/german_meta_family_corrected_per_plant_oomycete.csv)
  * [Fungi](https://github.com/tkarasov/controlled_metagenomics/blob/master/data/german_meta_family_corrected_per_plant_fungi.csv)

#### Laboratory infection trials
* Pst DC3000 infections
   * [Bacteria](https://github.com/tkarasov/controlled_metagenomics/blob/master/data/dc3000_meta_family_corrected_per_plant_bacteria.csv)
   * [Oomycetes](https://github.com/tkarasov/controlled_metagenomics/blob/master/data/dc3000_meta_family_corrected_per_plant_oomycete.csv)
   * [Fungi](https://github.com/tkarasov/controlled_metagenomics/blob/master/data/dc3000_meta_family_corrected_per_plant_fungi.csv)
 
* Hpa infections
   * [Bacteria](https://github.com/tkarasov/controlled_metagenomics/blob/master/data/hpa_meta_family_corrected_per_plant_bacteria.csv)
   * [Oomycetes](https://github.com/tkarasov/controlled_metagenomics/blob/master/data/hpa_meta_family_corrected_per_plant_oomycete.csv)
   * [Fungi](https://github.com/tkarasov/controlled_metagenomics/blob/master/data/hpa_meta_family_corrected_per_plant_fungi.csv)
