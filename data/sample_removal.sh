#!/usr/bin/sh

#some samples need to be removed if they did not have enough coverage. This is the script that goes through and removes them systematically

#dc3000: *control*
awk 'BEGIN { FS = "," } ; {if ($2 < 30000) print $0}' dc3000_count.txt


#hpa: *control* All good except control

#sweden: No control? All good

#germany: 
#These have too few reads: NextMet175, NextMet193, NextMet194, NextMet55
#These are the controls
#duplicated: NextMet24, NextMet50, NextMet66, NextMet70, NextMet83, NextMet85, NextMet87, NextMet96, NextMet119, NextMet131, NextMet142, NextMet144


for direc in `ls /ebio/abt6_projects9/metagenomic_controlled/data/raw_reads/german_samples`; do echo $direc; my_fil=`ls /ebio/abt6_projects9/metagenomic_controlled/data/raw_reads/german_samples/$direc | grep R2`; for rec in $my_fil; do num_rec=`zcat /ebio/abt6_projects9/metagenomic_controlled/data/raw_reads/german_samples/$direc/$rec | grep "^@" | wc -l`; echo $direc/rec, $num_rec >> german_count.txt; done; done
