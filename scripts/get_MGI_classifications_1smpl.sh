#!/bin/bash

# script: get_MGI_classifications_1smpl.sh
# extract classification info from MGI demultiplexed fastq
# add up data from all 4 lanes to a single summary file
# the resulting file can be compared with the output of get_NC_classifications_1smpl.sh
#
# Stéphane Plaisance - VIB-Nucleomics Core - 2023-02-09 v1.00
# visit our Git: https://github.com/Nucleomics-VIB

# L01_DemulX/V350134580_L01_435301ABE103GEX_2.fq.gz
# => 435301ABE103GEX
infq=$1

label=$(basename ${infq} | awk 'BEGIN{FS="_";OFS=""}{print $3}')
# echo ${label}

# initialize results
results="${label}-mgi_classification.txt"
cat /dev/null > ${results}

echo "# processing ${infq}"

bioawk -c fastx -v class="${label}" 'BEGIN{OFS="\t"}{print $name, class}' ${infq} >> ${results}
