#!/bin/bash

# script: get_NC_classifications_1smpl.sh
# extract classification info from NC demultiplexed fastq
# data from SplitDualBarcodes_NC.pl
# the resulting file can be compared with the output of get_MGI_classifications_1smpl.sh
#
# Stéphane Plaisance - VIB-Nucleomics Core - 2023-02-09 v1.00
# visit our Git: https://github.com/Nucleomics-VIB

# V350134580_L01_e1_f91_demux/V350134580_L01_435301ABE103GEX_2.fq.gz
# 35301ABE103GEX

infq=$1

# extract info
fname=$(basename ${infq})
#echo ${fname}

pfx=${fname%_2.fq.gz}
#echo ${pfx}

label=$(echo ${pfx} | awk 'BEGIN{FS="_";OFS=""}{print $3}')
# echo ${label}

# initialize results
results="${label}-nc_classification.txt"
cat /dev/null > ${results}

echo "# processing ${infq}"

bioawk -c fastx -v class="${label}" 'BEGIN{OFS="\t"}{print $name, class, $comment}' ${infq} >> ${results}
