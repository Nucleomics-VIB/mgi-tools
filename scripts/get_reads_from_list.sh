#!/bin/bash

# script: get_reads_from_list.sh
# extract a list of reads from fastq paired files
# Stéphane Plaisance - VIB-Nucleomics Core - 2023-02-09 v1.00
# visit our Git: https://github.com/Nucleomics-VIB

# depends on seqtk: https://github.com/lh3/seqtk

infq=$1
list=$2

for r in 1 2; do

# V350134580_L01_read_1.fq.gz
fastq=${infq/read_1/read${r}}

seqtk subseq infq list | bgzip -c > subset_${r}.fq.gz

done
