#!/bin/bash

# script: get_reads_from_list.sh
# extract a list of reads from fastq paired files
# Stéphane Plaisance - VIB-Nucleomics Core - 2023-02-22 v1.00
# visit our Git: https://github.com/Nucleomics-VIB

# depends on seqtk: https://github.com/lh3/seqtk present in a conda env

# L01/V350134580_L01_read_1.fq.gz: @V350134580L1C001R00100000912/1

myenv=seqtk
source /etc/profile.d/conda.sh
conda activate ${myenv} || \
  ( echo "# the conda environment ${myenv} was not found on this machine" ;
    echo "# please read the top part of the script!" \
    && exit 1 )

fastq1="/data/analyses/MGI_demux/V350134580/L01/V350134580_L01_read_1.fq.gz"
fastq2=${fastq1/_1/_2}

# run both for each unique list
for run in "nc_only" "mgi_only"; do

cp ${run}_reads.list ${run}_reads2.list
cat ${run}_reads.list | sed -e 's:/2:/1:' > ${run}_reads1.list

readlist1="/data/analyses/MGI_demux/compare/${run}_reads1.list"
readlist2="/data/analyses/MGI_demux/compare/${run}_reads2.list"

outfastq1=${run}_$(basename ${fastq1})
outfastq2=${run}_$(basename ${fastq2})

readlist1="/data/analyses/MGI_demux/compare/${run}_reads1.list"
readlist2="/data/analyses/MGI_demux/compare/${run}_reads2.list"

# extract read_1
seqtk subseq ${fastq1} ${readlist1} | bgzip -@9 -c > ${outfastq1} &

# extract read_2
seqtk subseq ${fastq2} ${readlist2} | bgzip -@9 -c > ${outfastq2} &

done
