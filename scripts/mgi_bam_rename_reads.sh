#!/bin/bash

# script mgi_bam_rename_reads.sh
# SP@NC; 2023-06-09; v1.0
# runid is used twice in absence of flowcell id 
# to obtain 7 fields as expected for Illumina

inbam=$1
outbam=${inbam%.bam}_ren.bam

{ samtools view -H ${inbam};
samtools view ${inbam} \
  | bioawk -c sam 'BEGIN{FS="\t"}{
    # V1000 02807 L3 C002 R026 632273
    dev=substr($qname,1,5);
    rid=substr($qname,6,5);
    lan=substr($qname,12,1);
    col=substr($qname,14,3);
    row=substr($qname,18,3);
    til=substr($qname,21,6);
    $qname=dev":"rid":"rid":"lan":"til":"row":"col;
    print $0}' - 
    } | samtools view -b > ${outbam}
