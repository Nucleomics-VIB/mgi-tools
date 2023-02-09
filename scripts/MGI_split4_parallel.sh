#!/bin/bash

# script: MGI_split4_parallel.sh
# split by barcode from 4 lanes of a MGI run
#
# Stephane Plaisance (VIB-NC) 2023/02/09; v1.0

# find all read1 files in path
find . -name "*_L0?_read_1.fq.gz" | sort > R1.list

# set parameters for SplitDualBarcodes_NC.pl

# -e 2 -f 91 -c Y -rc N 
# -r1 V350134580_1m_L01_read_1.fq.gz -r2 V350134580_1m_L01_read_2.fq.gz 
# -b bc.list -o new1m

# later passed from command
e=2
f=91
c=Y
rc=N
bclist="bc.list"

##########################
# run 4 lanes in parallel
##########################

parallel --plus -j 4 SplitDualBarcodes_NC.pl \
  -e ${e} \
  -f ${f} \
  -c ${c} \
  -rc ${rc} \
  -b ${bclist} \
  -r1 {} \
  -r2 {= s:read_1:read_2: =} \
  -o {= s:.*/::,  s:_read_1.fq.gz:: =}_e${e}_demux ::: $(cat R1.list)
  
######################################
# merge fastq from 4 lanes by barcode 
######################################

# create Lane folder array from R1.list
# each line is stored in ${a[0]} .. ${a[3]}
readarray -t a < R1.list
lanecnt=$(wc -l < R1.list)

# loop through all files in L01 (should also exist in the other 3 Lanes)
l1folder=$(basename ${a[0]%_read_1.fq.gz})_e${e}_demux
readarray -t read1list < <(find "${l1folder}" -name "*_1.fq.gz")

# concatenate same barcode in all Lanes to a single file
outfolder="merged_reads_e${e}"
mkdir -p ${outfolder}

for fq in "${read1list[@]}"; do
    pfx=$(basename ${fq%_1.fq.gz})
    pfx=${pfx/_L01/_merged}
    # empty if existing
    cat /dev/null > ${outfolder}/${pfx}_1.fq.gz
    cat /dev/null > ${outfolder}/${pfx}_2.fq.gz
    
    for l in $(seq 1 ${lanecnt}); do
        r1=${fq//L01/L0${l}}
        r2=${r1/_1.fq.gz/_2.fq.gz}
        # add to merge
        zcat ${r1} >> ${outfolder}/${pfx}_1.fq.gz
        zcat ${r2} >> ${outfolder}/${pfx}_2.fq.gz
    done
done

# switch to R for the merging of the summary files and plotting
