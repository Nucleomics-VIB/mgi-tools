#!/bin/bash

# script: SplitDualBarcodes_all_NC.sh
# split by barcode from 4 lanes of a MGI run
#
# REM:
# this script runs a modified version of the original 'SplitDualBarcodes.pl'
# (https://github.com/gateswell/SplitBarcode)
# this modified version can be found in our repo 
# (https://github.com/Nucleomics-VIB/mgi-tools)
#
# Requirements:
# SplitDualBarcodes_NC.pl
# gnu parallel
#
# Stephane Plaisance (VIB-NC) 2023/02/09; v1.0

version="1.00, 2023_02_09"

usage='# Usage: SplitDualBarcodes_all_NC.sh
# -b <path to sample + barcode(s) list file> (required)
# -f <position of the barcode in read2> (required)
# optional -i <path to the folder containing the Lane folders> (default to current folder)
# optional -o <path to save result folders and files> (default to current folder)
# optional -e <max nbr of mismatches with published barcode> (default:0))
# optional -c <compress the outout to fq.gz (Y/N)> (default Y)
# optional -rc <barcodes need be reverse-complemented (Y/N)> (default N)>
# script version '${version}'
# [-h for this help]'

while getopts "b:f:i:o:e:c:rc:h" opt; do
  case $opt in
    b)  optb=${OPTARG} ;;
    f)  optf=${OPTARG} ;;
    i)  opti=${OPTARG} ;;
    o)  opto=${OPTARG} ;;
    e)  opte=${OPTARG} ;;
    c)  optc=${OPTARG} ;;
    rc) optrc=${OPTARG} ;;
    h)  echo "${usage}" >&2; exit 0 ;;
    \?) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    *)  echo "this command requires arguments, try -h" >&2; exit 1 ;;
  esac
done

# disable buffering to get output during long process (loop)
#$|=1;

# test input are provided and present
# check if required arguments were provided
if [ -z "${optb+x}" ] || [ -z "${optf+x}" ]; then
   echo "# please provide mandatory arguments -b, and -f!"
   echo "${usage}"
   exit 1
fi

# test if input files exist
if [ ! -f "${optb}" ]
then
   echo "# please provide mandatory files designated by -b is present!"
   echo "${usage}"
   exit 1
fi

# test if -f position is defined and an integer
re='^[0-9]+$'
if [ -z "${optf+x}" ] || ! [[ ${optf} =~ $re ]]; 
   then echo "error: -f is not defined or not an integer"
   exit 1 
fi

echo "# start $(date)" >&2

# folders
infolder=${opti:-"."}
outfolder=${opto:-"."}

# build a list of all read1 files in path
find ${infolder} -name "*_L0?_read_1.fq.gz" | sort > /tmp/R1.list

# test if list is 1 to 4 lines or die
listlen=$(wc -l < /tmp/R1.list)
if [ -z "${listlen}" ] || [ "${listlen}" -gt 4 ]; then
  echo "# I should find between 1 and 4 lanes of data (read1 files) but instead found ${listlen}!"
  echo "# make sure the input folder contains only 1 up to 4 Lanes of a single run (not demultiplexed)"
  exit 1
fi

echo "# Demultiplexing ${listlen} Lanes in parallel" >&2

# set other parameters for SplitDualBarcodes_NC.pl
b=${optb}
f=${optf}
e=${opte:-0}
c=${optc:-"Y"}
rc=${optrc:-"N"}

##########################
# run 4 lanes in parallel
##########################

parallel --plus -j 4 SplitDualBarcodes_NC.pl \
  -r1 {} \
  -r2 {= s:read_1:read_2: =} \
  -b ${b} \
  -e ${e} \
  -f ${f} \
  -c ${c} \
  -rc ${rc} \
  -o ${outfolder}/{= s:.*/::,  s:_read_1.fq.gz:: =}_e${e}_f${f}_demux ::: $(cat /tmp/R1.list)
  
######################################
# merge fastq from 4 lanes by barcode 
######################################

echo "# Merging read1 and read2 files from ${listlen} Lanes by barcode" >&2

# create Lane folder array from R1.list
# each line is stored in ${a[0]} .. ${a[3]}
readarray -t a < /tmp/R1.list
# count lanes (1..4)
lanecnt=$(wc -l < /tmp/R1.list)

# loop through all files in L01 (same samples should also exist in the other 3 Lanes)
l1folder=${outfolder}/$(basename ${a[0]%_read_1.fq.gz})_e${e}_f${f}_demux
readarray -t read1list < <(find "${l1folder}" -name "*_1.fq.gz")

# concatenate same barcode in all Lanes to a single file (both for r1 and r2)
readfolder="${outfolder}/merged_reads_e${e}_f${f}"
mkdir -p ${readfolder}

for fq in "${read1list[@]}"; do
    pfx=$(basename ${fq%_1.fq.gz})
    pfx=${pfx/_L01/_merged}
    # empty if existing
    cat /dev/null > ${readfolder}/${pfx}_1.fq.gz
    cat /dev/null > ${readfolder}/${pfx}_2.fq.gz
    
    for l in $(seq 1 ${lanecnt}); do
        r1=${fq//L01/L0${l}}
        r2=${r1/_1.fq.gz/_2.fq.gz}
        # add to merge
        zcat ${r1} >> ${readfolder}/${pfx}_1.fq.gz
        zcat ${r2} >> ${readfolder}/${pfx}_2.fq.gz
    done
done

echo "# Merging done!" >&2
echo "# please run MGI_mergeStats4.R to merge summary files and produce a barcode frequency plot" >&2
echo
echo "# end $(date), using version ${version}" >&2
