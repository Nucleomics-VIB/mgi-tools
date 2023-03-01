#!/bin/bash

# script: NCDemuxMGI_parallel.sh
#
# aim: demultiplex a apired-end MGI run with 4 parallel jobs for speed-up
# uses variable values stored in run_config.yaml
# uses custom functions from /shared/custom_functions.sh
# uses gnu parallel

# Stéphane Plaisance - VIB-Nucleomics Core - 2023-03-01 v1.0

################
# initialisation
################

scriptname=$(basename "$0")
scriptversion='1.0_2023-03-01'

usage='# Usage: '${scriptname}' ...
#
# script version '${scriptversion}'
#
# [optional: -c <path to run_config.yaml (default to run_config.yaml in current folder)>]
# [optional: -l => debug: lists all variables and ends]
# [optional: -h <this help text>]'

# handle getops
while getopts "c:lh" opt; do
	case $opt in
		c) opt_config=${OPTARG};;
		l) opt_listactions=1;;
		h) echo "${usage}" >&2; exit 0;;
		\?) echo "# Invalid option: -${OPTARG}" >&2; exit 1;;
		*) echo "# this command requires arguments, try -h" >&2; exit 1;;
	esac
done
shift $(( OPTIND - 1 ))

# the quote symbols can be nice to have in complex calls
# use them as ${q} and ${d} where needed (hint: sqlite queries)
q="'"
d='"'

# source shared functions from the shared subfolder (required!)
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

[[ -n $(find ${SCRIPT_DIR}/shared -maxdepth 1 -type f -name '*_functions.sh') ]] \
	&& { for f in ${SCRIPT_DIR}/shared/*_functions.sh; do . $f; done; } \
	|| { echo "# *_functions.sh not found in ${SCRIPT_DIR}/shared/"; exit 1; }

# read yaml configuration and process it
run_config=${opt_config:-"run_config.yaml"}
process_config ${run_config}

# list defined variables and die
if [ -n "${opt_listactions}" ]
then
	echo "#-----------------------------------------------------------------------------"
	echo "# List of defined variables:"
	echo "#-----------------------------------------------------------------------------"
	parse_yaml ${run_config} "CONF_"
	echo "#-----------------------------------------------------------------------------"
	exit 0
fi

##############################################
# Run SplitBarcode_v2 in parallel on all 4 lanes
################################################

echo "# start $(date)" >&2

# folders
if [ ! -d "${CONF_outfolder}" ]; then
	echo "# creating ${CONF_outfolder} folder";
	mkdir -p ${CONF_outfolder}
fi

# save all future output to log files
log=${CONF_outfolder}/runlog_$(date +%s).txt
exec &> >(tee -i ${log})

# build a list of all read1 files in CONF_infolder
find ${CONF_infolder} -name "*_L0?_read_1.fq.gz" | sort > /tmp/R1.list

# test if list is 1 to 4 lines or die
listlen=$(wc -l < /tmp/R1.list)
if [ -z "${listlen}" ] || [ "${listlen}" -gt 4 ]; then
  echo "# I should find between 1 and 4 lanes of data (read1 files) but instead found ${listlen}!"
  echo "# make sure the input folder contains only 1 up to 4 Lanes of a single run (not demultiplexed)"
  exit 1
fi

echo "# Demultiplexing ${listlen} Lanes in parallel" >&2

# build b1 and b2 from run_config variables
b1=$(( ${CONF_read1_length}+${CONF_barcode1_atpos}-1 ))" "${CONF_barcode1_length}" "${CONF_barcode1_maxerr}
b2=$(( ${CONF_read1_length}+${CONF_barcode2_atpos}-1 ))" "${CONF_barcode2_length}" "${CONF_barcode2_maxerr}

# test if -r needed
if [ ${CONF_revcomp} = "Y" ]; then
  revcomp="-r"
else
  revcomp=""
fi

# test if barcode should be printed in the read name
if [ ${CONF_addumi} = "Y" ]; then
  addumi="--umi"
else
  addumi=""
fi

################################
# run up to 4 lanes in parallel
################################

parallel --plus -j 4 ${CONF_SplitBarcode_v2} \
  -1 {} \
  -2 {= s:read_1:read_2: =} \
  -B ${CONF_barcodelist} \
  -b ${b1} \
  -b ${b2} \
  ${revcomp} \
  -n ${CONF_maxthr} \
  -m ${CONF_maxmem} \
  ${addumi} \
  -o ${CONF_outfolder}/{= s:.*/::,  s:_read_1.fq.gz:: =}_demux ::: $(cat /tmp/R1.list)

echo "# end $(date), using version ${version}" >&2

exit 0
