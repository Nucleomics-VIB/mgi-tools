**MGI_mergeStats4.R**: a custom **R** script to merge results from 4-Lanes of MGI demultiplexing and plot the sample distribution of the merged data

<hr>

### Remarks

The code developed here applies to results from the script **NCDemuxMGI_parallel.sh** 

### Running

The script should be run from the folder where the 4-Lanes of demultiplexed data are located.
The resulting files are:

* BarcodeStat_per_lane.pdf
* BarcodeTotalCounts.csv
* BarcodeStat_density.pdf
* BarcodeStat_merged.csv

* SequenceStat_merged.csv
* SequenceStat_top100_undecoded.csv

NOTE: The script can also be run from within one of the Lane folders and will then only process a single BarcodeStat.txt file and produce text results and a plot


### Requirements

The script requires several **R** packages

* data.table
* dplyr
* reshape2
* ggplot2

REM: please report any misbehaviour so that I can improve this code.

_SP@NC last edit 2023-03-01_

<hr>

<h4>Please send comments and feedback to <a href="mailto:nucleomics.bioinformatics@vib.be">nucleomics.bioinformatics@vib.be</a></h4>

<hr>

![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)

Licensed under the **GNU General Public License v3.0** — see [`LICENSE`](../LICENSE) for the
full text, and [`NOTICE.md`](../NOTICE.md) for attribution and the relicensing history.
