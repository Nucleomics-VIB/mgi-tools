**MGI_mergeStats4.R**: a custom **R** script to merge results from 4-Lanes of MGI demultiplexing and plot the sample distribution of the merged data

<hr>

### Remarks

The code developed here applies to results from the script **NCDemuxMGI_parallel.sh** 

### Running

The script should be run from the folder where the 4-Lanes of demultiplexed data are located. It can also be run from within one of the Lane folders and will then only process a single BarcodeStat.txt file and produce text results and a plot


### Requirements

* The script requires several **R** packages

* data.table
* dplyr
* reshape2
* ggplot2

REM: please report any misbehaviour so that I can improve this code.

_SP@NC last edit 2023-03-01_

<hr>

<h4>Please send comments and feedback to <a href="mailto:nucleomics.bioinformatics@vib.be">nucleomics.bioinformatics@vib.be</a></h4>

<hr>

![Creative Commons License](http://i.creativecommons.org/l/by-sa/3.0/88x31.png?raw=true)

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/).
