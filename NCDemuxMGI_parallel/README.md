**NCDemuxMGI_parallel.sh**: a custom script to run the MGI off-device demultiplexer in parallel for 4 Lanes

<hr>

### Remarks

The code developed here applies to 10X data sequenced with paired reads of length 28 and 110 bases respectively.

The sequencing includes two barcodes of 10 bases each located at the end of read-2 (between position 90 and 110)

### Running

The script initializes and takes all necessary variable values from the **run_config.yaml** file (using custom functions sourced in the script). **Please review the run_config.yaml file and edit what needs be before running, this code is not performing any check on the provided values**

The sample to barcode list is provided as a two-column tab-separated text file in which column #1 reports sample labels and column#2 the merged string of barcode#1 and barcode#2 (both provided in Forward orientation).


### Requirements

* The script expects the MGI SplitBarcode executable (version 2) to be installed and running on the machine (its path should be edited in run_config.yaml)
* GNU parallel is used to execute 4 jobs, one for each lane

REM: please report any misbehaviour so that I can improve this code.

_SP@NC last edit 2023-02-13_

<hr>

<h4>Please send comments and feedback to <a href="mailto:nucleomics.bioinformatics@vib.be">nucleomics.bioinformatics@vib.be</a></h4>

<hr>

![Creative Commons License](http://i.creativecommons.org/l/by-sa/3.0/88x31.png?raw=true)

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/).
