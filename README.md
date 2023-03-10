# MGI-tools
[(Nucleomics-VIB)](https://github.com/Nucleomics-VIB)
![mgi-tools](pictures/MGI.png) - MGI-Tools
==========

*All tools presented below have only been tested by me and may contain bugs, please let me know if you find some. Each tool relies on dependencies normally listed at the top of the code (cpan for perl and cran for R will help you add them)*

*[[back-to-top](#top)]*  

<hr>

<h4>Please send comments and feedback to <a href="mailto:nucleomics.bioinformatics@vib.be">nucleomics.bioinformatics@vib.be</a></h4>

<hr>

![Creative Commons License](http://i.creativecommons.org/l/by-sa/3.0/88x31.png?raw=true)

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/).


# Demultiplexing a MGI run

## Usage example

### Install the following executables on your machine (and add them to your PATH)

* _SplitBarcode_NC.pl (not yet thoroughly tested!)_
* _SplitBarcodes_all_NC.sh (not yet produced, wil be added later)_

* **SplitDualBarcodes_NC.pl** 
* **SplitDualBarcodes_all_NC.sh**

* **MGI_mergeStats4.R** (depends on several R packages, please check its header)

### Move to the folder where the Lanes of undemultiplexed data have been saved

**Note:** You should have 1-4 Lanes in that folder and should not have another run or processed folders resulting from a former run

* Prepare a barcode_list.txt file with 2-columns (name,barcode) or three columns (name, barcode1, barcode2) depending on the library type
* Identify the position where the barcode is expected in read#2 (used with parameter -f)

* Run the demultiplexing command like:

```
SplitDualBarcodes_all_NC.sh -h
# Usage: SplitDualBarcodes_all_NC.sh
# -b <path to sample + barcode(s) list file> (required)
# -f <position of the barcode in read2> (required)
# optional -i <path to the folder containing the Lane folders> (default to current folder)
# optional -o <path to save result folders and files> (default to current folder)
# optional -e <max nbr of mismatches with published barcode> (default:0))
# optional -c <compress the outout to fq.gz (Y/N)> (default Y)
# optional -rc <barcodes need be reverse-complemented (Y/N)> (default N)
# optional -t <threads for merging> (default=8)
# script version 1.00, 2023_02_09
# [-h for this help]
```

```
SplitDualBarcodes_all_NC.sh -b bc.list -f 91 -o /opt/biotools/tmp -t 8 -e 2
```

The first part of the run will demultiplex each of the four Lanes of data to a folder of barcode read pairs (using a parallel process for each Lane) and two Lane summary stat files. 

**Note:** _Each parallel process is using several threads for archiving/deachiving, and for processing the data._

When this is done, a second part of the script will merge the barcodes pairs into a single read pair for each barcode from the provided list (and a pair of unbarcoded reads). This process is running in parallel to save time, adapt the thread number to fit your server.

* When the run has completed, run the R script:

```
MGI_mergeStats4.R
```

It will take care of:

  * merging the **BarcodeStat.txt** and **TagStat.txt** files from each of the 1-4 Lanes into single stat files. The merged stats are saved to csv files for downstream use.
  * also save a list of the top 100 next barcode pairs based on read counts (barcode combinations without a name in the barcode list)
  * plotting the barcode frequencies from the first merged stat file (BarcodeStat.txt) to a PDF file. 

REM: please report any misbehaviour so that I can improve this code.

_SP@NC last edit 2023-02-13_
