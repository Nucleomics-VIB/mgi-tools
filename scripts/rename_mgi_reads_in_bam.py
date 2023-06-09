#!/usr/bin/env python

# script rename_mgi_reads_in_bam.py
# SP@NC 2023-06-02, v1.0

import sys
import os
import io
import re
import pysam

# Check if an argument is provided
if len(sys.argv) < 2:
    print("Error: No input BAM file provided.")
    sys.exit(1)

# Get the input file path from the command-line argument
input_bam_file = sys.argv[1]

def rename_read_names(input_bam, output_bam, regex_pattern):
    # Open the input BAM file
    input_samfile = pysam.AlignmentFile(input_bam, "rb")

    # Create a new BAM file for writing
    output_samfile = pysam.AlignmentFile(output_bam, "wb", header=input_samfile.header)

    # Iterate over the reads in the input BAM file
    for read in input_samfile.fetch(until_eof=True):
        # Extract the read name
        read_name = read.query_name

        # Apply the regex pattern to extract substrings
        match = re.search(regex_pattern, read_name)

        if match:
            # Capture the substrings using regex capture groups
            captured_substrings = match.groups()

            # Perform any required modifications on the captured substrings
            new_read_name = f'V{match.group(1)}:{match.group(2)}:{match.group(2)}:{match.group(3)}:{match.group(6)}:{match.group(4)}:{match.group(5)}'

            # Update the read name
            read.query_name = new_read_name

        # Write the modified read to the output BAM file
        output_samfile.write(read)

    # Close the BAM files
    input_samfile.close()
    output_samfile.close()

# Running part
output_bam_file = "output.bam"
regex_pattern = r'V([0-9]{4})([0-9]{5})L([1-4])C([0-9]{3})R([0-9]{3})([0-9]+)'

rename_read_names(input_bam_file, output_bam_file, regex_pattern)
