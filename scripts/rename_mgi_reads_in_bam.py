#!/usr/bin/env python

# script rename_mgi_reads_in_bam.py
# SP@NC 2023-06-02, v1.0

# !!! at work, not validated yet

import sys
import io
import re
import pysam

def rename_read_names(input_bam, output_bam, regex_pattern, replace_pattern):
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
            # For example, you can use replace_pattern to rename the read name
            new_read_name = replace_pattern.format(*captured_substrings)

            # Update the read name
            read.query_name = new_read_name

        # Write the modified read to the output BAM file
        output_samfile.write(read)

    # Close the BAM files
    input_samfile.close()
    output_samfile.close()

# Example usage
input_bam_file = "input.bam"
output_bam_file = "output.bam"
regex_pattern = r'@V([0-9]{4})([0-9]{5})L([1-4])C([0-9]{3})R([0-9]{3})([0-9]+)'
replace_pattern = r"\1:\2:\3:\6:\4:\5"

rename_read_names(input_bam_file, output_bam_file, regex_pattern, replace_pattern)
