#!/usr/bin/env python

# script rename_mgi_reads.py
# SP@NC 2023-06-02, v1.0

import sys
import re
import gzip
import io

def rename_mgi_reads(file_path):
    # MGI pattern
    # V1000 02807 L3 C002 R026 632273
    mgi_pattern = r'@V([0-9]{4})([0-9]{5})L([1-4])C([0-9]{3})R([0-9]{3})([0-9]+)'

    # Process the file as a stream
    if file_path == '-':
        # Read from stdin
        data = io.TextIOWrapper(sys.stdin.buffer, encoding='utf-8')
    else:
        # Read from file
        if file_path.endswith('.gz'):
            # Read from gzipped file
            data = gzip.open(file_path, 'rt', encoding='utf-8')
        else:
            # Read from plain file
            data = open(file_path, 'r')

    with data as f:
        for line_num, line in enumerate(f, start=0):
            line = line.strip()

            # Modify the first line and print it
            if line_num % 4 == 0:
                match = re.match(mgi_pattern, line)
                if match:
                    illumina_name = f'@{match.group(1)}:{match.group(2)}:{match.group(3)}:{match.group(6)}:{match.group(4)}:{match.group(5)}'
                    print(illumina_name)
                else:
                    print(line)

            # Print other lines as-is
            else:
                print(line)

###############
# running part
###############

if len(sys.argv) < 2:
    print("Please provide the input file name as a command-line argument or pipe the data using '-' as the file name.")
    sys.exit(1)

file_path = sys.argv[1]
rename_mgi_reads(file_path)