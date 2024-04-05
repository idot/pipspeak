#!/bin/bash

# from bc1.txt
# removed ^GTCA , AACC$
# .*(\w{8}$) $1 in VSCode to truncate to last 8 characters
#
# perl -pe 's/^GTCA|AACC$//g; s/.*(\w{8}$)/$1/' bc1.txt > bc1_custom_primer.txt
# perl -pe 's/^GTCA|AACC$//g;' bc1.txt > bc1_no_fixed.txt
# the shortest was 8, the longest was 12 i.e. 4 characters prefix, default offset is 5


#./generate_test_fastq.py

# runs pipspeak on the generated files
#../target/debug/pipspeak --help

#Usage: pipspeak [OPTIONS] --r1 <R1> --r2 <R2> --config <CONFIG>

#Options:
#  -i, --r1 <R1>            Input file for R1
#  -I, --r2 <R2>            Input file for R2
#  -p, --prefix <PREFIX>    Output file prefix (output files will be named <prefix>_R[12].fq.gz) [default: pipspeak]
#  -t, --threads <THREADS>  Number of threads to use in gzip compression (0 = all threads) [default: 1]
#  -s, --offset <OFFSET>    The amount of nucleotides away from the start of R1 to accept a barcode [default: 5]
#  -c, --config <CONFIG>    The yaml config file describing the file paths of the 4 barcodes and the spacers
#  -u, --umi-len <UMI_LEN>  The length of the UMI [default: 12]
#  -x, --exact              Use exact matching instead of one mismatch
#  -l, --linkers            Include linkers in the output
#  -q, --quiet              Do not write anything to stderr
#  -h, --help               Print help
#  -V, --version            Print version

# Running pipspeak with the specified options
../target/debug/pipspeak -c vbcf_20240216_config.yaml --r1 R1.fastq.gz --r2 r2.fastq.gz --umi-len 8 --prefix test --offset 5
