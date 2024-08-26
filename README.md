# pipspeak

a CLI tool to whitelist filter pipseq reads and convert them to a 10X-style format.

## Overview

This tool is used to filter PIPSeq reads against their respective barcode whitelists
and then output fastq file formats in the style of 10X reads.

This parses the PIPSeq format, identifies the cell barcodes, and writes out a new
file to resemble the 10X sequence construct to be used with other tools that
have not yet adopted the PIPSeq format.

This will also output a whitelist of all the cell barcodes found to be supplied to
downstream mapping tools.

### PIPSeq v3 Sequence Construct

The PIPSeq sequence constructs are organized in the following way

``` txt
                                        ┌─'illumina_p5:29'
                                        ├─I2.fastq.gz────── ──'index5:8'
                                        ├─'truseq_read1:33'
                                        │                   ┌─'cb1:8'
                                        │                   ├─'linker1:3'
                                        │                   ├─'cb2:6'
                                        │                   ├─'linker2:3'
                                        ├─R1.fastq.gz───────┤
                                        │                   ├─'cb3:6'
─────────────────── ──RNA───────────────┤                   ├─'linker3:5'
                                        │                   ├─'cb4:8'
                                        │                   └─'umi:12'
                                        │                   ┌─'cDNA:98'
                                        ├─R2.fastq.gz───────┤
                                        │                   └─'ligationT:1'
                                        │                   ┌─'ME2:19'
                                        ├─nextera_read2─────┤
                                        │                   └─'s7:15'
                                        ├─I1.fastq.gz────── ──'index7:8'
                                        └─'illumina_p7:24'
```

And so the resulting R1 and R2 files boil down to:

``` txt
# R1
[barcode]ATG[barcode]GAG[barcode]TCGAG[barcode][UMI]

# R2
[cDNA]
```

The cell barcodes come from 4 different whitelists.
The ultimate cell-barcode is one of a combination of each of those lists.

### 10X File Format

The 10X sequence construct is organized in the following way

``` txt
                                        ┌─'illumina_p5:29'
                                        ├─'truseq_read1:10'
                                        │                   ┌─'barcode:16'
                                        ├─R1.fastq.gz───────┤
                                        │                   └─'umi:12'
─────────────────── ──RNA───────────────┤
                                        ├─R2.fastq.gz────── ──'cDNA:98'
                                        ├─'truseq_read2:34'
                                        ├─I1.fastq.gz────── ──'index7:8'
                                        └─'illumina_p7:24'
```

And so the resulting R1 and R2 files boil down to:

``` txt
# R1
[barcode][UMI]

# R2
[cDNA]
```

## Usage

This is a single command CLI tool.
It requires just the R1 and R2 filepaths alongside a configuration yaml
which provides the filepaths for the barcodes and describes the spacers.

For the v3 barcodes you can use the configuration and barcode files in
this github repository under `data/`.

``` bash
pipspeak -c data/config_v3.yaml \
    -i data/example_v3/example_R1.fq.gz \
    -I data/example_v3/example_R1.fq.gz
```

### Outputs

This program will output 3 files per run:

1. `<args.prefix>_R1.fq.gz`: A fastq with the `[barcode][UMI]` construct for all reads passing the whitelist.
2. `<args.prefix>_R2.fq.gz`: An unaltered fastq of the R2 for all reads passing the whitelist.
3. `<args.prefix>_whitelist.txt`: a whitelist of all the barcodes found in the dataset.
4. `<args.prefix>_log.yaml`: A log file containing the filtering statistics of the run.

### Configuration

The configuration yaml is very barebones and looks like the following.
It provides the file paths for the barcodes, and then sets the spacer
sequences.

``` yaml
barcodes:
  bc1: "data/barcodes_v3/fb_v3_bc1.tsv"
  bc2: "data/barcodes_v3/fb_v3_bc2.tsv"
  bc3: "data/barcodes_v3/fb_v3_bc3.tsv"
  bc4: "data/barcodes_v3/fb_v3_bc4.tsv"
spacers:
  s1: "ATG"
  s2: "GAG"
  s3: "TCGAG"
```
```
target/debug/pipspeak --loglevel debug -c data/config_v3.yaml   -i data/example_v3/example_R1.fq.gz  -I data/example_v3/example_R1.fq.gz

  Processed 250 reads, 198 passed filters (79.2000%)                                                                                                  
  parameters:
  offset: 5
  umi_len: 12
  exact_matching: false
  write_linkers: false
  pipspeak_version: 0.1.9
file_io:
  readpath_r1: data/example_v3/example_R1.fq.gz
  readpath_r2: data/example_v3/example_R1.fq.gz
  writepath_r1: pipspeak_R1.fq.gz
  writepath_r2: pipspeak_R2.fq.gz
  whitelist_path: pipspeak_whitelist.txt
statistics:
  total_reads: 250
  passing_reads: 198
  fraction_passing: 0.792
  whitelist_size: 198
  num_filtered_1: 41
  num_filtered_2: 6
  num_filtered_3: 3
  num_filtered_4: 2
  num_filtered_umi: 0
timing:
  timestamp: 2024-08-25 17:40:11.885440 +02:00
  elapsed_time: 0.007741375

```
