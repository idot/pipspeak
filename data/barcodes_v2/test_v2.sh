#!/bin/bash

../../target/debug/pipspeak --config data/barcodes_v2/config_v2.yaml \
                            --prefix SRR19180490_pipspeak \
                            --r1 data/barcodes_v2/SRR19180490_1.head.fastq.gz \
                            --r2 data/barcodes_v2/SRR19180490_2.head.fastq.gz  \
                            --umi-len 12