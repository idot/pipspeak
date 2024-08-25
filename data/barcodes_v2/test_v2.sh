#!/bin/bash

#because the first base in spaceer2 is ambigous, I 

../../target/debug/pipspeak --config config_v2.yaml \
                            --prefix SRR19180490_pipspeak \
                            --r1 SRR19180490_1.head.fastq.gz \
                            --r2 SRR19180490_2.head.fastq.gz  \
                            --umi-len 12