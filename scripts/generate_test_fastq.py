#!/usr/bin/env python3

import gzip
import random
from typing import List, Dict

MAXLEN = 90
MAXSEQ = 100

def read_barcodes(file: str) -> List[str]:
    with open(file, 'r') as f:
        return [line.strip() for line in f]

def mutate(sequence: str, num_mutations: int, mutation_char: str) -> str:
    sequence_list = list(sequence)
    for _ in range(num_mutations):
        pos = random.randint(0, len(sequence_list) - 1)
        sequence_list[pos] = mutation_char
    return ''.join(sequence_list)

def generate_polya(num_mutations: int) -> str:
    return 'T' * num_mutations + 'A' * (90 - num_mutations)

def generate_umi() -> str:
    return ''.join(random.choice('ACGT') for _ in range(8))

def get_mutations(mutation_map: Dict[str, str], file: str, num_mutations: int) -> str:
    return mutation_map[file] * num_mutations

def generate_read_r1(bc1: str, bc2: str, bc3: str) -> str:
    read = 'GTCA' + bc1 + 'AACC' + bc2 + 'ACAG' + bc3 + 'CCTA' + 'TTCGAG' + generate_umi() + generate_polya(0)
    return (read + 'A' * MAXLEN)[:MAXLEN]

def generate_read_r2(num_mutations_bc1: int, num_mutations_bc2: int, num_mutations_bc3: int) -> str:
    read = get_mutations(mutation_map, 'bc1.txt', num_mutations_bc1) + get_mutations(mutation_map, 'bc2.txt', num_mutations_bc2) + get_mutations(mutation_map, 'bc3.txt', num_mutations_bc3)
    return (read + 'A' * MAXLEN)[:MAXLEN]

def generate_quality() -> str:
    return 'I' * MAXLEN


def mut() -> int:
    return random.choices([0, 1, 2, 3], weights=[70, 10, 10, 10], k=1)[0]

bc1s = read_barcodes('bc1.txt')
bc2s = read_barcodes('bc2.txt')
bc3s = read_barcodes('bc3.txt')

mutation_map = {'bc1.txt': 'C', 'bc2.txt': 'G', 'bc3.txt': 'T'}



with gzip.open('R1.fastq.gz', 'wt') as r1, gzip.open('R2.fastq.gz', 'wt') as r2:
    for i in range(MAXSEQ):
        num_mutations_bc1 = mut()
        num_mutations_bc2 = mut()
        num_mutations_bc3 = mut()
        bc1 = mutate(bc1s[random.randint(0, len(bc1s) - 1)], num_mutations_bc1, mutation_map['bc1.txt'])
        bc2 = mutate(bc2s[random.randint(0, len(bc2s) - 1)], num_mutations_bc2, mutation_map['bc2.txt'])
        bc3 = mutate(bc3s[random.randint(0, len(bc3s) - 1)], num_mutations_bc3, mutation_map['bc3.txt'])
        read_r1 = generate_read_r1(bc1, bc2, bc3)
        read_r2 = generate_read_r2(num_mutations_bc1, num_mutations_bc2, num_mutations_bc3)
        quality = generate_quality()
        r1.write(f'@SEQ_ID_{i}\n{read_r1}\n+\n{quality}\n')
        r2.write(f'@SEQ_ID_{i}\n{read_r2}\n+\n{quality}\n')