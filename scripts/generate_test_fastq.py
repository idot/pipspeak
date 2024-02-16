import random
from typing import List

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

def generate_read(bc1: str, bc2: str, bc3: str, num_mutations: int) -> str:
    return 'GTCA' + bc1 + 'AACC' + bc2 + 'ACAG' + bc3 + 'CCTA' + 'TTCGAG' + generate_umi() + generate_polya(num_mutations)

def generate_quality() -> str:
    return 'I' * 90

bc1s = read_barcodes('bc1.txt')
bc2s = read_barcodes('bc2.txt')
bc3s = read_barcodes('bc3.txt')

mutation_map = {'bc1.txt': 'C', 'bc2.txt': 'G', 'bc3.txt': 'T'}

with open('R1.fastq', 'w') as r1, open('R2.fastq', 'w') as r2:
    for i in range(96):
        num_mutations = random.randint(0, 3)
        bc1 = mutate(bc1s[random.randint(0, len(bc1s) - 1)], num_mutations, mutation_map['bc1.txt'])
        bc2 = mutate(bc2s[random.randint(0, len(bc2s) - 1)], num_mutations, mutation_map['bc2.txt'])
        bc3 = mutate(bc3s[random.randint(0, len(bc3s) - 1)], num_mutations, mutation_map['bc3.txt'])
        read = generate_read(bc1, bc2, bc3, num_mutations)
        quality = generate_quality()
        r1.write(f'@SEQ_ID_{i}\n{read}\n+\n{quality}\n')
        r2.write(f'@SEQ_ID_{i}\n{read}\n+\n{quality}\n')