#!/usr/bin/env python3
"""Filters a GFF3 to retain only genes with IPR assignments."""

import sys

IN_TSV = sys.argv[1]
IN_GFF3 = sys.argv[2]

MRNA_SET = set()

with open(IN_TSV, 'r') as in_handle:
    for line in in_handle:
        line = line.rstrip('\n').split('\t')
        if len(line) > 11:
            mrna = line[0]
            ipr = line[11]
            if ipr and ipr[:3] == 'IPR':
                MRNA_SET.add(mrna)

with open(IN_GFF3, 'r') as in_handle:
    GFF3_BUFFER = in_handle.readlines()

GENE_SET = set()

for line in GFF3_BUFFER:
    if line[0] != '#':
        line = line.rstrip('\n').split('\t')
        if line[2] == 'mRNA':
            line[8] = line[8].split(';')
            mrna_id = line[8][0][3:]
            gene_id = line[8][1][7:]
            if mrna_id in MRNA_SET:
                GENE_SET.add(gene_id)

PRINT_LINE = True

for line in GFF3_BUFFER:
    if line[0] != '#' and line.rstrip().split('\t')[2] == 'gene':
        gene_id = line.rstrip().split('\t')[8].split(';')[0][3:]
        if gene_id in GENE_SET:
            PRINT_LINE = True
        else:
            PRINT_LINE = False
    if PRINT_LINE:
        sys.stdout.write(line)
