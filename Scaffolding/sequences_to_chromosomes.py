#!/usr/bin/env python3
"""Transforms a GFF3 to match scaffolded sequences."""

import sys

from Bio import SeqIO

IN_PAD = int(sys.argv[1])
IN_MAP = sys.argv[2]
IN_FASTA = sys.argv[3]
IN_GFF3 = sys.argv[4]

CHROM_DICT = {}
SEQ_DICT = {}

with open(IN_MAP, 'r') as in_handle:
    for line in in_handle:
        if line[0] != '#':
            line = line.rstrip('\n').split('\t')
            chrom, seqs = line[0], line[1:]
            CHROM_DICT[chrom] = [s[1:] for s in seqs]
            for seq in seqs:
                strand, seq = seq[0], seq[1:]
                SEQ_DICT[seq] = [strand]

for rec in SeqIO.parse(IN_FASTA, 'fasta'):
    if rec.id in SEQ_DICT:
        SEQ_DICT[rec.id].append(len(rec.seq))

GFF3_DICT = {}

with open(IN_GFF3, 'r') as in_handle:
    for line in in_handle:
        if line[0] != '#':
            line = line.rstrip('\n').split('\t')
            seq = line[0]
            if seq in SEQ_DICT:
                if seq not in GFF3_DICT:
                    GFF3_DICT[seq] = []
                line[3] = int(line[3])
                line[4] = int(line[4])
                if SEQ_DICT[seq][0] == '-':
                    seq_offset = SEQ_DICT[seq][1]+1
                    new_start = abs(line[4]-seq_offset)
                    new_end = abs(line[3]-seq_offset)
                    line[3] = new_start
                    line[4] = new_end
                    if line[6] == '+':
                        line[6] = '-'
                    elif line[6] == '-':
                        line[6] = '+'
                GFF3_DICT[seq].append(line)

for chrom in CHROM_DICT:
    seq_offset = 0
    for seq in CHROM_DICT[chrom]:
        if seq in GFF3_DICT:
            for line in GFF3_DICT[seq]:
                line[0] = chrom
                line[3] = str(line[3]+seq_offset)
                line[4] = str(line[4]+seq_offset)
                print('\t'.join(line))
        seq_offset += SEQ_DICT[seq][1]+IN_PAD
