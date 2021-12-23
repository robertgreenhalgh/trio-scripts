#!/usr/bin/env python3
"""Extracts and manually orients sequences from a FASTA file."""

import sys

from Bio import SeqIO

IN_FASTA = sys.argv[1]
IN_ORDER = sys.argv[2:]

FASTA_DICT = {}

for rec in SeqIO.parse(IN_FASTA, 'fasta'):
    if rec.id in IN_ORDER:
        FASTA_DICT[rec.id] = rec
    elif '*{}'.format(rec.id) in IN_ORDER:
        rec.id = '*{}'.format(rec.id)
        rec.seq = rec.seq.reverse_complement()
        FASTA_DICT[rec.id] = rec

if all(rec in FASTA_DICT for rec in IN_ORDER):
    for rec in IN_ORDER:
        SeqIO.write(FASTA_DICT[rec], sys.stdout, 'fasta')
else:
    for rec in [r for r in IN_ORDER if r not in FASTA_DICT]:
        print('Not in FASTA: {}'.format(rec))
