#!/usr/bin/env python3
"""Joins records in a multi-FASTA file into a single sequence."""

import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

IN_FASTA = sys.argv[1]
OUT_NAME = sys.argv[2]
OUT_PAD = int(sys.argv[3]) if len(sys.argv) > 3 else 0

OUT_SEQ = []
OUT_DESC = []

for rec in SeqIO.parse(IN_FASTA, 'fasta'):
    OUT_DESC.append(rec.id)
    OUT_SEQ.append(str(rec.seq))

OUT_REC = SeqRecord(
    Seq((OUT_PAD*'N').join(OUT_SEQ)),
    id = OUT_NAME,
    description = ' '.join(OUT_DESC)
)

SeqIO.write(OUT_REC, sys.stdout, 'fasta')
