#!/usr/bin/env python3
"""Modifies a BRAKER GFF3 to work with GenomeTools and the
   GFF3toolkit."""

import sys

from Bio import SeqIO

IN_FASTA = sys.argv[1]
IN_GFF3 = sys.argv[2]

sys.stdout.write('##gff-version 3\n')

for rec in SeqIO.parse(IN_FASTA, 'fasta'):
    rec.id = rec.id.split(' ')[0]
    sys.stdout.write('##sequence-region   {} 1 {}\n'.format(rec.id,
                                                            len(rec.seq)))

GENE_DICT = {}
MRNA_DICT = {}

with open(IN_GFF3, 'r') as in_handle:
    for line in in_handle:
        line = line.rstrip().rstrip(';').split('\t')
        gff_src = line[1]
        gff_type = line[2]
        gff_beg = int(line[3])
        gff_end = int(line[4])
        gff_attr_dict = {k:v for k, v in [a.split('=') for a in
                                          line[8].split(';')]}
        gff_id = gff_attr_dict['ID']
        gff_par = gff_attr_dict['Parent'].split('.')[0]
        if gff_type == 'mRNA':
            line[8] = '{};Name={}'.format(line[8], gff_id)
        sys.stdout.write('{}\n'.format('\t'.join(line))) 
        if gff_src == 'AUGUSTUS':
            if gff_par not in GENE_DICT:
                GENE_DICT[gff_par] = [
                    line[0], line[1], 'gene', gff_beg, gff_end, '.', line[6],
                    '.', 'ID={};Name={}'.format(gff_par, gff_par)
                ]
            else:
                if GENE_DICT[gff_par][3] > gff_beg:
                    GENE_DICT[gff_par][3] = gff_beg
                if GENE_DICT[gff_par][4] < gff_end:
                    GENE_DICT[gff_par][4] = gff_end
            if gff_type == 'CDS':
            	line[7] = '.'
            	sys.stdout.write(
                    '{}\n'.format('\t'.join(line).replace('CDS', 'exon'))
                )
        elif gff_src == 'GeneMark.hmm':
            if gff_par not in GENE_DICT:
                GENE_DICT[gff_par] = [
                    line[0], line[1], 'gene', gff_beg, gff_end, '.', line[6],
                    '.', 'ID={};Name={}'.format(gff_par[:-2], gff_par[:-2])
                ]
                MRNA_DICT[gff_par] = [
                    line[0], line[1], 'mRNA', gff_beg, gff_end, '.', line[6],
                    '.', 'ID={};Parent={};Name={}'.format(
                        gff_par, gff_par[:-2], gff_par
                    )
                ]
            else:
                if GENE_DICT[gff_par][3] > gff_beg:
                    GENE_DICT[gff_par][3] = gff_beg
                    MRNA_DICT[gff_par][3] = gff_beg
                if GENE_DICT[gff_par][4] < gff_end:
                    GENE_DICT[gff_par][4] = gff_end
                    MRNA_DICT[gff_par][4] = gff_end

for gene in GENE_DICT:
    sys.stdout.write('{}\n'.format('\t'.join(str(f) for f in GENE_DICT[gene])))
for mrna in MRNA_DICT:
    sys.stdout.write('{}\n'.format('\t'.join(str(f) for f in MRNA_DICT[mrna])))
