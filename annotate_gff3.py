#!/usr/bin/env python3
"""Appends InterProScan and BLAST hit information to a GFF3."""

import sys

from Bio import SeqIO

IN_GFF = sys.argv[1]
IN_FASTA = sys.argv[2]
IN_IPR = sys.argv[3]
IN_DESC = sys.argv[4]
IN_BLAST = sys.argv[5]
OUT_PRE = sys.argv[6]

with open(IN_GFF, 'r') as in_handle:
    GFF3_BUFFER = in_handle.readlines()

MRNA_TO_GENE_DICT = {}
GENE_SRC_DICT = {}

for line in GFF3_BUFFER:
    if line[0] != '#':
        line = line.rstrip('\n').split('\t')
        gff_type = line[2]
        if gff_type == 'mRNA':
            gff_src = line[1]
            gff_attrs = {k:v for k, v in [a.split('=') for a in
                                          line[8].split(';')]}
            mrna_id = gff_attrs['ID']
            gene_id = gff_attrs['Parent']
            MRNA_TO_GENE_DICT[mrna_id] = gene_id
            if gene_id not in GENE_SRC_DICT:
                GENE_SRC_DICT[gene_id] = set()
            GENE_SRC_DICT[gene_id].add(gff_src)

ISO_LEN_DICT = {}

for rec in SeqIO.parse(IN_FASTA, 'fasta'):
    mrna_id = rec.id.split(' ')[0]
    if mrna_id in MRNA_TO_GENE_DICT:
        gene_id = MRNA_TO_GENE_DICT[mrna_id]
        seq_len = len(rec.seq)
        if gene_id not in ISO_LEN_DICT:
            ISO_LEN_DICT[gene_id] = {}
        if seq_len not in ISO_LEN_DICT[gene_id]:
            ISO_LEN_DICT[gene_id][seq_len] = []
        ISO_LEN_DICT[gene_id][seq_len].append(mrna_id)

ISO_NUM_DICT = {}

for gene_id in ISO_LEN_DICT:
    iso_num = 0
    for seq_len in sorted(ISO_LEN_DICT[gene_id], reverse=True):
        for mrna_id in sorted(ISO_LEN_DICT[gene_id][seq_len]):
            iso_num += 1
            ISO_NUM_DICT[mrna_id] = iso_num

IPR_DICT = {}
GO_DICT = {}

with open(IN_IPR, 'r') as in_handle:
    for line in in_handle:
        line = line.rstrip('\n').split('\t')
        if len(line) > 11:
            mrna = line[0]
            ipr = line[11]
            if mrna in MRNA_TO_GENE_DICT:
                gene_id = MRNA_TO_GENE_DICT[mrna]
            else:
                gene_id = None
            if ipr and ipr[:3] == 'IPR' and gene_id:
                if gene_id not in IPR_DICT:
                    IPR_DICT[gene_id] = set()
                IPR_DICT[gene_id].add(ipr)
            if len(line) > 13:
                gos = line[13]
                if gos and gos[:2] == 'GO' and gene_id:
                    if gene_id not in GO_DICT:
                        GO_DICT[gene_id] = set()
                    for go in gos.split('|'):
                        GO_DICT[gene_id].add(go)

DESC_DICT = {}

with open(IN_DESC, 'r') as in_handle:
    for line in in_handle:
        line = line.rstrip('\n').split('\t')
        if len(line) == 3:
            prot, symb, desc = line
        else:
            prot, desc = line
            symb = 'Unknown'
        DESC_DICT[prot] = (symb, desc)

BLAST_DICT = {}

with open(IN_BLAST, 'r') as in_handle:
    for line in in_handle:
        line = line.rstrip('\n').split('\t')
        mrna_id = line[0]
        mrna_hit = line[1]
        e_val = float(line[10])
        gene_id = MRNA_TO_GENE_DICT[mrna_id]
        if mrna_hit in DESC_DICT:
            mrna_desc = DESC_DICT[mrna_hit]
            if gene_id not in BLAST_DICT:
                BLAST_DICT[gene_id] = (e_val, mrna_desc)
            else:
                if (
                    e_val < BLAST_DICT[gene_id][0] or
                    (mrna_desc[0] != 'Unknown' and
                     BLAST_DICT[gene_id][1][0] == 'Unknown') or
                    (
                        mrna_desc[0] and mrna_desc[0][:3] != 'LOC' and
                        BLAST_DICT[gene_id][1][0] and
                        BLAST_DICT[gene_id][1][0][:3] == 'LOC'
                    )
                ):
                    BLAST_DICT[gene_id] = (e_val, mrna_desc)

GENE_NAME_DICT = {}

NEW_MRNA_ID = None

for line in GFF3_BUFFER:
    line = line.rstrip('\n')
    if line[0] != '#':
        gff_fields = line.split('\t')
        gff_type = gff_fields[2]
        gff_attrs = {k:v for k, v in [a.split('=') for a in
                                      gff_fields[8].split(';')]}
        if gff_type == 'gene':
            gene_id = gff_attrs['ID']
            gff_fields[1] = ';'.join(sorted(GENE_SRC_DICT[gene_id]))
            if gene_id in BLAST_DICT:
                blast_symb, blast_hit = BLAST_DICT[gene_id][1]
            else:
                blast_symb, blast_hit = 'Unknown', 'None'
            if blast_symb not in GENE_NAME_DICT:
                GENE_NAME_DICT[blast_symb] = 0
            GENE_NAME_DICT[blast_symb] += 1
            new_gene_id = '{}_{}_{}'.format(
                OUT_PRE, blast_symb.replace('.', '-').upper(),
                str(GENE_NAME_DICT[blast_symb]).zfill(4)
            )
            line = '{}\tID={};Name={};description={}'.format(
                '\t'.join(gff_fields[:8]), new_gene_id, blast_symb, blast_hit
            )
            if gene_id in GO_DICT:
                gos = ','.join(sorted(GO_DICT[gene_id]))
                line = '{};go={}'.format(line, gos)
            if gene_id in IPR_DICT:
                iprs = ','.join(sorted(IPR_DICT[gene_id]))
                line = '{};ipr={}'.format(line, iprs)
        elif gff_type == 'mRNA':
            mrna_id = gff_attrs['ID']
            mrna_iso = ISO_NUM_DICT[mrna_id]
            NEW_MRNA_ID = '{}.{}'.format(new_gene_id, mrna_iso)
            new_mrna_name = '{}.{}'.format(blast_symb, mrna_iso)
            line = '{}\tID={};Parent={};Name={}'.format(
                '\t'.join(gff_fields[:8]), NEW_MRNA_ID, new_gene_id,
                new_mrna_name
            )
        else:
            if NEW_MRNA_ID:
                line = line.replace(mrna_id, NEW_MRNA_ID)
    print(line)
