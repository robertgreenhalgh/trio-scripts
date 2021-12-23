#!/usr/bin/env python3
"""Determines the number of SNPs and small indels detected by MUMmer."""

import sys

IN_SNPS = sys.argv[1]

SNP_NUM = 0
INDEL_NUM = 0

P1_POSS = set()
P2_POSS = set()

with open(IN_SNPS) as in_handle:
    for line in in_handle:
        line = line.rstrip().split()
        if len(line) == 15:
            p1_pos = int(line[0])
            p1_nuc = line[1]
            p2_nuc = line[2]
            p2_pos = int(line[3])
            if p1_nuc != '.' and p2_nuc != '.':
                SNP_NUM += 1
            else:
                if p1_pos not in P1_POSS and p2_pos not in P2_POSS:
                    INDEL_NUM += 1
            P1_POSS.add(p1_pos)
            P2_POSS.add(p2_pos)

print('SNPs: {}'.format(SNP_NUM))
print('Indels: {}'.format(INDEL_NUM))
print('Total: {}'.format(SNP_NUM+INDEL_NUM))
