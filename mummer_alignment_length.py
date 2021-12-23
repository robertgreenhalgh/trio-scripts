#!/usr/bin/env python3
"""Determines the length of aligned MUMmer sequences."""

import sys

IN_FILTER = sys.argv[1]

RAW_COORDS = []

with open(IN_FILTER, 'r') as in_handle:
    for line in in_handle:
        line = line.rstrip().split(' ')
        if len(line) == 7:
            start = int(line[2])
            end = int(line[3])
            RAW_COORDS.append((start, end))

FILT_COORDS = []

PREV_START = None
PREV_END = None
for coords in sorted(RAW_COORDS):
    start = min(coords)
    end = max(coords)
    if PREV_START == None:
        PREV_START = start
        PREV_END = end
    else:
        if start <= PREV_END:
            PREV_START = min(PREV_START, start)
            PREV_END = max(PREV_END, end)
        else:
            FILT_COORDS.append((PREV_START, PREV_END))
            PREV_START = start
            PREV_END = end
FILT_COORDS.append((PREV_START, PREV_END))

ALN_LEN = 0

for coord in FILT_COORDS:
    ALN_LEN += (coord[1]-coord[0])+1

print(ALN_LEN)
