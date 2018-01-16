#!/usr/bin/env python

import intervaltree
import argparse
import Tools
import sys

ap = argparse.ArgumentParser(description="Given a SAM file of locally aligned contigs, determine which contigs best overlap the genome segmented by contig alignment positions. ")
ap.add_argument("sambed", help="Alignments of the tiled assembly, translated to bed.")
ap.add_argument("tiling", help="Output file.")
ap.add_argument("--minContigLength", help="Minimum contig length.", default=0,type=int)
args = ap.parse_args()

samFile = open(args.sambed)
outFile = open(args.tiling, 'w')
chromIntervals = {}
chromLengths = {}
lineNumber = 0
chroms = []
for line in samFile:
    vals = line.split()
    if len(vals) > 0 and vals[0] not in chromIntervals:
        chromIntervals[vals[0]] = intervaltree.IntervalTree()

    if len(vals) > 0:
        start = int(vals[1])
        end   = int(vals[2])
        if end > start:
            chromIntervals[vals[0]].addi(int(vals[1]), int(vals[2]), vals[3])

for chrom in chromIntervals:
    intvs = chromIntervals[chrom]
    chromPos = {}
    for intv in intvs.items():
        chromPos[intv[0]] = 1
        chromPos[intv[1]] = 1
    chromPos = chromPos.keys()
    chromPos.sort()
    i = 0
    j = i
    nCondense = 0

    if (len(chromPos) < 2):
        continue

    for i in range(0,len(chromPos)-1):
        midPoint = (chromPos[i] + chromPos[i+1])/2
        if (chromPos[i+1] - chromPos[i] < 2):
            continue
        ovpMidPoint = intvs.search(midPoint)
        opt = None
        optDist = 0
        if (len(ovpMidPoint) == 0):
            continue
        for intv in ovpMidPoint:
            dist = abs(((intv[1]+intv[0])/2) - midPoint)
            if (opt == None or dist < optDist):
                optDist = dist
                opt = intv[2]
        outFile.write(chrom + "\t" + str(chromPos[i]) + "\t" + str(chromPos[i+1]) + "\t" + opt + "\t" + str(optDist) + "\n")
    
    




