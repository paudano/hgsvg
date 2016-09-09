#!/usr/bin/env python

import sys
import argparse

ap = argparse.ArgumentParser(description="Generate a file of reigons for freebayes.")
ap.add_argument("fai", help="Reference fai file.")
ap.add_argument("--regionSize", help="Region size (100000).", default=100000, type=int)
ap.add_argument("--startPos", help="Starting position, should be either 0 or 1 (1).", type=int, default=1)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--no-decoy", help="Do not print a region for a decoy contig", action='store_true',dest='noDecoy')
args = ap.parse_args()

faiFile = open(args.fai)

outFile = open(args.out,'w')

for line in faiFile:
    v = line.split()
    cl = int(v[1])
    regions = ""
    if (args.noDecoy and v[0].find("decoy") >= 0):
        continue
    for startPos in range(0,cl,args.regionSize):
        start = max(startPos, args.startPos)
        end   = min(startPos + args.regionSize, cl)
        regions = regions + "{}:{}-{}\n".format(v[0], start, end)
    outFile.write(regions)



