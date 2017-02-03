#!/usr/bin/env python



import argparse

ap = argparse.ArgumentParser(description="Remove homopolymer/indel artifacts from gap bed")
ap.add_argument("--bed", help="Gap bed file (with header)", required=True)
ap.add_argument("--maxLength", help="maximum length to filter.", default=0, type=int)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()


gapBedFile = open(args.bed)
outFile = open(args.out,'w')
headerLine = gapBedFile.readline()
outFile.write(headerLine)
headerVals = headerLine[1:].split()
header = {headerVals[i] : i for i in range(0,len(headerVals)) }

for line in gapBedFile:
    vals = line.split()
    seq = vals[header["svSeq"]]
    seq = seq.upper()
    counts = [seq.count(n) for n in ["A", "C", "T", "G"]]
    if args.maxLength == 0 or len(seq) >= args.maxLength:
        outFile.write(line)
        continue

    if max(counts) != len(seq):
        outFile.write(line)

