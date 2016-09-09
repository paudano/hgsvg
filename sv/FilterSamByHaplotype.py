#!/usr/bin/env python

import argparse

ap = argparse.ArgumentParser(description="Split samfiles by haplotype in title of contig")
ap.add_argument("fofn", help="File of samfile names")
ap.add_argument("--h0", help="H1 samfile.", default="alignments.h0.sam")
ap.add_argument("--h1", help="H0 samfile.", default="alignments.h1.sam")

args = ap.parse_args()

hapFofn = open(args.fofn)
hap0 = open(args.h0,'w')
hap1 = open(args.h1,'w')

headerFile = open("/net/eichler/vol24/projects/structural_variation/nobackups/projects/HGSVG/analysis/header.sam")

# get the header
for line in headerFile:
    if line[0] == '@':
        hap0.write(line)
        hap1.write(line)        
    else:
        break
        

for hapFileName in hapFofn:
    
    hapFile = open(hapFileName.strip())
    for line in hapFile:
        v = line.split()
        hap = v[-1]
        if hap == "HA:i:1":
            hap0.write(line)
        elif hap == "HA:i:2":
            hap1.write(line)

