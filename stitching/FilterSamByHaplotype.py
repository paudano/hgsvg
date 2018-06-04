#!/usr/bin/env python
import sys
import argparse

ap = argparse.ArgumentParser(description="Split samfiles by haplotype in title of contig")
ap.add_argument("fofn", help="File of samfile names")
ap.add_argument("--h0", help="H1 samfile.", default="alignments.h0.sam")
ap.add_argument("--h1", help="H0 samfile.", default="alignments.h1.sam")
ap.add_argument("--header", help="Use this header", default="/net/eichler/vol5/home/mchaisso/projects/HGSVG/hgsvg/sv/header.sam")
args = ap.parse_args()

hapFofn = open(args.fofn)
hap0 = open(args.h0,'w')
hap1 = open(args.h1,'w')

headerFile = open(args.header)

# get the header
for line in headerFile:
    if line[0] == '@':
        hap0.write(line)
        hap1.write(line)        
    else:
        break
        

for hapFileName in hapFofn:
    try:    
        hapFile = open(hapFileName.strip())
    except:	
        continue


    for line in hapFile:
        v = line.split()
        hap = v[-1]
        if v[0] == "*":
            continue
        sys.stderr.write(hapFileName)
        if v[0][-4:] == "/0/1":
            hap0.write(line)
            continue
        if v[0][-4:] == "/0/2":
            hap1.write(line)
            continue
        if hap == "HA:i:1":
            hap0.write(line)
        elif hap == "HA:i:2":
            hap1.write(line)

