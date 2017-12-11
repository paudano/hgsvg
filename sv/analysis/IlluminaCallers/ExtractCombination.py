#!/usr/bin/env python
import sys
import argparse

ap = argparse.ArgumentParser(description="Select rows with a particular combination")
ap.add_argument("--comb", help="Combination")
ap.add_argument("--minn", help="Minimum number of callers to print.", default=1,type=int)
ap.add_argument("--svBed", help="Input sv file")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()

bedFile = open(args.svBed)

hv = bedFile.readline()[1:].split("\t")
h = { hv[i] : i for i in range(0,len(hv)) }


combs = args.comb.split(",")
sys.stdout.write("#" + "\t".join(hv))
for line in bedFile:
    vals = line.split()
    callers = vals[h["CALLER"]]
    nCallers=0
    for c in combs:
        if c in callers:
            nCallers +=1
    
    if nCallers >= args.minn:
        sys.stdout.write(line)

    


