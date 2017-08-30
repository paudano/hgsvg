#!/usr/bin/env python
import sys
import numpy
import argparse

ap = argparse.ArgumentParser(description="Add variance to table of bin net gain")
ap.add_argument("netgain", help="netgain file")
ap.add_argument("--max-zero", help="maximum number of zeros to add", dest="maxZero", type=int, default=0)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()

netGainFile = open(args.netgain)
outFile = open(args.out,'w')
for line in netGainFile:
    vals = line.split()
    counts = [int(v) for v in vals[3:]]
    nZero = sum([v == 0 for v in counts])
    if nZero > args.maxZero:
        sd = 0
    else:
        sd = numpy.sqrt(numpy.var(counts))
    outFile.write("\t".join(vals) + "\t" + "{:2.2f}".format(sd) + "\n")
outFile.close()
        
    

