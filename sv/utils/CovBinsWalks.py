#!/usr/bin/env python

import argparse

ap = argparse.ArgumentParser(description="")
ap.add_argument("bins", help="Bins file, in bed format with 3 additional columns: coverage ins-count del-count")
ap.add_argument("--op", help="Operation (ins or del)", choices=["ins", "del"], required=True)
ap.add_argument("--out", help="Output file ", required=True)

args = ap.parse_args()
runOut = open(args.out,'w')
if (args.op == "ins"):
    col = 4
else:
    col = 5



bins = open(args.bins)
prev = None
run  = 0
prevLines = []
for line in bins:
    vals = line.rstrip().split()
    count = int(vals[col])
    if count == 0:
        if len(prevLines) > 0:
            prevLines[0][2] = prevLines[-1][2]
            prevLines[0][3] = str(run)
            runOut.write("\t".join(prevLines[0]) + "\n")
        runOut.write("\t".join(vals[0:3]+[vals[col]]) + "\n")
        prevLines = []
        run = 0
    else:
        run += count
        prevLines.append(vals[0:3]+[vals[col]])

if len(prevLines) > 0:
    runOut.write("\n".join("\t".join(p) for p in prevLines) + "\n")

        
        
    
