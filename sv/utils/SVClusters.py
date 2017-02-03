#!/usr/bin/env python
import sys
import argparse
ap = argparse.ArgumentParser(description="Write clusters of SVs that are within N kbp")
ap.add_argument("--svs", help="SVs,with header line.", required=True)
ap.add_argument("--window", help="window to search for svs", type=int,required=True)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")

args = ap.parse_args()

svFile = open(args.sv)
headerVals = svFile.readline()[1:].split()
header = {headervals[i] : i for i in range(0,len(headerVals)) }

def GetSVBounds(vals):
    return [vals[header["chrom"]], int(vals[header["tStart"]]), int(vals[header["tEnd"]])]

sv = [GetSVBounds(line.split()) for line in svFile ]
for
