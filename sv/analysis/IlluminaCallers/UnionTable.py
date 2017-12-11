#!/usr/bin/env python
import sys


import argparse

ap = argparse.ArgumentParser(description="Make a table  specific to the method union")
ap.add_argument("--noHeader", help="no header.", action='store_true', default=False)
args = ap.parse_args()

combs = ["PacBio", "Illumina", "BioNano", "PacBio,Illumina", "PacBio,BioNano","Illumina,BioNano", "All"]
counts = { c: 0 for c in combs }

sizes = { c: 0 for c in combs }

headerLine = sys.stdin.readline()
hv= headerLine[1:].split()
h= { hv[i] : i for i in range(0,len(hv)) }
if "sumsvLen" in h:
    svLenIdx = "sumsvLen"
else:
    svLenIdx = "svLen"
for line in sys.stdin:
    vals = line.split()
    u = vals[h["union"]]
    counts[u] += 1
    if svLenIdx in h:
        sizes[u] += int(vals[h[svLenIdx]])

for c in combs:
    outstr = ""
    if args.noHeader is False:
        outstr += c + "\t"
    outstr += str(counts[c])
    if "svLen" in h:
        avgLen = 0
        if counts[c] > 0:
            avgLen = float(sizes[c])/counts[c]
        outstr += "\t{:2.2f}".format(avgLen)
    outstr += "\n"
    sys.stdout.write(outstr)



           
