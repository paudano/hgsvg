#!/usr/bin/env python

import argparse

ap = argparse.ArgumentParser(description="Output from gap bed to bigBed compatible format")
ap.add_argument("inFile", help="input gaps .")
ap.add_argument("--outFile", help="Output file", default="/dev/stdout")
ap.add_argument("--fai", help="Use fai to trim intervals before endds of contigs.",default=None)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()
#!/usr/bin/env python
import sys
fai=None
if args.fai is not None:
    faiFile=open(args.fai)
    fai = {l.split()[0]: int(l.split()[1]) for l in faiFile}

    
inFile = open(args.inFile)
outFile = open(args.outFile,'w')
inFile.readline()
for line in inFile:
    vals = line.split()
    if fai is not None:
        if vals[0] in fai:
            vals[2] = str(min(fai[vals[0]], int(vals[2])))
    outLine=[""]*9
    outLine[0:3] = vals[0:3]
    outLine[3] = str(vals[4])
    outLine[4] = "1000"
    outLine[5] = "+"
    outLine[6] = vals[1]
    outLine[7] = vals[2]
    if vals[3] == "deletion":
        outLine[8] = "255,0,0"
    else:
        outLine[8] = "0,0,255"
    outFile.write("\t".join(outLine)+  "\n")

