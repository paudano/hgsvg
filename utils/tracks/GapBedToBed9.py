#!/usr/bin/env python
import sys
inFile = open(sys.argv[1])
outFile= open(sys.argv[2],'w')

inFile.readline()
for line in inFile:
    vals = line.split()
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

