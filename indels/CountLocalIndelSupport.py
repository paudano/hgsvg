#!/usr/bin/env python
import sys
import argparse

ap = argparse.ArgumentParser(description="Count indels consistent with length in sloppy file")
ap.add_argument("--isect", help="Input", default="/dev/stdin")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()


isect = open(args.isect)
hl    = isect.readline()[1:]
hv = hl.split()

h  = { hv[i] : i for i in range(0,len(hv)) }

prevRgn=None
nSup = 0
prevVals = None
sys.stdout.write("\t".join(hv[3:9] + ["locsup"]) + "\n")
for line in isect:

    vals = line.split()
    rgn = "_".join(vals[0:3])
    if rgn != prevRgn:
        if nSup > 0:
            indel = "\t".join(prevVals[3:9]) + "\t" + str(nSup)
            sys.stdout.write(indel + "\n")
        nSup = 0
    #if rgn == prevRgn:
    svLen = int(vals[h["svLen"]])
    qsvLen = int(vals[h["qsvlen"]])
    if svLen > 0.5*qsvLen and svLen < 2.0*qsvLen:
        nSup +=1

    prevRgn = rgn
    prevVals = vals
    

    
