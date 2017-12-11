#!/usr/bin/env python

import argparse
import sys

ap = argparse.ArgumentParser(description="Add int filt key based on merged nonredundant file and a key file")
ap.add_argument("--mnr", help="Merged nonredundant file.")
ap.add_argument("--values", help="Value file, the PBINT key will be used to look up value.")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()

mrFile = open(args.mnr)
outFile = open(args.out, 'w')
mrhv = mrFile.readline()[1:].split()
mrh = {mrhv[i]: i for i in range(0,len(mrhv)) }


valueFile = open(args.values)
vfv = valueFile.readline().split()
vfh = { vfv[i] : i for i in range(0,len(vfv)) }

values = {}
for line in valueFile:
    vals = line.split()
    # key is the pacbio, value is the original Illumina call, so there
    # is a 1-1 correspondence between the two.  If there are multiple
    # Illumina calls that match t
    values[vals[vfh["key"]]] = vals[vfh["value"]]


outFile.write("intkey\n")
nKey = 0
default=0
for line in mrFile:
    vals = line.split()
    key = "_".join(vals[0:3])
    if key in values:
        value = values[key]
        nKey+=1
    else:
        datasets = vals[mrh["union"]]
        
        if "Illumina" in datasets:
            value = "_".join(vals[0:3])
            default+=1
        else:
            value = "NO_MATCH"
    outFile.write(value + "\n")
sys.stderr.write("Matched " + str(nKey) + " " + str(default) + "\n")
        
            
