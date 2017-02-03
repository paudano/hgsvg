#!/usr/bin/env python

import argparse

ap = argparse.ArgumentParser(description="Given a gaps.bed file and a coverage file of breakpoint support, annotate coverage of breakpoints")
ap.add_argument("--gaps", help="gaps.bed file.", required=True)
ap.add_argument("--bpSupport", help="Support of breakpoints", required=True)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()
outFile = open(args.out,'w')
bpsFile = open(args.bpSupport)
bps = {}
for line in bpsFile:
    vals = line.split()
    keyList = vals[3].split("/")
    keyStr  = "/".join(keyList[0:3])
    bpAnnotation = keyList[3]
#    print keyList
#    print keyStr

    if keyStr not in bps:
        bps[keyStr] = []
    bps[keyStr].append((int(vals[4]), bpAnnotation))
#    print bps[keyStr]
    

def GetSupport(bpSupport):
    #
    # Some slightly archane rules for bp support
    #
    if len(bpSupport) == 2:
        return min(bpSupport[0][0], bpSupport[0][1])
    else:
        # If insertion (LR bp), both sides must be aligned.
        if bpSupport[0][1] in "LR":
            return 0
        else:
            # Deletion bp, return support at bp
            return bpSupport[0][0]
    
        
    
gapsFile = open(args.gaps)
for line in gapsFile:
    vals = line.split()
    key = "/".join(vals[7:10])
    if key in bps:
        support= GetSupport(bps[key])
    else:
        support= 0
#    print str(key)+ "\t" +str(  bps[key])
    outFile.write(str(support) + "\n")
