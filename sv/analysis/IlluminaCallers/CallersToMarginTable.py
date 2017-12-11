#!/usr/bin/env python

import argparse

ap = argparse.ArgumentParser(description="Convert a list of SVs and pass/fail, output table of callers")
ap.add_argument("--svbed", help="SV bed file.", default="/dev/stdin")
ap.add_argument("--callerKey", help="Key column for callers", default="CALLER")
ap.add_argument("--callers", help="Specify callers, otherwise determine from input.", nargs="+", default=[])
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()


svBedFile = open(args.svbed)
outFile = open(args.out, 'w')

svBed = []
for line in svBedFile:
    if line[0] == "#":
        headerVals = line[1:].split()
        header = {headerVals[i] : i for i in range(0,len(headerVals)) }
    else:
        svBed.append(line.split())
callers = []
if len(args.callers) == 0:
    callerDict={}
    for i in range(0,len(svBed)):
        for c in svBed[i][header[args.callerKey]].split(","):
            callerDict[c] = True
    if "NA" in callerDict:
        del callerDict["NA"]
    callers = callerDict.keys()
else:
    callers = args.callers

outFile.write("#chrom\ttStart\ttEnd\t" + "\t".join(callers))
if "orth_filter" in header:
    outFile.write("\torth_filter")
outFile.write("\n")

for i in range(0,len(svBed)):
    cdb = { c: "0" for c in callers }
    svCallers = svBed[i][header[args.callerKey]].split(",")
    for svCaller in svCallers:
        cdb[svCaller] = "1"
    cString = "\t".join(cdb[c] for c in callers)
    outFile.write("\t".join(svBed[i][0:3]) + "\t" + cString)
    if "orth_filter" in header:
        outFile.write("\t" + svBed[i][header["orth_filter"]])
    outFile.write("\n")
    


     
            
        

