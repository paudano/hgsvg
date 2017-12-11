#!/usr/bin/env python
import sys


import argparse

ap = argparse.ArgumentParser(description="")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--pli", help="Pli.")
ap.add_argument("--svBed", help="SV bed file")

args = ap.parse_args()
outFile = open(args.out,'w')
pliFile = open(args.pli)
pliHeader = pliFile.readline().split()
pliIdx = { pliHeader[i] : i for i in range(0,len(pliHeader)) }
plih = {}
for line in pliFile:
    vals = line.split()
    
    g = vals[pliIdx["gene"]]
    p = vals[pliIdx["pLI"]]
    plih[g]= float(p)

svBedFile = open(args.svBed)
svBedHeaderV = svBedFile.readline().split()
svh = { svBedHeaderV[i] :i for i in range(0,len(svBedHeaderV)) }

outFile.write("PLI\n")
for line in svBedFile:
    vals = line.split()
    geneName = vals[svh["geneName"]]
    pli = "NA-PLI-"+geneName
    if geneName in plih:
        pli = "{:2.3}".format(float(plih[geneName]))
    outFile.write(pli + "\n")





