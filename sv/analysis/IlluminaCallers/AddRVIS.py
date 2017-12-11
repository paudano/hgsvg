#!/usr/bin/env python


import argparse

ap = argparse.ArgumentParser(description="")
ap.add_argument("--svBed", help="Bed  file of svs")
ap.add_argument("--rvis", help="RVIS file.")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()

outFile = open(args.out,'w')
outFile.write("rvis\n")
rvisFile = open(args.rvis)
rhv = rvisFile.readline().split()
rh = {rhv[i] : i for i in range(0,len(rhv)) }

allRvis = {}
for line in rvisFile:
    vals = line.split()
    gene = vals[rh["GENE"]]
    rvis = vals[rh["%ALL_0.1%"]]
    allRvis[gene] = rvis

svBedFile = open(args.svBed)
svv = svBedFile.readline().split()
svh = {svv[i] : i for i in range(0,len(svv))}

for line in svBedFile:
    vals = line.split()
    geneName = vals[svh["geneName"]]
    rvis = "NA-"+geneName
    if geneName in allRvis:
        rvis = allRvis[geneName]
    outFile.write(rvis + "\n")

