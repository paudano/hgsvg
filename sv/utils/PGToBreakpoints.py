#!/usr/bin/env python


import argparse

ap = argparse.ArgumentParser(description="Print breakpoints near SV events")
ap.add_argument("--gaps", help="Gaps.bed file, transformed to query coordinates.", required=True)
ap.add_argument("--fai", help="Fasta index.", required=True)
ap.add_argument("--window", help="Search this window around breakpoints.", type=int, default=25)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()
gapsFile = open(args.gaps)
outFile = open(args.out,'w')

faiFile = open(args.fai)
fai = {line.split()[0]: int(line.split()[1]) for line in faiFile}

for line  in gapsFile:
    vals = line.split()
    ctg = vals[0]
    if ctg not in fai:
        continue
    if vals[3] == "insertion":
        bp1=[vals[0], str(max(0,int(vals[1])-args.window)), str(min(fai[ctg], int(vals[1])+args.window))]
        bp2=[vals[0], str(max(0,int(vals[2])-args.window)), str(min(fai[ctg], int(vals[2])+args.window))]
        outFile.write("\t".join(bp1) + "\t" + "/".join(vals[0:3] + ["L"]) + "\n")
        outFile.write("\t".join(bp2) + "\t" + "/".join(vals[0:3] + ["R"]) + "\n")
    if vals[3] == "deletion":
        bp = [vals[0], str(max(0,int(vals[1])-args.window)), str(min(fai[ctg], int(vals[1])+args.window))]
        outFile.write("\t".join(bp) + "\t" + "/".join(vals[0:3] + ["D"]) + "\n")        
