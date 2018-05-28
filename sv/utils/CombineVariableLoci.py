#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord

ap = argparse.ArgumentParser(description="Take table of cluster SV loci and their alt sequences, and produce a variant bed file.")

ap.add_argument("--regions", help="Bed file of regions.")
ap.add_argument("--lifted", help="Lifted regions", nargs=2)
ap.add_argument("--sequences", help="Lifted sequences", nargs=2)
ap.add_argument("--outFile", help="Out file", default="/dev/stdout")


args=ap.parse_args()
outFile = open(args.outFile,'w')

regionsFile = open(args.regions)
def ParseBed3(line):
    vals=line.rstrip().split()
    return [vals[0], int(vals[1]), int(vals[2])]


def ParseLifted(line):
    vals=line.rstrip().split()
    if vals[0] != "NA":
        vals[1] = int(vals[1])
        vals[2] = int(vals[2])
    return vals

regions = [ ParseBed3(line) for line in regionsFile ]

lb = [ open(f) for f in args.lifted ]
liftedTab = [ [ParseLifted(line) for line in lb[i] ] for i in range(0,2)]
liftedFasta = [ open(f) for f in args.sequences ]

idx=[0,1]
def GetSVLen(tab):
    if (tab[0] == "NA"):
        return None
    else:
        return tab[2]-tab[1]
parsers=[    SeqIO.parse(liftedFasta[j],"fasta") for j in idx ]


def GetDefinedLifted(tab, i, j):
    res = []
#    pdb.set_trace()
    for r in idx:
        if tab[r][i][j] != "NA":
            res.append(str(tab[r][i][j]))
    return ",".join(res)


header=["#chrom", "tStart","tEnd", "hap", "svType", "svLen", "svSeq", "qName", "qStart","qEnd", "region", "nAlt", "nRef"]
outFile.write("\t".join(header) + "\n")

for i in range(0,len(regions)):
    seqs = [next(p) for p in parsers]
    hap=None
    if liftedTab[0][i][0] == "NA" and liftedTab[1][i][0] == "NA":
        continue
    hap="HOM"
    if liftedTab[0][i][0] == "NA":
        hap="HAP1"
    elif liftedTab[1][i][0] == "NA":
        hap="HAP0"
    if hap is None:
        continue
    refLen = regions[i][2] - regions[i][1]
    
    svLens = [GetSVLen(liftedTab[j][i]) for j in idx]
    deltas = []
    for svLen in svLens:
        if svLen is not None:
            deltas.append(str(svLen - refLen))
            
#    import pdb
#    pdb.set_trace()
    region   = "{}:{}-{}".format(regions[i][0],regions[i][1], regions[i][2])
    svRecord = [str(s) for s in regions[i]] + [hap, "locus", ",".join(deltas), "|".join([str(s.seq) for s in seqs])] + [GetDefinedLifted(liftedTab, i, j) for j in [0,1,2]] + [region, "0","0"]
    outFile.write("\t".join(svRecord) + "\n")    
    





    
