#!/usr/bin/env python

import sys
import argparse

ap = argparse.ArgumentParser(description="Given a gaps.bed file, splice variants into the reference, remap reads to this, and count gaps. ")
ap.add_argument("--gaps", help="Gaps file (can be insertions and deletions).", required=True)
ap.add_argument("--window", help="Collect all gaps within this window.", type=int, default=1000)
ap.add_argument("--mssm", help="Parse mssm lines.", default=False, action='store_true')
ap.add_argument("--falcon", help="Parse in falcon format", default=False, action='store_true')
ap.add_argument("--out", help="Clusters file.", default='/dev/stdout')
ap.add_argument("--minCalls", help="Write out calls if they have at least this many events.", type=int, default=0)

args = ap.parse_args()
clusterOut = open(args.out,'w')

def ParseFalconLine(vals):
    return [vals[0], int(vals[1]), int(vals[2]), "deletion"]

def ParseMSSMLine(vals):

    if len(vals) > 0:
        return [vals[0], int(vals[1]), int(vals[2]), vals[3], vals[7], vals[8]]
    else:
        return None

def ParseGapsLine(vals):
    if len(vals) > 0:
        return [vals[header["chrom"]], int(vals[header["tStart"]]), int(vals[header["tEnd"]]), vals[header["svType"]], int(vals[header["svLen"]])] + vals[5:]
    else:
        return None


def GetEnd(v):
    if args.mssm:
        return v[1] + len(v[5])
    elif args.falcon:
        return v[2]
    else:
        if v[4] == "deletion":
            return int(v[2])
        else:
            return int(v[1])+1

def GetLength(v):
    return int(v[4])

def GetStart(v):
    return int(v[1])

def GetChrom(v):
    return v[0]

def GetType(v):
    return v[3]

gapsFile = open(args.gaps)
headerLine = gapsFile.readline()[1:].split()
header = {headerLine[i]: i for i in range(0,len(headerLine))}

gapLines = gapsFile.readlines()

if args.mssm:
    gaps = [ParseMSSMLine(line.split()) for line in gapLines]
elif args.falcon:
    gaps = [ParseFalconLine(line.split()) for line in gapLines]
else:
    gaps = [ParseGapsLine(line.split()) for line in gapLines]

svClusters = []
i=0
while i < len(gaps):
    svClusters.append([])
    svClusters[-1].append(gaps[i])

    while i < len(gaps)-1 and GetChrom(gaps[i]) == GetChrom(gaps[i+1]) and GetStart(gaps[i+1])-GetEnd(gaps[i]) < args.window:
        i+=1
        svClusters[-1].append(gaps[i])
    i+=1


def GetRegion(cluster):
    return "{}:{}-{}/{}".format(GetChrom(cluster), GetStart(cluster), GetStart(cluster) + GetLength(cluster),GetType(cluster))

for clust in svClusters:
    if len(clust) >= args.minCalls:
        start = min([int(c[1]) for c in clust])
        end = max([int(c[2]) for c in clust])
        clusterOut.write(GetChrom(clust[0]) + "\t{}\t{}\t{}\t".format(start,end,len(clust)) + ";".join([GetRegion(c) for c in clust]) + "\n")

clusterOut.close()

