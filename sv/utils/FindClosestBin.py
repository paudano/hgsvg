#!/usr/bin/env python

import sys


import argparse

ap = argparse.ArgumentParser(description="")
ap.add_argument("bins", help="Binsfile that has been intersected with an insertion or deletion file with -loj, and merged")
ap.add_argument("--minBin", help="Minimum bin to consider supported", type=int,default=5)
ap.add_argument("--maxDist",help="Maximum distance to search", type=int,default=15)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()

outFile = open(args.out,'w')

binsFile = open(args.bins)



distance={}
prevChrom = None
bins = []
for line in binsFile:
    v = line.split()
    if prevChrom is not None and v[0] != prevChrom:

        nIsect = 0
        for i in range(0,len(bins)):
            if bins[i][4] != ".":
                pv = 0
                p = i
                while p > 0 and i-p < args.maxDist and int(bins[p][3]) < args.minBin:
                    p-=1
                if p > 0 and i-p < args.maxDist:
                    pv=bins[p][3]
                    
                n = i
                nv = 0
                while n < len(bins) and n-i < args.maxDist and int(bins[n][3]) < args.minBin:
                    n+=1
                if  n < len(bins) and n-i < args.maxDist:
                    nv = bins[n][3]
                dist = min(i-p, n-i)
                sup = pv
                if dist==n-i:
                    sup=nv
                
                distance[bins[i][4]+"\t"+bins[i][5]+"\t"+ bins[i][6]] = (dist, sup)
                nIsect+=1
                
        bins = []
    else:
        bins.append(v)
    prevChrom = v[0]
        

        


for v in distance:
    outFile.write(v+ "\t" + str(distance[v][0]) + "\t" + str(distance[v][1]) + "\n")
