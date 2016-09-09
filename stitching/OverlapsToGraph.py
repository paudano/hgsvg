#!/usr/bin/env python
import sys
import networkx as nx
import argparse

ap = argparse.ArgumentParser(description="Convert overlaps to a graph, and remove transitively encoded edges.")
ap.add_argument("overlaps", help="Original overlaps file.")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()

#chr13:18849420-18929420/0       0       0       chr13:18829420-18909420/1       0       0       12222   12651   12953   13385

class Overlap:
    
    def __init__(self, line):
        v = line.split()
        
        self.a = v[0]
        self.aContig = (int(v[1]), int(v[2]))
        self.b = v[3]
        self.bContig = (int(v[4]), int(v[5]))
        if v[8] == "-1" or v[9] == "None":
            self.aStart = 0
            self.aEnd   = 0
            self.bStart = 0
            self.bEnd   = 0
        else:
            self.aStart = int(v[6])
            self.aEnd   = int(v[7])
            self.bStart = int(v[8])
            self.bEnd   = int(v[9])

    def OverlapLength(self):
        if self.aEnd != self.aStart and self.bEnd != self.bStart:
            ovp = self.aEnd - self.bStart
            if ovp > 0:
                return ovp
        return 0
overlapFile = open(args.overlaps)
overlaps = [ Overlap(line) for line in overlapFile ]
g = nx.DiGraph()
ovpIndex = 0
for ovp in overlaps:
    if ovp.OverlapLength() > 0:
        g.add_edge(ovp.a, ovp.b, index=ovpIndex)
    ovpIndex+=1

nx.write_gml(g, args.out)
    
    
for g in g.nodes():
    
