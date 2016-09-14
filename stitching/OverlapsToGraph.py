#!/usr/bin/env python
import sys
import networkx as nx
import argparse

import Overlap

ap = argparse.ArgumentParser(description="Convert overlaps to a graph, and remove transitively encoded edges.")
ap.add_argument("overlaps", help="Original overlaps file.")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()

#chr13:18849420-18929420/0       0       0       chr13:18829420-18909420/1       0       0       12222   12651   12953   13385
overlapFile = open(args.overlaps)
overlaps = [ Overlap.Overlap(line) for line in overlapFile ]
g = nx.DiGraph()
ovpIndex = 0
for ovp in overlaps:
    if ovp.OverlapLength() > 0 and ovp.Contained() == False and ovp.Extends(wiggle=1000) == True:
        g.add_edge(ovp.a, ovp.b, index=ovpIndex)
    ovpIndex+=1

nx.write_gml(g, args.out)
