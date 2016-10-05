#!/usr/bin/env python
import sys
inBed = open(sys.argv[1])
outBed = open(sys.argv[2],'w')
for line in inBed:
    v = line.split()
#   chr10   41705044        41776677        stitch.007.3.146459     0       3203    75940   254     0.989083        70851   613     1273    16
    strand = "+"
    if v[4]  == "1":
        strand = "-"
    outBed.write("\t".join(v[0:4]+ ["1000", strand]) + "\n")


