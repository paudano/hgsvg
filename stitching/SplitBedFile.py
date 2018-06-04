#!/usr/bin/env python


import argparse
import math
ap = argparse.ArgumentParser(description="Split a bed file with some overlap between splits")
ap.add_argument("--bed", help="Input bed file.",required=True)
ap.add_argument("--n", help="Numbe of lines per file.",required=True,type=int)
ap.add_argument("--overlap", help="Overlap between files",required=True,type=int)
ap.add_argument("--base", help="Base name for output file",required=True)

args = ap.parse_args()

bedFile = open(args.bed)
bed = bedFile.readlines()
linesPerBin=math.ceil(float(len(bed))/args.n)

for i in range(0, args.n):
    start=int(min(i*linesPerBin,len(bed)))
    end=int(min(len(bed), (i+1)*linesPerBin))
    outFile = open(args.base + "." + str(i) + ".bed",'w')
    if start < end:
        outFile.write(''.join(bed[start:end]))
    outFile.close()
    
