#!/usr/bin/env python


import argparse

ap = argparse.ArgumentParser(description="Split a bed file with some overlap between splits")
ap.add_argument("--bed", help="Input bed file.",required=True)
ap.add_argument("--n", help="Numbe of lines per file.",required=True,type=int)
ap.add_argument("--overlap", help="Overlap between files",required=True,type=int)
ap.add_argument("--base", help="Base name for output file",required=True)

args = ap.parse_args()

bedFile = open(args.bed)
bed = bedFile.readlines()
if len(bed) == 0:
    outFile= open(args.base + ".0.bed",'w')
    outFile.close()
else:
    for i in range(0, len(bed), args.n):
        start = max(0,i-args.overlap)
        end = min(i+args.n,len(bed))
        outFile = open(args.base + "." + str(i) + ".bed",'w')
        outFile.write(''.join(bed[start:end]))
        outFile.close()
    
