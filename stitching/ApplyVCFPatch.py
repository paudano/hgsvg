#!/usr/bin/env python

import pysam
import sys
import argparse
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq
ap = argparse.ArgumentParser(description="")
ap.add_argument("--genome", help="Input genome.")
ap.add_argument("--vcf", help="Variant vcf.")
ap.add_argument("--minab", help="Minimum alternate allele balance.", type=float, default=0.85)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()


vcfFile = open(args.vcf)
refFile = pysam.FastaFile(args.genome)
outFile = open(args.out,'w')
vcf = {}
for line in vcfFile:
    if line[0] == "#":
        continue
    vals = line.split()
    tags = vals[7].split(";")
    ro=0
    ao=0
    for t in tags:
        if "RO=" in t[0:3]:
            ro = int(t.split(",")[0][3:])
        if "AO=" in t[0:3]:
            ao = int(t.split(",")[0][3:])

    if (ao+ro > 0 and float(ao)/(ao+ro) >= args.minab and "," not in vals[3] and "," not in vals[4]):
        if vals[0] not in vcf:
            vcf[vals[0]] = []
        vcf[vals[0]].append(vals)

chroms = refFile.references



for chrom in chroms:
    chromSeq = refFile.fetch(chrom)
    
    if chrom not in vcf:
        rec = SeqRecord.SeqRecord(Seq.Seq(chromSeq),id=chrom,name="",description="")
        SeqIO.write(rec, outFile, "fasta")
        continue
    sys.stderr.write(chrom + "\n")
    var = vcf[chrom]
    varPos = [0] + [int(v[1])-1 for v in var]
    varRefLen = [0] + [len(v[3]) for v in var]
    refSeqs = [ chromSeq[varPos[i-1]+varRefLen[i-1]: varPos[i] ] + var[i-1][4] for i in range(1,len(varPos))] + [chromSeq[varPos[-1]+varRefLen[-1]:] ]
    newSeq = "".join(refSeqs)
    rec = SeqRecord.SeqRecord(Seq.Seq(newSeq),id=chrom,name="",description="")    
    SeqIO.write(rec, outFile, "fasta")

    
    
    
    
                        
    
    

    
        


