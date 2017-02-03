#!/usr/bin/env python


import argparse

ap = argparse.ArgumentParser(description="")
ap.add_argument("--vcf", help="Read variants as vcf", default=None)
ap.add_argument("--bed", help="Read variants as bed", default=None)
ap.add_argument("--bin", help="Bin size", default=10000,type=int)
ap.add_argument("--op", help="Override operation", default=None)
ap.add_argument("--radiff", help="Compute sv size as ref-alt", default=False, action='store_true')

ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()



outFile = open(args.out,'w')



if args.vcf is not None:
    vcfFile = open(args.vcf)
    for line in vcfFile:
        if line[0] == '#':
            continue
        vals = line.split()
        pos = int(vals[1])
        chrom = vals[0]
        info = vals[7].split(";")
        svlen = 0
        op  = None
        seq = "T" # fake
        if args.radiff is False:
            for i in info:
                kv = i.split("=")
                
                if len(kv) == 2:
                    if kv[0] == "SVLEN":
                        svlen = int(float(kv[1]))
                    elif kv[0] == "SVTYPE":
                        if "INS" in kv[1]:
                            op = "INS"
                        elif "DEL" in kv[1]:
                            op = "DEL"
                    elif kv[0] == "SEQ":
                        seq = kv[1]
                            
            if op is not None:
                if op == "INS":
                    outFile.write(chrom + "\t" + str(pos) + "\t" + str(pos + abs(svlen)) + "\tinsertion\t" + str(len(seq)) + "\t" + seq + "\n")

                else:
                    outFile.write(chrom + "\t" + str(pos) + "\t" + str(pos + abs(svlen)) + "\tdeletion\t" + str(len(seq)) + "\t" + seq + "\n")


        else:
            svlen = len(vals[4]) - len(vals[3])
            bins[chrom][pos/args.bin] += svlen
                        
