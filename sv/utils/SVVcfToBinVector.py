#!/usr/bin/env python


import argparse

ap = argparse.ArgumentParser(description="")
ap.add_argument("--vcf", help="Read variants as vcf", default=None)
ap.add_argument("--bed", help="Read variants as bed", default=None)
ap.add_argument("--bin", help="Bin size", default=10000,type=int)
ap.add_argument("--op", help="Override operation", default=None)
ap.add_argument("--radiff", help="Compute sv size as ref-alt", default=False, action='store_true')
ap.add_argument("--genome", help="Genome size table to use (e.g. fai file)", default=None, required=True)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()


faiFile = open(args.genome)
outFile = open(args.out,'w')

fai = { l.split()[0]: int(l.split()[1]) for l in faiFile}
bins = { c: [0]*(fai[c]/args.bin+1) for c in fai }

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
        if args.radiff is False:
            for i in info:
                kv = i.split("=")
                
                if len(kv) == 2:
                    if kv[0] == "SVLEN":
                        svlen = int(kv[1])
                    elif kv[0] == "SVTYPE":
                        if "INS" in kv[1]:
                            op = "INS"
                        elif "DEL" in kv[1]:
                            op = "DEL"
            if op is not None:
                if op == "INS":
                    bins[chrom][pos/args.bin] += abs(svlen)
                else:
                    bins[chrom][pos/args.bin] -= svlen

        else:
            svlen = len(vals[4]) - len(vals[3])
            bins[chrom][pos/args.bin] += svlen
                        

if args.bed is not None:
    bedFile = open(args.bed)
    for line in bedFile:
        vals = line.split()
        if vals[0] not in bins:
            continue
	if args.op is not None:
	    bins[vals[0]][int(vals[1])/args.bin] += int(vals[2]) - int(vals[1])
	else:
	    if vals[3] == "insertion":
                bins[vals[0]][int(vals[1])/args.bin] += int(vals[2]) - int(vals[1])
            if vals[3] == "deletion":
                bins[vals[0]][int(vals[1])/args.bin] -= int(vals[2]) - int(vals[1])

for c in bins:
    cstr = "\n".join(["{}\t{}\t{}\t{}".format(c, i*args.bin, (i+1)*args.bin, bins[c][i]) for i in range(0,len(bins[c]))])
    outFile.write(cstr + "\n")
    
    


