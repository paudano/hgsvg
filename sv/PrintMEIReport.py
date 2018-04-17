#!/usr/bin/env python

import sys
import numpy as np
import os
import argparse
ap = argparse.ArgumentParser(description="Print a mei summary.  This should be ran from the directory where AnnotationPipeline.mak was ran.")
ap.add_argument("--vcf", help="VCF with annotations.")
ap.add_argument("--out", help="out", default="/dev/stdout")
ap.add_argument("--sep", help="Separating character.", default="\t")

args = ap.parse_args()
out = open(args.out,'w')

col='f4'

def GetValue(kvs, k):
    for kv in kvs:
        kvp = kv.split("=")
        if kvp[0] == k:
            return kvp[1]


def ParseVCF(vcfName):
    vcf = open(vcfName)
    dels = []
    ins = []
    loc = []
    for line in vcf:
        if line[0] == "#":
            continue

        vals = line.split()
        op = vals[4]
        kvs = vals[7].split(";")
        if "INS" in op or "DEL" in op:
            svlen = abs(int(GetValue(kvs,"SVLEN")))
        
            svtype = GetValue(kvs,"REPEAT_TYPE")
            if svtype is None:
                svtype = GetValue(kvs, "SVCLASS")
            
            if "INS" in op:
                ins.append([svtype,int(svlen),0])
            else:
                dels.append([svtype,abs(int(svlen)),0])
        elif "LOC" in op:
            svlen=GetValue(kvs,"SVLEN")
            vals=[int(i) for i in svlen.split(",")]
            avgSVLen = sum(vals)/2
            if (avgSVLen > 0):
                ins.append(["locus", avgSVLen, 0])
            else:
                dels.append(["locus", abs(avgSVLen),0])
    sys.stderr.write(str(len(ins)) + "\t" + str(len(dels)) + "\t" + str(len(ins)) + "\n")
    return ins, dels


def GetReport(sv_types, ins, dels, minLength=0, maxLength=10000000):

    stats = []
    insLen = []
    for idx in range(0,len(ins)):
        i=ins[idx]
        if i[0] in sv_types or sv_types[0] == "any":
            insLen.append(i[1])
            ins[idx][2]+=1
    delLen = []
    for idx in range(0,len(dels)):
        d=dels[idx]
        if d[0] in sv_types or sv_types[0] == "any":
            delLen.append(d[1])
            dels[idx][2]+=1
    stats.append([len(insLen), abs(np.mean(insLen)), np.std(insLen), np.sum(insLen)])
    stats.append([len(delLen), abs(np.mean(delLen)), np.std(delLen), np.sum(delLen)])
    return stats

def PrintReport(name, s, out, sep="\t"):
    out.write(sep.join(s for s in [name, str(s[0][0]),"{:2.2f}".format(s[0][1]),"{:2.2f}".format(s[0][2]),str(s[0][3]),"{}".format(s[1][0]), "{:2.2f}".format(s[1][1]),"{:2.2f}".format(s[1][2]),"{}".format(s[1][3])]) + "\n")


datasets= ["AluY", "L1HS", "Alu.Mosaic", "Alu.STR","SVA","L1P","AluS","HERV", "MER","HSAT", "LTR", "Beta", "STR", "Complex", "ALR", "Singletons", "L1", "TRF", "NONE" , "PartialTandem", "BSR/Beta", "TandemRepeat", "locus"]

ins,dels = ParseVCF(args.vcf)
out.write("\tInsertion\t\t\t\tDeletion\t\t\t\t\n")
out.write("Repeat type\tCount\tMean\ts.d.\tTotal\tCount\tMean\ts.d.\tTotal\n")
s = GetReport(["Complex" ], ins, dels)
PrintReport("Not-repeat", s, out,sep=args.sep)

s = GetReport(["AluY"], ins, dels)
PrintReport("AluY", s, out,sep=args.sep)

s = GetReport(["AluS"], ins, dels)
PrintReport("AluS", s, out,sep=args.sep)

s = GetReport(["L1HS"], ins, dels)
PrintReport("L1Hs", s, out,sep=args.sep)

s = GetReport(["L1"], ins, dels)
PrintReport("L1", s, out,sep=args.sep)

s = GetReport(["L1P"], ins, dels)
PrintReport("L1P", s, out,sep=args.sep)

s = GetReport(["SVA"], ins, dels)
PrintReport("SVA", s, out,sep=args.sep)

s = GetReport(["HERV"], ins, dels)
PrintReport("HERV", s, out,sep=args.sep)

s = GetReport(["Alu.Mosaic", "Alu.STR"], ins, dels)
PrintReport("Alu-mosaic", s, out,sep=args.sep)

s = GetReport(["STR"], ins, dels)
PrintReport("STR", s, out,sep=args.sep)

s = GetReport(["MER","HSAT", "LTR", "Beta", "ALR", "BSR/Beta", "BSR"], ins, dels)
PrintReport("Satellite", s, out,sep=args.sep)

s=GetReport(["Singletons"], ins, dels)
PrintReport("Singletons", s, out,sep=args.sep)

s=GetReport(["TandemRepeat", "PartialTandem"], ins, dels)
PrintReport("Tandem_Repeats", s,out,sep=args.sep)

s=GetReport(["locus"], ins, dels)
PrintReport("loci", s,out,sep=args.sep)

s = GetReport(["any"], ins, dels)
PrintReport("Total", s, out,sep=args.sep)

sys.stderr.write("dups\n")
for i in range(0,len(dels)):
    if dels[i][2] > 2:
        sys.stderr.write(str(dels[i]) + "\n")

for i in range(0,len(ins)):
    if ins[i][2] > 2:
        sys.stderr.write(str(ins[i]) + "\n")
        
