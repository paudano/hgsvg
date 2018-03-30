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
    for line in vcf:
        if line[0] == "#":
            continue

        vals = line.split()
        op = vals[4]
        kvs = vals[7].split(";")
        svlen = abs(int(GetValue(kvs,"SVLEN")))
        svtype = GetValue(kvs,"REPEAT_TYPE")
        if svtype is None:
            svtype = GetValue(kvs, "SVCLASS")
            
        if "INS" in op:
            ins.append((svtype,int(svlen)))
        else:
            dels.append((svtype,abs(int(svlen))))

    return ins, dels   
        
def GetReport(sv_types, ins, dels, minLength=0, maxLength=10000000):
    dirs = ["insertion", "deletion"]
    stats = []
    insLen = []
    for i in ins:
        if i[0] in sv_types:
            insLen.append(i[1])
    delLen = []
    for d in dels:
        if d[0] in sv_types:
            delLen.append(d[1])
    stats.append([len(insLen), np.mean(insLen), np.std(insLen), np.sum(insLen)])
    stats.append([len(delLen), np.mean(delLen), np.std(delLen), np.sum(delLen)])
    return stats

def PrintReport(name, s, out, sep="\t"):
    out.write(sep.join(s for s in [name, str(s[0][0]),"{:2.2f}".format(s[0][1]),"{:2.2f}".format(s[0][2]),str(s[0][3]),"{}".format(s[1][0]), "{:2.2f}".format(s[1][1]),"{:2.2f}".format(s[1][2]),"{}".format(s[1][3])]) + "\n")


datasets= ["AluY", "L1HS", "Alu.Mosaic", "Alu.STR","SVA","L1P","AluS","HERV", "MER","HSAT", "LTR", "Beta", "STR", "Complex", "ALR", "Singletons", "L1", "TRF", "NONE" , "PartialTandem", "BSR/Beta", "TandemRepeat"]

ins,dels = ParseVCF(args.vcf)
out.write("\tInsertion\t\t\t\tDeletion\t\t\t\t\n")
out.write("Repeat type\tCount\tMean\ts.d.\tTotal\tCount\tMean\ts.d.\tTotal\n")
s = GetReport(["NONE" ], ins, dels)
PrintReport("Complex", s, out,sep=args.sep)

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

s = GetReport(["MER","HSAT", "LTR", "Beta", "ALR", "BSR/Beta"], ins, dels)
PrintReport("Satellite", s, out,sep=args.sep)

s=GetReport(["Singletons"], ins, dels)
PrintReport("Singletons", s, out,sep=args.sep)

s=GetReport(["TandemRepeat", "PartialTandem"], ins, dels)
PrintReport("Tandem_Repeats", s,out,sep=args.sep)

s = GetReport(datasets, ins, dels)
PrintReport("Total", s, out,sep=args.sep)

