#!/usr/bin/env python

import argparse
import re
import sys

ap = argparse.ArgumentParser(description="Transform a variant VCF into a bed file that can be processed by standard pipelines")
ap.add_argument("--vcf", help="VCF file.", required=True)
ap.add_argument("--bionano", help="Convert bionano vcf", default=False, action='store_true')
ap.add_argument("--operation", help="Specifically extract insertions or deletions.",default=None)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()


# First extract which additoinal fields to write.
vcfFile = open(args.vcf)
outFile = open(args.out,'w')
madeHeader = False

fieldVal=re.compile("##INFO=<ID=([^,]+),.*")
fields = []

def Tuple(kv):
    if "=" in kv:
        kvp = kv.split("=")
        return (kvp[0], kvp[1])
    else:
        return (0,0)

splitSVAnn= re.compile("(.*)(SVANN=.*)(SVREP.*)")

def GetOp(v):

    if v == "<INS>" or v == "<ins>":
        return "insertion"
    else:
        return "deletion"
    
for line in vcfFile:
    if line[0] == "#":
        if line[0:6] == "##INFO":
            fieldMatch = fieldVal.match(line)
            if fieldMatch is not None and len(fieldMatch.groups()) > 0:
                fields.append(fieldMatch.groups()[0])
    else:

        vals = line.split()
        
        # I messed up by using ;'s in SVANN,
        
        
        g=vals[7].split(";")
        info = [Tuple(v) for v in g]
        infokv = {i[0]: i[1] for i in info}

        if 0 in infokv:
            del(infokv[0])

        if args.bionano is False:
            svLen = 0
            seq ="NONE"
            
            if "SEQ" in infokv:
                end = int(vals[1]) + len(infokv["SEQ"])
                seq = infokv["SEQ"]
                svLen = len(infokv["SEQ"])
            elif "END" in infokv:
                end = int(infokv["END"])
                svLen = int(infokv["END"]) - int(vals[1])
            svType = "NONE"
            if "SVTYPE" in infokv:
                svType = infokv["SVTYPE"]
            elif "MERGE_TYPE" in infokv:
                svType = infokv["MERGE_TYPE"]
            
            lineVals = [vals[0], str(int(vals[1])-1), str(end), svType, str(svLen), seq]
        else:
            lineVals = [vals[0], vals[1], infokv["END"], GetOp(vals[4]), str(abs(int(infokv["SVLEN"] ))), vals[5]]
        if args.operation is not None and GetOp(vals[4]) != args.operation:
            continue
        
        if madeHeader == False:
            #
            # header was missing some keys, maybe
            #
            for key in infokv.keys():
                if key not in fields:
                    print key
                    fields.append(key)
            #
            # We don't want these fields.
            #
            
            for key in ["SEQ", "SVLEN", "SAMPLES", "DP", "CONTIG_DEPTH", "CONTIG_SUPPORT", "OLD_VARIANT"]:
                if key in fields:
                    fields.remove(key)
            i = 0

            while i < len(fields):
                if fields[i] not in infokv:
                    del fields[i]
                else:
                    i+=1
            order = { fields[i] : i for i in range(0,len(fields)) }

            if args.bionano is False:
                header = ["chrom", "tStart", "tEnd", "svType", "svLen", "svSeq"] + fields
            else:
                header = ["chrom", "tStart", "tEnd", "svType", "svLen", "svQual"] + fields
            
            outFile.write("#" + "\t".join(header) + "\n")
            madeHeader = True

                    
        lineVals += [ infokv[k] for k in fields ]            
        outFile.write("\t".join(lineVals) + "\n")

                
                
            
