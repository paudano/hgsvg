#!/usr/bin/env python

import argparse
import re
import sys

ap = argparse.ArgumentParser(description="Transform a variant VCF into a bed file that can be processed by standard pipelines")
ap.add_argument("--vcf", help="VCF file.", required=True)
ap.add_argument("--bionano", help="Convert bionano vcf", default=False, action='store_true')
ap.add_argument("--operation", help="Specifically extract insertions or deletions.",default=None)
ap.add_argument("--indv", help="Write for this individual", default=None)
ap.add_argument("--filter", help="Add filter value.",default=False,action='store_true')
ap.add_argument("--gtfields", help="Add these genotype fields.", default=[], nargs="+")
ap.add_argument("--fields", help="keep these fields", default=[], nargs="+")
ap.add_argument("--groupCommas", help="Group INFO values separated by commas into the same key", action="store_true",default=False)
ap.add_argument("--ignore-seqlen", help="Ignore length of the sequence when setting length of event",
                action="store_true",
                dest="ignoreSeqlen",
                default=False)
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
    infoFields = []
    if line[0] == "#":
        if line[0:6] == "##INFO":
            fieldMatch = fieldVal.match(line)
            if fieldMatch is not None and len(fieldMatch.groups()) > 0:
                infoFields.append(fieldMatch.groups()[0])
        if line[0:6] == "#CHROM" and args.indv is not None:
            vals    = line.split()
            samples = vals[9:]
            sampleIndex= None
            for s in range(0,len(samples)):
                if samples[s] == args.indv:
                    sampleIndex = s
                    break
            if sampleIndex is None:
                print "ERROR. Could not find sample " + args.indv
                sys.exit(1)
            
                
            
            
    else:

        vals = line.split()
        
        # I messed up by using ;'s in SVANN,
        
        infoStr = vals[7]
        info = []
        for infoGroup in infoStr.split(";"):
            if args.groupCommas is False:
                for infoVal in infoGroup.split(","):
                    if ":" in infoVal:
                        infoVal = infoVal.split(":")[1]
                    if "=" in infoVal:
                        tup = Tuple(infoVal)
                        info.append(tup)
            else:
                if "=" in infoGroup:
                    tup = Tuple(infoGroup)
                    info.append(tup)
                
        infokv = {i[0]: i[1] for i in info}

        gtfieldValues=  []
        if len(args.gtfields) != 0:
            gtKeys = vals[8].split(":")
            gtValues = vals[9].split(":")
            for gtKey in args.gtfields:
                val = "."
                for ki in range(0,len(gtKeys)):
                    if gtKey == gtKeys[ki]:
                        val = gtValues[ki]
                        break
                gtfieldValues.append(val)

        if 0 in infokv:
            del(infokv[0])

        if args.bionano is False:
            svLen = 0
            seq ="NONE"
            
            if "SEQ" in infokv:
                seq = infokv["SEQ"]
            if "SEQ" in infokv and args.ignoreSeqlen is False:
                end = int(vals[1]) + len(seq)
                svLen = len(seq)
            elif "END" in infokv:
                end = int(infokv["END"])
                svLen = int(infokv["END"]) - int(vals[1])
            elif "SVLEN" in infokv:
                end = int(vals[1]) + abs(int(infokv["SVLEN"]))
                svLen = infokv["SVLEN"]
                infokv["END"] = str(end)
            svType = "NONE"
            if "SVTYPE" in infokv:
                svType = infokv["SVTYPE"]
            elif "MERGE_TYPE" in infokv:
                svType = infokv["MERGE_TYPE"]
            svOp = GetOp(vals[4])
            if svOp == "deletion":
                start = vals[1]
            else:
                start = str(int(vals[1])-1)
            
            lineVals = [vals[0], start, str(end), svType, str(svLen), seq]
            if len(args.fields) > 0:
                lineVals += [infokv[f] for f in args.fields]
                
        else:
            lineVals = [vals[0], vals[1], infokv["END"], GetOp(vals[4]), str(abs(int(infokv["SVLEN"] ))), vals[5]]
        if args.operation is not None and GetOp(vals[4]) != args.operation:
            continue
        
        if madeHeader == False:
            #
            # header was missing some keys, maybe
            #
            
            if args.bionano is False:
                header = ["chrom", "tStart", "tEnd", "svType", "svLen", "svSeq"]
                if len(args.fields) > 0:
                    header+= args.fields
            else:
                header = ["chrom", "tStart", "tEnd", "svType", "svLen", "svQual"] 

            if args.indv is not None:
                header += ["gt"]

            if args.filter is True:
                header += ["FILTER"]

            if len(args.gtfields) > 0:
                header += args.gtfields
            
            outFile.write("#" + "\t".join(header) + "\n")
            madeHeader = True

        if args.indv is not None:
            lineVals += [vals[9+sampleIndex]]
        if args.filter is True:
            lineVals += [vals[6]]
        if len(args.gtfields) > 0:
            lineVals += gtfieldValues
        outFile.write("\t".join(lineVals) + "\n")

                
                
            
