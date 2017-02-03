#!/usr/bin/env python
import sys
outFile = open(sys.argv[1],'w')
dirName = sys.argv[2]
events = [ ["AluY", "AluY.simple.bed"],
           ["AluS", "AluS.simple.bed"],
           ["STR",  "STR.bed"],
           ["L1HS",    "L1HS.simple.bed"],
           ["Alu.Mosaic",     "Alu.Mosaic.bed"],
           ["Alu.STR",     "Alu.STR.bed"],
           ["ALR",     "ALR.bed"],
           ["SVA",     "SVA.simple.bed"],
           ["HERV",    "HERV.simple.bed"],
           ["L1P",    "L1P.bed"],
           ["BSR/Beta",   "Beta.bed"],
           ["HSAT",   "HSAT.bed"],
           ["MER",    "MER.bed"],
           ["L1",     "L1.bed"],
           ["LTR",    "LTR.bed"],
           ["Singletons", "Singletons.bed"],
           ["TandemRepeat", "TRF.bed"],
           ["NONE", "NONE.bed"],
           ["Complex", "Complex.bed"],
           ["PartialTandem", "partial_annotation.bed"]]

first=True
for e in events:
    with open(dirName + "/" + e[1]) as f:
        for l in f:
            if l[0] == "#":
                if first:
                    outFile.write(l.strip() + "\tsvClass\n")
                    first = False
                continue
            outFile.write(l.rstrip() + "\t" + e[0] + "\n")

    
    

