#!/usr/bin/env python

import sys

bedFile = open(sys.argv[1])
bed6File = open(sys.argv[2],'w')

headerVals = bedFile.readline()[1:].split()

header = { headerVals[i] : i for i in range(0,len(headerVals)) }
for line in bedFile:
    vals = line.split()
    bed6File.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(vals[header["chrom"]],
                                                     vals[header["tStart"]],
                                                     vals[header["tEnd"]],
                                                     vals[header["svType"]],
                                                     1000,
                                                     "+"))

                                                     
