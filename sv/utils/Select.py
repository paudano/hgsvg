#!/usr/bin/env python
import argparse
import datetime
import numpy as np
import pandas as pd
import pysam
import sys



ap = argparse.ArgumentParser(description="Select columns from a table")
ap.add_argument("--table", help="input table", default="/dev/stdin")
ap.add_argument("--cols", help="Column names", required=True, nargs="+")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()
inFile = open(args.table)
outFile = open(args.out,'w')
tab = pd.read_table(inFile, keep_default_na=False)#, header=None, usecols=columns, names=names)
for col in args.cols:
    if col not in tab:
        sys.stderr.write("ERROR! Could not find " + col + " in " + args.table + "\n")
        sys.exit(1)
sub = tab[args.cols]
headerList = args.cols[:]
headerList[0] = "#"+ headerList[0]
sub.columns = headerList
sub.to_csv(outFile, sep="\t", index=False,header=True)
