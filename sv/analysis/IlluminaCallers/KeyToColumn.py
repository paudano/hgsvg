#!/usr/bin/env python

import sys

import argparse

ap = argparse.ArgumentParser(description="")
ap.add_argument("--bed", help="Original bed file.")
ap.add_argument("--keys", help="Keys, one line per key.")
ap.add_argument("--header", help="Use this header", default=None)
ap.add_argument("--reverse", help="Write reverse keys in order here", default=None)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()

keyFile = open(args.keys)
bedFile = open(args.bed)
outFile = open(args.out,'w')
keys = {}
counts = {}
reverse = []
idx = 0
for line in keyFile:
    vals = line.split()
    key = vals[0]
    if len(key.split("_")) > 3:
        key = "_".join(key.split()[0:3])
    if key not in keys:
        keys[key] = [idx]
    else:
        keys[key].append(idx)

    counts[key] = 0
    idx +=1
sys.stderr.write("Stored " + str(len(keys)) + " keys of " + str(idx) + "\n")
reverse = ["NONE"] * idx
if args.header is not None:
    outFile.write(args.header + "\n")
# skip the header line
bedFile.readline()
nMatches = 0
for line in bedFile:
    vals = line.split()
    key = "_".join(vals[0:3])
    if key in keys:
        outFile.write(key +"\n")
        for idx in keys[key]:
            reverse[idx] = key
        counts[key] +=1
        nMatches +=1
    else:
        outFile.write("NO\n")
outFile.close()
sys.stderr.write("Matched " + str(nMatches) + "\n")
n =0
if args.reverse is not None:
    rf = open(args.reverse, 'w')
    rf.write("pbkey\n")
    rf.write("\n".join(reverse) + "\n")
    rf.close()
