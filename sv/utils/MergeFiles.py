#!/usr/bin/env python


import argparse
import sys

ap = argparse.ArgumentParser(description="Merge files with headers so they have the same number of columns")
ap.add_argument("--files", help="Input files.", nargs="+")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--default", help="Default value", default="0")
ap.add_argument("--append-dups", help="Append duplicate columns.", dest="appendDups", default=False, action="store_true")

args = ap.parse_args()

out = open(args.out,'w')
files = [open(fn) for fn in args.files]

headers = [ files[i].readline()[1:].split() for i in range(0,len(files))]
seen = {}
if args.appendDups:
    print "appendign dups"
    for i in range(0,len(headers)):
        for hi in range(0,len(headers[i])):
            if headers[i][hi] in seen:
                headers[i][hi] = headers[i][hi] + "_" + str(i)
            seen[headers[i][hi]] = True

masterHeader = headers[0][:]
for header in headers[1:]:
    for v in header:
        if v not in masterHeader:
            masterHeader.append(v)

out.write("#"+ "\t".join(masterHeader)+ "\n")
for fi in range(0, len(files)):
    #
    # Determine the order to write the values
    #
    cols = {}

    for hi in range(0,len(headers[fi])):
        cols[headers[fi][hi]] = hi

    for l in files[fi]:
        v = l.split()
        line = []
#        if len(v) != len(cols):
#            print "col mismatch in " + args.files[fi]
#            sys.exit(0)
#            continue
        for col in masterHeader:
            if col in cols.keys():
                if cols[col] >= len(v):
#                    import pdb
#                    pdb.set_trace()
                    #sys.stderr.write("Invalid merge on " + col + " from " + str(cols) + "\n")
                    #sys.exit(1)
                    line.append(".")
                    
                else:
                    line.append( v[cols[col]])
            else:
                line.append(args.default)
        out.write("\t".join(line) + "\n")
        
        
    

