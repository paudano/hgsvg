#!/usr/bin/env python

import argparse
from multiprocessing import Process, Lock, Semaphore, Pool
import sys
import subprocess

ap = argparse.ArgumentParser(description="Given an alignmetn bed file, run overlaps between contigs that are near each other")
ap.add_argument("bed", help="Alignment bed file.")
ap.add_argument("asm", help="Assemblies file")
ap.add_argument("--chrom", help="Process this chromosome only", default =None)
ap.add_argument("--path", help="Path to HGSVG stitching scripts", default="/net/eichler/vol5/home/mchaisso/projects/HGSVG/hgsvg/stitching")
ap.add_argument("--ahead", help="Allow overlap betewen alignemnts ending this far ahead of current.", type=int, default=20000)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--nproc", help="Number of processors", default=1,type=int)
ap.add_argument("--tmpdir", help="Pass along tempdir", default=".")
ap.add_argument("--blasr", help="Blasr command to run", default=None)
args = ap.parse_args()

def ParseBedLine(line):
    v = line.split()
    return [v[0], int(v[1]), int(v[2]), v[3]]

bedFile = open(args.bed)
outFile = open(args.out,'w')
alignments= [ParseBedLine(line) for line in bedFile]

        
commands = []
def Run(command):
    res = subprocess.check_output(command.split())
    return res.rstrip()
if args.blasr is not None:
    blasrCommand = " --blasr " + args.blasr
else:
    blasrCommand = ""
for i in range(0,len(alignments)-1):
    if args.chrom is not None and alignments[i][0] != args.chrom:
        continue
    j = i+1
    
    while j < len(alignments) and alignments[j][0] == alignments[i][0]  and  alignments[j][1] < alignments[i][2]+args.ahead:
        commands.append( "{}/OverlapTwoContigs.py --a {} --b {} --asm {} --tmpdir {} {}".format(args.path, alignments[i][3], alignments[j][3], args.asm, args.tmpdir, blasrCommand))
	sys.stderr.write(commands[-1] + "\n")
        j+=1
    i+=1

#
# Now process alignments
#
pool = Pool(args.nproc)
ovp = pool.map(Run, commands)
pool.close()
outFile.write("\n".join(ovp))
