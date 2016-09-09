#!/usr/bin/env python

import tempfile
import subprocess
import argparse
import pysam

ap = argparse.ArgumentParser(description="")
ap.add_argument("--a", help="First file")
ap.add_argument("--b", help="Second file")
ap.add_argument("--asm", help="Assemblies")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--tmpdir", help="Use this directory for temp", default=".")
args = ap.parse_args()


def GetSeq(name, ref, tempFile):
    faidx = "samtools faidx {} {}".format(ref, name)
    p = subprocess.Popen(ct.split(), stdout=tempFile)
    p.wait()

aTempFile = tempfile.NamedTemporaryFile(suffix=".fasta", dir=args.tmpdir)
bTempFile = tempfile.NamedTemporaryFile(suffix=".fasta", dir=args.tmpdir)
samTempFile = tempfile.NamedTemporaryFile(suffix=".sam", dir=args.tmpdir)
GetSeq(args.a, args.asm, aTempFile.name)
GetSeq(args.b, args.asm, bTempFile.name)

command = "blasr {} {} -bestn 1 -sam -out {}".format(aTempFile.name, bTempFile.name, samTempFile.name)
subprocess.call(command.split())


#
# parse the sam file to get the cooridnates of the overlap
#
aln = pysam.AlignmentFile(samTempFile.name, "r")





