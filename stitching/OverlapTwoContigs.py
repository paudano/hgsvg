#!/usr/bin/env python

import tempfile
import subprocess
import argparse
import pysam
import re
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq

ap = argparse.ArgumentParser(description="")
ap.add_argument("--a", help="First file")
ap.add_argument("--b", help="Second file")
ap.add_argument("--asm", help="Assemblies")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--tmpdir", help="Use this directory for temp", default=".")
ap.add_argument("--keep", help="Keep output file.",default=True,action='store_false')
ap.add_argument("--blasr",  help="Specify alternative blasr", default="/net/eichler/vol5/home/mchaisso/projects/blasr/cpp/alignment/bin/blasr")
args = ap.parse_args()


def GetSeq(name, ref, tempFile):
    faidx = "samtools faidx {} {}".format(ref, name)
    p = subprocess.Popen(faidx.split(), stdout=tempFile)
    p.wait()

aTempFile = tempfile.NamedTemporaryFile(suffix=".fasta", dir=args.tmpdir, delete=args.keep)
bTempFile = tempfile.NamedTemporaryFile(suffix=".fasta", dir=args.tmpdir, delete=args.keep)
samTempFile = tempfile.NamedTemporaryFile(suffix=".sam", dir=args.tmpdir, delete=args.keep)
GetSeq(args.a, args.asm, aTempFile)
GetSeq(args.b, args.asm, bTempFile)

seqA = SeqIO.read(open(aTempFile.name), "fasta")
seqB = SeqIO.read(open(bTempFile.name), "fasta")
nre=re.compile("^(N+)[^N]*(N+)$")


def GetN(seq):
    m= nre.match(seq.seq.tostring())
    if m is None:
        return (0,len(seq.seq))
    else:
        g =m.groups()
        return (len(g[0]), len(seq.seq) - len(g[1]))


command = "{} {} {} -bestn 1 -sam -out {} -clipping soft -maxMatch 25 -sdpTupleSize 13 -extend -maxExtendDropoff 50 ".format(args.blasr, aTempFile.name, bTempFile.name, samTempFile.name)
subprocess.call(command.split())

(aNPre, aNPost) = GetN(seqA)
(bNPre, bNPost) = GetN(seqB)

#
# parse the sam file to get the cooridnates of the overlap
#
aln = pysam.AlignmentFile(samTempFile.name, "r")
reads = [ r for r in aln.fetch() ]
if len(reads) == 0:
    coords = [0,0,0,0]
else:
    coords = [reads[0].qstart, reads[0].qend, reads[0].reference_start, reads[0].reference_end]

print "{}\t{}\t{}\t".format(args.a, aNPre, aNPost) + "{}\t{}\t{}\t".format(args.b, bNPre, bNPost) +  "\t".join([str(i) for i in coords])
