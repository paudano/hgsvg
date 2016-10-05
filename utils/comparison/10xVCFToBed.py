#!/usr/bin/env python

import sys
inFile = open(sys.argv[1])
outFile = open(sys.argv[2],'w')
#chr1    789462  call_2099_1     G       ]chr1:224014604]G       3       LOWQ;SEG_DUP    CIPOS=-8,8;SVTYPE=BND;SVTYPE2=BND;PS=560256;MATEID=call_2099_2;EVENT=call_2099;HAP_ALLELIC_FRAC=0.5;ALLELIC_FRAC=0.0142405063291;PAIRS=45;SPLIT=0;WildCov=.;MolTtl=.;MolTtlNoR=.;MolDel=.;MolDelNoR=.;MolWild=.;MolWildNoR=.;PVAL=.;BGSize=.;BGImparityPval=.;BGTtlRCnt=.;BGHP1RCnt=.;BGHP2RCnt=.;BGBaseCov=.;BGPhaseFrac=.;NRead=.;SOURCE=.    GT      0/1
#chr1    224236782       call_5034       T       <UNK>   2       LOWQ    CIEND=-18,18;CIPOS=-18,18;SVTYPE=UNK;SVLEN=45057;PS=224026116;HAP_ALLELIC_FRAC=0.117647058824;ALLELIC_FRAC=0.00689655172414;PAIRS=2;SPLIT=0;WildCov=.;MolTtl=.;MolTtlNoR=.;MolDel=.;MolDelNoR=.;MolWild=.;MolWildNoR=.;PVAL=.;BGSize=.;BGImparityPval=.;BGTtlRCnt=.;BGHP1RCnt=.;BGHP2RCnt=.;BGBaseCov=.;BGPhaseFrac=.;NRead=.;SOURCE=.;END=224281839    GT      0|1

def GetEnd(info):
    v = info.split(";")
    for p in v:
	kv = p.split("=")
	if len(kv) == 2:
            [key,value] = p.split("=")
            if key == "END":
	        return int(value)
    return None

for line in inFile:
    if line[0] == "#":
        continue
    if line[0:3] == "LEN":
        continue
#    if "END" in line:
#	import pdb
#	pdb.set_trace()
    v=line.split()
    start = int(v[1])
    info  = v[7]
    end = GetEnd(info)

    if end is not None:
        outFile.write("{}\t{}\t{}\t{}\t{}\n".format(v[0], start, end, end-start, v[4]))

