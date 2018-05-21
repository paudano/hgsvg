#!/usr/bin/env python
import sys
tr1f=open(sys.argv[1])
tr2f=open(sys.argv[2])
tr1ann=open(sys.argv[3],'w')
tr2ann=open(sys.argv[4],'w')

tr1l=tr1f.readlines()
tr2l=tr2f.readlines()

def GetTR(l):
    tr = {}
    for line in l:
        vals = line.split()
        if vals[0] == "svType":
            continue
        tr["_".join(vals[0:3])] = int(vals[4])
    return tr

tr1=GetTR(tr1l)
tr2=GetTR(tr2l)


def WriteZyg(lines,h1,h2,out):

    for l in lines:
        vals=l.split()
        k="_".join(vals[0:3])
        zyg="HOM"
        diff=h1[k]
        if k in h2:
            diff = max(h1[k],h2[k])-min(h1[k],h2[k])

            if diff > 50:
                zyg="HET"
        else:
            zyg="HET"
                
        out.write(str(zyg)+"\t"+str(diff)+"\n")


WriteZyg(tr1l,tr1,tr2,tr1ann)
WriteZyg(tr2l,tr2,tr1,tr2ann)



