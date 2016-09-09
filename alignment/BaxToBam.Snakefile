shell.prefix("source config.sh;")

configfile: "config.json"

import os
import glob
import csv



f=open(config["fofn"])
filenames=[l.rstrip() for l in f.readlines()]
datasets = {}
for f in filenames:
    base = os.path.basename(f).split(".")[0]
    if base not in datasets:
        datasets[base] = []
    datasets[base].append(f)


rule all:
    input:
        expand("{base}.subreads.bam", base=datasets.keys())

BAX2BAM="/net/eichler/vol18/zevk/great_apes/iso_seq/cc2_analysis/pitchfork/deployment/bin/bax2bam"

rule bax2bam:
    input:
        BAX=lambda wildcards: expand("{files}", files=datasets[wildcards.base])
    output:
        "{base}.subreads.bam"
    params:
        sge_opts="-l mfree=5G -l h_rt=24:00:00 -q eichler-short.q"
    shell:
        "{BAX2BAM} {input.BAX}; rm -f {base}.scraps.bam; rm -f {base}.scraps.bam.pbi"


