configfile: "config.json"

import os
bamFile = open(config['bams'])
bams = [ l.strip() for l in bamFile ]
base = { b : os.path.basename(b) for b in bams }
bamMap = { os.path.basename(b): b for b in bams }

rule all: 
	input:
		expand("{b}.vcf", b=bamMap.keys())
def BamToBase(wildcards):
	return bamMap[wildcards.base]

rule makevcf:
	input:
		bam=BamToBase,
		ref=config['ref'],
	params:
		coverage=config['coverage'],
	        sge_opts="-l mfree=8G -l h_rt=24:00:00 -q eichler-short.q"
	output:
		"{base}.vcf"
	
	shell:
		"samtools mpileup -q0 -Q 0 {input.bam} | /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/Phasing/MpileupToFreq.py  /dev/stdin | /net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/Phasing/PrintHetFreq.py {params.coverage} --maxCount 60000 | ~/projects/AssemblyByPhasing/scripts/abp/FreqToSimpleVCF.py --freq /dev/stdin  --out {output} --ref {input.ref}"
