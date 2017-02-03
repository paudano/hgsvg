import os
import tempfile

#
# A little complicated to find the temp dir
#
SSD_TMP_DIR = "/data/scratch/ssd"
if "TMPDIR" in os.environ:
    TMPDIR = os.environ['TMPDIR']
elif "TMPDIR" in config:
    TMPDIR = config['TMPDIR']
elif os.path.exists(SSD_TMP_DIR):
    TMPDIR = SSD_TMP_DIR
else:
    TMPDIR = tempfile.gettempdir()

configfile: "config.json"


SLOP_FOR_SV_SEQUENCE_POSITIONS = 5000


faiFile = open(config['ref']+".fai")
chroms = [l.split()[0].rstrip() for l in faiFile]

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
cwd=os.getcwd()

shell.prefix(". {SNAKEMAKE_DIR}/config.sh; ")

haps=["h0","h1"]

rule all:
    input:
        asmBed      = expand("{base}.{hap}.bam.bed",base=config['alnBase'], hap=haps),
        contigBed   = expand("overlaps/overlaps.{hap}.{chrom}.ctg0.bed",hap=haps,chrom=chroms),
        asmFasta    = expand("{base}.{hap}.bam.fasta", base=config['alnBase'], hap=haps),
        asmOverlaps = expand("overlaps/overlap.{hap}.{chrom}.txt", hap=haps, chrom=chroms),
        asmGraphs   = expand("overlaps/overlap.{hap}.{chrom}.txt.gml", hap=haps, chrom=chroms),
        asmPaths    = expand("overlaps/overlap.{hap}.{chrom}.txt.path", hap=haps, chrom=chroms),
        asmContigs  = expand("contigs/patched.{hap}.{chrom}.fasta", hap=haps, chrom=chroms),
	chrFasta    = expand("contigs.{hap}.fasta", hap=haps),
	chrAln      = expand("contigs.{hap}.fasta.sam", hap=haps),
        chrBed      = expand("contigs.{hap}.fasta.sam.bed", hap=haps),
	chrBed6     = expand("contigs.{hap}.fasta.sam.bed6", hap=haps),
	chrBB       = expand("contigs.{hap}.fasta.sam.bb", hap=haps),
        annotation  = "stitching_hap_gaps/diploid/insertions.bed"
        


rule MakeAnnotation:
    input:
        asmSam=expand("contigs.{hap}.fasta.sam",hap=haps)
    output:
        annotation="stitching_hap_gaps/diploid/insertions.bed"
    params:
        sge_opts="-l mfree=1G -pe serial 4 -l h_rt=24:00:00"
    shell:
        "make -f " + SNAKEMAKE_DIR + "/DiploidAnnotation.mak H0SAM=contigs.h0.fasta.sam H1SAM=contigs.h1.fasta.sam DIR=stitching_hap_gaps -j 2"
    
rule MakeChrAsmFasta:
    input:
        asmContig=expand("contigs/patched.{{hap}}.{chrom}.fasta", chrom=chroms)
    output:
        asmFasta="contigs.{hap}.fasta"
    params:
        sge_opts="-l mfree=1G -pe serial 1 -l h_rt=01:00:00"
    shell:
        "cat {input.asmContig} > {output.asmFasta}"

rule MakeAsmAln:
    input:
        asmContigs=expand("contigs/patched.{{hap}}.{chrom}.fasta.sam", chrom=chroms),
        aln="alignments.{hap}.bam"
    output:
        asmSam="contigs.{hap}.fasta.sam"
    params:
        sge_opts="-l mfree=1G -pe serial 1 -l h_rt=01:00:00"
    shell:
        "samtools view -H {input.aln} > {output.asmSam}; grep -h -v \"^@\" {input.asmContigs} >> {output.asmSam}"
    
    
rule MakeContigAsmAln:
    input:
        asmFasta="contigs/patched.{hap}.{chrom}.fasta"
    output:
        asmSam=temp("contigs/patched.{hap}.{chrom}.fasta.sam")
    params:
        sge_opts="-l mfree=4G -pe serial 4 -l h_rt=04:00:00",
    	ref=config['ref']
    shell:
        SNAKEMAKE_DIR+ "/MapContigs.py --contigs {input.asmFasta} --ref {params.ref} --tmpdir " + TMPDIR + " --blasr blasr --out {output.asmSam} --nproc 4"

rule MakeChrAsmBed:
    input:
        asmSam="contigs.{hap}.fasta.sam"
    output:
        asmBed="contigs.{hap}.fasta.sam.bed"
    params:
        sge_opts="-l mfree=1G -pe serial 1 -l h_rt=01:00:00",
        hgsvg=SNAKEMAKE_DIR+ "/.."
    shell:
        "/net/eichler/vol5/home/mchaisso/projects/mcutils/bin/samToBed {input.asmSam} --useXS --reportIdentity | bedtools sort > {output.asmBed}"

rule MakeChrAsmBed6:
    input:
        asmBed="contigs.{hap}.fasta.sam.bed"
    output:
        asmBed6="contigs.{hap}.fasta.sam.bed6"
    params:
        sge_opts="-l mfree=1G -pe serial 1 -l h_rt=01:00:00",
	hgsvg=SNAKEMAKE_DIR+ "/.."
    shell:
        "{params.hgsvg}/utils/tracks/SamBedToBed6.py {input.asmBed} {output.asmBed6} "

rule MakeAsmBB:
    input:
        asmBed="contigs.{hap}.fasta.sam.bed6"
    output:
        asmBB="contigs.{hap}.fasta.sam.bb"
    params:
        sge_opts="-l h_rt=01:00:00 -l mfree=1G -pe serial 1",
        ref=config['ref']
    shell:
        "module load ucsc  && bedToBigBed {input.asmBed} {params.ref}.fai {output.asmBB} -type=bed6"


rule MakeAsmContigs:
    input:
        asmPath="overlaps/overlap.{hap}.{chrom}.txt.path",
        asmOverlap="overlaps/overlap.{hap}.{chrom}.txt",
        asmFasta="alignments.{hap}.bam.fasta"
    output:
        asmContig="contigs/patched.{hap}.{chrom}.fasta"
    params:
        sge_opts=" -l h_rt=01:00:00 -l mfree=4G -pe serial 1"
    shell:
        "~/projects/HGSVG/hgsvg/stitching/PatchPaths.py {input.asmOverlap} {input.asmFasta} {input.asmPath} {output.asmContig}"

rule MakeAsmPaths:
    input:
        asmOverlap="overlaps/overlap.{hap}.{chrom}.txt",
        asmOverlapGraph="overlaps/overlap.{hap}.{chrom}.txt.gml"
    output:
        asmPath="overlaps/overlap.{hap}.{chrom}.txt.path"
    params:
        sge_opts="-l h_rt=06:00:00 -l mfree=3G  -pe serial 1"
    shell:
        "~/projects/HGSVG/hgsvg/stitching/OverlapGraphToPaths.py {input.asmOverlap} {input.asmOverlapGraph} {output.asmPath} "
        
rule MakeAsmGraphs:
    input:
        asmOverlap="overlaps/overlap.{hap}.{chrom}.txt"
    output:
        asmOverlapGraph="overlaps/overlap.{hap}.{chrom}.txt.gml"
    params:
        sge_opts="-l h_rt=06:00:00 -l mfree=2G  -pe serial 1"
    shell:
        "~/projects/HGSVG/hgsvg/stitching/OverlapsToGraph.py {input.asmOverlap} --out {output.asmOverlapGraph}"


#subworkflow AsmOverlapsWorkflow:
#    snakefile: SNAKEMAKE_DIR +"/MakeAsmOverlaps.Snakefile"
#    workdir:cwd
#
rule MakeAsmOverlaps:
    input:
        bed="overlaps/overlaps.{hap}.{chrom}.ctg0.bed",
        asm="alignments.{hap}.bam.fasta"
#        asmOverlap="overlaps/overlap.{hap}.{chrom}.txt"
#        asmOverlap=AsmOverlapsWorkflow("overlaps/overlap.{hap}.{chrom}.txt")
    output:
        asmOverlapOut="overlaps/overlap.{hap}.{chrom}.txt"
    params:
        sge_opts="-l mfree=1G -pe serial 12 -l h_rt=72:00:00 -N ovp -e /dev/null -o /dev/null ",
        tmpdir=TMPDIR
    shell:
        "module load anaconda; mkdir -p overlaps/sub_snakemake.{wildcards.hap}.{wildcards.chrom}/overlaps ; cd overlaps/sub_snakemake.{wildcards.hap}.{wildcards.chrom}/overlaps; ln -s ../../../{input.bed} .; cd ..; ln -s ../../{input.asm}* .; ln -s ../../config.json .;  snakemake -cwd -p -s " + SNAKEMAKE_DIR +"/MakeAsmOverlaps.Snakefile -j 20 --cluster \"qsub {params.sge_opts} \" {output.asmOverlapOut} --config " + " ".join(["{}={}".format(k,config[k]) for k in config.keys()] + ["bed={input.bed}"]) + "; mv -f {output.asmOverlapOut} ../;"

                                                                                                                                                          


rule MakeContigBed:
    input:
        asmBed = "alignments.{hap}.bam.bed"
    output:
        contigBed = expand("overlaps/overlaps.{{hap}}.{chrom}.ctg0.bed",chrom=chroms)
    params:
        sge_opts="-l h_rt=06:00:00 -l mfree=2G  -pe serial 1"
    shell:
        """for c in {} ; do  
        grep \"^$c\t\" {{input.asmBed}} > overlaps/overlaps.{{wildcards.hap}}.$c.bed;
        grep "/0	"  overlaps/overlaps.{{wildcards.hap}}.$c.bed > overlaps/overlaps.{{wildcards.hap}}.$c.ctg0.bed
        done || true """.format(" ".join(chroms))
    
rule MakeAsmBed:
    input:
        asmBam = "alignments.{hap}.bam"
    output:
        asmBed = "alignments.{hap}.bam.bed"
    params:
        sge_opts="-l h_rt=06:00:00 -l mfree=2G  -pe serial 1"
    shell:
        "samtools view {input.asmBam} | samToBed /dev/stdin --reportIdentity > {output.asmBed}"

rule MakeAsmFasta:
    input:
        asmBam = "alignments.{hap}.bam"
    output:
        asmFasta = "alignments.{hap}.bam.fasta"
    params:
        sge_opts="-l h_rt=06:00:00 -l mfree=2G  -pe serial 1"
    shell:
        "samtools view {input.asmBam} | awk '{{ print \">\"$1; print $10;}}' | fold | sed '/^$/d' > {output.asmFasta}; samtools faidx {output.asmFasta}"

