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

shell.prefix(". {SNAKEMAKE_DIR}/config.sh; ")

haps=[0,1]

rule all:
    input:
        asmBed      = expand("{base}.{hap}.bam.bed",base=config['alnBase'], hap=haps),
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
        "make -f " + SNAKEMAKE_DIR + "/DiploidAnnotation.mak H0SAM=contigs.0.fasta.sam H1SAM=contigs.1.fasta.sam DIR=stitching_hap_gaps -j 2"
    
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
        asmFasta="contigs.{hap}.fasta"
    output:
        asmSam="contigs.{hap}.fasta.sam"
    params:
        sge_opts="-l mfree=15G -pe serial 4 -l h_rt=04:00:00",
    	ref=config['ref']
    shell:
        "blasr {input.asmFasta} {params.ref} -alignContigs -sam -minMapQV 30 -out {output.asmSam} -piecewise -nproc 4"

rule MakeChrAsmBed:
    input:
        asmSam="contigs.{hap}.fasta.sam"
    output:
        asmBed="contigs.{hap}.fasta.sam.bed"
    params:
        sge_opts="-l mfree=1G -pe serial 1 -l h_rt=01:00:00",
        hgsvg=SNAKEMAKE_DIR+ "/.."
    shell:
        "/net/eichler/vol5/home/mchaisso/projects/mcutils/bin/samToBed {input.asmSam} --reportIdentity | bedtools sort > {output.asmBed}"

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
        "module load uccc  && bedToBigBed {input.asmBed} {params.ref}.fai {output.asmBB} -type=bed6"


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



rule MakeAsmOverlaps:
    input:
        asmBed = "alignments.h{hap}.bam.bed",
        asm = "alignments.{hap}.bam.fasta"
    output:
        asmOverlap="overlaps/overlap.{hap}.{chrom}.txt"
    params:
        sge_opts="-l mfree=1G -pe serial 12 -l h_rt=72:00:00 -N ovp"
    shell:
        "mkdir -p " + TMPDIR + "; mkdir -p overlaps; /net/eichler/vol5/home/mchaisso/projects/HGSVG/hgsvg/stitching/OverlapContigsOrderedByBed.py {input.asmBed} {input.asm} --chrom {wildcards.chrom} --out {output.asmOverlap} --nproc 12 --tmpdir " + TMPDIR + " --blasr " + config['blasr']
    
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
        asmBam = "alignments.h{hap}.bam"
    output:
        asmFasta = "alignments.{hap}.bam.fasta"
    params:
        sge_opts="-l h_rt=06:00:00 -l mfree=2G  -pe serial 1"
    shell:
        "samtools view {input.asmBam} | awk '{{ print \">\"$1; print $10;}}' | fold | sed '/^$/d' > {output.asmFasta}; samtools faidx {output.asmFasta}"


rule convert_psl_to_chain:
    input: "contigs_to_merged_assemblies.psl"
    output: "contigs_to_merged_assemblies.chain"
    params: sge_opts="-l h_rt=00:05:00"
    shell: "pslToChain {input} {output}"

rule merge_local_assemblies_into_chromosomes:
    input: "local_assemblies_for_genotyping.fasta", "tiling_path_in_contigs.bed"
    output: fasta="merged_local_assemblies_for_genotyping.fasta", psl="contigs_to_merged_assemblies.psl"
    params: sge_opts="-l h_rt=00:15:00"
    shell: "python /net/eichler/vol4/home/jlhudd/fasta_tools/merge_assemblies_by_tiling_path.py --gap_size=5000 {input} {output.fasta} _merged > {output.psl}"

rule get_local_assemblies_for_genotyping:
    input: config["local_assembly_alignments"], "contigs_in_final_tiling_path.txt"
    output: "local_assemblies_for_genotyping.fasta", "local_assemblies_for_genotyping.fasta.fai"
    params: sge_opts=""
    shell: "python /net/eichler/vol4/home/jlhudd/fasta_tools/filter_assemblies_by_name.py {input} > {output[0]}; samtools faidx {output[0]}"

rule get_contigs_in_final_tiling_path:
    input: "tiling_path_in_contigs.bed"
    output: "contigs_in_final_tiling_path.txt"
    params: sge_opts="-l h_rt=00:05:00"
    shell: "cut -f 4 {input} | sort | uniq > {output}"

# Convert tiling path through uncontained assemblies with SVs into corresponding
# space in contigs
rule convert_tiling_path_into_local_assembly_coordinates:
    input: local_assemblies=config["local_assembly_alignments"], tiling_path="padded_tiling_path_of_contigs_with_variants.bed"
    output: "tiling_path_in_contigs.bed"
    params: sge_opts="-l h_rt=02:00:00"
    shell: "python /net/eichler/vol4/home/jlhudd/src/smrtsv/scripts/tiling_path_in_reference_to_contigs.py {input.local_assemblies} {input.tiling_path} | sort -k 1,1 -k 2,2n > {output}"
