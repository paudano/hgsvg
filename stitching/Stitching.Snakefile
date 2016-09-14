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
        asmPaths    = expand("overlaps/overlap.{hap}.{chrom}.txt.path", hap=haps, chrom=chroms)

rule MakeAsmPaths:
    input:
        asmOverlap="overlaps/overlap.{hap}.{chrom}.txt",
        asmOverlapGraph="overlaps/overlap.{hap}.{chrom}.txt.gml"
    output:
        asmPath="overlaps/overlap.{hap}.{chrom}.txt.path"
    params:
        sge_opts="-l h_rt=06:00:00 -l mfree=2G  -pe serial 1"
    shell:
        "~/projects/HGSVG/hgsvg/stitching/OverlapGraphToPaths.py {input.asmOverlap} {input.asmOverlapGraph} {output.asmPath}"
        
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
        sge_opts="-l mfree=1G -pe serial 12 -l h_rt=24:00:00 -N ovp"
    shell:
        "mkdir -p " + TMPDIR + "; mkdir -p overlaps; /net/eichler/vol5/home/mchaisso/projects/HGSVG/scripts/StitchingContigs/OverlapContigsOrderedByBed.py {input.asmBed} {input.asm} --chrom {wildcards.chrom} --out {output.asmOverlap} --nproc 12 --tmpdir " + TMPDIR + " --blasr " + config['blasr']
    
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
    shell: "python ~jlhudd/fasta_tools/merge_assemblies_by_tiling_path.py --gap_size=5000 {input} {output.fasta} _merged > {output.psl}"

rule get_local_assemblies_for_genotyping:
    input: config["local_assembly_alignments"], "contigs_in_final_tiling_path.txt"
    output: "local_assemblies_for_genotyping.fasta", "local_assemblies_for_genotyping.fasta.fai"
    params: sge_opts=""
    shell: "python ~jlhudd/fasta_tools/filter_assemblies_by_name.py {input} > {output[0]}; samtools faidx {output[0]}"

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
    shell: "python ~jlhudd/src/smrtsv/scripts/tiling_path_in_reference_to_contigs.py {input.local_assemblies} {input.tiling_path} | sort -k 1,1 -k 2,2n > {output}"
