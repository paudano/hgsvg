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


SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
shell.prefix(". {SNAKEMAKE_DIR}/config.sh; ")

cwd=os.getcwd()

svTypes=["DEL", "INS"]
pbSVTypeMap={"DEL": "deletion", "INS": "insertion"}
bnSVTypeMap={"DEL": "del", "INS": "ins"}

callHash={}
callFile = open(config["calls"])
callFile.readline()
for line in callFile:
    vals = line.split()
    callHash[vals[5]] = True

callers = callHash.keys()

samples=config["individuals"]



rule all:
    input:
        methodConfig=expand("{method}.{sample}/sv_support.json", method=callers, sample=samples),
        allCalls=expand("{method}.{sample}/calls.bed", method=callers, sample=samples),
        allVCF=expand("{method}.{sample}/calls.vcf", method=callers, sample=samples),
#        allMergged=expand("{method}.{sample}/merged_ortho.{svtype}.bed", method=callers, sample=samples, svtype=svTypes)
        




rule ExtractCalls:
    input:
        calls=config["calls"]
    output:
        sampleCalls="{method}.{sample}/calls.bed"
    params:
        sge_opts=config["sge_small"],
    shell:"""
echo -e "#chrom\ttStart\ttEnd\tsvType\tgenotype\tcaller\tsample\tEND\tsvLen\tqName\tqStart\tqEnd\tsvSeq\tMERGE_TYPE\tNUM_CALLER\tCALLER" > {output.sampleCalls}
grep {wildcards.method} {input.calls} | \
grep  {wildcards.sample} | \
awk '{{ print $0"\t"$3"\t"$3-$2"\tnone\t0\t0\tA\tunique\t1\t{wildcards.method}";}}' |\
 sed "s/\(DEL[^\\t]*\\t\)/deletion\\t/" | \
 sed "s/\(INS[^\\t]*\\t\)\|\(DUP[^\\t]*\t\)/insertion\\t/"  >> {output.sampleCalls}
"""

rule MakeVCF:
    input:
        sampleCalls="{method}.{sample}/calls.bed"
    output:
        sampleVcf="{method}.{sample}/calls.vcf"
    params:
        sge_opts=config["sge_small"],
        sd=SNAKEMAKE_DIR,       
        ref=config["ref"],
    shell:"""
{params.sd}/../../utils/variants_bed_to_vcf.py --bed {input.sampleCalls} --vcf {output.sampleVcf} --type sv    --addci 5 --fields  MERGE_TYPE MERGE_TYPE  NUM_CALLER  NUM_CALLER CALLER CALLER END END  --reference {params.ref}
"""
    
rule MakeConfig:
    input:
        calls=config["calls"]
    output:
        config="{method}.{sample}/sv_support.json"
    params:
        sge_opts=config["sge_small"],
    run:
        configStr="""
{{
                "svvcf" : "calls.vcf",
                "bionano-vcf" : "/net/eichler/vol24/projects/structural_variation/nobackups/projects/HGSVG/analysis/{}/BioNano/{}.bionano.vcf.gz",
                "pbfinal": "/net/eichler/vol24/projects/structural_variation/nobackups/projects/HGSVG/analysis/{}/uw_mssm_merging/sv_calls.bed",
                "localasm" : [ "/net/eichler/vol24/projects/structural_variation/nobackups/projects/HGSVG/analysis/{}/alignments.h0.bam",
                   "/net/eichler/vol24/projects/structural_variation/nobackups/projects/HGSVG/analysis/{}/alignments.h1.bam"],
                "local": ["/net/eichler/vol24/projects/structural_variation/nobackups/projects/HGSVG/analysis/{}/hap_gaps/hap0/gaps.bed",
                                                        "/net/eichler/vol24/projects/structural_variation/nobackups/projects/HGSVG/analysis/{}/hap_gaps/hap1/gaps.bed"],
                "stitch": ["/net/eichler/vol24/projects/structural_variation/nobackups/projects/HGSVG/analysis/{}/SVQC/hap0/gaps.bed.support",
                                                         "/net/eichler/vol24/projects/structural_variation/nobackups/projects/HGSVG/analysis/{}/SVQC/hap1/gaps.bed.support"],
                "ref": "/net/eichler/vol2/eee_shared/assemblies/GRCh38/GRCh38.fasta",
                "window": 1000,
                "fofn" : "/net/eichler/vol24/projects/structural_variation/nobackups/projects/HGSVG/Alignments/{}/{}.alignments.fofn",
                "sge_small" : " -pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
                "sample" : "{}",
                "bn_read_overlap": "/net/eichler/vol24/projects/structural_variation/nobackups/projects/HGSVG/analysis/IlluminaCombined/{}/bn_read_overlap",
                "bn_read_cov": "/net/eichler/vol24/projects/structural_variation/nobackups/projects/HGSVG/analysis/IlluminaCombined/{}/bn_read_cov",
                "method": "{}"
        
     
}}
""".format(wildcards.sample, wildcards.sample, wildcards.sample, wildcards.sample, wildcards.sample, wildcards.sample, wildcards.sample, wildcards.sample, wildcards.sample, wildcards.sample,wildcards.sample, wildcards.sample, wildcards.sample, wildcards.sample, wildcards.method)
        outFile = open(output.config, 'w')
        outFile.write(configStr)

