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

configfile: "bins.json"


SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
cwd=os.getcwd()


rule all:
    input:
        ljGaps=expand("{gaps}.lj", gaps=config["gaps"]),
        wideGaps=expand("{gaps}.wide", gaps=config["gaps"]),
        netGain=expand("{gaps}.netGain", gaps=config["gaps"]),        
        allBins="bin_expansion/regions_with_vars.bed",
        allNetGain="bin_expansion/netgain.tab",
        allNetGainSD="bin_expansion/netgain.tab.sd",
        topNetGainSD="bin_expansion/netgain.tab.sd.top",
        topNetGainSDBed="bin_expansion/netgain.tab.sd.bed",
        topNetGainFasta=expand("bin_expansion/{sample}.{hap}.fasta", sample=config["samples"],hap=config["haps"]),
        expanded="bin_expansion/Expanded.pdf",
        dotplots="bin_expansion/Dotplots.pdf",
        geneGains="bin_expansion/netgain.tab.sd.genes.bed",

samples=config["samples"]
haps=config["haps"]
sampleStr = "\t".join([samples[s] + "_" +haps[h] for s in range(0,len(samples)) for h in range(0,len(haps)) ])

rule MakeNetGainTopBed:
    input:
        top="bin_expansion/netgain.tab.sd.top",
    output:
        bed="bin_expansion/netgain.tab.sd.bed",
    params:
        ref=config["ref"],
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        sampleHeader=sampleStr
    shell:"""
echo -e "#chrom\tstart\tend\t"{params.sampleHeader}"\tsd" > {output.bed}
cat {input.top} | tr "\:" "\\t" | sed "s/\\-/\\t/" >> {output.bed}
"""

rule IntersectBinsWithExons:
    input:
        gains="bin_expansion/netgain.tab.sd.bed",
        exons=config["genes"],
    output:
        geneGains="bin_expansion/netgain.tab.sd.genes.bed",
    params:
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        ref=config["ref"],
    shell:"""
module load bedtools/2.25.0
bedtools slop -header -i {input.gains} -g {params.ref}.fai -b 1000 | \
 bedtools intersect -header -a stdin -b {input.exons} -wa -wb | \
 bedtools groupby -header -g 1-3 -c 10 -o first -full | \
 bioawk -c hdr '{{ if (substr($0,0,1) == "#") {{ print $0"\texonChrom\texonStart\texonEnd\tgene";}} else {{ print $0; }} }}'| \
 cut -f 1-14 > {output.geneGains}
"""
    
rule MakeLJ:
    input:
        gaps="{gaps}"
    output:
        lj="{gaps}.lj"
    params:
        ref=config["ref"],
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        sd=SNAKEMAKE_DIR
    shell:"""
cat {input.gaps} | {params.sd}/../../utils/LeftJustifyBed.py > {output.lj}
"""

rule MakeWide:
    input:
        lj="{gaps}.lj"
    output:
        wide="{gaps}.wide"
    params:
        ref=config["ref"],
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        sd=SNAKEMAKE_DIR,
        slop=config["slop"]
    shell:"""
 bedtools slop -b {params.slop} -g {params.ref}.fai -i {input.lj} | bedtools sort  |bedtools merge > {output.wide}
"""

rule MergeWide:
    input:
        wide=expand("{gaps}.wide",gaps=config["gaps"])
    output:
        bins="bin_expansion/regions_with_vars.bed"
    params:
        ref=config["ref"],
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        sd=SNAKEMAKE_DIR,
        slop=config["slop"]
    shell:"""
cat {input.wide} | bedtools sort | bedtools merge > {output.bins}
"""
    
rule GetBinNetGain:
    input:
        lj="{bins}.lj",
        combined="bin_expansion/regions_with_vars.bed"
    output:
        netGain="{bins}.netGain"
    params:
        ref=config["ref"],
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        sd=SNAKEMAKE_DIR,
    shell:"""
bedtools intersect -a {input.combined} -b {input.lj} -loj | ~mchaisso/projects/HGSVG/hgsvg/sv/analysis/bin_analysis/FlipLOJ.py | bedtools groupby -g 1-3 -c 8 -o  sum -full | cut -f 14 > {output.netGain}
"""

rule CombineNetGain:
    input:
        netGain=expand("{gaps}.netGain", gaps=config["gaps"]),
        combined="bin_expansion/regions_with_vars.bed"
    output:
        allNetGain="bin_expansion/netgain.tab",
    params:
        ref=config["ref"],
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        sd=SNAKEMAKE_DIR,
    shell:"""
paste {input.combined} {input.netGain} > {output.allNetGain}
"""
rule AddSD:
    input:
        allNetGain="bin_expansion/netgain.tab",
    output:
        allNetGainSD="bin_expansion/netgain.tab.sd",
    params:
        ref=config["ref"],
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
    shell:"""
module unload python/3.5.2
module load python/2.7.3
module load numpy
~mchaisso/projects/HGSVG/hgsvg/sv/analysis/bin_analysis/AddSD.py {input.allNetGain} --out {output.allNetGainSD}
"""

rule SelectTopSD:
    input:
        allNetGainSD="bin_expansion/netgain.tab.sd",
    output:
        topNetGainSD="bin_expansion/netgain.tab.sd.top",
    params:
        ref=config["ref"],
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=1G",
        minSD=config["minSD"]
    shell:"""
cat {input.allNetGainSD} | awk '{{ if ($10  > {params.minSD}) print $1":"$2"-"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10;}}' > {output.topNetGainSD}
"""

# If you wanted  to make the grid dotplots of a number of regions, the
# entry point is here.  You would use a different file for
# nentgain.tab.sd.top, with the first column a bunch of regions. 


rule ExtractFasta:
    input:
        topNetGainSD="bin_expansion/netgain.tab.sd.top",
    output:
        topNetGainFasta="bin_expansion/{sample}.{hap}.fasta",
    params:
        ref=config["ref"],
        sge_opts="-pe serial 1 -l h_rt=1:00:00 -l mfree=8G",
    shell:"""
~mchaisso/projects/mcutils/src/samSubseq -s {wildcards.sample}/contigs.{wildcards.hap}.fasta.sam -R {input.topNetGainSD} -f {output.topNetGainFasta} -n
"""

        
    
rule MakeExpandRegionMSAs:
    input:
        topNetGainFasta=expand("bin_expansion/{sample}.{hap}.fasta",sample=config["samples"],hap=config["haps"]),
        topNetGainSD="bin_expansion/netgain.tab.sd.top",
    output:
        split=dynamic("bin_expansion/{rgn}/msa.fasta"),
    params:
        sge_opts="-pe serial 1 -l mfree=1G -l h_rt=1:00:00 -cwd -e /dev/null -o /dev/null  ",
        sd=SNAKEMAKE_DIR,
        samples=config["samples"],
        haps=config["haps"]
    shell:"""
module unload python/3.5.2
module load python/2.7.3
for rgn in `cut -f 1 {input.topNetGainSD}`; do
    echo $rgn
    rgnName=`echo $rgn | tr ":-" "__"`
    mkdir -p bin_expansion/$rgnName
    for sample in {params.samples}; do
        for hap in {params.haps}; do
            samtools faidx bin_expansion/$sample.$hap.fasta $rgn | ~mchaisso/projects/HGSVG/hgsvg/sv/analysis/bin_analysis/rename.py $sample.$hap >> bin_expansion/$rgnName/msa.fasta
        done
    done
done
"""

rule PlotRegions:
    input:
        msa="bin_expansion/{rgn}/msa.fasta",
    output:
        plot="bin_expansion/{rgn}/threshold300.pdf"
    params:
        sge_opts="-pe serial 1 -l mfree=4G -l h_rt=1:00:00 -cwd -e /dev/null -o /dev/null  ",
        sd=SNAKEMAKE_DIR,
        samples=config["samples"],
        haps=config["haps"]
    shell:"""
cd bin_expansion/{wildcards.rgn}
module unload python/3.5.2
module load python/2.7.3
module load numpy/1.7.0
module load biopython
~mchaiss/projects/HGSVG/hgsvg/sv/analysis/bin_analysis/FixShort.py msa.fasta msa.fix.fasta
miropeats -s 300 -onlyinter msa.fix.fasta
~mchaisso/projects/HGSVG/hgsvg/sv/analysis/FixMiropeatsPS.py threshold300 threshold300.ps
ps2pdf threshold300.ps
"""


rule MakePairwiseDotplots:
    input:
        msa="bin_expansion/{rgn}/msa.fasta",
    output:
        pdf="bin_expansion/{rgn}/msa.pdf",
    params:
        sge_opts="-pe serial 1 -l mfree=4G -l h_rt=1:00:00 -cwd -e /dev/null -o /dev/null  ",
    shell:"""
module unload python/3.5.2
module load python/2.7.3
module load R/latest
cd bin_expansion/{wildcards.rgn}
/net/eichler/vol5/home/mchaisso/projects/HGSVG/hgsvg/sv/analysis/bin_analysis/MakePairwiseDotplots.py --msa msa.fasta --out msa.pdf --region {wildcards.rgn}
    convert -density 150 -background white msa.pdf msa.png
"""
    
    
rule MakeExpandRegionsList:
    input:
        split=dynamic("bin_expansion/{rgn}/threshold300.pdf"),
    output:
        expanded="bin_expansion/Expanded.pdf"
    params:
        sge_opts="-pe serial 1 -l mfree=1G -l h_rt=1:00:00 -cwd -e /dev/null -o /dev/null  ",
        sd=SNAKEMAKE_DIR,
    shell:"""
pdfunite {input.split} {output.expanded}
"""
    
        
rule CombineDotplots:
    input:
        split=dynamic("bin_expansion/{rgn}/msa.pdf"),
    output:
        comb="bin_expansion/Dotplots.pdf"
    params:
        sge_opts="-pe serial 1 -l mfree=1G -l h_rt=1:00:00 -cwd -e /dev/null -o /dev/null  ",
        sd=SNAKEMAKE_DIR,
    shell:"""
pdfunite {input.split} {output.comb}
"""
    
