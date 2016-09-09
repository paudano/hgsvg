PBS=/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts
HGSVGS=/net/eichler/vol5/home/mchaisso/projects/HGSVG/scripts
BLASR=/net/eichler/vol5/home/mchaisso/software/blasr_2/cpp/alignment/bin/blasr
REF=/net/eichler/vol2/eee_shared/assemblies/GRCh38/GRCh38.fasta

all: hap_gaps/hap0/gaps.bed \
  hap_gaps/hap1/gaps.bed \
  hap_gaps/hap0/insertions.bed \
  hap_gaps/hap1/insertions.bed \
  hap_gaps/diploid/insertions.bed \
  hap_gaps/diploid/deletions.bed \
  alignments.h0.bam \
  alignments.h1.bam \
  hap_gaps/diploid/deletion/L1.bed \
  hap_gaps/diploid/insertion/L1.bed 

samfiles.fofn: samfiles/chr16.12240000-12300000.sam
	ls samfiles | awk '{ print "samfiles/"$$1;}' > $@

alignments.h0.sam: samfiles.fofn
	$(HGSVGS)/FilterSamByHaplotype.py samfiles.fofn

alignments.h0.bam: alignments.h0.sam
	mkdir -p /var/tmp/mchaisso
	samtools view -bS alignments.h0.sam | samtools sort -T /var/tmp/mchaisso/aln.0.$$PPID -o $@
	samtools index $@

alignments.h1.bam: alignments.h1.sam
	mkdir -p /var/tmp/mchaisso
	samtools view -bS alignments.h1.sam | samtools sort -T /var/tmp/mchaisso/aln.1.$$PPID -o $@
	samtools index $@

hap_gaps/hap0/gaps.bed: alignments.h0.bam
	mkdir -p hap_gaps/hap0
	samtools view $^ | $(PBS)/PrintGaps.py $(REF) /dev/stdin --minAlignmentLength 30000 --flankIdentity 0.7 --flank 1000 --condense 20 --outFile $@

hap_gaps/hap1/gaps.bed: alignments.h1.bam
	mkdir -p hap_gaps/hap1
	samtools view $^ | $(PBS)/PrintGaps.py $(REF) /dev/stdin --minAlignmentLength 30000 --flankIdentity 0.7 --flank 1000 --condense 20 --outFile $@

hap_gaps/hap0/insertions.bed: hap_gaps/hap0/gaps.bed
	cd hap_gaps/hap0 && make -f $(PBS)/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=$(PWD)/alignments.h0.bam

hap_gaps/hap1/insertions.bed: hap_gaps/hap1/gaps.bed
	cd hap_gaps/hap1 && make -f $(PBS)/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=$(PWD)/alignments.h1.bam

hap_gaps/diploid/deletions.bed: hap_gaps/hap0/deletions.bed hap_gaps/hap1/deletions.bed
	mkdir -p hap_gaps/diploid
	bedtools intersect -a hap_gaps/hap0/deletions.bed -b hap_gaps/hap1/deletions.bed -r -f 0.5 -wao | \
    awk '{ if ($$11 != ".") print; }' | \
    cut -f 1-10 | \
    bedtools groupby -c 4 -o first -full > hap_gaps/diploid/deletions.0.r.bed
	bedtools intersect -b hap_gaps/hap0/deletions.bed -a hap_gaps/hap1/deletions.bed -r -f 0.5 -wao | \
    awk '{ if ($$11 != ".") print; }' | \
    cut -f 1-10 | \
    bedtools groupby -c 4 -o first -full > hap_gaps/diploid/deletions.1.r.bed
	bedtools intersect -v -a hap_gaps/hap0/deletions.bed -b hap_gaps/diploid/deletions.0.r.bed -r -f 1.0 > hap_gaps/diploid/deletions.h0.bed
	bedtools intersect -v -a hap_gaps/hap1/deletions.bed -b hap_gaps/diploid/deletions.1.r.bed -r -f 1.0 > hap_gaps/diploid/deletions.h1.bed
	cp hap_gaps/diploid/deletions.1.r.bed hap_gaps/diploid/deletions.hom.bed
	cat hap_gaps/diploid/deletions.h0.bed hap_gaps/diploid/deletions.h1.bed hap_gaps/diploid/deletions.hom.bed > $@
#	-rm -f hap_gaps/diploid/deletions.1.r.bed 

hap_gaps/diploid/insertions.bed: hap_gaps/hap0/insertions.bed hap_gaps/hap1/insertions.bed
	mkdir -p hap_gaps/diploid
	bedtools intersect -a hap_gaps/hap0/insertions.bed -b hap_gaps/hap1/insertions.bed -r -f 0.5 -wao | \
    awk '{ if ($$11 != ".") print; }' | \
    cut -f 1-10 | \
    bedtools groupby -c 4 -o first -full > hap_gaps/diploid/insertions.0.r.bed
	bedtools intersect -b hap_gaps/hap0/insertions.bed -a hap_gaps/hap1/insertions.bed -r -f 0.5 -wao | \
    awk '{ if ($$11 != ".") print; }' | \
    cut -f 1-10 | \
    bedtools groupby -c 4 -o first -full > hap_gaps/diploid/insertions.1.r.bed
	bedtools intersect -v -a hap_gaps/hap0/insertions.bed -b hap_gaps/diploid/insertions.0.r.bed -r -f 1.0 > hap_gaps/diploid/insertions.h0.bed
	bedtools intersect -v -a hap_gaps/hap1/insertions.bed -b hap_gaps/diploid/insertions.1.r.bed -r -f 1.0 > hap_gaps/diploid/insertions.h1.bed
	cp hap_gaps/diploid/insertions.0.r.bed hap_gaps/diploid/insertions.hom.bed
	cat hap_gaps/diploid/insertions.h0.bed hap_gaps/diploid/insertions.h1.bed hap_gaps/diploid/insertions.hom.bed > $@
#	-rm -f hap_gaps/diploid/insertions.1.r.bed 


hap_gaps/diploid/insertion/L1.bed: hap_gaps/diploid/insertions.bed
	-mkdir hap_gaps/diploid/insertion
	-mkdir hap_gaps/diploid/insertion/rm
	touch hap_gaps/diploid/aln.sam
	touch hap_gaps/diploid/gaps.bed
	# Fake prerequisites for insertions.bed
	cd hap_gaps/diploid && make -t -f $(PBS)/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam insertions.bed
	cd hap_gaps/diploid && make -f $(PBS)/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam insertion/L1.bed

hap_gaps/diploid/deletion/L1.bed: hap_gaps/diploid/deletions.bed
	-mkdir deletion
	-mkdir deletion/rm
	touch hap_gaps/diploid/gaps.bed
	touch hap_gaps/diploid/aln.sam
	# Fake prerequisites for insertions.bed
	cd hap_gaps/diploid && make -t -f $(PBS)/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam deletions.bed
	cd hap_gaps/diploid && make -f $(PBS)/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam deletion/L1.bed

