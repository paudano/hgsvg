PBS=/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts
HGSVG=/net/eichler/vol5/home/mchaisso/projects/HGSVG/hgsvg
BLASR=/net/eichler/vol5/home/mchaisso/software/blasr_2/cpp/alignment/bin/blasr
REF=/net/eichler/vol2/eee_shared/assemblies/GRCh38/GRCh38.fasta

H0SAM?=$(PWD)/alignments.h0.sam
H1SAM?=$(PWD)/alignments.h1.sam
DIR?=hap_gaps
NFIELDS?=10
include $(PARAMS)
all:$(DIR)/hap0/gaps.bed \
  $(DIR)/hap1/gaps.bed \
  $(DIR)/hap0/insertions.bed \
  $(DIR)/hap1/insertions.bed \
  $(DIR)/diploid/insertions.bed \
  $(DIR)/diploid/deletions.bed \
  $(PWD)/alignments.h0.bam \
  $(PWD)/alignments.h1.bam \
  $(DIR)/diploid/deletions.annotated.bed \
  $(DIR)/diploid/insertions.annotated.bed \
  $(DIR)/diploid/sv_calls.vcf \
  $(DIR)/diploid/deletions.hap.bed \
  $(DIR)/diploid/insertions.hap.bed \
  $(DIR)/diploid/deletions.hap.bb \
  $(DIR)/diploid/insertions.hap.bb \
  $(DIR)/diploid/deletion/NONE.bed \
  $(DIR)/diploid/insertion/NONE.bed 


samfies.fofn: samfiles/chr22.18953772-19013772.sam
	ls samfiles/*| awk '{ print "samfiles/"$1; }' > $@


$(PWD)/alignments.h0.sam: samfiles.fofn
	$(HGSVG)/sv/CombineAssemblies.py --alignments samfiles.fofn --header $(HGSVG)/sv/header.sam

$(PWD)/alignments.h0.bam: $(H0SAM)
	-mkdir -p $TMPDIR/mchaisso
	samtools view -bS $(H0SAM) | samtools sort -T $TMPDIR/aln.0.$$PPID -o $@
	samtools index $@

$(PWD)/alignments.h1.bam: $(H0SAM)
	-mkdir -p /var/tmp/mchaisso
	samtools view -bS $(H1SAM) | samtools sort -T $TMPDIR/mchaisso/aln.1.$$PPID -o $@
	samtools index $@

$(DIR)/hap0/gaps.bed: $(PWD)/alignments.h0.bam
	-mkdir -p $(DIR)/hap0
	samtools view $^ | $(PBS)/PrintGaps.py $(REF) /dev/stdin --minAlignmentLength 30000 --ignoreHP 3 --minDist 1000 --condense 20 --outFile $@

$(DIR)/hap1/gaps.bed: $(PWD)/alignments.h1.bam
	mkdir -p $(DIR)/hap1
	samtools view $^ | $(PBS)/PrintGaps.py $(REF) /dev/stdin --minAlignmentLength 30000  --ignoreHP 3 --minDist 1000 --condense 20 --outFile $@

$(DIR)/hap0/insertions.bed: $(DIR)/hap0/gaps.bed $(PWD)/alignments.h0.bam
	cd $(DIR)/hap0 && make -f $(HGSVG)/sv/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=$(PWD)/alignments.h0.bam

$(DIR)/hap1/insertions.bed: $(DIR)/hap1/gaps.bed $(PWD)/alignments.h1.bam
	cd $(DIR)/hap1 && make -f $(HGSVG)/sv/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=$(PWD)/alignments.h1.bam

$(DIR)/hap1/deletions.bed: $(DIR)/hap1/gaps.bed  $(PWD)/alignments.h1.bam
	cd $(DIR)/hap1 && make -f $(HGSVG)/sv/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=$(PWD)/alignments.h1.bam deletions.bed

$(DIR)/hap0/deletions.bed: $(DIR)/hap0/gaps.bed  $(PWD)/alignments.h0.bam
	cd $(DIR)/hap0 && make -f $(HGSVG)/sv/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=$(PWD)/alignments.h0.bam deletions.bed

$(DIR)/diploid/deletions.bed: $(DIR)/hap0/deletions.bed $(DIR)/hap1/deletions.bed
	mkdir -p $(DIR)/diploid
	$(HGSVG)/sv/utils/MergeHaplotypes.sh $(DIR)/hap0/deletions.bed $(DIR)/hap1/deletions.bed $@ "svType svLen svSeq qName qStart qEnd"

$(DIR)/diploid/insertions.bed: $(DIR)/hap0/insertions.bed $(DIR)/hap1/insertions.bed
	mkdir -p $(DIR)/diploid
	$(HGSVG)/sv/utils/MergeHaplotypes.sh $(DIR)/hap0/insertions.bed $(DIR)/hap1/insertions.bed $@ "svType svLen svSeq qName qStart qEnd"

$(DIR)/diploid/insertion/L1.bed: $(DIR)/diploid/insertions.bed
	-mkdir $(DIR)/diploid/insertion
	-mkdir $(DIR)/diploid/insertion/rm
	touch $(DIR)/diploid/aln.sam
	touch $(DIR)/diploid/gaps.bed
	cd $(DIR)/diploid && make -t -f $(PBS)/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam insertions.bed
	cd $(DIR)/diploid && make -f $(PBS)/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam insertion/L1.bed

$(DIR)/diploid/deletion/L1.bed: $(DIR)/diploid/deletions.bed
	-mkdir deletion
	-mkdir deletion/rm
	touch $(DIR)/diploid/gaps.bed
	touch $(DIR)/diploid/aln.sam
	cd $(DIR)/diploid && make -t -f $(PBS)/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam deletions.bed
	cd $(DIR)/diploid && make -f $(PBS)/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam deletion/L1.bed



$(DIR)/diploid/insertion/NONE.bed: $(DIR)/diploid/insertions.bed
	-mkdir $(DIR)/diploid/insertion
	-mkdir $(DIR)/diploid/insertion/rm
	touch $(DIR)/diploid/aln.sam
	touch $(DIR)/diploid/gaps.bed
	# Fake prerequisites for insertions.bed
	cd $(DIR)/diploid && make -t -f $(HGSVG)/stitching/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam insertions.bed
#	cd $(DIR)/diploid && make -f $(HGSVG)/stitching/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam insertion/L1.bed
	cd $(DIR)/diploid && make -f $(HGSVG)/stitching/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam insertion/NONE.bed

$(DIR)/diploid/deletion/NONE.bed: $(DIR)/diploid/deletions.bed
	-mkdir deletion
	-mkdir deletion/rm
	touch $(DIR)/diploid/gaps.bed
	touch $(DIR)/diploid/aln.sam
	# Fake prerequisites for insertions.bed
	cd $(DIR)/diploid && make -t -f $(HGSVG)/stitching/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam deletions.bed
#	cd $(DIR)/diploid && make -f $(HGSVG)/stitching/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam deletion/L1.bed
	cd $(DIR)/diploid && make -f $(HGSVG)/stitching/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam deletion/NONE.bed


$(DIR)/diploid/insertions.hap.bb: $(DIR)/diploid/insertions.hap.bed
	module load ucsc && bedToBigBed $^ $(REF).fai $@ -type=bed9

$(DIR)/diploid/deletions.hap.bb: $(DIR)/diploid/deletions.hap.bed
	module load ucsc && bedToBigBed $^ $(REF).fai $@ -type=bed9

$(DIR)/diploid/insertions.bb: $(DIR)/diploid/insertions.bed
	cut -f 1-3 $^ | bedtools sort  >  "$^".sorted
	module load ucsc; bedToBigBed "$^".sorted $(REF).fai $@

$(DIR)/diploid/deletions.bb: $(DIR)/diploid/deletions.bed
	cut -f 1-3 $^ | bedtools sort  > "$^".sorted
	module load ucsc; bedToBigBed "$^".sorted $(REF).fai $@

$(DIR)/diploid/deletions.annotated.bed: $(DIR)/diploid/deletion/NONE.bed
	$(HGSVG)/sv/utils/CombineAnnotatedEvents.sh $@ $(DIR)/diploid/deletion/

$(DIR)/diploid/insertions.annotated.bed: $(DIR)/diploid/insertion/NONE.bed
	$(HGSVG)/sv/utils/CombineAnnotatedEvents.sh $@ $(DIR)/diploid/insertion/

$(DIR)/diploid/sv_calls.bed: $(DIR)/diploid/insertions.annotated.bed $(DIR)/diploid/deletions.annotated.bed
	$(HGSVG)/sv/utils/MergeFiles.py --files $(DIR)/diploid/insertions.annotated.bed $(DIR)/diploid/deletions.annotated.bed | \
  bedtools sort -header > $@

$(DIR)/diploid/sv_calls.vcf:   $(DIR)/diploid/sv_calls.bed
	module load pandas && $(HGSVG)/sv/utils/variants_bed_to_vcf.py --bed $(DIR)/diploid/sv_calls.bed --ref $(REF) --sample $(SAMPLE) --type sv --vcf $@ 

$(DIR)/diploid/insertions.hap.bed: $(DIR)/diploid/insertions.annotated.bed
	awk '{ col="0,0,255"; if ($$13 == "HAP0") {col="255,0,0";} if ($$13 == "HAP1") {col="0,255,0";} print $$1"\t"$$2"\t"$$3"\t"$$14"\t1000\t+\t"$$2"\t"$$3"\t"col}' <  $^ | sort -k1,1 -k2,2n > $@

$(DIR)/diploid/deletions.hap.bed: $(DIR)/diploid/deletions.annotated.bed
	awk '{ col="0,0,255"; if ($$13 == "HAP0") {col="255,0,0";} if ($$13 == "HAP1") {col="0,255,0";} print $$1"\t"$$2"\t"$$3"\t"$$14"\t1000\t+\t"$$2"\t"$$3"\t"col}' <  $^ | sort -k1,1 -k2,2n > $@

