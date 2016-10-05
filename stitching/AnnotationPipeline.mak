PBS=/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts
BLASR=/net/eichler/vol5/home/mchaisso/software/blasr_2/cpp/alignment/bin/blasr
REF=/net/eichler/vol2/eee_shared/assemblies/GRCh38/GRCh38.fasta
MASKER?=repeatmasker
ALIGNMENTS?=alignments.sam
HGSVC=$(HOME)/projects/HGSVG/hgsvg/
GAPS?=gaps.bed

.PRECIOUS: $(ALIGNMENTS) $(GAPS)

all: $(GAPS) insertions.bed \
		deletions.bed \
	insertion/insertions.fasta\
	deletion/deletions.fasta\
	deletions.labeled.bed \
  insertions.labeled.bed \
  insertions.bb \
  deletions.bb

#	deletion/full/deletions.partial_masked.bed.trf\
#	insertion/full/insertions.partial_masked.bed.trf
#	insertion/rm/insertions.fasta.out \
#	deletion/rm/deletions.fasta.out \
#	deletion/deletions.annotated.bed \
#	insertion/insertions.annotated.bed \
#	insertion/L1.bed \
#	deletion/L1.bed \
# insertion/NONE.bed\
#  deletion/NONE.bed\

help:
	@echo "Usage: make -f AnnotationPipeline.mak alignments=assembly.sam [INPUT_FASTA=<assemblies_file_name>] "


insertions.bed: $(GAPS)
	grep insertion $(GAPS) | awk '{ if ( $$11 < 0.5) print;}' | bedtools sort > $@

deletions.bed: $(GAPS)
	grep deletion $(GAPS) | awk '{ if ( $$11 < 0.5) print;}' | bedtools sort > $@


#
# First step to annotation is to create a fasta file.
#

insertion/insertions.fasta: insertions.bed
	-mkdir -p insertion
	$(PBS)/GapBedToFasta.py --seqidx 5 insertions.bed insertion/insertions.fasta

deletion/deletions.fasta: deletions.bed
	-mkdir -p deletion
	$(PBS)/GapBedToFasta.py --seqidx 5 deletions.bed deletion/deletions.fasta

#
# Next, the repeat masked files are generated.
#

insertion/rm/insertions.fasta.out:  insertion/insertions.fasta
	mkdir -p insertion/rm
	cd insertion && $(HGSVG)/sv/utils/MaskFasta.sh $(MASKER) insertions.fasta rm

deletion/rm/deletions.fasta.out:  deletion/deletions.fasta
	mkdir -p deletion/rm
	cd deletion && $(HGSVG)/sv/utils/MaskFasta.sh $(MASKER) deletions.fasta rm


MIN_MASK=0.70
deletion/deletions.annotated.bed: deletion/rm/deletions.fasta.out
	$(PBS)/AnnotateGapBed.py deletions.bed $@ deletion/rm/deletions.fasta.out deletion/rm/deletions.fasta.masked

#	-grep NONE deletion/deletions.annotated.bed > deletion/NONE.bed
#	grep -v NONE deletion/deletions.annotated.bed >  deletion/deletions.annotated.all.bed
	cat deletion/deletions.annotated.bed | $(PBS)/PrintFractionMasked.py $(MIN_MASK) lt  > deletion/deletions.partial_masked.bed
	cat deletion/deletions.annotated.bed | $(PBS)/PrintFractionMasked.py $(MIN_MASK) gte > deletion/deletions.full_masked.bed


insertion/insertions.annotated.bed: insertion/rm/insertions.fasta.out
	$(PBS)/AnnotateGapBed.py insertions.bed $@ insertion/rm/insertions.fasta.out insertion/rm/insertions.fasta.masked
	cat insertion/insertions.annotated.bed | $(PBS)/PrintFractionMasked.py $(MIN_MASK) lt > insertion/insertions.partial_masked.bed
	cat insertion/insertions.annotated.bed | $(PBS)/PrintFractionMasked.py $(MIN_MASK) gte > insertion/insertions.full_masked.bed

#
# The resulting repeat masked files should have mostly sequence annotated as repeat inside of it.
#

insertion/L1.bed: insertion/insertions.annotated.bed
	$(HGSVG)/sv/utils/PrintUniqueEvents.sh insertion/insertions.full_masked.bed insertion
	rm -rf insertion/[0-9]*
	mkdir -p insertion/full
	mv insertion/ins* insertion/full

# Partially masked sequences may be fully tandem repeats, annotate this now.

insertion/full/insertions.partial_masked.bed.trf: insertion/insertions.annotated.bed
	$(PBS)/RunTRF.sh insertion/full/insertions.partial_masked.bed

deletion/full/deletions.partial_masked.bed.trf: deletion/deletions.annotated.bed
	$(PBS)/RunTRF.sh deletion/full/deletions.partial_masked.bed


# 
# For remaining annotations, first detect all sequences that are mostly tandem repeat
#   
insertion/NONE.bed: insertion/L1.bed insertion/full/insertions.partial_masked.bed.trf
	cat insertion/full/insertions.partial_masked.trf.bed |$(HGSVC)/sv/utils/FilterPartialMasked.py 0.70 gte | cut -f 1-14 > insertion/TRF.bed
	cat insertion/full/insertions.partial_masked.trf.bed | $(HGSVC)/sv/utils/FilterPartialMasked.py 0.70 lt | cut -f 1-14 >  insertion/full/not_trf.bed
	grep NONE insertion/full/not_trf.bed > insertion/NONE.bed
	grep -v NONE insertion/full/not_trf.bed > insertion/partial_annotation.bed

deletion/L1.bed: deletion/deletions.annotated.bed
	$(HGSVG)/sv/utils/PrintUniqueEvents.sh deletion/deletions.full_masked.bed deletion
	rm -rf deletion/[0-9]*
	mkdir -p deletion/full
	mv deletion/del* deletion/full

deletion/NONE.bed: deletion/L1.bed deletion/full/deletions.partial_masked.bed.trf
	cat deletion/full/deletions.partial_masked.trf.bed | $(HGSVC)/sv/utils/FilterPartialMasked.py 0.70 gte  | cut -f 1-14 > deletion/TRF.bed
	cat deletion/full/deletions.partial_masked.trf.bed | $(HGSVC)/sv/utils/FilterPartialMasked.py 0.70 lt   | cut -f 1-14 > deletion/full/not_trf.bed
	grep NONE deletion/full/not_trf.bed > deletion/NONE.bed
	grep -v NONE deletion/full/not_trf.bed > deletion/partial_annotation.bed

insertions.bb: insertions.bed
	cut -f 1-3 $^ | egrep -v "Un|GL|KI" | $(HGSVC)/stitching/TruncateInsertions.py $(REF).fai |  bedtools sort > insertions.lim.bed
	bedToBigBed insertions.lim.bed $(REF).fai $@

deletions.bb: deletions.bed
	cut -f 1-3 $^ | egrep -v "Un|GL|KI" | $(HGSVC)/stitching/TruncateInsertions.py $(REF).fai |  bedtools sort > deletions.lim.bed
	bedToBigBed deletions.lim.bed $(REF).fai $@

insertions.labeled.bed: deletion/NONE.bed

clean:
	rm -rf insertion
	rm -rf deletion
	rm -f insertions.bed
	rm -f deletions.bed
