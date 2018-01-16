PBS=/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts
BLASR=/net/eichler/vol5/home/mchaisso/software/blasr_2/cpp/alignment/bin/blasr
REF=/net/eichler/vol2/eee_shared/assemblies/GRCh38/GRCh38.fasta
HGSVG=/net/eichler/vol5/home/mchaisso/projects/HGSVG/hgsvg
MASKER?=repeatmasker
ALIGNMENTS?=alignments.sam

GAPS?=gaps.bed

.PRECIOUS: $(ALIGNMENTS) $(GAPS)

all: $(GAPS) \
   $(ALIGNMENTS).bed \
  insertions.bed \
	deletions.bed \
	insertion/insertions.fasta\
	deletion/deletions.fasta\
	insertion/rm/insertions.fasta.out \
	deletion/rm/deletions.fasta.out \
	deletion/deletions.annotated.bed \
	insertion/insertions.annotated.bed \
  tiling_contigs.tab \
  insertions.all.bed \
  deletions.all.bed \
  insertions.bb \
  deletions.bb \
	insertion/L1.bed \
	deletion/L1.bed \
  insertion/NONE.bed\
  deletion/NONE.bed\
	deletion/deletions.partial_masked.bed.trf\
	insertion/insertions.partial_masked.bed.trf

help:
	@echo "Usage: make -f AnnotationPipeline.mak alignments=assembly.sam [INPUT_FASTA=<assemblies_file_name>] "


#
# Next print the gapped sequences, split into insertion and deletion, and simplify redundant (duplicate) gaps.
#
#$(GAPS): $(alignments)
#	$(PBS)/PrintGaps.py $(REF)  $(alignments) --condense 20 --tsd 20 --outFile $(GAPS) --minLength 50  --qpos 
#

#$(GAPS): $(GOR_GAPS)
#	cat $(GOR_GAPS) | awk 'BEGIN {OFS="\t";} { print $$4,$$5,$$6,$$7,$$8,$$9,$$10,$$11,$$1,$$2,$$3;}' > $@
#

$(ALIGNMENTS).bed: $(ALIGNMENTS)
	samtools view -h $(ALIGNMENTS) | samToBed /dev/stdin  --reportIdentity | bedtools sort > $@

tiling_contigs.tab: $(ALIGNMENTS).bed
	./TilingPath.py $^ $@

insertions.all.bed: $(GAPS)
	egrep "^#|insertion" $(GAPS) | bedtools sort -header  > $@

deletions.all.bed: $(GAPS)
	egrep "^#|deletion" $(GAPS) | bedtools sort -header > $@

insertions.bed: insertions.all.bed tiling_contigs.tab
	$(PBS)/SmallIndelAnalysis/FilterGapsByTilingPath.py insertions.all.bed tiling_contigs.tab --out $@

deletions.bed: deletions.all.bed tiling_contigs.tab
	$(PBS)/SmallIndelAnalysis/FilterGapsByTilingPath.py deletions.all.bed tiling_contigs.tab --out $@




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
	-mkdir -p insertion/rm
	cd insertion && $(HGSVG)/sv/utils/MaskFasta.sh $(MASKER) insertions.fasta rm

deletion/rm/deletions.fasta.out:  deletion/deletions.fasta
	-mkdir -p deletion/rm
	cd deletion && $(HGSVG)/sv/utils/MaskFasta.sh $(MASKER) deletions.fasta rm


MIN_MASK=0.70
deletion/deletions.annotated.bed: deletion/rm/deletions.fasta.out
	$(PBS)/AnnotateGapBed.py deletions.bed $@ deletion/rm/deletions.fasta.out deletion/rm/deletions.fasta.masked


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
#	mv insertion/ins* insertion/full

# Partially masked sequences may be fully tandem repeats, annotate this now.

insertion/insertions.partial_masked.bed.trf: insertion/insertions.annotated.bed
	$(PBS)/RunTRF.sh insertion/insertions.partial_masked.bed

deletion/deletions.partial_masked.bed.trf: deletion/deletions.annotated.bed
	$(PBS)/RunTRF.sh deletion/deletions.partial_masked.bed


# 
# For remaining annotations, first detect all sequences that are mostly tandem repeat
#   
insertion/NONE.bed: insertion/insertions.partial_masked.bed.trf
	cat insertion/insertions.partial_masked.trf.bed | awk '{ if ($$14 >= 0.70) print;}' > insertion/TRF.bed
	cat insertion/insertions.partial_masked.trf.bed | awk '{ if ($$14 < 0.70) print;}' > insertion/not_trf.bed
	grep NONE insertion/not_trf.bed > insertion/NONE.bed
	grep -v NONE insertion/not_trf.bed > insertion/partial_annotation.bed

deletion/L1.bed: deletion/deletions.annotated.bed
	$(HGSVG)/sv/utils/PrintUniqueEvents.sh deletion/deletions.full_masked.bed deletion
	rm -rf deletion/[0-9]*
	mkdir -p deletion


deletion/NONE.bed: deletion/deletions.partial_masked.bed.trf
	cat deletion/deletions.partial_masked.trf.bed | awk '{ if ($$14 >= 0.70) print;}' > deletion/TRF.bed
	cat deletion/deletions.partial_masked.trf.bed | awk '{ if ($$14 < 0.70) print;}' > deletion/not_trf.bed
	grep NONE deletion/not_trf.bed > deletion/NONE.bed
	grep -v NONE deletion/not_trf.bed > deletion/partial_annotation.bed

insertions.bb: insertions.bed
	cut -f 1-3 $^ | egrep -v "Un|GL|KI" | bedtools sort > insertions.lim.bed
	bedToBigBed insertions.lim.bed $(REF).fai $@

deletions.bb: deletions.bed
	cut -f 1-3 $^ | egrep -v "Un|GL|KI" | bedtools sort > deletions.lim.bed
	bedToBigBed deletions.lim.bed $(REF).fai $@



clean:
	rm -rf insertion
	rm -rf deletion
	rm -f insertions.bed
	rm -f deletions.bed
