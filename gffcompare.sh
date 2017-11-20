#--------------------------------------------------
# input=./IsoSeqAnalysis/gmap_sam/qv.hq.all3.fa.rename.sam.mu.collapsed.gff
# gffread $input -o ${input}3
# gffcompare -r ../../annotation/maker.gff3.MAKER.rename -R -M ${input}3
#--------------------------------------------------

input=./IsoSeqAnalysis/gmap_sam/ES02.qv.hq.all3.fa.rename.sam.mu.collapsed.gff3
gffcompare -r ../../annotation/maker.gff3.MAKER.rename -R -M ${input} -o ES02
input=./IsoSeqAnalysis/gmap_sam/ES04.qv.hq.all3.fa.rename.sam.mu.collapsed.gff3
gffcompare -r ../../annotation/maker.gff3.MAKER.rename -R -M ${input} -o ES04
#--------------------------------------------------
# awk '$3=="transcript"' gffcmp.annotated.gtf  | grep -o -P 'gene_name.*?;' | sort | uniq -c | sort -k1,1nr > gffcmp.annotated.gtf.isoform.count2
#--------------------------------------------------
