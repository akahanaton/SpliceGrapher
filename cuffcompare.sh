#!/usr/bin/env bash
set –euo pipefail
#--------------------------------------------------
# cuffcompare -G -r ./eonSpe/GCF_20170320_genomic.gtf -s ./eonSpe/GCF_20170320_genomic.fa -V ./IsoSeqAnalysis/gmap_sam/ES02.qv.hq.all3.fa.rename.sam.mu.collapsed.gff ./IsoSeqAnalysis/gmap_sam/ES04.qv.hq.all3.fa.rename.sam.mu.collapsed.gff -o ref
#--------------------------------------------------
cuffcompare -G -M -N -r ./IsoSeqAnalysis/gmap_sam/ES02.qv.hq.all3.fa.rename.sam.mu.collapsed.gff -s ./eonSpe/GCF_20170320_genomic.fa -V  ./IsoSeqAnalysis/gmap_sam/ES04.qv.hq.all3.fa.rename.sam.mu.collapsed.gff -o ES0204
cuffcompare -G -M -N -r ./IsoSeqAnalysis/gmap_sam/ES04.qv.hq.all3.fa.rename.sam.mu.collapsed.gff -s ./eonSpe/GCF_20170320_genomic.fa -V  ./IsoSeqAnalysis/gmap_sam/ES02.qv.hq.all3.fa.rename.sam.mu.collapsed.gff -o ES0402
