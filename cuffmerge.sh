#--------------------------------------------------
# for list in eonSpe.*.txt
#--------------------------------------------------
cuffmerge -o cuffmerge -g ./eonSpe/GCF_20170320_genomic.gff -s ./eonSpe/GCF_20170320_genomic.fa -p 40 gtf.list
