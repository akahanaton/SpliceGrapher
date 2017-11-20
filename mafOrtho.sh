#--------------------------------------------------
# awk '$3~/gene/' Homo_sapiens.GRCh38.84.chr.gff3 > Homo_sapiens.GRCh38.84.chr.gff3.gene
#--------------------------------------------------
#--------------------------------------------------
# awk 'OFS="\t"{match($9, "Name=[^;]*;", a); print "chr"$1,$4,$5,a[0]}' Homo_sapiens.GRCh38.84.chr.gff3.gene | sed -e 's/Name=//' -e 's/;$//' | awk '{print $0 > $1".bed" }'
#--------------------------------------------------

rm -rf maker.gff3.MAKER.rename.*matched

for bed in *.bed
do
    echo $bed
    chrName=`basename $bed .bed`
    echo `cut -f4 $bed | sort | uniq -d`
    #--------------------------------------------------
    # ~/software/src/kent/mafsInRegion $bed -outDir $chrName.maf ./eonSpe2.hg38.maf
    #--------------------------------------------------
    qsub -j n -o maker.gff3.MAKER.rename.$chrName.matched -e maker.gff3.MAKER.rename.$chrName.unmatched -N $chrName python mafOrtho.py Homo_sapiens.GRCh38.84.chr.gff3.gene maker.gff3.MAKER.rename $chrName.maf hg38 eonSpe2
done
