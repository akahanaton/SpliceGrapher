#--------------------------------------------------
# import util/log util/exception util/tryCatch
#--------------------------------------------------
#--------------------------------------------------
# import util/type #data structure
#--------------------------------------------------

#--------------------------------------------------
# ./gff_split_per_chr_per_gene.py Homo_sapiens.GRCh38.84.chr.gff3 Homo_sapiens.GRCh38.84.chr.gff3.split
#--------------------------------------------------

plot_root='all_plot'
if [[ ! -e $plot_root ]]; then mkdir $plot_root; fi
if [[ ! -e sam_file ]]; then mkdir sam_file; fi
if [[ ! -e reads_sam_sh ]]; then mkdir reads_sam_sh; fi

#--------------------------------------------------
# cat maker.gff3.MAKER.rename.chr*.matched |
#--------------------------------------------------
for file in maker.gff3.MAKER.rename.chr*.matched
do
    chrom=`echo $file| cut -d '.' -f5`
    chromNum=`echo $chrom | sed -e 's/://' -e 's/chr//'`
    if [[ ! -s sam_file/$chrom.sh ]]; then
        echo "awk '\$3==\"$chromNum\"' tophat2/tophat2.JW_HEK293T_Mock_0h/accepted_hits.sam tophat2/tophat2.JW_HEK293T_HeV_24h/accepted_hits.sam > sam_file/$chrom.sam" > sam_file/$chrom.sh
        qsub -o /dev/null -N awk bash sam_file/$chrom.sh
    fi
done

while [[ `qstat -u gmswenm | grep awk | wc -l` -gt 0  ]]
do
    echo `qstat -u gmswenm | grep awk | wc -l`
    sleep 30
done

#--------------------------------------------------
# head -1 maker.gff3.MAKER.rename.chr*.matched | grep -v '=' | grep -v '^$' |
#--------------------------------------------------
cat ./Homo_sapiens.GRCh38.84.chr.gff3.gene.pro |
while read line
do
    chrom=`echo $line|awk '{print "chr"$1}'`
    chromNum=`echo $line|awk '{print $1}'`
    id=`echo $line | grep -P -o 'Name=.*?;' | sed -e 's/Name=//' -e 's/;//' -e 's/\//-/'`
    outdir="$plot_root/$chrom/$id"

    echo $id

    if [[ ! -e $outdir ]]; then mkdir -p $outdir ; fi

    if [[ ! -s $outdir/graph.gff ]]; then
        awk -F '\t' 'OFS="\t"{if($3 ~ /transcript/){$3="mRNA"}; $2="human"; print}' ./Homo_sapiens.GRCh38.84.chr.gff3.split/$chromNum/$id.gff > $outdir/model.gff
        gene_model_to_splicegraph.py -m $outdir/model.gff -a -o ${outdir}/graph.gff
    fi

    #--------------------------------------------------
    # begin=`awk '$3=="gene" {print $4}' $outdir/model.gff`
    # end=`awk '$3=="gene" {print $5}' $outdir/model.gff`
    # awk -v begin=$begin -v end=$end '$16>=begin && $17 <= end' 000090F.psl  | sort -k 10,10 -k 16,16n | awk '$1 > 500' > $outdir/ests.psl
    #--------------------------------------------------
    #--------------------------------------------------
    # ests_to_splicegraph.py $outdir/ests.psl -i 40 -m ${outdir}/model.gff
    #--------------------------------------------------
    # mv ${outdir}*_ests.gff $outdir/ests.gff
    #--------------------------------------------------

    if [[ ! -e $outdir/reads.sam ]]; then
        begin=`awk '$3=="gene" {print $4}' $outdir/model.gff`
        end=`awk '$3=="gene" {print $5}' $outdir/model.gff`
        echo  "awk '\$4>=$begin && \$4 <= $end' sam_file/$chrom.sam  | sort -k 3,3 -k 4,4n > $outdir/reads.sam" > reads_sam_sh/$id.sh
    fi

    #--------------------------------------------------
    # pred_params="-d ${outdir}/reads.sam -s ${outdir}/ests.gff -o ${outdir}/predicted.gff -J 1 -M 1 -T 20 -v"
    # predict_splicegraph.py ${outdir}/graph.gff ${pred_params}
    #--------------------------------------------------
done

sed "s/DirName/reads_sam_sh/" ~/software/myProgram/qsub.bash.snake > reads.sam.snake
snakemake -q -j 20 -k --cluster "qsub -o ./logs -j y" -s reads.sam.snake

while [[ `qstat -u gmswenm | grep awk | wc -l` -gt 0  ]]
do
    echo `qstat -u gmswenm | grep awk | wc -l`
    sleep 30
done

#--------------------------------------------------
# head -1 maker.gff3.MAKER.rename.chr*.matched | grep -v '=' | grep -v '^$' |
#--------------------------------------------------
cat ./Homo_sapiens.GRCh38.84.chr.gff3.gene.pro |
while read line
do
    chrom=`echo $line|awk '{print "chr"$1}'`
    chromNum=`echo $line|awk '{print $1}'`
    id=`echo $line | grep -P -o 'Name=.*?;' | sed -e 's/Name=//' -e 's/;//' -e 's/\//-/'`
    outdir="$plot_root/$chrom/$id"
    idPlot=`echo $line | grep -P -o 'ID=.*?;' | sed -e 's/ID=//' -e 's/;//' -e 's/\//-/'`
    echo $idPlot

    #--------------------------------------------------
    # view_params="-d ${id}/reads.sam -m ${id}/model.gff -o ${id}.pdf -s ${id}/predicted.gff -G ${id}/graph.gff -J 2 -c -L"
    #--------------------------------------------------
    view_params=''
    if [[ -e $outdir/pacbio.gff ]]; then
        view_params="-d ${outdir}/reads.sam -m ${outdir}/model.gff -o ${outdir}/$id.pdf -s ${outdir}/pacbio.gff -G ${outdir}/graph.gff -J 2 -c -L"
    else
        view_params="-d ${outdir}/reads.sam -m ${outdir}/model.gff -o ${outdir}/$id.pdf -G ${outdir}/graph.gff -J 2 -c -L"
    fi

    if [[ ! -s $outdir/$id.pdf ]]; then
        if [[ -s $outdir/model.gff && -s $outdir/graph.gff ]]; then
            view_splicegraph_multiplot.py ${idPlot} ${view_params} -H 7 -W 10 -J 1 -c -L
        fi
    fi
done
