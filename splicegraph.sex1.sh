#--------------------------------------------------
# import util/log util/exception util/tryCatch
#--------------------------------------------------

reads_sh_dir=reads_sam_sh
if [[ ! -e all_plot ]]; then mkdir all_plot; fi
if [[ ! -e $reads_sh_dir ]]; then mkdir $reads_sh_dir; fi
if [[ ! -e ${reads_sh_dir}.done ]]; then mkdir ${reads_sh_dir}.done; fi # store the completed sh files

if [[ -e logs ]]; then find logs/ -delete ; fi
mkdir logs

#--------------------------------------------------
# cat ./maker.gff3.MAKER.rename.gene | head -200  |
#--------------------------------------------------

DIRCOUNT=$(ls -1Af $reads_sh_dir |wc -l)
if [ $DIRCOUNT -eq 0  ]; then
    cat ./maker.gff3.MAKER.rename.gene |
    while read line
    do
        idFull=`echo $line | grep -P -o 'ID=.*?;' | sed -e 's/ID=//' -e 's/;//' -e 's/\//-/'`
        chrom=`echo $line|awk '{print $1}'`
        id=`echo $idFull | cut -d ':' -f 1`
        outdir="all_plot/$chrom/$id"

        if [[ ! -e $outdir  ]]; then mkdir -p $outdir; fi

        if [[ ! -s $outdir/graph.gff ]]; then
            grep $idFull maker.gff3.MAKER.rename > $outdir/model.gff
            #--------------------------------------------------
            # ./gff_separate_parent.py $outdir/model.gff.raw > $outdir/model.gff
            #--------------------------------------------------
            gene_model_to_splicegraph.py -m $outdir/model.gff -a -o ${outdir}/graph.gff
        fi

        #--------------------------------------------------
        # if [[ ! -s $outdir/pacbio.gff.male.full ]]; then
        #     if [[ -e $outdir/pacbio.gff.raw ]]; then rm $outdir/pacbio.gff.raw; fi
        #     for pacbioID in `grep $id ../pacbio/gffcmp.annotated.gtf | grep $chrom | cut -f9 | grep -o -P 'gene_id.*?;' | sed -e 's/gene_id //' -e 's/"//g' -e 's/;//' | sort | uniq`
        #     do
        #         grep "${pacbioID}\." ./csus.all3.fa.rename.sam.allmapping.collapsed.gff3 >> $outdir/pacbio.gff.raw
        #     done
        #     if [[ -s $outdir/pacbio.gff.raw ]]; then
        #         if [[ -e $outdir/pacbio.gff.male.tmp ]]; then rm $outdir/pacbio.gff.male.tmp; fi
        #         if [[ -e $outdir/pacbio.gff.female.tmp ]]; then rm $outdir/pacbio.gff.female.tmp; fi
        #         for transcriptID in `awk '$3=="transcript"' $outdir/pacbio.gff.raw | grep -o -P 'ID=.*?;' | sed -e 's/ID=//' -e 's/;//'`
        #         do
        #             gender=`grep -F "$transcriptID" ./IsoSeqAnalysis/gmap_sam/csus.all3.fa.rename.sam.allmapping.collapsed.group.txt.2rec.isoform.biased | awk '{if($3>0) {print $2} else{print $4}}'` 
        #             if [[ $gender == "ES02" ]]; then
        #                 grep -F $transcriptID $outdir/pacbio.gff.raw >> $outdir/pacbio.gff.male.tmp
        #             elif [[ $gender == "ES04" ]]; then
        #                 grep -F $transcriptID $outdir/pacbio.gff.raw >> $outdir/pacbio.gff.female.tmp
        #             fi
        #         done
        #         geneid=`head -1 $outdir/pacbio.gff.raw | awk '{print $9}' | grep -o -P 'geneID=.*' | sed 's/geneID/ID/'`
        #         strand=`head -1 $outdir/pacbio.gff.raw | awk '{print $7}'`
        #         head -1 $outdir/model.gff | awk -v geneid=$geneid -v str=$strand 'OFS="\t"{print $1,"Pacbio",$3,$4,$5,$6,str,$8,geneid}' > $outdir/pacbio.gff.full
        #         if [[ -s $outdir/pacbio.gff.male.tmp ]]; then
        #             cat $outdir/pacbio.gff.full $outdir/pacbio.gff.male.tmp > $outdir/pacbio.gff.male.full
        #             gene_model_to_splicegraph.py -m $outdir/pacbio.gff.male.full -a -o $outdir/pacbio.gff.male
        #         fi
        #         if [[ -s $outdir/pacbio.gff.female.tmp ]]; then
        #             cat $outdir/pacbio.gff.full $outdir/pacbio.gff.female.tmp > $outdir/pacbio.gff.female.full
        #             gene_model_to_splicegraph.py -m $outdir/pacbio.gff.female.full -a -o $outdir/pacbio.gff.female
        #         fi
        #         cat $outdir/pacbio.gff.raw >> $outdir/pacbio.gff.full
        #         gene_model_to_splicegraph.py -m $outdir/pacbio.gff.full -a -o $outdir/pacbio.gff
        #     fi
        # fi
        #--------------------------------------------------


        rm $outdir/*sam*
        echo "#!/bin/bash" > $reads_sh_dir/$id.sh
        for dir in `ls -d [[:upper:]]*_sam`
        do
            begin=`awk '$3=="gene" {print $4}' $outdir/model.gff`
            end=`awk '$3=="gene" {print $5}' $outdir/model.gff`
            echo  "awk '\$4>=$begin && \$4 <= $end' $dir/$chrom.sam  | sort -k 3,3 -k 4,4n > $outdir/`echo $dir | sed 's/_/./g'`" >> $reads_sh_dir/$id.sh
        done

        #--------------------------------------------------
        # pred_params="-d ${outdir}/reads.sam -s ${outdir}/ests.gff -o ${outdir}/predicted.gff -J 1 -M 1 -T 20 -v"
        # predict_splicegraph.py ${outdir}/graph.gff ${pred_params}
        #--------------------------------------------------
    done
else
    sed "s/DirName/$reads_sh_dir/" ~/software/myProgram/qsub.bash.snake > reads.sam.snake
    snakemake -q -j 100 -k --cluster "qsub -o ./logs -j y" -s reads.sam.snake
fi
