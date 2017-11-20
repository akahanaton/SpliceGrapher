plot_sh_dir=plot_sh_dir
plot_out_dir=sex_biased_pdf
modelDir=all_plot.20171111
if [[ ! -e $plot_sh_dir ]]; then mkdir $plot_sh_dir; fi
if [[ ! -e ${plot_sh_dir}.done ]]; then mkdir ${plot_sh_dir}.done; fi
if [[ ! -e $plot_out_dir ]]; then mkdir $plot_out_dir; fi
if [[ -e logs ]]; then find logs/ -delete ; fi
mkdir logs

cat ./maker.gff3.MAKER.rename.clean.gene |
while read line
do
    break
    idFull=`echo $line | grep -P -o 'ID=.*?;' | sed -e 's/ID=//' -e 's/;//' -e 's/\//-/'`
    chrom=`echo $line|awk '{print $1}'`
    id=`echo $idFull | cut -d ':' -f 1`
    outdir="$modelDir/$chrom/$id"

    view_params="-H 30 -W 15 -J 2 -c -L -o $outdir/$id.pdf"


    pteGff=pteAle/GCF_000325575.1_ASM32557v1_genomic.gff.gene.chr.gff.gene

    idCombine=$idFull
    model="$outdir/model.gff"
    idFull=`head -1 $model | grep -P -o 'ID=.*?;' | sed -e 's/ID=//' -e 's/;//' -e 's/\//-/'`
    graph="1,$outdir/graph.gff"
    depth=`find $outdir/ -name "*.sam" -size +0 | paste -sd ',' | sed -e 's/,/,1,/g' -e 's/^/1,/'`
    if [[ $depth == "1," ]]; then
        depth=""
    fi
    splice=`find $outdir/ -name "pacbio.*" -size +0 | grep -E -v 'tmp|raw|full' | paste -sd ',' | sed -e 's/,/,1,/' -e 's/^/1,/'`
    pteLine=`grep -F "=${id};" $pteGff|head -1`
    if [[ $pteLine != "" ]]; then
        pteIDPlot=`echo $pteLine | grep -P -o 'ID=.*?;' | sed -e 's/ID=//' -e 's/;//' -e 's/\//-/'`
        pteID=`echo $pteLine | grep -P -o 'Name=.*?;' | sed -e 's/Name=//' -e 's/;//' -e 's/\//-/'`
        pteChrom=`echo $pteLine|awk '{print $1}'`
        idCombine="$idCombine,$pteIDPlot"
        model=$model",pteAle/all_plot/$pteChrom/$id/model.gff"
        graph=$graph",2,pteAle/all_plot/$pteChrom/$id/graph.gff"
        depth=$depth`find pteAle/all_plot/$pteChrom/$id/ -name "sam.*" -size +0 | paste -sd ',' | sed -e 's/,/,2,/g' -e 's/^/,2,/'`
    fi

    #--------------------------------------------------
    # if [[ -e ${outdir}/$id.pdf ]]; then rm $outdir/$id.pdf; fi
    #--------------------------------------------------

    if [[ $splice != "1," ]]; then
        view_params=$view_params" -s $splice"
    fi

    if [[ $depth == *"2," ]]; then
        if [[ $depth == "1,"* ]]; then
            depth=`echo $depth | sed 's/2,//'`
            echo "view_splicegraph_multiplot.py $idCombine ${view_params} -m $model -G $graph -d $depth" > $plot_sh_dir/$id.sh
        else
            echo "view_splicegraph_multiplot.py $idCombine ${view_params} -m $model -G $graph" > $plot_sh_dir/$id.sh
        fi
    else
        echo "view_splicegraph_multiplot.py $idCombine ${view_params} -m $model -G $graph -d $depth" > $plot_sh_dir/$id.sh
    fi
    if [[ $splice == *"male"* ]]; then
        echo "cp ${outdir}/$id.pdf $plot_out_dir" >> $plot_sh_dir/$id.sh
    fi

    #--------------------------------------------------
    # view_splicegraph_multiplot.py 'CCDC50:maker-000000F-snap-gene-9.85,gene:ENSG00000152492' -m 'model.gff,/wlflab/Ming/genome_assembly/EonycterisSpelaea/results/alternativesplice/hg38/all_plot/chr3/CCDC50/model.gff' -G '1,graph.gff,2,/wlflab/Ming/genome_assembly/EonycterisSpelaea/results/alternativesplice/hg38/all_plot/chr3/CCDC50/graph.gff' -H 30 -W 15 -J 2 -c -L -d '1,male.reads.sam,1,female.reads.sam,2,/wlflab/Ming/genome_assembly/EonycterisSpelaea/results/alternativesplice/hg38/all_plot/chr3/CCDC50/reads.sam' -o CCDC50.pdf -s '1,pacbio.gff,1,pacbio.gff.male,1,pacbio.gff.female,2,/wlflab/Ming/genome_assembly/EonycterisSpelaea/results/alternativesplice/hg38/all_plot/chr3/CCDC50/graph.gff'
    #--------------------------------------------------

done

sed "s/DirName/$plot_sh_dir/" ~/software/myProgram/qsub.bash.snake > plot.snake
snakemake -q -j 100 -k --cluster "qsub -o ./logs -j y" -s plot.snake
