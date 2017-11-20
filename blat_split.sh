#!/usr/bin/env bash

if [ $# -lt 1  ]; then
    echo "Usage : $0 input.fa genome.fa/cut.dir chunks.num output.prefix"
    exit;
fi

inputFa=$1
dbFile=$2
chunksNum=$3
prefix=$4
inputdir=`dirname $inputFa`

#--------------------------------------------------
# source "/home/gmswenm/software/bash-oo-framework/lib/oo-bootstrap.sh"
# import util/log util/exception util/tryCatch
#--------------------------------------------------

if [ -e logs ]; then
    find logs -delete
fi
mkdir logs

if [ ! -e blat_psl ]; then
    mkdir blat_psl
fi

#--------------------------------------------------
# try{
#--------------------------------------------------
    if [[ ! -e $prefix.psl ]] || [[ ! -s $prefix.psl ]]; then
        if [[ ! -e splite.done ]]; then
            ~/wlflab/software/src/maker/bin/fasta_tool --chunks ${chunksNum} ${inputFa}
            touch splite.done
        fi
        if [[ -d $dbFile ]]; then
            echo $dbFile
        else
            for part in `ls $inputdir/*.fasta | grep part`
            do
                echo $part
                fa=`basename $part`
                qsub -N blat -o logs blat -t=dna -q=dna -noHead -minScore=15 -minIdentity=50 -extendThroughN ${dbFile} $part ./blat_psl/$prefix.$fa.psl
            done
        fi
    fi
#--------------------------------------------------
# } catch {
#     echo "Caught Exception:$(UI.Color.Red) $__BACKTRACE_COMMAND__ $(UI.Color.Default)"
#     echo "File: $__BACKTRACE_SOURCE__, Line: $__BACKTRACE_LINE__"
#     ## printing a caught exception couldn't be simpler, as it's stored in "${__EXCEPTION__[@]}"
#     Exception::PrintException "${__EXCEPTION__[@]}"
# }
#--------------------------------------------------

while [ `qstat -u $USER | grep 'blat' | wc -l` -gt 0  ]
do
    sleep 60
done

if [[ ! -e blat.done ]]; then
    cat ./blat_psl/$prefix.*part*.psl > $prefix.psl
    touch blat.done
fi

#--------------------------------------------------
# /home/gmswenm/wlflab/software/src/scipio-1.4/scipio.1.4.1.pl --blat_output=all.pep.psl $dbFile $inputFa > all.pep.yaml
# cat all.pep.yaml | /home/gmswenm/wlflab/software/src/scipio-1.4/yaml2gff.1.4.pl > all.pep.scipiogff
# /home/gmswenm/wlflab/software/src/augustus.2.5.5/scripts/scipiogff2gff.pl --in=all.pep.scipiogff -out=all.pep.gff
# /home/gmswenm/wlflab/software/src/augustus.2.5.5/scripts/gff2gbSmallDNA.pl all.pep.gff $dbFile 1000 all.pep.raw.gb
# /home/gmswenm/wlflab/software/src/augustus.2.5.5/bin/etraining --species=bat --stopCodonExcludedFromCDS=true all.pep.raw.gb 2>all.pep.train.err
# cat all.pep.train.err | perl -pe 's/.*in sequence (\S+): .*/$1/' > all.pep.badgenes.lst
# /home/gmswenm/wlflab/software/src/augustus.2.5.5/scripts/filterGenes.pl all.pep.badgenes.lst all.pep.raw.gb > all.pep.gb
#--------------------------------------------------

#--------------------------------------------------
# rule scipio:
#     input: allPsl='all.pep.psl', dbFile=rules.blat.input.dbFile, inputFa=rules.fasta_tool.input.faFile
#     output: yaml='all.pep.yaml'
#
# rule scipiogff:
#     input: yaml='all.pep.yaml'
#     output: scipiogff='all.pep.scipiogff'
#
# rule gff:
#     input: scipiogff='all.pep.scipiogff'
#     output: gff='all.pep.gff'
#
# rule rawGb:
#     input: gff='all.pep.gff',dbFile='batsGenome/eonSpe/p_ctg.fa'
#     output: rawGb='all.pep.raw.gb'
#
# rule etraining:
#     input: rawGb='all.pep.raw.gb'
#     output: trainErr='all.pep.train.err'
#
# rule badgenes:
#     input: trainErr='all.pep.train.err'
#     output: badList='all.pep.badgenes.lst'
#
# rule gb:
#     input: badList='all.pep.badgenes.lst',rawGb='all.pep.raw.gb'
#     output: gb='all.pep.gb'
#--------------------------------------------------
