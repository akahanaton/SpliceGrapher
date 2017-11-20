#!/usr/bin/env bash
set –euo pipefail
for sam in ./tophat2/tophat2.*ES-0*/accepted_hits.sam
do
    echo $sam
    gender=""
    if  `echo $sam | grep -q 'ES-02'` ; then
        gender="male"
    else
        gender="female"
    fi
    echo $gender
    tissue=`echo $sam | grep -P -o '[A-Z][a-z]*?_' | sed 's/_//'`
    dir=${tissue}_${gender}_sam
    echo $dir
    if [[ ! -e $dir ]]; then
        mkdir $dir;
        awk -v dir=$dir '{if(! /^@/ ) {print > dir"/"$3".sam"}}' $sam
    fi
done
