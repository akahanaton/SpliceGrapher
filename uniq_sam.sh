#!/usr/bin/env bash
for sam in *_sam/000000F.sam
do
    echo $sam
    uniq_field.py $sam 1 | grep 'NH:i:1' > $sam.uniq
done
