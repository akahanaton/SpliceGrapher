#!/usr/bin/env bash

snakemake -q -j 100 -k --cluster "qsub -o ./logs -j oe -d ./" -s ./convert.snake
