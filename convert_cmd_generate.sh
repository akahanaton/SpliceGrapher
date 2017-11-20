#!/usr/bin/env bash
for pdf in sex_biased_pdf/*;
do
    shFile=convert_sh/`basename $pdf .pdf`.sh
    echo "cd \$PBS_O_WORKDIR" > $shFile
    echo "convert $pdf `echo $pdf | sed 's/pdf/png/g'`" >> $shFile
done
