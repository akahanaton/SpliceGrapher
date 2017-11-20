#!/usr/bin/env bash
if [ $# -lt 1  ]; then
    echo "Usage : $0 DirName"
    exit;
fi
#--------------------------------------------------
# set –euo pipefail
#--------------------------------------------------
shDir=$1
for file in `ls -f $shDir/`; do if [[ -e $shDir/$file.done  ]]; then mv ${shDir}/$file* ${shDir}.done/ ; fi; done
