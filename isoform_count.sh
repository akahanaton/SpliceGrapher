#!/usr/bin/env bash
if [[ $# -lt 1  ]]; then
    echo "Usage : $0  input.gff"
    exit;
fi
gff=$1

set –euo pipefail

awk '$3=="mRNA"' $gff | head -10



