if [[ ! -e ./Homo_sapiens.GRCh38.84.chr.gtf ]]; then
    gffread ./Homo_sapiens.GRCh38.84.chr.gff3 -T -o Homo_sapiens.GRCh38.84.chr.gtf
fi
python3 SUPPA/suppa.py generateEvents -i Homo_sapiens.GRCh38.84.chr.gtf -o Homo_sapiens.GRCh38.84.chr.gtf -b V -e {SE,SS,MX,RI,FL} #-m DEBUG
