#!/bin/env python
from __future__ import print_function
from Bio import SeqIO
from Bio.Seq import Seq
from BCBio import GFF
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from collections import defaultdict
import argparse
import sys
import pprint
import re
import copy
import os
from subprocess import call

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def nested_dict(n, type):
    if n == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: nested_dict(n-1, type))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Concat sequence from UCSC aligment')
    parser.add_argument('inGffFile', type=str, help='input gff file')
    parser.add_argument('outputDir', type=str, help='output gff file')
    #--------------------------------------------------
    # parser.add_argument('-c', dest='cleanCodon', type=str, default='Y', choices=['Y','N'], help='clean stop codon: Y/N')
    # parser.add_argument('-r', dest='cleanByRef', type=str, default=None, help='aligment format')
    #--------------------------------------------------

    args = parser.parse_args()

    if not os.path.isdir(args.outputDir):
        os.makedirs(args.outputDir)

    inHandle = open(args.inGffFile)
    for chrName in GFF.parse(inHandle):
        chrDir = args.outputDir + "/" + chrName.id
        print(chrDir)
        chrDoneFile = chrDir + '.done'
        if os.path.exists(chrDoneFile):
            continue
        if not os.path.isdir(chrDir):
            os.makedirs(chrDir)
        for i in range(0,len(chrName.features)): # gene
            #--------------------------------------------------
            # print(chrName.features[i].type)
            #--------------------------------------------------
            if chrName.features[i].type != 'chromosome':
                if 'Name' in chrName.features[i].qualifiers:
                    outHandle = open(chrDir + "/" + chrName.features[i].qualifiers['Name'][0].replace('/','-') + ".gff", 'w')
                    rec=SeqRecord(Seq(''),chrName.id)
                    rec.features = [chrName.features[i]]
                    #--------------------------------------------------
                    # chrName.features[i].seq=Seq('')
                    # chrName.features[i].annotations=dict()
                    #--------------------------------------------------
                    GFF.write([rec], outHandle)
                    outHandle.close()
        call(["touch", chrDoneFile])
    inHandle.close()
