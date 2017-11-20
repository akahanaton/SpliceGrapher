#!/usr/bin/env python
from __future__ import print_function
from Bio import SeqIO
from Bio.Seq import Seq
from BCBio import GFF
from collections import defaultdict
import argparse
import sys
import pprint
import re
import copy

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def nested_dict(n, type):
    if n == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: nested_dict(n-1, type))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='seperate exon records with multiple parents annotation')
    parser.add_argument('inGffFile', type=str, help='input gff file')
    #--------------------------------------------------
    # parser.add_argument('-c', dest='cleanCodon', type=str, default='Y', choices=['Y','N'], help='clean stop codon: Y/N')
    # parser.add_argument('-r', dest='cleanByRef', type=str, default=None, help='aligment format')
    #--------------------------------------------------
    args = parser.parse_args()
    parent=re.compile(r"Parent=.*?;|Parent=.*$")
    with open(args.inGffFile) as f:
        for line in f:
            arr = line.rstrip().split("\t")
            if arr[2] == "exon":
                parents=parent.search(arr[8]).group().replace("Parent=","").split(',')
                if len(parents) > 1:
                    for p in parents:
                        print( "%s\t%s" % ("\t".join(arr[0:8]), re.sub(r"Parent=.*?;|Parent=.*$", "Parent=%s" % p, arr[8])))
                else:
                    print("%s" % line.rstrip())
            else:
                print("%s" % line.rstrip())
