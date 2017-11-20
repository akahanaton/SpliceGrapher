#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

import psyco_full
import bx.align.maf

from Bio import SeqIO
from BCBio import GFF
from collections import defaultdict
import argparse
import sys
import pprint
import GTF                      # https://gist.github.com/slowkow/8101481
from docopt import docopt
import pandas as pd
import gzip
import time
from contextlib import contextmanager
import re
from collections import namedtuple
import os.path

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def nested_dict(n, type):
    if n == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: nested_dict(n-1, type))

@contextmanager
def log(message):
    """Log a timestamp, a message, and the elapsed time to stderr."""
    start = time.time()
    sys.stderr.write("{} # {}\n".format(time.asctime(), message))
    yield
    elapsed = int(time.time() - start + 0.5)
    sys.stderr.write("{} # done in {} s\n".format(time.asctime(), elapsed))
    sys.stderr.flush()

def main(args):
    with log("Reading compare Gencode annotation file: {}".format(args.compGffFile)):
        gc = GTF.dictionary(args.compGffFile,"ID")
    compGeneInfo = gc['gene']
    #--------------------------------------------------
    # gene['Name'] = gene['ID'].map(lambda x: re.sub(r':maker.*','',x))
    #--------------------------------------------------

    with log("Reading reference Gencode annotation file: {}".format(args.refGffFile)):
        gc = GTF.dataframe(args.refGffFile)

    # Select just genes of protein coding genes, and columns that we want to use.
    idx = (gc.feature == 'gene')
    gene = gc.ix[idx, ['seqname','start','end','Name']]
    #--------------------------------------------------
    # print(gene)
    #--------------------------------------------------
    # Convert columns to proper types.
    gene.start = gene.start.astype(int)
    gene.end = gene.end.astype(int)

    for geneID in gene['Name']:
        if geneID in compGeneInfo:
            # gene annotated in both species, read coordinates projecting information in maf
            mafFile = args.mafPath + "/" + geneID + ".maf"
            if not os.path.exists(mafFile):
                continue
            with log("Reading the Maf file: {}".format(mafFile)):
                with open(mafFile) as maf:
                    out_files = dict()
                    geneCoords = dict()
                    for block in bx.align.maf.Reader(maf):
                        ref_comp = block.components[0]
                        refSpecies, refChrom = ref_comp.src.split('.')[:2]
                        if refSpecies not in geneCoords:
                            geneCoords[refSpecies] = nested_dict(2,str)
                            geneCoords[refSpecies]['refInfo']['start'] = ref_comp.forward_strand_start
                            geneCoords[refSpecies]['refInfo']['end'] = ref_comp.forward_strand_end
                            geneCoords[refSpecies]['refInfo']['chr'] = refChrom
                        for comp in block.components[1:]:
                            comp_species, compChrom = comp.src.split('.')[:2]
                            if comp_species not in geneCoords:
                                geneCoords[comp_species] = nested_dict(2,str)
                                geneCoords[comp_species][compChrom]['start'] = comp.start
                                geneCoords[comp_species][compChrom]['end'] = int(comp.end)
                            if compChrom not in geneCoords[comp_species]:
                                geneCoords[comp_species][compChrom]['start'] = comp.start
                                geneCoords[comp_species][compChrom]['end'] = int(comp.end)
                            if comp_species not in out_files:
                                bedfile = "%s/%s.%s.bed" % (args.mafPath, geneID, comp_species )
                                f = open( bedfile , "w" )
                                out_files[comp_species] = f
                            pid = block_pid( ref_comp, comp )

                            if pid:
                                #--------------------------------------------------
                                # print("%s\t%s" % (comp.end, geneCoords[comp_species][compChrom]))
                                #--------------------------------------------------
                                if geneCoords[refSpecies]['refInfo']['start'] > ref_comp.forward_strand_start:
                                    geneCoords[refSpecies]['refInfo']['start'] = ref_comp.forward_strand_start
                                if geneCoords[refSpecies]['refInfo']['end'] < ref_comp.forward_strand_end:
                                    geneCoords[refSpecies]['refInfo']['end'] = ref_comp.forward_strand_end
                                if geneCoords[comp_species][compChrom]['start'] > comp.start:
                                    geneCoords[comp_species][compChrom]['start'] = comp.start
                                if geneCoords[comp_species][compChrom]['end'] <= int(comp.end):
                                    geneCoords[comp_species][compChrom]['end'] = int(comp.end)
                                out_files[comp_species].write( "%s\t%d\t%d\t%s:%d-%d,%s\t%f\n" %
                                                ( refChrom, ref_comp.forward_strand_start, ref_comp.forward_strand_end, \
                                                compChrom, comp.start, comp.end, comp.strand, pid ) )

                    for f in out_files.values():
                        f.close()
            if args.compSpecies in geneCoords:
                for chrom in geneCoords[args.compSpecies]:
                    if chrom in compGeneInfo[geneID]:
                        annoStart = int(compGeneInfo[geneID][chrom]['start'])
                        annoEnd = int(compGeneInfo[geneID][chrom]['end'])
                        compStart = int(geneCoords[args.compSpecies][chrom]['start'])
                        compEnd = int(geneCoords[args.compSpecies][chrom]['end'])
                        if compEnd > annoEnd and compStart < annoEnd or compEnd > annoStart and compStart < annoEnd :
                            print("Matched\t%s\t%s\t%s\tanno: %d - %d\tmapped: %d - %d\t%s\t%s\t%s: %s - %s" % \
                                    (geneID, args.compSpecies, chrom, annoStart,annoEnd, compStart, compEnd, compGeneInfo[geneID][chrom]['ID'], \
                                    args.refSpecies,geneCoords[args.refSpecies]['refInfo']['chr'],  geneCoords[args.refSpecies]['refInfo']['start'], geneCoords[args.refSpecies]['refInfo']['end'] ))
                        else:
                            eprint("unMatched\t%s\t%s\t%s\tanno: %d - %d\tmapped: %d - %d" % \
                                    (geneID, args.compSpecies, chrom, annoStart, annoEnd, compStart, compEnd))
                    else:
                        eprint("Error Chrom\t%s\t%s\t%s\tmapped: %s - %s" % \
                                (geneID, args.compSpecies, chrom, geneCoords[args.compSpecies][chrom]['start'],geneCoords[args.compSpecies][chrom]['end']))

def block_pid( comp1, comp2 ):
    match = 0
    total = 0
    t1 = comp1.text.lower()
    t2 = comp2.text.lower()
    for i in range( 0, len(t1) ):
        a, b = t1[i], t2[i]
        if a == '-' or b == '-':
            continue
        elif a == b:
            match += 1
        total += 1
    if total == 0: return None
    return ( match / total )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get the ortholog region according to maf alignment')
    parser.add_argument('refGffFile', type=str, help='reference gff file')
    parser.add_argument('compGffFile', type=str, help='compare gff file')
    parser.add_argument('mafPath', type=str, help='maf dir path')
    parser.add_argument('refSpecies', type=str, help='name of reference species')
    parser.add_argument('compSpecies', type=str, help='name of compare species')
    parser.add_argument('-f', dest='fastaFile', type=str, help='fasta file matching the gff file')
    #--------------------------------------------------
    # parser.add_argument('-c', dest='cleanCodon', type=str, default='Y', choices=['Y','N'], help='clean stop codon: Y/N')
    # parser.add_argument('-r', dest='cleanByRef', type=str, default=None, help='aligment format')
    #--------------------------------------------------
    args = parser.parse_args()
    main(args)
