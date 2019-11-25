#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This script use a GFF file and a FNA file to produce DNA gene sequences FFN file.
usage:
    python gff2faa.py annotation.gff sequences.fna output.ffn
date: 2019-11-25
'''
__author__ = 'Yan Li'
__version__= '1.0.0'
__email__= 'yan.li2@outlook.com'

import sys
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re

def main():
    
    in_annotation = sys.argv[1]
    in_seq = sys.argv[2]
    out_file = sys.argv[3]
    #print('Start')
    gff2ffn(in_annotation, in_seq, out_file)

    return

def gff2ffn(in_annotation, in_seq, out_file):

    gff = pd.read_csv(in_annotation, sep='\t', header=None)
    gff.columns = ['chr', 'source', 'type', 'start', 'end', 5, 'strand', 7, 'annotation']
    #print(gff.describe())
    #print(gff.head)

    proteins = gff.to_dict('record')    #dataframe to dict
    #print(proteins)

    with open(in_seq, 'r') as fi_seq:
        chrs = SeqIO.to_dict(SeqIO.parse(fi_seq, "fasta"))
    #print(chrs)

    out_list = [cut_seq(x, chrs) for x in proteins]
    SeqIO.write(out_list, out_file, "fasta")

    return

def cut_seq(gff_record, chrs):

    chr_to_cut = chrs[gff_record['chr']]
    gene_seq = chr_to_cut.seq[gff_record['start']-1:gff_record['end']]
    gene_name = re.search('Name=(.+?);',gff_record['annotation'])[1]    # '?' use non-greedy mode
  
    if gff_record['strand'] == '-':
        gene_seq = gene_seq.reverse_complement()
    
    gene_record = SeqRecord(gene_seq, id=gene_name, name=gene_name, description=gff_record['type']+'|'+gff_record['annotation'])

    return gene_record

if __name__ == "__main__":
    main()