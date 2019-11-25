#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This script use a GFF file and a FNA file to produce protein sequences FAA file.
usage:
    python gff2faa.py annotation.gff sequences.fna output.faa
date: 2019-11-25
'''
__author__ = 'Yan Li'
__version__= '1.0.0'
__email__= 'yan.li2@outlook.com'

import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Bio.Data.CodonTable import TranslationError
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
    gff_proteins = gff[gff['type']=='CDS']
    #print(gff_proteins.describe())
    #print(gff_proteins.head)

    proteins = gff_proteins.to_dict('record')    #dataframe to dict
    #print(proteins)

    with open(in_seq, 'r') as fi_seq:
        chrs = SeqIO.to_dict(SeqIO.parse(fi_seq, "fasta"))
    #print(chrs)

    out_list = [cut_seq(x, chrs) for x in proteins]
    SeqIO.write(out_list, out_file, "fasta")

    return

def cut_seq(gff_record, chrs):

    chr_to_cut = chrs[gff_record['chr']]
    gene_name = re.search('Name=(.+?);',gff_record['annotation'])[1]    # '?' use non-greedy mode
    gene_seq = chr_to_cut.seq[gff_record['start']-1:gff_record['end']]

    if gff_record['strand'] == '-':
        gene_seq = gene_seq.reverse_complement()

    '''
    if len(gene_seq)%3 != 0:
        print("\nWarning: {} length is not a multiple of 3".format(gene_name))

    protein_seq = gene_seq.translate()
    protein_seq = protein_seq.rstrip('*')
    '''
    try:
        protein_seq = gene_seq.translate(table="Bacterial", cds=True)
    except TranslationError as e:
        print("\n{} is not a CDS".format(gene_name))
        print(e)
        protein_seq = Seq('NNN')
    #print(len(gene_seq))
    
    gene_record = SeqRecord(protein_seq, id=gene_name, name=gene_name, description=gff_record['annotation'])

    return gene_record

if __name__ == "__main__":
    main()