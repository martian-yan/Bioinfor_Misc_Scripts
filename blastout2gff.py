#!/usr/bin/python3
# -*- coding: utf-8 -*-
'''
A script make GFF file from blast output.
A key feature is decide the +/- strand
date: 2019-08-05
'''
__author__ = 'Yan Li'
__version__= '1.0.0'
__email__= 'yan.li2@outlook.com'

import sys

def main(blast_outfile, gff_file):

    with open(blast_outfile, 'r') as fi:
        records = fi.readlines()[1:] # ignore the header line

    with open(gff_file, 'w') as fo:
        fo.write('#GFF file created from blast result {}'.format(blast_outfile))
        fo.write('\n')
        
        for record in records:
            record = record.split('\t')
            gff_feature = [record[0], 'blast', 'CDS', record[6], record[7], '.', '+', '.', 'ID=Gene_{0};Name={1};Identity={2};Length={3}'.format(record[6],record[1],record[2],record[3])]
            if record[8] > record[9]: # the start point > end point
                gff_feature[6] = '-'
            gff_feature = '\t'.join(gff_feature)
            fo.write(gff_feature)
            fo.write('\n')
    return

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])