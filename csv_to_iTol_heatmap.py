#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
A script to write a iTol heatmap file from columns of a csv file.
The options should be customized

Usage:
python csv_to_iTol_heatmap.py your_csv_file.csv
date: 2019-03-05
'''
__author__ = 'Yan Li'
__version__= '1.0.0'
__email__= 'yan.li2@outlook.com'

import pandas as pd
import sys
import re

options = {
    "dataset_label": "AMR genes",
    "color": "#CD0102",
    "margin": 0,
    "strip_width": 120,
    "show_internal": 0,
    "show_tree": 1,
    "color_min": "#e6e6e6",
    "color_max": "#CD0102",
    "use_mid_color": 0,
    #"color_mid": "#ffff00" 
}

def main():
    
    df = pd.read_csv(sys.argv[1], index_col=0,sep='\t')
    columns = df.columns.tolist()
    df = df[columns]
    write_hm(options, df)
    return

def write_hm(options, df):
    filename = re.sub('\.[a-z]+','',sys.argv[1])

    with open('heatmap_'+filename+'.txt', 'w') as f:
        f.write("DATASET_HEATMAP\n")
        f.write("SEPARATOR TAB\n")
        
        for key in options.keys():
            f.write(key.upper()+'\t'+str(options[key])+'\n')
                        
        labels = "FIELD_LABELS"
        for column in df.columns:
            labels = labels+'\t'+column
        f.write(labels+'\n')

        f.write("DATA\n")

    df.to_csv('heatmap_'+filename+'.txt', header=None, sep='\t', mode='a')

    return

if __name__ == "__main__":
    main()