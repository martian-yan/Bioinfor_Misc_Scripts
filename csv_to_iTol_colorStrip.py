#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
A script to write a iTol colorstrip file from columns of a csv file.
Colour palettable is from http://jiffyclub.github.io/palettable/, currently using "tableau_10"

Usage:
python csv_to_iTol_colorStrip.py your_csv_file.csv
date: 2019-07-12
'''
__author__ = 'Yan Li'
__version__= '1.0.0'
__email__= 'yan.li2@outlook.com'


__author__ = 'Yan Li'
__version__= '1.0.0'
__email__= 'yan.li2@outlook.com'

import pandas as pd
import sys
import re
from palettable.cartocolors import qualitative
from palettable import tableau
from collections import Iterable

def main():
    
    MAXIMUM_QUALITIES = 10 # any int no more than 30
    df = pd.read_csv(sys.argv[1], index_col=0)
    columns = df.columns.tolist()
    print(columns)
    barcodes = df.index.values.tolist()

    #colors = qualitative.Vivid_10.hex_colors + qualitative.Bold_10.hex_colors + qualitative.Prism_10.hex_colors
    colors = tableau.Tableau_10.hex_colors
    
    n = 0
    for each_column in columns:
        
        values = df[each_column].tolist()
        values = list(map(lambda a:str(a), values)) # turn all element into str in case 'zip' can't parse int

        # sort values by count
        sorted_values = []
        for count, elem in sorted(((values.count(e), e) for e in set(values)), reverse=True):
            #print(elem)
            sorted_values.append(elem)

        # replace values by "others" if more than X different values
        count = len(sorted_values)
        if count > MAXIMUM_QUALITIES:
            values = list(map(lambda x:"others" if x in sorted_values[MAXIMUM_QUALITIES-1:] else x, values))
            sorted_values = sorted_values[:MAXIMUM_QUALITIES-1]
            sorted_values.append("others")
            count = MAXIMUM_QUALITIES
        
        # create a list of correspondent colors of each value
        print(type(sorted_values), len(sorted_values), isinstance(sorted_values, Iterable))
        color_dict = dict(zip(sorted_values, colors[:count]))
        color_values = [color_dict.get(n, n) for n in values]

        legend_shapes = '\t'.join(['1']*count)
        legend_colors = '\t'.join(colors[:count])
        legend_labels = '\t'.join(sorted_values)

        options = {
            "dataset_label": each_column,
            # legend
            "legend_title": each_column,
            "legend_shapes": legend_shapes,
            "legend_colors": legend_colors,
            "legend_labels": legend_labels, 
            # dataset color
            "color": colors[n],
            # others
            "color_branches":1,
            "strip_width":20,
            "margin":20,
            "border_width":2,
            "border_color":"#ffffff",
            "show_internal":0
        }

        if n < MAXIMUM_QUALITIES:
            n += 1
        else:
            n = MAXIMUM_QUALITIES

        with open('colorStrip_{}.txt'.format(each_column), 'w') as f:
            f.write("DATASET_COLORSTRIP\n")
            f.write("SEPARATOR TAB\n")
            
            for key in options.keys():
                f.write(key.upper()+'\t'+str(options[key])+'\n')
                print(key.upper()+'\t'+str(options[key])+'\n')

            f.write("DATA\n")
            for barcode, color in zip(barcodes, color_values):
                f.write("{0}\t{1}\n".format(barcode, color))
                print("{0}\t{1}".format(barcode, color))   
    return

if __name__ == "__main__":
    main()