#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Visualize plasmids, color the genes according to the source: AMR gene, IS gene or Plasmids probes.
Usage:
python visualize_plasmids.py genebank_file
date: 2021-02-07
'''
__author__ = 'Yan Li'
__email__ = 'yan.li2@outlook.com'

import matplotlib.colors as mcolors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

import sys

genbank_file = sys.argv[1]
svg_file = genbank_file.replace("genbank", "svg")

record = SeqIO.read(genbank_file, 'genbank')
gd_diagram = GenomeDiagram.Diagram(record.id)
gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
gd_feature_set = gd_track_for_features.new_set()
print(record)

for feature in record.features:
    print(feature)
    color = mcolors.TABLEAU_COLORS['tab:gray']
    if "ISfinder" in feature.qualifiers['note']:
        color = mcolors.TABLEAU_COLORS['tab:blue']
    elif ("resfinder"  in feature.qualifiers['note']) or ("card" in feature.qualifiers['note']):
        color = mcolors.TABLEAU_COLORS['tab:red']
    elif "plasmidfinder" in feature.qualifiers['note']:
        color = mcolors.TABLEAU_COLORS['tab:orange']

    angle = 45
    if feature.location.strand == 1:
        angle = 45
    elif feature.location.strand == -1:
        angle = 135

    gd_feature_set.add_feature(feature,
                               sigil="BIGARROW",
                               arrowhead_length=0.25,
                               color=color,
                               label=True,
                               label_size=14,
                               label_angle=45,
                               label_position='middle',
                               label_strand=1)

gd_diagram.draw(format="linear",
                pagesize=(5*cm, 30*cm),
                yt=0.4,
                fragments=1,
                start=0,
                end=len(record))

gd_diagram.write(svg_file, "SVG")