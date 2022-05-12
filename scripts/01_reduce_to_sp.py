#!/usr/bin/env python3

"""
Author: Megan Brenes Guallar
Contact: m.brenesguallar@genetwister.nl

Description:    this script reduces the input GFF containing structural
                genome annotation features to only those features related to a
                signal peptide. The input GFF file is obtained from the
                combination of the signal peptide GFF file with genomic
                coordinates (output from process_signalp.py) and the structural
                genome annotation GFF file.
Usage:  python3 01_reduce_to_sp.py <sp_orf_gff> <outpath>
        Where:
            <sp_orf_gff>    is a GFF file with structural annotations, including
                            signal peptide predictions
            <outpath>       is the path to the output GFF

Output:
    gff3 file (reduced to features containing signal peptides)

This script is a part of the effector detection pipeline for F. oxysporum and is
based on the FoEC pipeline by van Dam et al (2016).
"""
# Imports

import os
import re
from sys import argv

# Functions
def get_ids(sp_orf_gff):
    """Gets feature IDs for signal peptides/mRNAs/genes with a signal peptide

    sp_orf_gff  --  input, path to GFF file
    clear_ids   --  output, list, signal peptide related feature IDs
    """
    # {gene_id: {mrna_id: [features, with, mrna, parent]}}
    gene_dict = {}
    # Signal peptide feature IDs
    sp_ids = []
    # All IDs somehow related to signal peptides
    clear_ids = []
    with open(sp_orf_gff, 'r') as in_gff:
        for line in in_gff:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                line_type = line[2]
                id_match = re.search('ID=(.*?);', str(line[8]), re.IGNORECASE)
                if id_match:
                    line_id = id_match.group(1)
                match = re.search('Parent=(.*?);', str(line[8]), re.IGNORECASE)
                if line_type == 'gene':
                    gene_dict[line_id] = {}
                    cur_gene = line_id
                elif line_type == 'mRNA':
                    if match:
                        mrna_parent = match.group(1)
                        gene_dict[mrna_parent][line_id] = []
                        cur_mrna = line_id
                else:
                    gene_dict[cur_gene][cur_mrna].append(line_id)
                    if line_type == 'signal_peptide':
                        sp_ids.append(line_id)
    # Determine which IDs to keep
    for gene_key, mrna_dict in gene_dict.items():
        for mrna_id, children in mrna_dict.items():
            for sp_id in sp_ids:
                if sp_id in children:
                    clear_ids.append(gene_key)
                    clear_ids.append(mrna_id)
                    for child in children:
                        clear_ids.append(child)

    return clear_ids

def write_gff(sp_orf_gff, clear_ids, outpath):
    """Write gff lines for signal peptide, adapted to genomic coordinates

    sp_orf_gff  --  input, path to GFF file
    clear_ids   --  input, list, signal peptide related feature IDs
    outpath     --  input, path to reduced GFF file
    """
    gff_lines = []
    with open(sp_orf_gff, 'r') as in_gff:
        for line in in_gff:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                id_match = re.search('ID=(.*?);', str(line[8]), re.IGNORECASE)
                if id_match:
                    line_id = id_match.group(1)
                    if line_id in clear_ids:
                        gff_lines.append(line)
    # Write reduced GFF
    with open(outpath, 'w') as sp_gff_out:
        sp_gff_out.write('##gff-version 3\n')
        for gff_entry in gff_lines:
            sp_gff_out.write('\t'.join(gff_entry) + '\n')

    return

def main():
    """Functions to be executed in main.

    """
    sp_orf_gff = argv[1]
    outpath = argv[2]
    clear_ids = get_ids(sp_orf_gff)
    write_gff(sp_orf_gff, clear_ids, outpath)

# Main
if __name__ == "__main__":
    main()