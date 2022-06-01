#!/usr/bin/env python3

"""
Author: Megan Brenes Guallar
Contact: m.brenesguallar@genetwister.nl

Description:    This is a script to filter a presence absence variation (PAV)
                table. The filter consists in calculating the column average
                (hits per putative effector) and removing columns with an
                average above a set mean_threshold (default = 15). These columns are
                likely transposable elements and not putative effectors.
Usage:  python3 03_filter_pav.py <pav> <mean_threshold> <individual_threshold>
        <pav_out> <peff_outfile> <peff_csv>
        Where:
            <pav>                   is the unfiltered PAV table
            <mean_threshold>        is max column mean allowed (number of
                                    hits for a putative effector in a genome).
                                    Columns exceeding this mean will be dropped.
            <individual_threshold>  is the max value allowed. If a value exceeds
                                    it, the putative effector group will be
                                    dropped
            <pav_out>               is the PAV csv outfile
            <peff_outfile>          is the path to output file of putative
                                    effector IDs
            <peff_csv>              is the path to output putative effector csv
                                    file for visualization

Output:
    '01_presence_absence.csv', csv with PAV of putative effectors in genomes

This script is a part of the effector detection pipeline for F. oxysporum and is
based on the FoEC pipeline by van Dam et al (2016).
"""
# Imports

import os
import pandas as pd
from sys import argv

# Functions
def filt(pav, mean_threshold, individual_threshold, outfile):
    """Filters PAV columns out

    pav                     --  input, presence absence variations for putative
                                effectors (cols) per genome (rows). Numbers
                                represent putative effector hits to genomes
    mean_threshold          --  input, int, max column mean allowed (number of
                                hits for a putative effector in a genome).
                                Columns exceeding this will be dropped.
    individual_threshold    --  input, int, max value allowed. If a value
                                exceeds it, the putative effector group will be
                                dropped
    outfile                 --  input, path to output file with reduced pav csv
    col_list                --  output, list of strs, putative effector labels
    """
    pav_tab = pd.read_csv(pav, sep = '\t', header = 0, index_col = 0)
    pav_mean = pav_tab.mean() >= float(mean_threshold)
    pav_tab = pav_tab[pav_tab.columns[~pav_mean]]
    pav_tab = pav_tab.loc[:, ~(pav_tab > float(individual_threshold)).any()]
    pav_tab.to_csv(outfile, sep = '\t')
    col_list = list(pav_tab.columns)

    return col_list

def get_eff_numbers(col_list, peff_outfile, peff_csv):
    """Gets number IDs from putative effectors.

    col_list        --  input, list of strs, putative effector labels
    peff_outfile    --  input, path to output file of putative effector IDs
    peff_csv        --  input, path to output putative effector csv file for
                        visualization
    """
    with open(peff_outfile, 'w') as pout, open(peff_csv, 'w') as csv:
        csv.write('Putative effectors\n')
        for peff in col_list:
            # if peff.startswith('p_effector_'):
            #     e_id = peff.split('_')[-1]
            # else:
            #     e_id = peff
            e_id = peff.split('_')[-1]
            pout.write(e_id + '\n')
            csv.write('%s\n' % peff)

def main():
    """Functions to be executed in main.

    """
    pav = argv[1]
    mean_threshold = argv[2]
    individual_threshold = argv[3]
    outfile = argv[4]
    peff_outfile = argv[5]
    peff_csv = argv[6]
    col_names = filt(pav, mean_threshold, individual_threshold, outfile)
    get_eff_numbers(col_names, peff_outfile, peff_csv)

# Main
if __name__ == "__main__":
    main()