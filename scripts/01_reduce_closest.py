#!/usr/bin/env python3

"""
Author: Megan Brenes Guallar
Contact: m.brenesguallar@genetwister.nl

Description:    this script reduces the tsv file of closest mimps by taking into
                account incomplete mimp directionality. If an incomplete mimp
                points upstream, it cannot be considered the 'closest mimp' to
                an ORF which is downstream of the mimp.
Usage:  python3 01_reduce_closest.py <closest> <outpath> <sample>
        Where:
            <closest>   is a tsv with the mimps (both complete and incomplete)
                        closest to each ORF
            <outpath>   is the path for the output BED file
            <sample>    is the name of the sample (used for logging)
Output:
        tsv with mimps closest to a given ORF.
        Columns:
            orf region name, orf start, orf stop, orf name, orf strand,
            mimp region name, mimp start, mimp stop, mimp name, mimp strand

This script is a part of the effector detection pipeline for F. oxysporum and is
based on the FoEC pipeline by van Dam et al (2016).
"""
# Imports

import pandas as pd
from datetime import datetime
from sys import argv

# Functions

def true_closest(closest, outpath):
    """ Reduce closest mimp tsv based on incomplete mimp directionality

    closest --  input, str, path to tsv file with closest mimps
    outpath --  input, str, path to reduced tsv file with closest mimps
    outpath --  output, str, path to reduced tsv file with closest mimps
    """
    df = pd.read_csv(closest, sep="\t", header = None)
    # Remove rows which have a right pointing incomplete mimp matched to an orf
    # which is upstream or a left pointing incomplete mimp match to an orf which
    # is downstream
    df = df.loc[~(((df[11] > 0) & (df[10].str.contains('_r_')) | (df[11] < 0) &
    (df[10].str.contains('_l_'))))]
    # Add a new column with absolute values (distance)
    df[12] = df[11].abs()
    # Sort by orf position and then mimp distance
    df = df.sort_values(by = [1, 12], ascending = [True, True])
    df = df.drop_duplicates(3, keep = 'first').drop([6, 12], axis = 1)
    # Write to tsv
    df.to_csv(outpath, sep = '\t', index = False, header = False)

    return outpath

def write_log(outpath, sample):
    """Writes log for ORFs per mimp/TIR

    outpath --  input, str, path to reduced tsv file with closest mimps
    sample  --  input, sample name
    """
    orf_dict = {}
    # Count lines (number of ORFs)
    i = 0
    with open(outpath, 'r') as incsv:
        for line in incsv:
            line = line.strip().split('\t')
            if line[9] not in orf_dict:
                orf_dict[line[9]] = []
            orf_dict[line[9]].append(line[3])
            i += 1
    print('// Sample: "%s"' % sample)
    print('// Total number of ORFs: %i' % i)
    print('// Transposable elements with downstream ORFs: %i' % len(orf_dict))
    print('*' * 80)
    for te_key, orf_value in orf_dict.items():
        print('// Transposable element: %s' % te_key)
        print('// ORFs found in downstream regions: %i' % len(orf_value))
        orf_value = sorted(orf_value, key = lambda orf_name: float(orf_name[4:]))
        for orf in orf_value:
            print('\t%s' % orf)

def main():
    """Functions to be executed in main.

    """
    start_time = datetime.now()
    closest = argv[1]
    outpath = argv[2]
    sample = argv[3]
    true_closest(closest, outpath)
    print('*' * 80)
    print(' ' * 31 + 'reduce_closest.py' + ' ' * 31)
    print('*' * 80)
    write_log(outpath, sample)
    total_time = datetime.now() - start_time
    print('// Time needed: %s' % total_time)

# Main
if __name__ == "__main__":
    main()