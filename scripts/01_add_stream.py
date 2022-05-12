#!/usr/bin/env python3

"""
Author: Megan Brenes Guallar
Contact: m.brenesguallar@genetwister.nl

Description:    this script adds whether a region is up or downstream from the
                mimp to a BED file containing mimp flanking regions
Usage:  python3 01_add_stream.py <region> <outpath>
        Where:
            <region>    is a BED file
            <outpath>   is the path for the output BED file
Output:
        BED file with mimp flanking regions

This script is a part of the effector detection pipeline for F. oxysporum and is
based on the FoEC pipeline by van Dam et al (2016).
"""
# Imports

from sys import argv

# Functions

def add_sense(bedfile, outpath):
    """ Add stream to name

    bedfile --  input, str, path to BED file with downstream region
    outpath --  input, str, path to output BED file with sense column (col 6)
    outpath --  output, str, path to BED file with complete TIRs
    """
    prev_line = ''
    i = 1
    with open(bedfile, 'r') as bed_in, open(outpath, 'w') as outfile:
        for line in bed_in:
            line = line.strip().split('\t')
            label = line[3].split('_')
            if label[1] == 'com':
                if i % 2 == 0:
                    line[3] = '{}_ds'.format(line[3])
                    stream = '+'
                else:
                    line[3] = '{}_us'.format(line[3])
                    stream = '-'
            elif label[1] == 'inc':
                if label[2] == 'l':
                    line[3] = 'mimp_inc_{}_us'.format(line[3][-1])
                    stream = '-'
                elif label[2] == 'r':
                    line[3] = 'mimp_inc_{}_ds'.format(line[3][-1])
                    stream = '+'
            i += 1
            line.append('0')
            line.append(stream)
            outfile.write('{}\n'.format('\t'.join(line)))

    return outpath

def main():
    """Functions to be executed in main.

    """
    bedfile = argv[1]
    outpath = argv[2]
    add_sense(bedfile, outpath)

# Main
if __name__ == "__main__":
    main()