#!/usr/bin/env python3

"""
Author: Megan Brenes Guallar
Contact: m.brenesguallar@genetwister.nl

Description:    this script detects complete and incomplete mimps in a given
                genome based on a BED file containing mimp TIR coordinates
                (output from mimpTIR_finder.py). Incomplete mimps also show if
                they are pointing left (l) or right (r).
Usage:  python3 01_complete_mimps.py <tir_bed> <outpath> <outpath_left>
        <outpath_right> <gff_out>
        Where:
            <tir_bed>       is a BED file containing mimp TIR coordinates
            <outpath>       is the path for the output BED file with complete
                            mimps
            <outpath_left>  is the path for the output BED file with left
                            pointing TIRs (TIRrc)
            <outpath_right> is the path for the output BED file with right
                            pointing TIRs (TIR)
            <gff_out>       is the path for the output GFF (without header) file
                            with mimp and TIR info
Output:
        BED file with (in)complete mimp sequence coordinates
        GFF lines with complete mimps and TIRs

This script is a part of the effector detection pipeline for F. oxysporum and is
based on the FoEC pipeline by van Dam et al (2016).
"""
# Imports

import os
from sys import argv

# Functions

def check_TIRs(tir_bed, outpath):
    """ Checks if TIRs in a sorted BED file form pairs (within 400 bp).

    tir_bed         --  input, str, path to BED file with TIRs
    outpath         --  input, str, path to output BED file with complete TIRs
    complete        --  output, set, labels of TIRs forming complete mimps
    outpath         --  output, str, path to BED file with complete TIRs
    mimp_tir_dict   --  output, dict with mimps and corresponding TIR pairs
    """
    prev_line = ''
    complete = set()
    mimp_tir_dict = {}
    tir_dict = {}
    mimp_dict = {}
    i = 1
    n = 1
    it = 0
    with open(tir_bed, 'r') as bed_in, open(outpath, 'w') as outfile:
        for line in bed_in:
            line = line.strip().split('\t')
            tir_dict[line[-1]] = [int(line[1]), int(line[2])]
            tir_type = line[3].split('_')[0]
            if prev_line == '':
                prev_line = line
                continue
            else:
                prev_tir_type = prev_line[3].split('_')[0]
            # A pair is formed by TIRrc + TIR, and they must be in the same
            # genomic region (line[0])
            if prev_tir_type == 'TIRrc' and tir_type != prev_tir_type and \
            line[0] == prev_line[0]:
                # The max distance between TIRs in a pair is 400 bp
                dist = int(prev_line[2]) + 400
                if dist >= int(line[2]):
                    label = 'mimp_com_%s' % i
                    mimp = '{}\t{}\t{}\t{}\n'.format(line[0], prev_line[1], \
                    line[2], label)
                    outfile.write(mimp)
                    i += 1
                    complete.add(prev_line[-1])
                    complete.add(line[-1])
                    mimp_tir_dict[label] = [prev_line[3], line[3]]
                    mimp_dict[label] = [line[0], int(prev_line[1]), int(line[2])]
            it += 1
            prev_line = line

    return complete, tir_dict, mimp_dict, mimp_tir_dict

def add_incomplete(tir_bed, out_l, out_r, complete):
    """Adds incomplete mimps (only 1 TIR found) to BED file

    tir_bed     --  input, str, path to BED file with TIRs (incomplete)
    out_l       --  input, str, path to output BED file with incomplete TIRs
                    pointing left (TIRrc)
    out_r       --  input, str, path to output BED file with incomplete TIRs
                    pointing right
    complete    --  input, set, labels of TIRs forming complete mimps
    incomplete  --  output, dict, single TIRs with region label as value
    out_[l|r]   --  output, str, path to BED file with TIRs
    """
    i = 1
    incomplete = {}
    with open(tir_bed, 'r') as bed_in, open(out_l, 'w') as left, open(\
    out_r, 'w') as right:
        for line in bed_in:
            line = line.strip().split('\t')
            tir_type = line[3].split('_')[0]
            if line[-1] not in complete:
                incomplete[line[-1]] = line[0]
                # Indicate whether you will need left flank (TIRrc points left)
                # or right flank (TIR points right)
                if tir_type == 'TIRrc':
                    tir_lab = 'l'
                    mimp = '{}\t{}\t{}\t{}\n'.format(line[0], line[1], \
                    line[2], 'mimp_inc_%s_%s' % (tir_lab, i))
                    i += 1
                    left.write(mimp)
                else:
                    tir_lab = 'r'
                    mimp = '{}\t{}\t{}\t{}\n'.format(line[0], line[1], \
                    line[2], 'mimp_inc_%s_%s' % (tir_lab, i))
                    i += 1
                    right.write(mimp)

    return incomplete

def write_gff(tir_dict, mimp_dict, mimp_tir_dict, incomplete, gff_out):
    """Writes GFF3 file with mimp and TIR information

    tir_dict        --  input, dict with TIRs and coords
    mimp_dict       --  input, dict with mimps and coords
    mimp_tir_dict   --  input, dict with mimps and corresponding TIR pairs
    incomplete      --  input, dict, single TIRs with region label as value
    gff_out         --  input, path to output GFF file
    """
    with open(gff_out, 'w') as outfile:
        # Write GFF lines for complete mimps and corresponding TIRs
        for key, value in mimp_dict.items():
            mimp_id = '{}-mimp:{}'.format(value[0], key.split('_')[-1])
            # mimp 1-based coordinates
            mimp_start = int(value[1]) + 1
            mimp_stop = int(value[2])
            mimp_attributes = 'ID={};Name={}'.format(mimp_id, key)
            m_info = '{}\tFOEC\ttransposable_element\t{}\t{}\t.\t.\t.\t{}\n'\
            .format(value[0], mimp_start, mimp_stop, mimp_attributes)
            # TIR name
            tir_a = mimp_tir_dict[key][0]
            tir_b = mimp_tir_dict[key][1]
            # TIR 1-based coordinates
            tir_a_start = int(tir_dict[tir_a][0]) + 1
            tir_a_stop = int(tir_dict[tir_a][1])
            tir_b_start = int(tir_dict[tir_b][0]) + 1
            tir_b_stop = int(tir_dict[tir_b][1])
            tir_attr_a = 'ID={}-tir:{}a;Name={};Parent={}'.format(value[0], \
            tir_a.split('_')[-1], tir_a, mimp_id)
            tir_attr_b = 'ID={}-tir:{}b;Name={};Parent={}'.format(value[0], \
            tir_b.split('_')[-1], tir_b, mimp_id)
            tir_info_a = '{}\tFOEC\tterminal_inverted_repeat\t{}\t{}\t.\t-\t.\t{}\n'\
            .format(value[0], tir_a_start, tir_a_stop, tir_attr_a)
            tir_info_b = '{}\tFOEC\tterminal_inverted_repeat\t{}\t{}\t.\t+\t.\t{}\n'\
            .format(value[0], tir_b_start, tir_b_stop, tir_attr_b)
            outfile.write(m_info)
            outfile.write(tir_info_a)
            outfile.write(tir_info_b)
        # Write GFF lines for single TIRs
        for tir, reg_label in incomplete.items():
            tir_start = int(tir_dict[tir][0]) + 1
            tir_stop = int(tir_dict[tir][1])
            tir_lab = tir.split('_')
            if tir_lab[0].endswith('c'):
                tir_side = 'l'
                tir_sense = '-'
            else:
                tir_side = 'r'
                tir_sense = '+'
            tir_attr = 'ID={}-single_tir:{}{}-;Name={}'.format(reg_label, \
            tir_lab[-1],  tir_side, tir)
            tir_info = '{}\tFOEC\tterminal_inverted_repeat\t{}\t{}\t.\t{}\t.\t{}\n'\
            .format(reg_label, tir_start, tir_stop, tir_sense, tir_attr)
            outfile.write(tir_info)

def main():
    """Functions to be executed in main.

    """
    tir_bed = argv[1]
    outpath = argv[2]
    outpath_left = argv[3]
    outpath_right = argv[4]
    gff_out = argv[5]
    complete, tir_dict, mimp_dict, mimp_tir_dict = check_TIRs(tir_bed, outpath)
    incomplete = add_incomplete(tir_bed, outpath_left, outpath_right, complete)
    write_gff(tir_dict, mimp_dict, mimp_tir_dict, incomplete, gff_out)

# Main
if __name__ == "__main__":
    main()
