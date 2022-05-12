#!/usr/bin/env python3

"""
Author: Megan Brenes Guallar
Contact: m.brenesguallar@genetwister.nl

Description:    this script processes the gff3 output from SignalP v[4|5] by
                correcting the seqid, start and stop columns to refer to the
                original genomic positions.
Usage:  python3 01_process_signalp.py <sp_gff> <orf_gff> <sp_summary> <outpath>
        <sample>
        Where:
            <sp_gff>        is the SignalP output gff3 file
            <orf_gff>       is the path to the ORF gff3 (orf_finder.py output)
            <sp_summary>    is the path to the SignalP summary file
            <outpath>       is the path for the output gff file
            <sample>        is the name of the sample (used for logging)
Output:
        gff3 file with signal peptides refering to ORF parent nodes

This script is a part of the effector detection pipeline for F. oxysporum and is
based on the FoEC pipeline by van Dam et al (2016).
"""
# Imports

import os
import re
from sys import argv

# Functions
def write_gff(sp_gff, orf_gff, sp_sum, outpath, version):
    """Write gff lines for signal peptide, adapted to genomic coordinates

    sp_gff  --  input, path to input SignalP GFF file
    orf_gff --  input, path to ORF GFF file
    sp_sum  --  input, path to SignalP summary file
    outpath --  input, path to signal peptide genomic GFF file
    version --  input, int, SignalP version [4|5]
    sp_orfs --  output, list of ORF names with signal peptides
    """
    # Labels of ORFs with signal peptides
    sp_orfs = {}
    sp_cds = {}
    final_gff = []
    sp_version = ''
    with open(sp_gff, 'r') as gff_in:
        for line in gff_in:
            if not line.startswith('#'):
                sp_lines = line.strip().split('\t')
                # SignalP seqid
                sp_orf = sp_lines[0]
                # Start, stop, SignalP score
                sp_orfs[sp_orf] = [sp_lines[3], sp_lines[4], sp_lines[5]]
                # Initiate signal peptide CDS dict with empty lists
                sp_cds[sp_orf] = []
                sp_version = sp_lines[1]
    with open(orf_gff, 'r') as orf_in:
        for orf_line in orf_in:
            if not orf_line.startswith('##'):
                orf_line = orf_line.strip().split('\t')
                orf_label = orf_line[8].split(';')[0].split('=')[-1]
                match = re.search('Parent=(.*?);', str(orf_line[8]), re.IGNORECASE)
                if match:
                    m = match.group(1)
                    m_type = orf_line[2]
                    if m in sp_orfs and m_type == 'CDS':
                        sp_cds[m].append(orf_line)
    for rna_key, cds_entry in sp_cds.items():
        sp_values = sp_orfs[rna_key]
        sp_0 = ((int(sp_values[0]) - 1) * 3)
        sp_1 = ((int(sp_values[1]) - 1) * 3)
        strand = cds_entry[0][6]
        if strand == '+':
            cds_start = int(cds_entry[0][3])
            cds_stop = int(cds_entry[0][4])
            sp_entry_start = cds_start + sp_0
            sp_entry_stop = cds_start + sp_1 + 2
        elif strand == '-':
            cds_start = int(cds_entry[-1][3])
            cds_stop = int(cds_entry[-1][4])
            sp_entry_start = cds_stop - sp_0
            sp_entry_stop = cds_stop - sp_1 - 2
        sp_start = min(sp_entry_start, sp_entry_stop)
        sp_stop = max(sp_entry_start, sp_entry_stop)
        region = cds_entry[0][0]
        sp_id = '{}:signal_peptide'.format(rna_key)
        sp_attributes = 'ID={};Parent={}'.format(sp_id, rna_key)
        gff_line = '{}\t{}\tsignal_peptide\t{}\t{}\t{}\t{}\t.\t{}'\
            .format(region, sp_version, sp_start, sp_stop, sp_values[2], strand, \
            sp_attributes)
        final_gff.append(gff_line.split('\t'))
    # Write GFF file with signal peptides and genomic coordinates
    with open(outpath, 'w') as sp_gff_out:
        sp_gff_out.write('##gff-version 3\n')
        for gff_entry in final_gff:
            sp_gff_out.write('\t'.join(gff_entry) + '\n')
    # Print logs
    not_sp = 0
    yes_sp = 0
    with open(sp_sum, 'r') as sp_log:
        for sum_line in sp_log:
            sum_line = sum_line.strip()
            if not sum_line.startswith('#'):
                sum_list = sum_line.split()
                if version == 5:
                    if sum_list[1] == 'OTHER':
                        not_sp += 1
                    elif sum_list[1].startswith('SP'):
                        yes_sp += 1
                elif version == 4:
                    if sum_list[9] == 'N':
                        not_sp += 1
                    elif sum_list[9] == 'Y':
                        yes_sp += 1
    print('// Total sequences analyzed: %i' % (not_sp + yes_sp))
    print('\t// Sequences with predicted signal peptides: %i' % yes_sp)
    print('\t// Sequences without predicted signal peptides: %i' % not_sp)
    print('*' * 80)

    return

def main():
    """Functions to be executed in main.

    """
    print('*' * 80)
    print(' ' * 31 + 'process_signalp.py' + ' ' * 31)
    print('*' * 80)
    sp_gff = argv[1]
    orf_gff = argv[2]
    sp_summary = argv[3]
    outpath = argv[4]
    sample = argv[5]
    version = int(argv[6])
    print('// Sample: "%s"' % sample)
    write_gff(sp_gff, orf_gff, sp_summary, outpath, version)

# Main
if __name__ == "__main__":
    main()