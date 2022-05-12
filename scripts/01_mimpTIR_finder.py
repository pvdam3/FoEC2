#!/usr/bin/env python3

"""
Author: Megan Brenes Guallar
Contact: m.brenesguallar@genetwister.nl

Description:    this script detects mimp TIRs in a given genome
Usage:  python3 01_mimpTIR_finder.py <genome> <outpath> <sample>
        Where:
            <genome>    is a genome FASTA file (DNA)
            <outpath>   is the path for the output BED file
            <sample>    is the name of the sample (used for logging)
Output:
    BED file with start and stop positions of found mimp TIRs.

This script is a part of the effector detection pipeline for F. oxysporum and is
based on the FoEC pipeline by van Dam et al (2016).
"""
# Imports

import re
from datetime import datetime
from sys import argv
from Bio import SeqIO

# Functions

def find_TIRs(genome, outpath):
    """ Finds TIR motifs present in a genome FASTA file and records coordinates.

    genome  --  input, FASTA file, F. oxysporum genome
    outpath --  input, str, path to output BED file
    outpath --  output, str, path to BED file with TIR coordinates
    """
    print('// Motifs:')
    TIR_motif = 'TT[TA]TTGC..CCCACTG..'
    TIR_motif_rc = '..CAGTGGG..GCAA[TA]AA'
    print('\t "%s"' % TIR_motif)
    print('\t "%s"' % TIR_motif_rc)
    # Count TIR and TIR rc occurrences (1-based)
    i = 1
    n = 1
    tot = 0
    motif_seqs = 0
    log_lines = []
    # Get records from FASTA file
    with open(outpath, 'w') as outfile:
        for seq_record in SeqIO.parse(genome, 'fasta'):
            motif = False
            tot += 1
            log_lines.append('// Analyzing \'%s\'... (%i bp long)' % \
            (seq_record.description, len(seq_record)))
            # Get TIR consensus sequence match
            match = re.finditer(TIR_motif, str(seq_record.seq), re.IGNORECASE)
            # Get TIR rc match
            match_rc = re.finditer(TIR_motif_rc, str(seq_record.seq),
            re.IGNORECASE)
            # match coords are 0-based, like BED file
            for m in match:
                outfile.write('{}\t{}\t{}\tTIR_{}\n'.format(\
                seq_record.id, m.start(), m.end(), i))
                i += 1
                motif = True
                log_lines.append('\t// Match found! TIR_%i' % i)
            for m_rc in match_rc:
                outfile.write('{}\t{}\t{}\tTIRrc_{}\n'.format(\
                seq_record.id, m_rc.start(), m_rc.end(), n))
                n += 1
                motif = True
                log_lines.append('\t// Reverse complement match found! TIRrc_%i' % n)
            if motif:
                motif_seqs += 1
    print('// Total sequences analyzed: %i' % tot)
    print('\t// Sequences containing a motif: %i' % motif_seqs)
    print('// Total motifs found: %i' % (i + n))
    print('\t// TIR motifs found: %i' % i)
    print('\t// Reverse complement TIR motifs found: %i' % n)
    print('*' * 80)
    for log_line in log_lines:
        print(log_line)

    return outpath

def main():
    """Functions to be executed in main.

    """
    start_time = datetime.now()
    genome = argv[1]
    outpath = argv[2]
    sample_name = argv[3]
    print('*' * 80)
    print(' ' * 31 + 'mimpTIR_finder.py' + ' ' * 31)
    print('*' * 80)
    print('// Sample: "%s"' % sample_name)
    find_TIRs(genome, outpath)
    total_time = datetime.now() - start_time
    print('// Time needed: %s' % total_time)

# Main
if __name__ == "__main__":
    main()
