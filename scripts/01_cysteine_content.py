#!/usr/bin/env python3

"""
Author: Megan Brenes Guallar
Contact: m.brenesguallar@genetwister.nl

Description:    This is a script to calculate the number of cysteines in a
                protein sequence and the total sequence lemgth.
Usage:  python3 01_cysteine_content.py <sequence> <cys_threshold> <min_len>
        <max_len> <out_seqs> <keep>
        Where:
            <sequence>      is the protein FASTA file
            <cys_threshold> is the minimum required number of cysteines in a
                            sequence to be considered a putative effector
            <min_len>       is the minimum protein length
            <max_len>       is the maximum protein length
            <out_seqs>      is the protein FASTA file with sequences that pass
                            the cysteine threshold
            <keep>          is the path to output text file with FASTA headers
                            from sequences that pass the filter criteria
Output:
    FASTA file (protein) with sequences that meet the length and cysteine
    requirements
    Text file with FASTA headers from sequences that meet the length and
    cysteine requirements

This script is a part of the effector detection pipeline for F. oxysporum and is
based on the FoEC pipeline by van Dam et al (2016).
"""
# Imports

import os
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sys import argv

# Functions
def parse_fasta(sequence):
    """Parses FASTA file and puts sequence into a dict

    sequence    --  input, protein sequence
    prot_dict   --  output, dict, {seq_name: seq}
    """
    prot_dict = {}
    with open(sequence, 'r') as prot_in:
        for line in prot_in:
            line = line.strip()
            if line.startswith('>'):
                fasta_id = line[1:]
                prot_dict[fasta_id] = []
            else:
                prot_dict[fasta_id].append(line)
    for fid, seqs in prot_dict.items():
        prot_dict[fid] = ''.join(seqs)

    return prot_dict

def get_cysteines(prot_dict, cys_thresh, min_len, max_len, out_seqs, keep):
    """Determines number of cysteines in a protein sequence.

    prot_dict   --  input, dict, {seq_name: seq}
    cys_thresh  --  input, int, threshold for number of cysteines in a
                    sequence to be considered a putative effector
    min_len     --  input, int, minimum protein length threshold
    max_len     --  input, int, maximum protein length threshold
    out_seqs    --  input, path to output FASTA file with sequences that pass
                    the filter criteria
    keep        --  input, path to output text file with FASTA headers from
                    sequences that pass the filter criteria
    """
    all_filters = {}
    putative_effectors = {}
    for header, seq in prot_dict.items():
        prot = ProteinAnalysis(seq)
        cys = prot.count_amino_acids()['C']
        seq_len = len(seq)
        all_filters[header] = {'cys': cys, 'len': seq_len}
    for entry, values in all_filters.items():
        if int(min_len) <= int(values['len']) <= int(max_len):
            if int(values['cys']) >= int(cys_thresh):
                putative_effectors[entry] = prot_dict[entry]

    with open(out_seqs, 'w') as outfile, open(keep, 'w') as keepfile:
        for key, value in putative_effectors.items():
            outfile.write('>%s\n%s\n' % (key, value))
            keepfile.write('%s\n' % key)

def main():
    """Functions to be executed in main.

    """
    sequence = argv[1]
    cys_thresh = argv[2]
    min_len = argv[3]
    max_len = argv[4]
    out_seqs = argv[5]
    keep = argv[6]
    prot_dict = parse_fasta(sequence)
    get_cysteines(prot_dict, cys_thresh, min_len, max_len, out_seqs, keep)
# Main
if __name__ == "__main__":
    main()