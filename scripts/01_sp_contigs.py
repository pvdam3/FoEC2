#!/usr/bin/env python3

"""
Author: Megan Brenes Guallar
Contact: m.brenesguallar@genetwister.nl

Description:    this script parses out FASTA sequences containing ORFs with
                predicted signal peptides (SPs) and returns them in a directory.
                Each file in the directory is a single sequence FASTA file.
Usage:  python3 01_sp_contigs.py <sp_gff> <genome> <fasta_dir>
        Where:
            <sp_gff>    is the path to a GFF file containing predicted SPs
            <genome>    is the path to a genome FASTA file
            <fasta_dir> is the path to the output directory

Output:
    Directory with individual FASTA files per sequence containing ORFs with
    predicted SPs.
This script is a part of the effector detection pipeline for F. oxysporum and is
based on the FoEC pipeline by van Dam et al (2016).
"""
# Imports

import os
from sys import argv

# Functions
def get_names(sp_gff):
    """Gets names of contigs containing ORFs with SPs

    sp_gff      --  input, path to GFF file with SPs
    contig_ids  --  output, list of contig IDs
    """
    contig_ids = []
    with open(sp_gff, 'r') as infile:
        for line in infile:
            line = line.strip()
            if not line.startswith('#'):
                contig = line.split()[0].split('-rna')[0]
                if contig not in contig_ids:
                    contig_ids.append(contig)

    return contig_ids

def get_sp_contigs(contig_ids, genome):
    """Gets FASTA sequences of contigs containing ORFs with SPs.

    contig_ids  --  input, list of contig IDs with SPs
    genome      --  input, path to genome
    fasta_dict  --  output, dict {name:[seq]}, sequences in a genome with SP
                    containing ORFs
    """
    fasta_dict = {}
    found = False
    with open(genome, 'r') as gen:
        for line in gen:
            line = line.strip()
            if not line:
                continue # Skip empty lines
            if line.startswith(">"):
                contig = line[1:]
                if line[1:] in contig_ids:
                    fasta_dict[contig] = []
                    found = True
                else:
                    found = False
            elif found == True:
                fasta_dict[contig].append(line)

    return fasta_dict

def write_fasta(fasta_dict, fasta_dir):
    """Writes output FASTA files

    fasta_dict  --  input, dict {name:[seq]}, sequences in a genome with SP
                    containing ORFs
    fasta_dir   --  input, path to dir with FASTA files
    """
    for key, value in fasta_dict.items():
        outname = '%s/%s.fasta' % (fasta_dir, key)
        with open(outname, 'w') as outfile:
            outfile.write('>%s\n' % key)
            seq = ''.join(value)
            outfile.write('%s\n' % seq)

def main():
    """Functions to be executed in main.

    """
    sp_gff = argv[1]
    genome = argv[2]
    fasta_dir = argv[3]
    contig_ids = get_names(sp_gff)
    fasta_dict = get_sp_contigs(contig_ids, genome)
    write_fasta(fasta_dict, fasta_dir)

# Main
if __name__ == "__main__":
    main()