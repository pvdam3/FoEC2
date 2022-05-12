#!/usr/bin/env python3

"""
Author: Megan Brenes Guallar
Contact: m.brenesguallar@genetwister.nl

Description:    this script
Usage:  python3 01_run_augustus.py <in_gff> <contig_dir> <outdir>
        Where:
            <in_gff>        is the path to the GFF file with SPs
            <contig_dir>    is the path to dir containing contig FASTA files
            <outdir>        is the path to the Augustus output directory

Output:
    Directory containing Augustus gene predictions

This script is a part of the effector detection pipeline for F. oxysporum and is
based on the FoEC pipeline by van Dam et al (2016).
"""
# Imports

import os
import subprocess
from sys import argv

# Functions
def parse_gff(in_gff):
    """Gets info about SP location in an ORF

    in_gff  --  input, path to GFF file with SPs
    sp_dict --  output, dict, {contig_name: [(sp_a_coord_1, sp_a_coord_2,
                strand_a, rna_id_a), (sp_b_coord_1, sp_b_coord_2, strand_b,
                rna_id_b)]}, SP info, strand and RNA ID
    """
    strand_dict = {'+':'forward', '-':'backward'}
    sp_dict = {}
    with open(in_gff, 'r') as gff:
        next(gff)
        for line in gff:
            line = line.strip().split()
            if line[2].startswith('sig'):
                rna_id = line[-1].split(';')[1].split('rna:')[-1]
                if line[0] not in sp_dict:
                    # Add contig name with SP info, strand and rna ID to dict
                    sp_dict[line[0]] = [(line[3], line[4], strand_dict[line[6]],
                    rna_id)]
                else:
                    sp_dict[line[0]].append((line[3], line[4],
                    strand_dict[line[6]], rna_id))

    return sp_dict

def augustus(sp_dict, contig_dir, outdir):
    """Gets info to write Augustus commands

    sp_dict     --  input, dict, SP info, strand and RNA ID
    contig_dir  --  input, path to dir containing contig FASTA files
    outdir      --  input, path to Augustus output directory
    aug_cmds    --  output, list of Augustus commands
    """
    aug_cmds = []
    for contig_file in os.listdir(contig_dir):
        if contig_file.endswith('.fasta'):
            fname = contig_file.split('.fasta')[0]
            contig_path = '%s/%s' % (contig_dir, contig_file)
            for key, value in sp_dict.items():
                for info in value:
                    strand = info[2]
                    if key == fname:
                        if strand == 'forward':
                            start = int(info[0]) - 1
                            end = int(info[1]) + 5000
                        elif strand == 'backward':
                            start = max((int(info[0]) - 5000), 0)
                            # Gives it more room to predict
                            end = int(info[1]) + 30
                        aug_cmd = 'augustus '\
                            '--species=fusarium '\
                            '--predictionStart={} '\
                            '--predictionEnd={} '\
                            '--strand={} '\
                            '--gff3=on '\
                            '--genemodel=complete '\
                            '--noInFrameStop=true '\
                            '{} > {}/augustus_{}.gff'.format(
                                start, end, info[2], contig_path, outdir,
                                info[3])
                        aug_cmds.append(aug_cmd)
    return aug_cmds

def run_cmds(aug_cmds):
    """Runs Augustus commands

    aug_cmds    --  output, list of Augustus commands
    """
    for cmd in aug_cmds:
        subprocess.run(cmd, shell = True)

def main():
    """Functions to be executed in main.

    """
    in_gff = argv[1]
    contig_dir = argv[2]
    outdir = argv[3]
    sp_dict = parse_gff(in_gff)
    aug_cmds = augustus(sp_dict, contig_dir, outdir)
    run_cmds(aug_cmds)


# Main
if __name__ == "__main__":
    main()