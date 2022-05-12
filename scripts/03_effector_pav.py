#!/usr/bin/env python3

"""
Author: Megan Brenes Guallar
Contact: m.brenesguallar@genetwister.nl

Description:    This is a script to create a presence absence variation (PAV)
                table to show which genomes have hits for detected putative
                effectors. If more than 1 hit is found, the total number of hits
                is shown.
Usage:  python3 03_effector_pav.py <genome_config> <hits> <pav_out>
        Where:
            <genome_config> is the 'genome_config.yaml' file created by the
                            pipeline
            <hits>          is a file that contains putative effector hits per
                            genome '0_genome_effector_hits.out'
            <pav_out>       is the PAV tsv outfile
Output:
    'presence_absence.tsv', tsv with PAV of putative effectors in genomes


This script is a part of the effector detection pipeline for F. oxysporum and is
based on the FoEC pipeline by van Dam et al (2016).
"""
# Imports

import os
from sys import argv

# Functions
def get_names(genome_config):
    """Creates a dict with genome files and their labels (taken from config)

    genome_config   --  input, path to genome config file
    genome_names    --  output, dict, {file: label}
    """
    genome_names = {}
    with open(genome_config, 'r') as infile:
        for line in infile:
            line = line.strip()
            if line.startswith('genomes:'):
                continue
            elif line.startswith('annotations:'):
                break
            else:
                line = line.split(': ')
                basename = line[-1].split('/')[-1]
                genome_names[basename] = line[0]

    return genome_names

def pav_dict(effector_hits, genome_names):
    """Creates PAV dict with putative effectors and genomes

    effector_hits   --  input, path to file with putative effector nhmmer hits
    genome_names    --  input, dict, {file: label}
    pav             --  output, dict, {effector: {genome: count}}
    """
    pav = {}
    with open(effector_hits, 'r') as hits_file:
        for line in hits_file:
            if line.startswith('genome'):
                continue
            else:
                line = line.strip().split()
                gen_label = genome_names[line[0]]
                eff_label = line[4]
                # If the p. effector is already in the PAV dict
                if eff_label in pav:
                    # If the genome is already in the p. effector dict, + 1 hit
                    if gen_label in pav[eff_label]:
                        pav[eff_label][gen_label] += 1
                    # If the genome is not in the p. effector dict, add it
                    else:
                        pav[eff_label][gen_label] = 1
                # If the p. effector is not already in the PAV dict, add it
                else:
                    pav[eff_label] = {gen_label: 1}

    return pav

def create_table(genome_names, pav, outfile):
    """Creates PAV table in tsv format

    genome_names    --  input, dict, {file: label}
    pav             --  input, dict, {effector: {genome: count}}
    outfile         --  input, path to output PAV tsv
    """
    effectors = pav.keys()
    genomes = list(genome_names.values())
    pav_list = []
    for genome in genomes:
        genome_list = []
        for effector in effectors:
            if genome in pav[effector]:
                genome_list.append(pav[effector][genome])
            else:
                genome_list.append(0)
        pav_list.append(genome_list)
    with open(outfile, 'w') as pav_out:
        pav_out.write('genome\t{}\n'.format('\t'.join(effectors)))
        for i, gen in enumerate(genomes):
            str_list = [str(x) for x in pav_list[i]]
            pav_out.write('{}\t{}\n'.format(gen, '\t'.join(str_list)))

def main():
    """Functions to be executed in main.

    """
    genome_config = argv[1]
    effector_hits = argv[2]
    outfile = argv[3]
    genome_names = get_names(genome_config)
    pav = pav_dict(effector_hits, genome_names)
    create_table(genome_names, pav, outfile)
# Main
if __name__ == "__main__":
    main()