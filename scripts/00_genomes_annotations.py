#!/usr/bin/env python3

"""
Author: Megan Brenes Guallar
Contact: m.brenesguallar@genetwister.nl

Description:    this script creates a YAML file with input genomes and their
                annotations. Genomes and annotations must have the same names
                (excluding extensions).
Usage:  python3 00_genomes_annotations.py <gen_dir> <ann_dir> <yaml>
        Where:
            <gen_dir>   directory containing genomes
            <ann_dir>   directory containing annotations
            <yaml>      path to output YAML config file
Output:
        YAML config file for genomes and annotations

This script is a part of the effector detection pipeline for F. oxysporum and is
based on the FoEC pipeline by van Dam et al (2016).
"""
# Imports

import pathlib
import textwrap
from sys import argv, exit

# Functions

def create_dict(gen_dir, ann_dir):
    """Creates a dictionary with genome and annotation file paths

    gen_dir         --  input, path to directory with genome files
    ann_dir         --  input, path to directory with annotation files
    gen_ann_dict    --  output, dict of lists, {genome_name:[genome], annotation}
    """
    fasta_exts = ('.fasta', '.fa', '.fna', '.faa')
    # Genomes in dir
    gen_num = 0
    # Annotations in dir
    ann_num = 0
    # Track annotations added
    anns_added = 0
    # Initialize dictionary which will contain root name as key and a list value
    # List contains genome path and annotation path if available
    gen_ann_dict = {}
    gen_paths = pathlib.Path(gen_dir).glob('**/*')
    # Cycle through genomes in directory
    for gen_filepath in gen_paths:
        abs_gen_path = gen_filepath.absolute()
        # Skip index files if present
        if str(abs_gen_path).endswith(fasta_exts):
            gen_num += 1
            gen_stem = abs_gen_path.stem
            gen_ann_dict[gen_stem] = [str(abs_gen_path)]
            if ' ' in gen_stem or '*' in gen_stem:
                print('*WARNING! Invalid file name "%s". Please do not use spaces or asterisks.' % gen_stem)
                exit()
    if ann_dir != 'none':
        ann_paths = pathlib.Path(ann_dir).glob('**/*')
        # Cycle through annotations in directory
        for ann_filepath in ann_paths:
            ann_num += 1
            abs_ann_path = ann_filepath.absolute()
            # Only take GFF files (.gff or .gff3)
            if str(abs_ann_path).endswith('.gff') or \
            str(abs_ann_path).endswith('.gff3'):
                ann_stem = abs_ann_path.stem
                # Check if stem already in dictionary
                if ann_stem in gen_ann_dict:
                    gen_ann_dict[ann_stem].append(str(abs_ann_path))
                    anns_added += 1
                # If not, no genome is associated with the annotation. Do not add.
                else:
                    print('*WARNING! No genome file found for \'{}\'.'.format(abs_ann_path))
            else:
                print('*WARNING! Invalid annotation file \'{}\'. Only GFF3 format accepted.'.\
                format(abs_ann_path))
    print('Genomes found in \'{}\': {}'.format(gen_dir, gen_num))
    print('Genomes added: %i' % len(gen_ann_dict.keys()))
    if ann_dir != 'none':
        print('Annotations found in \'{}\': {}'.format(ann_dir, ann_num))
        print('Annotations added: %i' % anns_added)

    return gen_ann_dict

def create_yaml(gen_ann_dict, ann_dir, outpath):
    """Writes YAML file with genome and annotation information

    gen_ann_dict    --  input, dict of lists, {genome_name:[genome], annotation}
    ann_dir         --  input, path to directory with annotation files
    outpath         --  input, path to output YAML
    """
    with open(outpath, 'w') as yaml:
        yaml.write('genomes:\n')
        for key, value in gen_ann_dict.items():
            yaml.write('    {}: {}\n'.format(key, value[0]))
        if ann_dir != 'none':
            yaml.write('annotations:\n')
            for key, value in gen_ann_dict.items():
                if len(value) == 2:
                    yaml.write('    {}: {}\n'.format(key, value[1]))

def main():
    """Functions to be executed in main.

    """
    print('*' * 80)
    print(' ' * 27 + '00_genomes_annotations.py' + ' ' * 27)
    print('*' * 80)
    description = \
    'This script creates two config files: genome_config.yaml and '\
    'visualization_config.csv. '\
    'The file genome_config.yaml contains all genome and annotation paths. These '\
    'paths are preceeded by labels which are used to refer to a given sample. i.e. '\
    '"GCF_XXXXX" can be refered to as "Fol4287". '\
    'The file visualization_config.csv contains a column with genome labels. More '\
    'columns can be added to describe the genomes (i.e. forma specialis, sample '\
    'location, year, etc.).'
    print(textwrap.fill(description, 80) + '\n')
    gen_dir = argv[1]
    ann_dir = argv[2]
    outpath = argv[3]
    gen_ann_dict = create_dict(gen_dir, ann_dir)
    create_yaml(gen_ann_dict, ann_dir, outpath)
    print('\nConfig files created!')
    print('\n*** GENOME CONFIG ***')
    outmessage = \
    'Labels can be changed in "config/genome_config.yaml". If '\
    'annotations are used, make sure that the label for the annotation is the '\
    'same as the label for the genome. ' \
    'TIP! If using find and replace, refactoring or any similar functionality ' \
    'to change labels, please include the colon (":") to avoid changing the file ' \
    'paths.'
    print('\n' + textwrap.fill(outmessage, 80))
    print('\n*** VISUALIZATION CONFIG ***')
    vis_message = \
    'The visualization config file "config/visualization_config.csv" can be used '\
    'to provide information about your genomes. This will help group your genomes '\
    'in downstream visualization steps. Edit the file to add information about '\
    'your genomes. For example:'
    print('\n' + textwrap.fill(vis_message, 80))
    print('genomes\t\tf.sp.\ngen_a.fasta\tcubense\ngen_b.fasta\tlilii')
    next_step = \
    'Continue by running the snakemake pipeline with:\n\tsnakemake --use-conda '\
    '--cores [N]\n--OR for macOS users--\n\tsnakemake --use-singularity '\
    '--use-conda --cores [N]\nTo see more options, use:\n\tsnakemake -h'
    print('\n' + next_step)

# Main
if __name__ == "__main__":
    main()