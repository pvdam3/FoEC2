#!/usr/bin/env python3

"""
Author: Megan Brenes Guallar
Contact: m.brenesguallar@genetwister.nl

Description:    This script parses the output from nhmmer (using the --tblout
                format).
Usage:  python3 03_parse_hmmer.py <hmm_dir> <profile_dir> <len_thresh> <hit_out>
        <eff>
        Where:
            <hmm_dir>       is a directory containing all nhmmer hit files
            <profile_dir>   is a directory containing all effector HMM profiles
            <len_thresh>    is the length threshold for a match [0-1]
            <hit_out>       is the hits output file
            <eff>           is a flag to search for a given list of effectors or
                            not  (0 = normal mode, 1 = search effectors mode)
Output:
    '00_genome_effector_hits.out', tsv (with a header) with effector hits in
    genomes

This script is a part of the effector detection pipeline for F. oxysporum and is
based on the FoEC pipeline by van Dam et al (2016).
"""
# Imports

import os
from datetime import datetime
from sys import argv

# Functions
def profile_length(profile_dir):
    """Parses out HMM profile length

    profile_length  --  input, path to directory with HMM profiles
    profile_dict    --  output, dict with profile name and length {name:len}
    """
    profile_dict = {}
    for hmm in os.listdir(profile_dir):
        full_hmm = '%s/%s' % (profile_dir, hmm)
        with open(full_hmm, 'r') as in_profile:
            for line in in_profile:
                line = line.strip()
                if line.startswith('NAME'):
                    hmm_name = line.split()[-1]
                    profile_dict[hmm_name] = ''
                if line.startswith('LENG'):
                    hmm_len = line.split()[-1]
                    profile_dict[hmm_name] = hmm_len

    return profile_dict

def parse_hmmer(hmm_out, profile_dict, len_thresh):
    """Parses nhmmer output

    hmm_out         --  input, nhmmer file
    profile_dict    --  input, dict with profile name and length {name:len}
    hit_lines       --  output, list of lists, hit lines meeting E value and
                        length thresholds
    """
    temp_lines = []
    hit_lines = []
    with open(hmm_out, 'r') as hmm_file:
        print('// Parsing file: %s' % hmm_out)
        found_hit = 0
        for line in hmm_file:
            line = line.strip()
            if not line.startswith('#'):
                # Last field may contain spaces!!
                line_info = line.split()[0:15]
                description = '_'.join(line.split()[15:])
                line_info.append(description)
                temp_lines.append(line_info)
                found_hit += 1
            if line.startswith('# Target file:'):
                target = line.split()[-1].split('/')[-1]
        if found_hit == 0:
            print('// No hits found.')
        else:
            print('// %i hits found.' % found_hit)
    for temp_line in temp_lines:
        if len(temp_line) == 16:
            temp_line.append(target)
            match_len = (abs(int(temp_line[5]) - int(temp_line[4]))) / int(profile_dict[temp_line[2]])
            print('// %s query coverage: %i' % (temp_line[0], match_len*100))
            if match_len >= float(len_thresh):
                hit_lines.append(temp_line)
        else:
            print('// Corrupted line! Hit not added to output:\n\t%s' % '\t'.join(temp_line))

    return hit_lines

def run_dir(hmmer_dir, profile_dict, len_thresh):
    """Runs hmmer parser for all files in a dir, creates list of output lines

    hmmer_dir       --  input, path to directory with nhmmer output
    profile_dict    --  input, dict with profile name and length {name:len}
    len_thresh      --  input, float, length threshold for a match [0-1]
    all_lines       --  output, list of lists, each list corresponds to a line
                        of a valid match from nhmmer output (valid = meeting
                        E-value and length threshold requirements)
    """
    all_lines = []
    for hmmer_out in os.listdir(hmmer_dir):
        hmmer_out = '%s/%s' % (hmmer_dir, hmmer_out)
        hit_lines = parse_hmmer(hmmer_out, profile_dict, len_thresh)
        for hit_line in hit_lines:
            all_lines.append(hit_line)

    return all_lines

def write_hits(all_lines, hit_out, eff):
    """Writes output file with summary of putative effector hits

    all_lines   --  input, list of lists, all hit lines meeting E value and
                    length thresholds
    hit_out     --  input, path to hit summary output file
    eff         --  input, binary, 0 = normal mode, 1 = search effectors mode
    """
    print('// Writing output file: %s' % hit_out)
    with open(hit_out, 'w') as outfile:
        outfile.write('genome\thit_name\thit_start\thit_stop\teffector\te_start\te_stop\tEVAL\tdescription\n')
        for hit in all_lines:
            if int(eff) == 0:
                eff_name = 'p_effector_' + hit[2].split('_')[-1]
            else:
                eff_name = hit[2]
            # genome, hit name, hit start, hit stop, effector name, effector start,
            # effector stop, E-value hit, description
            line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(\
            hit[16], hit[0], hit[6], hit[7], eff_name, hit[4], hit[5], hit[12], hit[15])
            outfile.write(line)

def main():
    """Functions to be executed in main.

    """
    start_time = datetime.now()
    print('*' * 80)
    print(' ' * 34 + 'parse_hmm.py' + ' ' * 34)
    print('*' * 80)
    hmm_dir = argv[1]
    profile_dir = argv[2]
    len_thresh = argv[3]
    hit_out = argv[4]
    eff = argv[5]
    profile_dict = profile_length(profile_dir)
    all_lines = run_dir(hmm_dir, profile_dict, len_thresh)
    write_hits(all_lines, hit_out, eff)
    total_time = datetime.now() - start_time
    print('// Time needed: %s' % total_time)
# Main
if __name__ == "__main__":
    main()