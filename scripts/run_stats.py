#!/usr/bin/env python3
"""Runs some stats on output directory generated by the FoEC2 pipeline.

python3 run_stats.py <output_directory>
"""

import os
import pandas as pd
import re
from sys import argv

def orf_detection(output_dir):
    info_list = []
    find_TIRs_dir = '%s/01.findeffectors/logs/find_TIRs' % output_dir
    run_signalp_log = '%s/01.findeffectors/logs/run_signalp' % output_dir
    clust_path = '%s/02.clustereffectors/clusters.out' % output_dir
    hmm_log = '%s/03.presenceabsence/logs/parse_hmm.log' % output_dir
    pav_tsv = '%s/03.presenceabsence/01_presence_absence.tsv' % output_dir
    for tir_log in os.listdir(find_TIRs_dir):
        m = re.search('find_TIRs_(.*).log', tir_log)
        genome = m.group(1)
        genome_dir = '%s/01.findeffectors/%s' % (output_dir, genome)
        sp_path = '%s/run_signalp_%s.log' % (run_signalp_log, genome)
        tir_path = '%s/%s' % (find_TIRs_dir, tir_log)
        with open(tir_path, 'r') as infile:
            for line in infile:
                line = line.strip()
                if line.startswith('// Total sequences'):
                    tot_seqs = line.split(': ')[-1]
                elif line.startswith('// Sequences containing'):
                    seqs_containing = line.split(': ')[-1]
                elif line.startswith('// Total motifs'):
                    tot_motifs = line.split(': ')[-1]
        with open(sp_path, 'r') as spfile:
            for line in spfile:
                line = line.strip()
                if line.startswith('Total proteins'):
                    after_cys = line.split(': ')[-1][:-1]
        with open(clust_path, 'r') as clustfile:
            all_clust = 0
            single_clust = 0
            for line in clustfile:
                if line != '\n':
                    all_clust += 1
                    if len(line.split()) == 1:
                        single_clust += 1
        with open(hmm_log, 'r') as hmmfile:
            hit_count = 0
            unique_hits = 0
            prev_line = False
            for line in hmmfile:
                line = line.strip()
                if line.startswith('// Parsing file'):
                    gen_compare = '_'.join(line.split('/')[-1].split('_')[2:]).split('.out')[0]
                    if gen_compare == genome:
                        prev_line = True
                elif line.endswith('hits found.') and prev_line == True:
                    num = line.split()[1]
                    if num != 'No':
                        hit_count += int(line.split()[1])
                        unique_hits += 1
                        prev_line = False
        with open(pav_tsv, 'r') as pav:
            header_line = pav.readline()
            num_peffs = len(header_line.split()) - 1
        for genout in os.listdir(genome_dir):
            genpath = '%s/%s' % (genome_dir, genout)
            if genpath.endswith('_00_complete_mimps.fasta'):
                mimp_count = 0
                with open(genpath, 'r') as mimpfile:
                    for line in mimpfile:
                        line = line.strip()
                        if line.startswith('>'):
                            mimp_count += 1
            elif genpath.endswith('_03_foec_orfs_genomic.bed'):
                orf_count = 0
                with open(genpath, 'r') as orffile:
                    for line in orffile:
                        if line != '\n':
                            orf_count += 1
            elif genpath.endswith('_08_putative_effectors_protein.fasta'):
                sp_count = 0
                with open(genpath, 'r') as spfile:
                    for line in spfile:
                        if line.startswith('>'):
                            sp_count += 1

        info_list.append([genome, tot_seqs, seqs_containing, tot_motifs, \
        mimp_count, orf_count, after_cys, sp_count, all_clust, single_clust, \
        hit_count, unique_hits, num_peffs])

    myCols = [
        'genome',
        'num_seqs',
        'num_tir_seqs',
        'num_tirs',
        'num_mimps',
        'num_orfs',
        'num_orfs_after_cysteine',
        'num_orfs_after_signalp',
        'num_clusters',
        'num_single_clusters',
        'nhmmer_total_hits',
        'nhmmer_cluster_hits',
        'final_peffectors'
    ]
    stats_df = pd.DataFrame(info_list, columns = myCols)
    stats_df.to_csv('./stats_report.csv', index = False)

def main():
    """Functions to be executed in main.

    """
    output_dir = argv[1]
    orf_detection(output_dir)

# Main
if __name__ == "__main__":
    main()