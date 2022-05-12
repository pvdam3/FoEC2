#!/usr/bin/env python3

"""
Author: Megan Brenes Guallar
Contact: m.brenesguallar@genetwister.nl

Description:    This script creates FASTA files for clusters determined by MCL.
                Files are named differently depending on if there is more than
                one sequence in the FASTA file ("multicluster") or not
                ("cluster"). This facilitates downstream processing.
Usage:  python3 02_get_cluster_fasta.py <protein_fasta> <mcl_out> <outdir>
        Where:
            <protein_fasta> is a protein FASTA file with all putative effectors
            <mcl_out>       is the output from MCL
            <outdir>        is the path to the output directory
Output:
    A directory with FASTA files, one per MCL cluster.

This script is a part of the effector detection pipeline for F. oxysporum and is
based on the FoEC pipeline by van Dam et al (2016).
"""
# Imports

import os
from sys import argv

# Functions
def parse_fasta(prots):
    """Parses FASTA file and puts sequences into a dict

    prots       --  input, protein sequences from all effectors
    prot_dict   --  output, dict, {seq_name: seq}
    """
    prot_dict = {}
    with open(prots, 'r') as prot_in:
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

def clust_dict(prot_dict, clust):
    """Creates a dict of effector clusters and their proteins sequences

    prot_dict   --  input, dict, {seq_name: seq}
    clust       --  input, mcl output
    clust_dict  --  output, dict, {cluster: [seq_name, seq]}
    """
    clust_dict = {}
    i = 0
    with open(clust, 'r') as clust_in:
        for line in clust_in:
            clust_line = line.strip().split('\t')
            i += 1
            clust_lab = 'cluster_%i' % i
            clust_dict[clust_lab] = []
            for prot_lab in clust_line:
                clust_dict[clust_lab].append([prot_lab, prot_dict[prot_lab]])

    return clust_dict

def write_fasta(clust_dict, outdir):
    """Writes FASTA files per cluster

    clust_dict  --  input, dict, {cluster: [seq_name, seq]}
    outdir      --  input, path to output directory
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for clust_id, clust_info in clust_dict.items():
        if len(clust_info) > 1:
            filename = outdir + '/multi%s.fasta' % clust_id
        else:
            filename = outdir + '/%s.fasta' % clust_id
        with open(filename, 'w') as out_clust:
            for seq in clust_info:
                out_clust.write('>%s\n' % seq[0])
                out_clust.write('%s\n' % seq[1])


def main():
    """Functions to be executed in main.

    """
    prots = argv[1]
    clust = argv[2]
    outdir = argv[3]
    prot_dict = parse_fasta(prots)
    c_dict = clust_dict(prot_dict, clust)
    write_fasta(c_dict, outdir)

# Main
if __name__ == "__main__":
    main()