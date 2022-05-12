#!/usr/bin/env python3

"""
Author: Megan Brenes Guallar
Contact: m.brenesguallar@genetwister.nl

Description:    this is a script to create a GFF file including information from
                two input GFFs: the pipeline GFF and the Augustus GFF. Each
                predicted ORF (with a signal peptide) is present. If Augustus
                created a gene model for it, the Augustus feature is used. If
                not, the simple model created by the pipeline is used.
Usage:  python3 01_contruct_gff.py <aug_dir> <in_gff> <out_gff>
        Where:
            <aug_dir>   is the directory containing Augustus output GFF files
            <in_gff>    is the pipeline GFF
            <out_gff>   is the path to the output GFF

Output:
    A GFF file where both pipeline and Augustus features are present.

This script is a part of the effector detection pipeline for F. oxysporum and is
based on the FoEC pipeline by van Dam et al (2016).
"""
# Imports

import os
from sys import argv

# Functions
def parse_original_gff(in_gff):
    """Creates a dict of GFF feature lines where ORF ID is the key

    in_gff      --  input, path to input GFF file with SP ORFs
    gff_dict    --  output, dict, {id:{gene:gene_feature, mrna:mrna_feautre,
                    exon:exon_feature, cds:cds_feature,signal_peptide:sp_feature}}
    """
    gff_dict = {}
    with open(in_gff, 'r') as infile:
        next(infile)
        for line in infile:
            line = line.strip().split()
            id_num = line[-1].split(';')[0].split(':')[1]
            category = line[2].lower()
            if id_num not in gff_dict:
                gff_dict[id_num] = {'gene':'', 'mrna':'', 'exon':'', 'cds':'',
                'signal_peptide':''}
            gff_dict[id_num][category] = line

    return gff_dict

def parse_augustus(aug_dir):
    """Parses Augustus output GFF

    aug_dir     --  input, directory with Augustus output
    aug_dict    --  output, dict with Augustus predictions. {gene_id: [[f_1],
                    [f_2], [f_3]]}
    """
    aug_dict = {}
    for gene_file in os.listdir(aug_dir):
        if gene_file.endswith('.gff'):
            rna_id = gene_file.split('_')[-1].split('.')[0]
            gene_path = '%s/%s' % (aug_dir, gene_file)
            if rna_id not in aug_dict:
                aug_dict[rna_id] = {}
            with open(gene_path, 'r') as infile:
                for line in infile:
                    line = line.strip()
                    if not line.startswith('#'):
                        feature = line.split('\t')
                        if feature[2] == 'gene':
                            gene_id = feature[-1].split('ID=')[-1]
                            aug_dict[rna_id][gene_id] = [feature]
                        elif gene_id in feature[-1].split('Parent=')[-1]:
                            aug_dict[rna_id][gene_id].append(feature)

    return aug_dict

def compare_gffs(gff_dict, aug_dict):
    """Compares pipeline GFF and Augustus GFF to get final features

    gff_dict    --  input, dict with pipeline predicted features (with SPs)
    aug_dict    --  input, dict with Augustus predicted features
    keep_feats  --  output, list of lists with accepted features
    """
    keep_feats = []
    for g_key, g_item in gff_dict.items():
        og_prediction = []
        # Double check that there's an Augustus prediction for the gene
        if g_key in aug_dict:
            # If there is an Augustus prediction
            if len(aug_dict[g_key]) != 0:
                for gene, features in aug_dict[g_key].items():
                    # When the Augustus prediction starts at a predicted ORF
                    # with a signal peptide
                    aug_strand = features[0][6]
                    if aug_strand == '+':
                        aug_start = features[0][3]
                        gff_start = g_item['gene'][3]
                    elif aug_strand == '-':
                        aug_start = features[0][4]
                        gff_start = g_item['gene'][4]
                    if aug_start == gff_start:
                        gene_id = '%s-gene:%s' % (features[0][0], g_key)
                        for feature in features:
                            new_id = feature[-1].replace(gene, gene_id)
                            new_feat = feature[:-1]
                            new_feat.append(new_id)
                            if new_feat[2] == 'transcript':
                                new_feat[2] = 'mRNA'
                            keep_feats.append(new_feat)
                            og_prediction.append(False)
                        # Add signal peptide prediction to Augustus features
                        sp_feat = g_item['signal_peptide'][:-1]
                        sp_desc = '{0}.sp;Parent={0}.t1'.format(gene_id)
                        sp_feat.append(sp_desc)
                        keep_feats.append(sp_feat)
                    # When the Augstus prediction doesn't match, use original
                    else:
                        og_prediction.append(True)
            # When Augustus doesn't have a prediction, use original
            else:
                og_prediction.append(True)
        # If og_prediction is ALL True (no Augustus prediction used)
        if all(og_prediction):
            for g_feats in g_item.values():
                keep_feats.append(g_feats)

    return keep_feats


def write_gff(keep_feats, out_gff):
    """Writes GFF rows (without header) with ORFs

    keep_feats  --  input, list of lists with accepted features
    out_gff     --  input, path to output GFF file
    """
    with open(out_gff, 'w') as outfile:
        outfile.write('##gff-version 3\n')
        for feature in keep_feats:
            out_feat = '\t'.join(feature)
            outfile.write('%s\n' % out_feat)

def main():
    """Functions to be executed in main.

    """
    aug_dir = argv[1]
    in_gff = argv[2]
    out_gff = argv[3]
    gff_dict = parse_original_gff(in_gff)
    aug_dict = parse_augustus(aug_dir)
    keep_feats = compare_gffs(gff_dict, aug_dict)
    write_gff(keep_feats, out_gff)

# Main
if __name__ == "__main__":
    main()
