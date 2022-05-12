#!/usr/bin/env python3

"""
Author: Megan Brenes Guallar
Contact: m.brenesguallar@genetwister.nl

Description:    this script detects possible ORFs in a given genome region from
                a BED file (in this pipeline, mimp flanking regions)
Usage:  python3 01_orf_finder.py <infasta> <outpath> <min_len> <max_len>
        <gen_bed> <gff_path> <mode>
        Where:
            <infasta>   is a genomic FASTA file with the up/down stream search
                        regions
            <outpath>   is the path for the output 6 frame translation of
                        <infasta>
            <min_len>   is the minimum length in bp that an ORF must have
            <max_len>   is the maximum length in bp that an ORF must have
            <gen_bed>   is the path for the output BED file with ORFs using
                        coordinates for the original genomic FASTA input
            <gff_path>  is the path for the file with GFF lines (no header yet)
            <mode>      is a flag [3|6] to determine the translation mode (three
                        frame or six frame tranlsation)
Output:
        FASTA file with translated search region
        BED file with genomic ORF coordinates
        FASTA file with genomic ORFs
        FASTA file with protein ORFs
        preliminary GFF file (no GFF3 header line) with gene, mRNA, exon and CDS

This script is a part of the effector detection pipeline for F. oxysporum and is
based on the FoEC pipeline by van Dam et al (2016).
"""
# Imports

from sys import argv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Functions

def translate(region_fasta, mode, outpath):
    """ Translate mimp surrounding region into 3 or 6 reading frames.

    region_fasta    --  input, FASTA file of surrounding mimp region
    mode            --  input, str, [3|6], 3 or 6 frame translation
    outpath         --  input, FASTA file of translated regions
    """
    translated_records = []
    for seq_record in SeqIO.parse(region_fasta, 'fasta'):
        search_ds = False
        search_us = False
        te_desc = seq_record.description.split(':')[0]
        # Only do 3 frame translation
        if mode == '3':
            # If multiple TEs for a single search region
            if ',' in te_desc:
                te_list = te_desc.split(',')
                directions = []
                # For each TE, check if pointing US or DS
                for te_elem in te_list:
                    te_type = te_elem.split('_')[1]
                    te_dir = te_elem.split('_')[-1]
                    directions.append(te_dir)
                # Set to see if both US and DS or only one
                set_directions = set(directions)
                # If only one direction is found
                if len(set_directions) == 1:
                    # Only search upstream region
                    if list(set_directions)[0] == 'us':
                        search_us = True
                    # Only search downstream region
                    else:
                        search_ds = True
                # If both directions are found, search both US and DS regions
                else:
                    search_us = True
                    search_ds = True
            # If only one TE for the search region
            else:
                te_label = te_desc.split('_')
                te_type = te_label[1]
                te_dir = te_label[-1]
                # Search upstream
                if te_dir == 'us':
                    search_us = True
                # Search downstream
                else:
                    search_ds = True
        # For 6 frame translation
        else:
            search_us = True
            search_ds = True

        for i in range(3):
            if search_ds:
                s = seq_record.seq[i:]
                while len(s) % 3 != 0:
                    s += 'N'
                translate_s = SeqRecord(seq = s.translate(table = 1), id = \
                seq_record.id, description = seq_record.description + \
                ' frame=' + str('+%i' % (i + 1)))
                translated_records.append(translate_s)
            if search_us:
                s_rc = seq_record.seq.reverse_complement()[i:]
                while len(s_rc) % 3 != 0:
                    s_rc += 'N'
                translate_rc = SeqRecord(seq = s_rc.translate(table = 1), id = \
                seq_record.id, description = seq_record.description + \
                ' frame=' + str('-%i' % (i + 1)))
                translated_records.append(translate_rc)
    SeqIO.write(translated_records, outpath, 'fasta')

def get_orfs(translated, min_len, max_len, bedout):
    """ Finds and records ORFs

    translated  --  input, FASTA file with translated mimp flanking regions
    min_len     --  input, int, minimum protein length in aa
    max_len     --  input, int, maximum protein length in aa
    bedout      --  input, path to output BED file with putative ORFs
    bedout      --  output, BED file with putative ORFs
    orf_dict    --  output, dict {orf_name:[seq id, start, stop, strand, frame]},
                    information about found orfs
    """
    orf_dict = {}
    orf_bed = []
    i = 1
    for seq_record in SeqIO.parse(translated, 'fasta'):
        frame = seq_record.description.split('frame=')[1]
        strand = frame[0]
        frame_num = int(frame[1])
        split_rec = seq_record.id.split(':')
        region_start = int(split_rec[3].split('-')[0])
        region_end = int(split_rec[3].split('-')[1].split('(')[0])
        searchreg = region_end - region_start
        sequence = seq_record.seq
        met_loc = sequence.find('M')
        stop_loc = met_loc + 1
        while met_loc >= 0 and stop_loc >= 0:
            stop_loc = sequence.find('*', (met_loc + 1))
            prot = sequence[met_loc:(stop_loc + 1)]
            # if min_len < len(prot) < max_len and \
            if min_len < len(prot) and prot.count('X') < len(prot) * 0.5:
                start_bp = (met_loc * 3) + frame_num
                stop_bp = (stop_loc * 3) + frame_num
                extract_start = start_bp - 1
                extract_end = stop_bp - 1
                if strand == '+':
                    genomic_a = region_start + extract_start
                    genomic_b = region_start + extract_end
                if strand == '-':
                    genomic_a = region_end - extract_start
                    genomic_b = region_end - extract_end
                orf_name = 'orf_%i'%i
                i += 1
                # Make sure it is written as start < end
                genomic_start = min(genomic_a, genomic_b)
                genomic_end = max(genomic_a, genomic_b)
                orf_dict[orf_name] = [seq_record.id, genomic_start, \
                genomic_end, strand, frame_num]
                # Outfiles
                orf_bed.append([split_rec[2], str(genomic_start), \
                str(genomic_end), orf_name, '0', strand])
            # Look for next M
            met_loc = sequence.find('M', met_loc + 1)
    with open(bedout, 'w') as outfile:
        for bed_line in orf_bed:
            outfile.write('\t'.join(bed_line) + '\n')

    return orf_dict

def write_gff(orf_dict, gff_out):
    """Writes GFF rows (without header) with ORFs

    orf_dict    --  input, dict {orf_name:[seq id, start, stop, strand, frame]},
                    information about found orfs
    gff_out     --  input, path to output GFF file
    """
    # Different levels needed to represent found ORFs in a GFF file
    # {level: [name, parent]}
    levels = {
        'gene' : ['gene', None],
        'mRNA' : ['rna', 'gene'],
        'exon' : ['exon', 'mRNA'],
        'CDS' : ['cds', 'mRNA']}
    with open(gff_out, 'w') as outfile:
        outfile.write('##gff-version 3\n')
        for key, value in orf_dict.items():
            # orf 1-based coordinates
            orf_start = int(value[1]) + 1
            orf_stop = int(value[2])
            seqid = value[0].split(':')[-2]
            lvl_ids = {}
            for level, level_list in levels.items():
                lvl_label = '{}:{}'.format(level_list[0], key.split('_')[-1])
                lvl_id = '{}-{}'.format(seqid, lvl_label)
                lvl_ids[level] = lvl_id
                lvl_attr = 'ID={};Name=foec_{}'.format(lvl_id, '_'.join(\
                lvl_label.split(':')))
                if level_list[1] != None:
                    lvl_attr = '{};Parent={}'.format(lvl_attr, \
                    lvl_ids[level_list[1]])
                if level == 'CDS':
                    phase = '0'
                else:
                    phase = '.'
                lvl_info = '{}\tFOEC\t{}\t{}\t{}\t.\t{}\t{}\t{}\n'\
                .format(seqid, level, orf_start, orf_stop, value[3], \
                phase, lvl_attr)
                outfile.write(lvl_info)

def main():
    """Functions to be executed in main.

    """
    infasta = argv[1]
    outpath = argv[2]
    min_len = int(argv[3])
    max_len = int(argv[4])
    gen_bed = argv[5]
    gff_out = argv[6]
    mode = argv[7]
    translate(infasta, mode, outpath)
    orf_dict = get_orfs(outpath, min_len, max_len, gen_bed)
    write_gff(orf_dict, gff_out)

# Main
if __name__ == "__main__":
    main()
