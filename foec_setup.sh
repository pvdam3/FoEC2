#!/bin/bash
# Create config files
settings_config='config/config.yaml'
genome_config='config/genome_config.yaml'
vis_config='config/visualization_config.csv'
# Get opts
annotation_dir='none'
make_config='false'
effectors=()
print_usage() {
    printf "
Welcome to the setup script for the Fusarium oxysporum Effector Clustering
pipeline: FoEC2! This script should be run to create the config files needed to
run the pipeline.

Usage $0:
    -a  <path>
        annotation directory (optional). Directory containing annotations for
        genome files. Must have the same root name as files present in genome
        directory and be in GFF3 format.
    -e  <path>
        effector file (optional). A FASTA file containing curated effector 
        nucleotide sequences. These effectors will be searched for in the input
        genomes.
    -g  <path>
        genome directory (required). Directory which only contains genomes in
        FASTA format to run through the pipeline.
    -h
        help. Shows this message.
"
}

while getopts ':a:e:g:h?' flag; do
    case "${flag}" in
        a)
            annotation_dir=$OPTARG
            echo "Using annotation directory: $annotation_dir" >&2
            ;;
        e)
            effectors=("$OPTARG")
            echo "Using effector FASTA file directory: $OPTARG" >&2
            ;;
        g)
            genome_dir=$OPTARG
            echo "Using genome directory: $genome_dir" >&2
            ;;
        \? | h | *)
            print_usage
            exit 1
            ;;
        esac
    done
if [ ${#effectors[@]} -gt 0 ]; then
    if [ -d $effectors ]; then
        echo "ERROR! A directory was provided to '-e' instead of a file."
        exit 21
    elif [ -f $effectors ]; then
        sed -i "s@^effectors: .*@effectors: $effectors@g" $settings_config
    else
        echo "ERROR! Unrecognized input provided to '-e'. "
        exit 1
    fi
else
# Make sure 'none' is there if effector dir isn't specified
    sed -i "s@^effectors: .*@effectors: none@g" $settings_config
fi
if [ -n "$genome_dir" ]; then
    if [ -d $genome_dir ]; then
        python3 scripts/00_genomes_annotations.py $genome_dir $annotation_dir $genome_config
        grep -oP '.*:' $genome_config > $vis_config && sed -i 's/://g' $vis_config
        sed -i 's/ //g' $vis_config
    elif [ -f $genome_dir ]; then
        echo "ERROR! A file was provided to '-g' instead of a directory."
        exit 20
    fi
else
    echo "Missing required argument -g <genome_directory>. Non-zero exit status."
    exit 2
fi

