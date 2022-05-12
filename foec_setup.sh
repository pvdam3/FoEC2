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
        effector directory. (optional) A directory containing nucleotide FASTA
        file(s) with curated effector sequences (one file per sequence). These
        effectors will be searched for in the input genomes. If more than one
        sequence is found in the FASTA file, multiple files will be created
        using the header as a file name.
    -g  <path>
        genome directory (required). Directory containing genomes to run through
        the pipeline.
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
            echo "Added effector FASTA file directory: $OPTARG" >&2
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
    sed -i "s@^effectors: .*@effectors: $effectors@g" $settings_config
else
# Make sure 'none' is there if effector dir isn't specified
    sed -i "s@^effectors: .*@effectors: none@g" $settings_config
fi
if [ -n "$genome_dir" ]; then
python3 scripts/00_genomes_annotations.py $genome_dir $annotation_dir $genome_config
grep -oP '.*:' $genome_config > $vis_config && sed -i 's/://g' $vis_config
sed -i 's/ //g' $vis_config
else
echo "Missing required argument -g <genome_directory>. Non-zero exit status."
exit 2
fi

