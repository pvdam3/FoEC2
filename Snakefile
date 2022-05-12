configfile: 'config/genome_config.yaml'
configfile: 'config/config.yaml'

import subprocess

eff_in = config['effectors']
def get_input(eff_in):
    infiles = []
    if eff_in == 'none':
        infiles.append('output/logs/check_input_genomes.log')
        infiles.append('output/01.findeffectors/all_complete_mimps.fasta')
        infiles.append('output/01.findeffectors/all_putative_effectors_protein.fasta')
        for infile in expand('output/01.findeffectors/{sample}/{sample}_00_complete_mimps.fasta', sample = config['genomes']):
            infiles.append(infile)
        for indir in expand('output/01.findeffectors/{sample}/augustus_out/', sample = config['genomes']):
            infiles.append(indir)
        infiles.append('output/02.clustereffectors/clusters_fasta')
    infiles.append('output/03.presenceabsence/putative_effector_fastas')
    infiles.append('output/03.presenceabsence/putative_effector_msas')
    infiles.append('output/03.presenceabsence/01_presence_absence.tsv')

    return infiles

rule all:
    input:
        get_input(eff_in)

include: './rules/01_findeffectors.smk'
include: './rules/02_clustereffectors.smk'
include: './rules/03_presenceabsence.smk'
