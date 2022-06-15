rule mafft_msa:
# Runs MAFFT to generate MSAs of clusters with more than one sequence
    input:
        in_dir = 'output/02.clustereffectors/clusters_fasta/'
    output:
        out_dir = temp(directory('output/03.presenceabsence/temp_hmm_in'))
    log:
        'output/03.presenceabsence/logs/mafft_msa.log'
    message:
        'Running MAFFT to generate MSAs...'
    conda:
        '../envs/mafft.yml'
    threads:
        4
    shell:
        """
        mkdir -p {output}
        for clust_file in {input}/*; do
            base=$(basename $clust_file .fasta)
            if [[ $base == multi* ]]; then
                mafft --thread {threads} --nuc $clust_file > {output}/$base.afa 2> {log}
            else
                cp $clust_file {output}/$base.fasta
            fi
        done
        """

def hmm_input(eff_file):
    """Determines input directory for hmmer_profiles rule

    eff_file --  input, path to file with effectors (if provided)
    """
    if os.path.exists('p_effector_check'):
            os.system('rm p_effector_check/*')
            os.system('rmdir p_effector_check')
    if eff_file == 'none':
        return 'output/03.presenceabsence/temp_hmm_in'
    else:
        os.system('mkdir -p p_effector_check')
        out_dir = 'p_effector_check'
        sep_cmd = \
        """awk '/^>/ {close(F) ; F = "%s/"substr($0,2)".fasta"} {print >> F}' %s""" \
        % (out_dir, eff_file)
        os.system(sep_cmd)
        return out_dir

rule hmmer_profiles:
# Creates HMMER profiles from clusters of putative effectors
    input:
        hmm_input(config['effectors'])
    output:
        temp(directory('output/03.presenceabsence/temp_hmm'))
    conda:
        '../envs/hmmer.yml'
    log:
        'output/03.presenceabsence/logs/hmm_profiles.log'
    message:
        'Creating HMM profiles...'
    shell:
        """
        mkdir -p {output}
        if [[ "{input}" == *p_effector_check ]]; then
            for infile in {input}/*; do
            if [[ -f "$infile" ]]; then
                new_name=$(echo "$infile" | sed 's/ /_/g')
                if [[ "$infile" == "$new_name" ]]; then
                    :
                else
                    mv "$infile" "$new_name"
                fi
            fi
        done
        fi
        for msa_file in {input}/*; do
            if [[ $msa_file == *.afa ]]; then
                hmm_file=$(basename $msa_file .afa)
            else
                hmm_file=$(basename $msa_file .fasta)
            fi
            hmm_file=$hmm_file.hmm
            hmmbuild {output}/$hmm_file $msa_file >> {log}
        done
        """

rule hmmemit:
# Emits a consensus sequence from a profile HMM
    input:
        'output/03.presenceabsence/temp_hmm'
    output:
        temp(directory('output/03.presenceabsence/putative_effectors'))
    conda:
        '../envs/hmmer.yml'
    message:
        'Generating putative effector consensus sequence...'
    params:
        # Threshold for showing weakly conserved residues as lower case [0-1]
        minl = config['minl'],
        # Threshold for showing strongly conserved residues as upper case [0-1]
        minu = config['minu'],
        # Check if effectors provided
        eff = config['effectors']
    shell:
        """
        mkdir -p {output}
        for hmm_profile in {input}/*; do
            eff_file=$(basename $hmm_profile .hmm)
            if [ "{params.eff}" = "none" ]; then
                eff_lab='p_effector'
                eff_lab+=_$(echo "$eff_file" | grep -o -E "[0-9]+$")
            else
                eff_lab=$eff_file
            fi
            hmmemit -C --minl {params.minl} --minu {params.minu} -o {output}/$eff_lab.fasta $hmm_profile
        done
        """

rule nhmmer:
# Runs nhmmer, searches for profiles in genome
    input:
        genome_file = lambda wildcards: config['genomes'][wildcards.sample],
        hmm_dir = 'output/03.presenceabsence/temp_hmm'
    output:
        nhmmer = temp(directory('output/03.presenceabsence/nhmmer_{sample}'))
    conda:
        '../envs/hmmer.yml'
    log:
        'output/03.presenceabsence/logs/nhmmer/nhmmer_{sample}.log'
    message:
        'Running nhmmer...'
    params:
        e_val = config['eval']
    threads:
        2
    benchmark:
        'output/03.presenceabsence/benchmarks/nhmmer/{sample}_nhmmer.benchmark.txt'
    shell:
        """
        hitfastas=output/03.presenceabsence/hit_fastas/
        mkdir -p {output.nhmmer}
        for cluster_profile in {input.hmm_dir}/*; do
            hmm_file=$(basename $cluster_profile .hmm)
            clust_num=$hmm_file
            hmm_file+=_{wildcards.sample}.out
            clusthit=$hitfastas/$clust_num
            mkdir -p $clusthit
            align=$clusthit/{wildcards.sample}
            nhmmer --cpu {threads} -A $align.sto -E {params.e_val} --tblout {output.nhmmer}/$hmm_file $cluster_profile {input.genome_file} >> {log}
            if [ -s "$align.sto" ]; then
                esl-reformat fasta $align.sto > $align.fa
            fi
            rm $align.sto
        done
        """

def merge_nhmmer(in_string, outdir):
    cmds = []
    in_list = list(in_string)[0]
    outdir = list(outdir)[0]
    for indir in in_list:
        cmd = 'mv {0}/* {1} && rm {0}/.snake* && rmdir {0}'.format(indir, outdir)
        cmds.append(cmd)
    cmds = ' && '.join(cmds)
    return cmds

rule merge_dirs:
# Merges nhmmer dirs
    input:
        nhmmer = expand('output/03.presenceabsence/nhmmer_{sample}', sample = config['genomes'])
    output:
        temp(directory('output/03.presenceabsence/nhmmer'))
    run:
        shell('mkdir -p {output}')
        shell(merge_nhmmer({input}, {output}))

rule parse_hmm:
# Parses output from nhmmer
    input:
        nhmmer = 'output/03.presenceabsence/nhmmer',
        hmm_dir = 'output/03.presenceabsence/temp_hmm'
    output:
        'output/03.presenceabsence/00_genome_effector_hits.out'
    message:
        'Parsing nhmmer output...'
    log:
        'output/03.presenceabsence/logs/parse_hmm.log'
    params:
        length = config['length'],
        # Check if effectors provided
        eff = config['effectors']
    shell:
        """
        if [ "{params.eff}" = "none" ]; then
            eff=0
        else
            eff=1
        fi
        python3 scripts/03_parse_hmm.py {input.nhmmer} {input.hmm_dir} {params.length} {output} $eff > {log}
        """

rule pav_table:
# Creates effector PAV table
    input:
        'output/03.presenceabsence/00_genome_effector_hits.out'
    output:
        pav = temp('output/03.presenceabsence/presence_absence.tsv')
    message:
        'Creating PAV table...'
    params:
        genome_config = 'config/genome_config.yaml'
    shell:
        'python3 scripts/03_effector_pav.py '
        '{params.genome_config} {input} {output}'

rule filter_pav:
# Filters PAV table, removes some putative effectors (likely TEs)
    input:
        pav = 'output/03.presenceabsence/presence_absence.tsv'
    output:
        sv = 'output/03.presenceabsence/01_presence_absence.tsv',
        pout = temp('output/03.presenceabsence/putative_effector_ids.txt')
    params:
        mean_thresh = config['mean_thresh'],
        individual_thresh = config['individual_thresh'],
        csv = 'config/visualization_config_effectors.csv'
    message:
        'Filtering PAV table...'
    benchmark:
        'output/03.presenceabsence/benchmarks/part_03/part_03.benchmark.txt'
    conda:
        '../envs/py38-biopython178.yml'
    shell:
        'python3 scripts/03_filter_pav.py '
        '{input.pav} {params.mean_thresh} {params.individual_thresh} {output} {params.csv}'

rule get_msas:
# Keeps MSAs of putative effector clusters that pass selection criteria
    input:
        id_list = 'output/03.presenceabsence/putative_effector_ids.txt',
        peffs = 'output/03.presenceabsence/putative_effectors'
    output:
        msa_in = temp(directory('output/03.presenceabsence/msa_input')),
        msa_out = directory('output/03.presenceabsence/putative_effector_msas'),
        fastas = directory('output/03.presenceabsence/putative_effector_fastas')
    log:
        'output/03.presenceabsence/logs/mafft_msa_final.log'
    message:
        'Running MAFFT to generate MSAs...'
    conda:
        '../envs/mafft.yml'
    threads:
        4
    params:
        effdir = config['effectors']
    shell:
        """
        fastas=output/03.presenceabsence/hit_fastas
        mkdir -p {output.msa_in}
        mkdir -p {output.msa_out}
        mkdir -p {output.fastas}
        readarray -t parray < {input.id_list}
        for clust_dir in $fastas/*; do
            if [ -z "$(ls -A $clust_dir)" ]; then
                continue
            else
                if [ "{params.effdir}" == none ]; then
                    clust_num=$(echo "${{clust_dir%/}}" | grep -o -E "[0-9]+$")
                else
                    clust_num=$(basename ${{clust_dir}})
                fi
                cat $clust_dir/* > {output.msa_in}/$clust_num.fa
            fi
        done
        for peffid in "${{parray[@]}}"; do
            for in_fa in {output.msa_in}/*; do
                catfile=$(basename $in_fa .fa)
                if [[ $catfile == $peffid ]]; then
                    mafft --thread {threads} --nuc $in_fa > {output.msa_out}/p_effector_$catfile.afa 2> {log}
                fi
            done
            for in_fasta in {input.peffs}/*; do
                fbase=$(basename $in_fasta .fasta)
                fnum=${{fbase##*p_effector_}}
                if [[ $fnum == $peffid ]]; then
                    cp $in_fasta {output.fastas}
                fi
            done
        done
        rm -rf $fastas
        """
