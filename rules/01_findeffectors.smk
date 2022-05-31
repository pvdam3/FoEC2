rule check_input:
# Check input FASTA files (must be multiline)
    input:
        genome_file = lambda wildcards: config['genomes'][wildcards.sample]
    output:
        temp('output/logs/check_input_genomes_{sample}.log')
    shell:
        """
        num_lines=$(cat {input} | wc -l | awk '{{print $1;}}')
        num_seqs=$(grep -c '>' {input})
        num_seqs=$(($num_seqs * 2))
        if [[ "$num_seqs" -eq "num_lines" ]]; then
            echo 'WARNING! Single line FASTA file detected in: {input}' > {output}
        else
            echo 'Okay! Multi line FASTA file detected in: {input}' > {output}
        fi
        """

rule combine_checks:
# Combines genome input checks into a single log
    input:
        expand('output/logs/check_input_genomes_{sample}.log', sample = config['genomes'])
    output:
        'output/logs/check_input_genomes.log'
    shell:
        'echo "## Check for single line FASTA files." > {output} && '
        'echo "## Single line FASTA files can lead to incorrect putative effector sequences." >> {output} && '
        'echo "## Please convert single line FASTA files to a multi line format to ensure correct results." >> {output} && '
        'cat {input} >> {output}'

rule find_TIRs:
# Find TIRs based on motif
    input:
        genome_file = lambda wildcards: config['genomes'][wildcards.sample]
    output:
        bed = temp('output/01.findeffectors/{sample}/{sample}_TIRs.bed')
    message:
        'Finding mimp TIRs for {wildcards.sample}'
    log:
        'output/01.findeffectors/logs/find_TIRs/find_TIRs_{sample}.log'
    conda:
        '../envs/py38-biopython178.yml'
    shell:
        'python3 scripts/01_mimpTIR_finder.py '
        '{input.genome_file} {output} {wildcards.sample} > {log}'

rule sort_bed:
# Sort bed file
    input:
        bed = 'output/01.findeffectors/{sample}/{sample}_{label_sort}.bed'
    output:
        sorted_bed = temp('output/01.findeffectors/{sample}/{sample}_{label_sort}_sorted.bed')
    message:
        'Sorting "{input}"'
    shell:
        'sort -k1,1 -k2,2n {input} > {output}'

rule get_mimp_bed:
# Identify complete and incomplete mimps from TIRs and save to bed files
    input:
        sorted_bed = 'output/01.findeffectors/{sample}/{sample}_TIRs_sorted.bed'
    output:
        mimp_bed = temp('output/01.findeffectors/{sample}/{sample}_mimps_b.bed'),
        l_bed = temp('output/01.findeffectors/{sample}/{sample}_mimps_l.bed'),
        r_bed = temp('output/01.findeffectors/{sample}/{sample}_mimps_r.bed'),
        all_bed = temp('output/01.findeffectors/{sample}/{sample}_all_mimps.bed'),
        unsort = temp('output/01.findeffectors/{sample}/annotations/{sample}_unsorted_mimps.gff'),
        gff = 'output/01.findeffectors/{sample}/annotations/{sample}_mimps.gff'
    message:
        'Classifying (in)complete mimps in {wildcards.sample}'
    shell:
        'python3 scripts/01_complete_mimps.py {input} '
        '{output.mimp_bed} {output.l_bed} {output.r_bed} {output.unsort} && '
        'cat {output.mimp_bed} {output.l_bed} {output.r_bed} > {output.all_bed} && '
        'echo "##gff-version 3" >> {output.gff} && '
        'cat {output.unsort} | sort -k1,1V -k4,4n -k5,5rn >> {output.gff}'

rule bedtools_getfasta:
# Get FASTA file from bed file
    input:
        genome_file = lambda wildcards: config['genomes'][wildcards.sample],
        in_bed = 'output/01.findeffectors/{sample}/{sample}_{label}.bed'
    output:
        out_fasta = 'output/01.findeffectors/{sample}/{sample}_{label}.fasta'
    conda:
        '../envs/bedtools.yml'
    message:
        'Writing genomic FASTA file for "{input}"'
    shell:
        'bedtools getfasta '
        '-s '
        '-name '
        '-fi {input.genome_file} '
        '-bed {input.in_bed} '
        '-fo {output}'

rule combine_mimps:
# Combine all complete mimps into a single FASTA file
    input:
        expand('output/01.findeffectors/{sample}/{sample}_00_complete_mimps.fasta', \
        sample = config['genomes']),
    output:
        'output/01.findeffectors/all_complete_mimps.fasta'
    message:
        'Writing "{output}"'
    shell:
        'cat {input} > {output} && '
        'rm features_level*.json' # Remove JSON from AGAT when tool is done

rule combine_fasta:
# Combine all putative effector genomic/protein FASTAs into a single FASTA
    input:
        eff_pro = expand(\
        'output/01.findeffectors/{sample}/{sample}_08_putative_effectors_protein.fasta', \
        sample = config['genomes']),
        eff_gen = expand(\
        'output/01.findeffectors/{sample}/{sample}_09_putative_effectors_genomic.fasta', \
        sample = config['genomes'])
    output:
        eff_pro = 'output/01.findeffectors/all_putative_effectors_protein.fasta',
        eff_gen = 'output/01.findeffectors/all_putative_effectors_genomic.fasta'
    message:
        'Writing genomic and protein FASTA files: '
        '"{output.eff_pro}" and "{output.eff_gen}"'
    shell:
        'cat {input.eff_pro} > {output.eff_pro} && '
        'cat {input.eff_gen} > {output.eff_gen}'

rule combine_bed:
# Combines left, right and both bed files into a single bed
    input:
        both = 'output/01.findeffectors/{sample}/{sample}_mimps_flanked_b.bed',
        left = 'output/01.findeffectors/{sample}/{sample}_mimps_flanked_l.bed',
        right = 'output/01.findeffectors/{sample}/{sample}_mimps_flanked_r.bed'
    output:
        temp('output/01.findeffectors/{sample}/{sample}_all_mimps_flanked.bed')
    message:
        'Combining BED files for {wildcards.sample}'
    shell:
        'cat {input} > {output}'

rule samtools_faidx:
# Create index for FASTA file
    input:
        genome_file = '{genome}'
    output:
        fai = temp('{genome}.fai')
    message:
        'Creating FASTA index file for "{input}"'
    conda:
        '../envs/samtools.yml'
    shell:
        'samtools faidx {input} --fai-idx {output}'

def get_side(flank):
    """ Determine which side to use for Bedtools flank.

    flank       --  input, str, b|l|r, which flag is needed
    flank_side  --  output, str, params needed for Bedtools flank
    """
    if flank == 'b':
        flank_side = '-b'
    # Flank number determined by config
    elif flank == 'l':
        flank_side = '-r 0 -l'
    elif flank == 'r':
        flank_side = '-l 0 -r'

    return flank_side

rule bedtools_flank:
# Create a new bed file with flanking regions of input bed file
    input:
        in_fai = lambda wildcards: config['genomes'][wildcards.sample] + '.fai',
        mimp_bed = 'output/01.findeffectors/{sample}/{sample}_{label}_{flank}.bed'
    output:
        temp('output/01.findeffectors/{sample}/{sample}_{label}_flanked_{flank}.bed')
    wildcard_constraints:
        flank = '[l|r|b]'
    conda:
        '../envs/bedtools.yml'
    message:
        'Getting flanking regions for "{input.mimp_bed}"'
    params:
        size = config['search_dist'],
        flank_side = lambda wildcards: get_side(wildcards.flank)
    shell:
        'bedtools flank '
        '-i {input.mimp_bed} '
        '-g {input.in_fai} '
        '{params.flank_side} '
        '{params.size} > {output}'

rule bedtools_merge:
# Combines overlapping elements in a bed file
    input:
        'output/01.findeffectors/{sample}/{sample}_all_mimps_streamed_sorted.bed'
    output:
       temp('output/01.findeffectors/{sample}/{sample}_mimps_flanked_merged.bed')
    conda:
        '../envs/bedtools.yml'
    message:
        'Combining overlapping elements in "{input}"'
    params:
        mode = config['mode']
    shell:
        """
        if [[ "{params.mode}" = "3" ]]; then
            bedtools merge -s -c 4 -o collapse -i {input} > {output}
        elif [[ "{params.mode}" = "6" ]]; then
            bedtools merge -c 4 -o collapse -i {input} > {output}
        fi
        """

rule bedtools_closest:
# Finds closest features in A and B
    input:
        orfs = 'output/01.findeffectors/{sample}/{sample}_orfs_genomic_sorted.bed',
        both_bed = 'output/01.findeffectors/{sample}/{sample}_mimps_b.bed',
        l_bed = 'output/01.findeffectors/{sample}/{sample}_mimps_l.bed',
        r_bed = 'output/01.findeffectors/{sample}/{sample}_mimps_r.bed'
    output:
        temp('output/01.findeffectors/{sample}/{sample}_closest_mimps_expanded.tsv')
    conda:
        '../envs/bedtools.yml'
    message:
        'Finding closest mimps to putative ORFs for {wildcards.sample}'
    shell:
    # -D to report distances (+ or - depending on down or upstream)
    # -io to ignore overlap, only get *closest*
        'bedtools closest '
        '-a {input.orfs} '
        '-b {input.both_bed} {input.l_bed} {input.r_bed} > {output} -io -D ref'

rule reduce_closest:
# Reduces closest mimp output to only one mimp per ORF
    input:
        'output/01.findeffectors/{sample}/{sample}_closest_mimps_expanded.tsv'
    output:
        'output/01.findeffectors/{sample}/{sample}_closest_mimps.tsv'
    log:
        'output/01.findeffectors/logs/reduce_closest/reduce_closest_{sample}.log'
    message:
        'Determining closest mimps upstream to putative ORFs for {wildcards.sample}'
    conda:
        '../envs/py38-biopython178.yml'
    shell:
        'python3 scripts/01_reduce_closest.py {input} {output} {wildcards.sample} '
        '> {log}'

rule add_stream:
# Add if an element is up or downstream to names in bed file
    input:
        mimp_flanked_bed = 'output/01.findeffectors/{sample}/{sample}_{label}_flanked.bed'
    output:
        mimp_reg_bed = temp('output/01.findeffectors/{sample}/{sample}_{label}_streamed.bed')
    message:
        'Labeling regions up/downstream to mimps for {wildcards.sample}'
    shell:
        'python3 scripts/01_add_stream.py {input} {output}'

rule get_orfs:
# Identify ORFs in FASTA file and get genomic/protein sequences
    input:
        merged = 'output/01.findeffectors/{sample}/{sample}_mimps_flanked_merged.fasta'
    output:
        translation = 'output/01.findeffectors/{sample}/{sample}_translated_flanks.fasta',
        put_orfs = 'output/01.findeffectors/{sample}/{sample}_orfs_genomic.bed',
        orf_gff = temp('output/01.findeffectors/{sample}/annotations/{sample}_orfs.gff')
    message:
        'Finding putative ORF coordinates for {wildcards.sample}'
    params:
        min_len = config['min_len'],
        max_len = config['max_len'],
        mode = config['mode']
    conda:
        '../envs/py38-biopython178.yml'
    shell:
        'python3 scripts/01_orf_finder.py {input.merged} {output.translation} '
        '{params.min_len} {params.max_len} {output.put_orfs} {output.orf_gff} '
        '{params.mode}'

rule bedtools_intersect:
    input:
        user_gff = lambda wildcards: config['annotations'][wildcards.sample],
        merged_flanks = 'output/01.findeffectors/{sample}/{sample}_mimps_flanked_merged.bed'
    output:
        temp('output/01.findeffectors/{sample}/{sample}_reduced_user.gff')
    conda:
        '../envs/bedtools.yml'
    message:
        'Reducing input GFF "{input.user_gff}" to regions surronding detected mimps for {wildcards.sample}'
    shell:
        'bedtools intersect '
        '-a {input.user_gff} '
        '-b {input.merged_flanks} '
        '-u > {output}'

def gff_input(wildcards):
    """Adds user GFF file to mix if provided

    """
    inputs = ['output/01.findeffectors/{sample}/annotations/{sample}_orfs.gff']
    try:
        gff_file = config['annotations'][wildcards.sample].split('/')[-1]
        reduced_gff = 'output/01.findeffectors/{sample}/{sample}_reduced_user.gff'
        inputs.append(reduced_gff)
    except:
        pass

    return inputs

rule gffread_extract_genomic:
# Extract CDS sequences from GFF and output genomic FASTA
    input:
        genome_fasta = lambda wildcards: config['genomes'][wildcards.sample],
        gff = 'output/01.findeffectors/{sample}/annotations/{sample}_putative_{tag}.gff'
    output:
        genomic_fasta = temp('output/01.findeffectors/{sample}/{sample}_putative_{tag}_genomic.fasta')
    message:
        'Writing genomic FASTA for "{input.gff}"'
    conda:
        '../envs/gffread.yml'
    shell:
        'gffread '
        '-x {output.genomic_fasta} '
        '-g {input}'

rule gffread_extract_protein:
# Extract CDS sequences from GFF and output protein FASTA
    input:
        genome_fasta = lambda wildcards: config['genomes'][wildcards.sample],
        gff = 'output/01.findeffectors/{sample}/annotations/{sample}_putative_{tag}.gff'
    output:
        protein_fasta = temp('output/01.findeffectors/{sample}/{sample}_putative_{tag}_protein.fasta'),
    message:
        'Writing genomic and protein FASTA for "{input.gff}"'
    conda:
        '../envs/gffread.yml'
    shell:
        'gffread '
        '-y {output.protein_fasta} '
        '-g {input}'

rule final_filter:
# Filters out sequences by size, cysteine threshold
    input:
        prot = 'output/01.findeffectors/{sample}/{sample}_putative_tofilter_protein.fasta'
    output:
        prot = temp('output/01.findeffectors/{sample}/{sample}_putative_final_protein_filtered.fasta'),
        keep = temp('output/01.findeffectors/{sample}/{sample}_putative_final_protein_filtered_keeplist.txt')
    message:
        'Removing sequences which do not meet minimum size and cysteine content requirements...'
    params:
        cys_thresh = config['cysteine_threshold'],
        min_len = config['min_len'],
        max_len = config['max_len']
    conda:
        '../envs/py38-biopython178.yml'
    shell:
        'python3 scripts/01_cysteine_content.py '
        '{input.prot} {params.cys_thresh} {params.min_len} {params.max_len} '
        '{output.prot} {output.keep}'

rule run_signalp:
# Runs SignalP on ORFs
    input:
        orf_fasta = 'output/01.findeffectors/{sample}/{sample}_putative_orfs_protein.fasta'
    output:
        sp_gff = temp('output/01.findeffectors/{sample}/annotations/{sample}_sp.gff3'),
        sp_summary = temp('output/01.findeffectors/{sample}/annotations/{sample}_sp_summary.txt'),
        keep_list = temp('output/01.findeffectors/{sample}/annotations/{sample}_keep_list.txt')
    log:
        run_sp = 'output/01.findeffectors/logs/run_signalp/run_signalp_{sample}.log'
    message:
        'Running SignalP version {params.version} on "{input}"'
    params:
        cwd = os.getcwd(),
        signalp = config['progs']['signalp']['path'],
        batch = config['progs']['signalp']['batch'],
        version = config['progs']['signalp']['version']
    benchmark:
        'output/01.findeffectors/benchmarks/signalp/{sample}_signalp.benchmark.txt'
    threads:
        2
    run:
        if params.version == 5:
            shell('echo "Sample: {wildcards.sample}" > {log} && '
            '{params.signalp} '
            '-fasta {input} '
            '-format short '
            '-gff3 '
            '-org euk '
            '-prefix {wildcards.sample}_sp '
            '-batch {params.batch} >> {log} && '
            'mv {params.cwd}/{wildcards.sample}_sp.gff3 {output.sp_gff} && '
            'mv {params.cwd}/{wildcards.sample}_sp_summary.signalp5 {output.sp_summary} && '
            'awk -F "\\t" "/^##/ {{next}}; {{print \$1}}" {output.sp_gff} > {output.keep_list}'
            )
        elif params.version == 4:
            shell('{params.signalp} '
            '-f short -v -l {log}'
            '-n temp_spv4_{wildcards.sample}.gff '
            '{input} >> {output.sp_summary} && '
            'grep -v \'##\' temp_spv4_{wildcards.sample}.gff | '
            'awk \'{{OFS="\\t"}}{{$3="signal_peptide"; $7=$8=$9="."; print $0}}\' > '
            '{output.sp_gff} && '
            'sed -i \'1s/^/##gff-version 3\\n/\' {output.sp_gff} && '
            'awk -F "\\t" "/^##/ {{next}}; {{print \$1}}" {output.sp_gff} > {output.keep_list}'
            )

rule process_signalp:
# Processes SignalP output
    input:
        sp_gff = 'output/01.findeffectors/{sample}/annotations/{sample}_sp.gff3',
        orfs_gff = 'output/01.findeffectors/{sample}/annotations/{sample}_putative_orfs.gff',
        sp_summary = 'output/01.findeffectors/{sample}/annotations/{sample}_sp_summary.txt'
    output:
        sp_out = temp('output/01.findeffectors/{sample}/annotations/{sample}_genomic_sp.gff')
    log:
        'output/01.findeffectors/logs/process_signalp/process_signalp_{sample}.log'
    message:
        'Processing SignalP output for {wildcards.sample}'
    params:
        version = config['progs']['signalp']['version']
    shell:
        'python3 scripts/01_process_signalp.py '
        '{input.sp_gff} {input.orfs_gff} {input.sp_summary} {output} {wildcards.sample} {params.version} > {log} && '
        'cat {input.sp_summary} >> {log}'

rule agat_merge:
# Merges multiple GFF files (user input and pipeline orfs)
    input:
        foec_gff = 'output/01.findeffectors/{sample}/annotations/{sample}_orfs.gff',
        gffs = gff_input
    output:
        temp('output/01.findeffectors/{sample}/annotations/{sample}_putative_orfs_tosort.gff')
    log:
        'output/01.findeffectors/logs/agat_merge/agat_merge_{sample}.log'
    message:
        'Merging input GFF files for {wildcards.sample}'
    conda:
        '../envs/agat.yml'
    shell:
        """
        if [ ! -f features_level3.json ]; then
            agat_convert_sp_gxf2gxf.pl --expose
            sed -i 's/"sig_peptide":"exon",/"signal_peptide":"mrna",/g' features_level3.json
        fi
        in_gffs=({input.gffs})
        if [ "${{#in_gffs[@]}}" -gt 1 ]; then
            instr=( "${{in_gffs[@]:0:1}}" "-f" "${{in_gffs[@]:1}}" )
            agat_sp_merge_annotations.pl -f ${{instr[@]}} -o {output} > {log}
        else
            mv {input.foec_gff} {output}
            echo 'No GFF files provided for {wildcards.sample}' > {log}
        fi
        """

rule agat_filter_keep_list:
# Filters GFF based on keep list
    input:
        combined_orfs_gff = 'output/01.findeffectors/{sample}/annotations/{sample}_putative_orfs.gff',
        keep_list = 'output/01.findeffectors/{sample}/annotations/{sample}_keep_list.txt'
    output:
        temp('output/01.findeffectors/{sample}/annotations/{sample}_filtered.gff'),
        temp('output/01.findeffectors/{sample}/annotations/{sample}_filtered_report.txt')
    log:
        'output/01.findeffectors/logs/agat_filter/agat_filter_{sample}.log'
    message:
        'Filtering GFF file for proteins with a signal peptide for {wildcards.sample}'
    conda:
        '../envs/agat.yml'
    shell:
        'agat_sp_filter_feature_from_keep_list.pl '
        '--gff {input.combined_orfs_gff} '
        '--keep_list {input.keep_list} '
        '--output {output} > {log}'

rule agat_filter_keep_list_final:
# Filters final GFF based on keep list
    input:
        gff = 'output/01.findeffectors/{sample}/annotations/{sample}_putative_tofilter.gff',
        keep = 'output/01.findeffectors/{sample}/{sample}_putative_final_protein_filtered_keeplist.txt'
    output:
        temp('output/01.findeffectors/{sample}/annotations/{sample}_putative_finalfiltered.gff'),
    log:
        'output/01.findeffectors/logs/agat_filter/agat_filter_final_{sample}.log'
    message:
        'Filtering GFF file for {wildcards.sample} (final filter)'
    conda:
        '../envs/agat.yml'
    shell:
        'agat_sp_filter_feature_from_keep_list.pl '
        '--gff {input.gff} '
        '--keep_list {input.keep} '
        '--output {output} > {log}'

rule agat_sort:
# Sorts and fixes GFF files
    input:
        gff = 'output/01.findeffectors/{sample}/annotations/{sample}_{description}_tosort.gff'
    output:
        sorted_gff = temp('output/01.findeffectors/{sample}/annotations/{sample}_{description}.gff')
    log:
        sort_log = 'output/01.findeffectors/logs/agat_sort/agat_sort_{sample}_{description}.log'
    message:
        'Sorting GFF file "{input}"'
    conda:
        '../envs/agat.yml'
    params:
        cwd = os.getcwd(),
        log_name = '{sample}_{description}_tosort.agat.log'
    shell:
        'agat_convert_sp_gxf2gxf.pl '
        '-g {input} '
        '-o {output} > {log} && '
        'rm {params.cwd}/{params.log_name}'

rule combine_gffs:
# Combines SP predicted features (with genomic coords) and filtered GFF
    input:
        sp_orfs = 'output/01.findeffectors/{sample}/annotations/{sample}_filtered.gff',
        sps = 'output/01.findeffectors/{sample}/annotations/{sample}_genomic_sp.gff'
    output:
        combined = temp('output/01.findeffectors/{sample}/annotations/{sample}_sp_orfs_tosort.gff')
    message:
        'Combining signal peptide and putative ORF GFF files for {wildcards.sample}'
    shell:
        'awk "FNR==1 && NR!=1 {{next}}; {{print}}" {input} > {output}'

rule remove_nonsp:
# Removes non SP features
    input:
        all_orfs = 'output/01.findeffectors/{sample}/annotations/{sample}_sp_orfs.gff'
    output:
        reduced_sp = temp('output/01.findeffectors/{sample}/annotations/{sample}_putative_effectors_tosort.gff')
    message:
        'Removing features unrelated to signal peptides for {wildcards.sample}'
    shell:
        'python3 scripts/01_reduce_to_sp.py {input} {output} && '
        'rm {input}'

rule get_sp_contigs:
# Gets contigs containing ORFs with SPs for Augustus prediction
    input:
        sp_gff = 'output/01.findeffectors/{sample}/annotations/{sample}_sp.gff3',
        genome_file = lambda wildcards: config['genomes'][wildcards.sample]
    output:
        outfasta = temp(directory('output/01.findeffectors/{sample}/contigs_{sample}/'))
    shell:
        'mkdir -p {output} && '
        'python3 scripts/01_sp_contigs.py {input} {output}'

rule run_augustus:
# Run Augustus per SP containing ORF
    input:
        gff = 'output/01.findeffectors/{sample}/annotations/{sample}_putative_effectors.gff',
        contigs = 'output/01.findeffectors/{sample}/contigs_{sample}/'
    output:
        aug_dir = temp(directory('output/01.findeffectors/{sample}/augustus_out/'))
    conda:
        '../envs/augustus.yml'
    benchmark:
        'output/01.findeffectors/benchmarks/augustus/{sample}_augustus.benchmark.txt'
    shell:
        'mkdir -p {output} && '
        'python3 scripts/01_run_augustus.py {input} {output}'

rule merge_augustus:
# Merges Augustus GFF + pipeline GFF (only includes valid Augustus predictions)
    input:
        p_gff = 'output/01.findeffectors/{sample}/annotations/{sample}_putative_effectors.gff',
        aug_dir = 'output/01.findeffectors/{sample}/augustus_out/'
    output:
        construct_gff = temp('output/01.findeffectors/{sample}/annotations/{sample}_constructed.gff'),
        merge_loci = temp('output/01.findeffectors/{sample}/annotations/{sample}_mergedloci.gff'),
        out_gff = temp('output/01.findeffectors/{sample}/annotations/{sample}_putative_tofilter_tosort.gff')
    log:
        gxf = 'output/01.findeffectors/logs/merge_augustus/agat_sort_{sample}.log',
        li = 'output/01.findeffectors/logs/merge_augustus/agat_longest_isoform_{sample}.log'
    params:
        cwd = os.getcwd()
    conda:
        '../envs/agat.yml'
    shell:
        'python3 scripts/01_construct_gff.py '
        '{input.aug_dir} '
        '{input.p_gff} '
        '{output.construct_gff} && '
        'agat_convert_sp_gxf2gxf.pl '
        '-g {output.construct_gff} '
        '-ml '
        '-o {output.merge_loci} > {log.gxf} && '
        'agat_sp_keep_longest_isoform.pl '
        '-gff {output.merge_loci} '
        '-o {output.out_gff} && '
        'mv {params.cwd}/{wildcards.sample}_constructed.agat.log {log.li}'

rule rename_output:
    input:
        comp_mimp_fasta = 'output/01.findeffectors/{sample}/{sample}_mimps_b.fasta',
        mimp_merged_fasta = 'output/01.findeffectors/{sample}/{sample}_mimps_flanked_merged.fasta',
        mimp_merged_trans = 'output/01.findeffectors/{sample}/{sample}_translated_flanks.fasta',
        orfs_genomic_bed = 'output/01.findeffectors/{sample}/{sample}_orfs_genomic.bed',
        orfs_genomic_fasta = 'output/01.findeffectors/{sample}/{sample}_orfs_genomic.fasta',
        closest_bed = 'output/01.findeffectors/{sample}/{sample}_closest_mimps.tsv',
        pro_fasta_orf = 'output/01.findeffectors/{sample}/{sample}_putative_orfs_protein.fasta',
        gen_fasta_orf = 'output/01.findeffectors/{sample}/{sample}_putative_orfs_genomic.fasta',
        pro_fasta_eff = 'output/01.findeffectors/{sample}/{sample}_putative_final_protein_filtered.fasta',
        gen_fasta_eff = 'output/01.findeffectors/{sample}/{sample}_putative_finalfiltered_genomic.fasta',
        mimp_gff = 'output/01.findeffectors/{sample}/annotations/{sample}_mimps.gff',
        orf_gff = 'output/01.findeffectors/{sample}/annotations/{sample}_putative_orfs.gff',
        sp_gff = 'output/01.findeffectors/{sample}/annotations/{sample}_sp.gff3',
        final_gff = 'output/01.findeffectors/{sample}/annotations/{sample}_putative_finalfiltered.gff'
    output:
        comp_mimp_fasta = 'output/01.findeffectors/{sample}/{sample}_00_complete_mimps.fasta',
        mimp_merged_fasta = 'output/01.findeffectors/{sample}/{sample}_01_mimps_downstream.fasta',
        mimp_merged_trans = 'output/01.findeffectors/{sample}/{sample}_02_translated_mimps_downstream.fasta',
        orfs_genomic_bed = 'output/01.findeffectors/{sample}/{sample}_03_foec_orfs_genomic.bed',
        orfs_genomic_fasta = 'output/01.findeffectors/{sample}/{sample}_04_foec_orfs_genomic.fasta',
        closest_bed = 'output/01.findeffectors/{sample}/{sample}_05_closest_mimps.tsv',
        pro_fasta_orf = 'output/01.findeffectors/{sample}/{sample}_06_putative_orfs_protein.fasta',
        gen_fasta_orf = 'output/01.findeffectors/{sample}/{sample}_07_putative_orfs_genomic.fasta',
        pro_fasta_eff = 'output/01.findeffectors/{sample}/{sample}_08_putative_effectors_protein.fasta',
        gen_fasta_eff = 'output/01.findeffectors/{sample}/{sample}_09_putative_effectors_genomic.fasta',
        mimp_gff = 'output/01.findeffectors/{sample}/annotations/{sample}_00_mimps.gff',
        orf_gff = 'output/01.findeffectors/{sample}/annotations/{sample}_01_orfs.gff',
        sp_gff = 'output/01.findeffectors/{sample}/annotations/{sample}_02_signalp.gff3',
        final_gff = 'output/01.findeffectors/{sample}/annotations/{sample}_03_putative_effectors.gff'

    message:
        'Renaming output files to improve organization for {wildcards.sample}'
    benchmark:
        'output/01.findeffectors/benchmarks/part_01/{sample}_part_1.benchmark.txt'
    run:
        inputs = list({input})[0]
        outputs = list({output})[0]
        zipped = zip(inputs, outputs)
        for (infile, outfile) in zipped:
            cmd = 'mv {} {}'.format(infile, outfile)
            subprocess.run(cmd, shell = True)
