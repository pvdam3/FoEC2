rule diamond_makedb:
# Creates a Diamond database
    input:
        prot_effs = 'output/01.findeffectors/all_putative_effectors_protein.fasta'
    output:
        'output/02.clustereffectors/diamond.dmnd'
    conda:
        '../envs/diamond.yml'
    log:
        'output/02.clustereffectors/logs/diamond_makedb.log'
    message:
        'Creating Diamond database'
    params:
        cwd = os.getcwd()
    shell:
        'diamond makedb '
        '--in {input} '
        '-d diamond '
        '--log && '
        'mv {params.cwd}/diamond.dmnd {output} && '
        'mv {params.cwd}/diamond.log {log}'

rule diamond_blastp:
# Runs BLASTP with Diamond
    input:
        prot_effs = 'output/01.findeffectors/all_putative_effectors_protein.fasta',
        dmnd_db = 'output/02.clustereffectors/diamond.dmnd'
    output:
        dmnd_tsv = 'output/02.clustereffectors/diamond_out.tsv'
    conda:
        '../envs/diamond.yml'
    log:
        'output/02.clustereffectors/logs/diamond_blastp.log'
    message:
        'Running Diamond BLASTX'
    params:
        cwd = os.getcwd()
    shell:
        'diamond blastp '
        '-q {input.prot_effs} '
        '-d {input.dmnd_db} '
        '-o {output} '
        '--log && '
        'mv {params.cwd}/diamond.log {log}'

rule mcl:
# Creates clusters with MCL
    input:
        dmnd_tsv = 'output/02.clustereffectors/diamond_out.tsv'
    output:
        mcl_out = 'output/02.clustereffectors/clusters.out'
    conda:
        '../envs/mcl.yml'
    message:
        'Running MCL to create putative effector clusters'
    params:
        mcl_i = config['mcl_i']
    shell:
        'mcl {input} '
        '--abc '
        '-I {params.mcl_i} '
        '-o {output}'

rule fasta_clusters:
# Uses MCL output to create cluster FASTA files
    input:
        gen_effs = 'output/01.findeffectors/all_putative_effectors_genomic.fasta',
        mcl_out = 'output/02.clustereffectors/clusters.out'
    output:
        directory('output/02.clustereffectors/clusters_fasta')
    message:
        'Getting cluster information'
    benchmark:
        'output/02.clustereffectors/benchmarks/part_02/part_02.benchmark.txt'
    shell:
        'python3 scripts/02_get_cluster_fasta.py {input} {output}'