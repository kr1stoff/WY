rule bcftools_consensus:
    input:
        rules.bgzip_index.output
    output:
        '4.consensus/{sample}.consensus.fa'
    params:
        bcftools = config['path']['bcftools'],
        ref = config['database']['reference']
    log:
        '.log/{sample}.bcftools_consensus.log'
    benchmark:
        '.log/{sample}.bcftools_consensus.bm'
    shell:
        'cat {params.ref} | {params.bcftools} consensus {input} > {output} 2> {log}'
