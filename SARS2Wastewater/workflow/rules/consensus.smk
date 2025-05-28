rule bcftools_mpileup_call:
    input:
        '2.align/{sample}.trmsrt.bam'
    output:
        '4.consensus/{sample}.calls.vcf.gz'
    params:
        bcftools = config['path']['bcftools'],
        ref = config['database']['reference'],
        prm_mpileup = config['params']['bcftools_mpileup'],
        prm_call = config['params']['bcftools_call']
    log:
        '.log/{sample}.bcftools_mpileup_call.log'
    benchmark:
        '.log/{sample}.bcftools_mpileup_call.bm'
    shell:
        """
        {params.bcftools} mpileup {params.prm_mpileup} -f {params.ref} {input} | {params.bcftools} call {params.prm_call} -o {output} 2> {log}
        {params.bcftools} index {output}
        """

rule bcftools_consensus:
    input:
        '4.consensus/{sample}.calls.vcf.gz'
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
