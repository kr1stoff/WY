rule samtools_mpileup:
    input:
        '2.align/{sample}.trmsrt.bam'
    output:
        '3.muts/{sample}.mpileup'
    params:
        samtools = config['path']['samtools'],
        prm = config['params']['samtools_mpileup'],
        ref = config['database']['reference']
    log:
        '.log/{sample}.mpileup.log'
    benchmark:
        '.log/{sample}.mpileup.bm'
    shell:
        '{params.samtools} mpileup {params.prm} --reference {params.ref} -o {output} {input} &> {log}'

rule ivar_variants:
    input:
        '3.muts/{sample}.mpileup'
    output:
        '3.muts/{sample}.ivar.tsv'
    params:
        ivar = config['path']['ivar'],
        ref = config['database']['reference'],
        gff = config['database']['genome_gff'],
        prm = config['params']['ivar_variants']
    log:
        '.log/{sample}.ivar_variants.log'
    benchmark:
        '.log/{sample}.ivar_variants.bm'
    run:
        pfx = output[0].replace('.tsv', '')
        shell('cat {input} | {params.ivar} variants -p {pfx} -g {params.gff} -r {params.ref} {params.prm} &> {log}')
