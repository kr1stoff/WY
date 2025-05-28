rule trim:
    input:
        rules.minimap2.output
    output:
        '2.align/{sample}.trmst.bam'
    params:
        ivar = config['path']['ivar'],
        samtools = config['path']['samtools'],
        bed = config['database']['primer_bed'],
        prms = '-e -q 1'
    log:
        '.log/{sample}.trim.log'
    benchmark:
        '.log/{sample}.trim.bm'
    shell:
        """
        {params.ivar} trim {params.prms} -b {params.bed} -i {input} | {params.samtools} sort -@ 4 -o {output} - 2> {log}
        {params.samtools} index {output}
        """
