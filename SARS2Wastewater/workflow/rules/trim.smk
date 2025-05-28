rule ivar_trim:
    input:
        '2.align/{sample}.sort.bam'
    output:
        trmb = '2.align/{sample}.trim.bam',
        trmsb = '2.align/{sample}.trmsrt.bam'
    params:
        ivar = config['path']['ivar'],
        samtools = config['path']['samtools'],
        bed = config['database']['primer_bed']
    log:
        '.log/{sample}.ivar_trim.log'
    benchmark:
        '.log/{sample}.ivar_trim.bm'
    shell:
        """
        {params.ivar} trim -e -b {params.bed} -i {input} > {output.trmb} 2>> {log}
        {params.samtools} sort -@ 4 {output.trmb} -o {output.trmsb} &>> {log}
        {params.samtools} index {output.trmsb}
        """
