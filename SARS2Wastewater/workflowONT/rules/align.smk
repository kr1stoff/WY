rule minimap2:
    input:
        '1.qc/{sample}.cln.fastq'
    output:
        '2.align/{sample}.sort.bam'
    threads: config['threads']['low']
    params:
        minimap2 = config['path']['minimap2'],
        samtools = config['path']['samtools'],
        ref = config['database']['reference'],
        prms = config['params']['minimap2']
    log:
        '.log/{sample}.minimap2.log'
    benchmark:
        '.log/{sample}.minimap2.bm'
    shell:
        """
        {params.minimap2} {params.prms} -t {threads} {params.ref} {input} \
            | {params.samtools} view -bS -F 4 - \
            | {params.samtools} sort -@ 4 -o {output} - \
            2> {log}
        {params.samtools} index {output}
        """
