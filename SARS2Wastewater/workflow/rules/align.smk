rule bowtie2_pe:
    input: 
        '1.qc/{sample}.clean.1.fastq',
        '1.qc/{sample}.clean.2.fastq'
    output: 
        sam = '2.align/{sample}.sam'
    params: 
        bowtie2 = config['path']['bowtie2'],
        bti = config['database']['reference']
    threads:
        config['threads']['low']
    log:
        '.log/{sample}.bowtie2_pe.log'
    benchmark:
        '.log/{sample}.bowtie2_pe.bm'    
    shell: 
        """
        {params.bowtie2} --no-unal --threads {threads} -x {params.bti} -1 {input[0]} -2 {input[1]} -S {output.sam} &> {log}
        """

rule bowtie2_se:
    input: 
        '1.qc/{sample}.clean.fastq'
    output: 
        sam = '2.align/{sample}.sam'
    params: 
        bowtie2 = config['path']['bowtie2'],
        bti = config['database']['reference']
    threads:
        config['threads']['low']
    log:
        '.log/{sample}.bowtie2_se.log'
    benchmark:
        '.log/{sample}.bowtie2_se.bm'
    shell: 
        """
        {params.bowtie2} --no-unal --threads {threads} -x {params.bti} -U {input} -S {output.sam} &> {log}
        """

rule sam_to_sort_bam:
    input:
        '2.align/{sample}.sam'
    output:
        bam = '2.align/{sample}.bam',
        stbm = '2.align/{sample}.sort.bam'
    params: 
        samtools = config['path']['samtools']
    threads:
        4
    log:
        '.log/{sample}.sam_to_sort_bam.log'
    benchmark:
        '.log/{sample}.sam_to_sort_bam.bm'
    shell:
        """
        {params.samtools} view -@ {threads} -hbS {input} -o {output.bam} &> {log}
        {params.samtools} sort -@ {threads} {output.bam} -o {output.stbm} &>> {log}
        """
