rule pre_fastqc_pe:
    input: 
        '.rawdata/{sample}.1.fastq.gz',
        '.rawdata/{sample}.2.fastq.gz'
    output:
        directory('1.qc/{sample}/before')
    params:
        fastqc = config['path']['fastqc'],
        prm = '--extract'
    threads:
        config['threads']['low']
    log:
        '.log/{sample}.pre_fastqc_pe.log'
    benchmark:
        '.log/{sample}.pre_fastqc_pe.bm'
    shell:
        """
        mkdir -p {output}
        {params.fastqc} -t {threads} {params.prm} -o {output} {input[0]} {input[1]} &> {log}
        """

rule pre_fastqc_se:
    input:
        '.rawdata/{sample}.fastq.gz'
    output:
        directory('1.qc/{sample}/before')
    params:
        fastqc = config['path']['fastqc'],
        prm = '--extract'
    threads:
        config['threads']['low']
    log:
        '.log/{sample}.pre_fastqc_se.log'
    benchmark:
        '.log/{sample}.pre_fastqc_se.bm'
    shell:
        """
        mkdir -p {output}
        {params.fastqc} -t {threads} {params.prm} -o {output} {input} &> {log}
        """

rule pre_fastqc_pe_uncompress:
    input: 
        '.rawdata/{sample}.1.fastq',
        '.rawdata/{sample}.2.fastq'
    output:
        directory('1.qc/{sample}/before')
    params:
        fastqc = config['path']['fastqc'],
        prm = '--extract'
    threads:
        config['threads']['low']
    log:
        '.log/{sample}.pre_fastqc_pe_uncompress.log'
    benchmark:
        '.log/{sample}.pre_fastqc_pe_uncompress.bm'
    shell:
        """
        mkdir -p {output}
        {params.fastqc} -t {threads} {params.prm} -o {output} {input[0]} {input[1]} &> {log}
        """

rule pre_fastqc_se_uncompress:
    input:
        '.rawdata/{sample}.fastq'
    output:
        directory('1.qc/{sample}/before')
    params:
        fastqc = config['path']['fastqc'],
        prm = '--extract'
    threads:
        config['threads']['low']
    log:
        '.log/{sample}.pre_fastqc_se_uncompress.log'
    benchmark:
        '.log/{sample}.pre_fastqc_se_uncompress.bm'
    shell:
        """
        mkdir -p {output}
        {params.fastqc} -t {threads} {params.prm} -o {output} {input} &> {log}
        """

rule post_fastqc_pe:
    input: 
        '1.qc/{sample}.clean.1.fastq',
        '1.qc/{sample}.clean.2.fastq'
    output:
        directory('1.qc/{sample}/after')
    params:
        fastqc = config['path']['fastqc'],
        prm = '--extract'
    threads:
        config['threads']['low']
    log:
        '.log/{sample}.post_fastqc_pe.log'
    benchmark:
        '.log/{sample}.post_fastqc_pe.bm'
    shell:
        """
        mkdir -p {output}
        {params.fastqc} -t {threads} {params.prm} -o {output} {input[0]} {input[1]} &> {log}
        """

rule post_fastqc_se:
    input:
        '1.qc/{sample}.clean.fastq'
    output:
        directory('1.qc/{sample}/after')
    params:
        fastqc = config['path']['fastqc'],
        prm = '--extract'
    threads:
        config['threads']['low']
    log:
        '.log/{sample}.post_fastqc_se.log'
    benchmark:
        '.log/{sample}.post_fastqc_se.bm'
    shell:
        """
        mkdir -p {output}
        {params.fastqc} -t {threads} {params.prm} -o {output} {input} &> {log}
        """
