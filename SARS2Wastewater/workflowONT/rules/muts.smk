rule lofreq:
    input:
        rules.trim.output
    output:
        '3.muts/{sample}.raw.vcf'
    params:
        lofreq = config['path']['lofreq'],
        samtools = config['path']['samtools'],
        ref = config['database']['reference'],
        prms = config['params']['lofreq']
    log:
        '.log/{sample}.lofreq.log'
    benchmark:
        '.log/{sample}.lofreq.bm'
    shell:
        """
        {params.lofreq} call {params.prms} -f {params.ref} -o {output} {input} 2> {log}
        """

rule snpeff:
    input:
        rules.lofreq.output
    output:
        vcf = '3.muts/{sample}.snpeff.vcf',
        html = '3.muts/{sample}.snpeff.html',
        csv = '3.muts/{sample}.snpeff.csv'
    params:
        snpEff = config['path']['snpEff']
    log:
        '.log/{sample}.snpeff.log'
    benchmark:
        '.log/{sample}.snpeff.bm'
    shell:
        """
        {params.snpEff} NC_045512.2 -htmlStats {output.html} -csvStats {output.csv} {input} > {output.vcf} 2> {log}
        """

rule bgzip_index:
    input:
        rules.lofreq.output
    output:
        '3.muts/{sample}.raw.vcf.gz'
    params:
        bgzip = config['path']['bgzip'],
        bcftools = config['path']['bcftools']
    log:
        '.log/{sample}.bgzip_index.log'
    benchmark:
        '.log/{sample}.bgzip_index.bm'
    shell:
        """
        {params.bgzip} -c {input} > {output} 2> {log}
        {params.bcftools} index {output}
        """
