rule nanofilt_uncmprs:
    input:
        '.rawdata/{sample}.fastq'
    output:
        '1.qc/{sample}.cln.fastq'
    params:
        NanoFilt = config['path']['NanoFilt'],
        activate = config['path']['activate'],
        prms = config['params']['NanoFilt']
    log:
        '.log/{sample}.nanofilt_uncmprs.log'
    benchmark:
        '.log/{sample}.nanofilt_uncmprs.bm'
    run:
        env = str(Path(params.NanoFilt).parents[1])
        shell("""
        source {params.activate} {env}
        {params.NanoFilt} {params.prms} {input} > {output} 2> {log}
        """)

rule nanofilt:
    input:
        '.rawdata/{sample}.fastq.gz'
    output:
        '1.qc/{sample}.cln.fastq'
    params:
        NanoFilt = config['path']['NanoFilt'],
        activate = config['path']['activate'],
        prms = config['params']['NanoFilt']
    log:
        '.log/{sample}.nanofilt.log'
    benchmark:
        '.log/{sample}.nanofilt.bm'
    run:
        env = str(Path(params.NanoFilt).parents[1])
        shell("""
        source {params.activate} {env}
        zcat {input} | NanoFilt {params.prms} > {output} 2> {log}
        """)
