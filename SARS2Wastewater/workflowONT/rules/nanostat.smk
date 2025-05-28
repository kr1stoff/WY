from pathlib import Path

rule pre_nanostat:
    input:
        '.rawdata/{sample}.fastq.gz'
    output:
        '1.qc/{sample}.raw.nanostat'
    params:
        NanoStat = config['path']['NanoStat'],
        activate = config['path']['activate']
    log:
        '.log/{sample}.nanostat.log'
    benchmark:
        '.log/{sample}.nanostat.bm'
    run:
        env = str(Path(params.NanoStat).parents[1])
        shell("""
        source {params.activate} {env}
        {params.NanoStat} --fastq {input} > {output} 2> {log}
        """)

rule pre_nanostat_uncmprs:
    input:
        '.rawdata/{sample}.fastq'
    output:
        '1.qc/{sample}.raw.nanostat'
    params:
        NanoStat = config['path']['NanoStat'],
        activate = config['path']['activate']
    log:
        '.log/{sample}.nanostat.log'
    benchmark:
        '.log/{sample}.nanostat.bm'
    run:
        env = str(Path(params.NanoStat).parents[1])
        shell("""
        source {params.activate} {env}
        {params.NanoStat} --fastq {input} > {output} 2> {log}
        """)

rule post_nanostat:
    input:
        '1.qc/{sample}.cln.fastq'
    output:
        '1.qc/{sample}.cln.nanostat'
    params:
        NanoStat = config['path']['NanoStat'],
        activate = config['path']['activate']
    log:
        '.log/{sample}.post_nanostat.log'
    benchmark:
        '.log/{sample}.post_nanostat.bm'
    run:
        env = str(Path(params.NanoStat).parents[1])
        shell("""
        source {params.activate} {env}
        {params.NanoStat} --fastq {input} > {output} 2> {log}
        """)

rule collate_nanostat:
    input:
        raw = '1.qc/{sample}.raw.nanostat',
        cln = '1.qc/{sample}.cln.nanostat'
    output:
        '1.qc/{sample}.stat.txt'
    benchmark:
        '.log/{sample}.collate_nanostat.bm'
    script:
        '../scripts/collate_nanostat.py'
