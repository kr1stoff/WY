from pathlib import Path

rule freyja_variants:
    input:
        '2.align/{sample}.trmsrt.bam'
    output:
        var = '6.demix/{sample}.freyja.variants.tsv',
        depth = '6.demix/{sample}.freyja.depths.tsv'
    params:
        freyja = config['path']['freyja'],
        activate = config['path']['activate'],
        ref = config['database']['reference']
    log:
        '.log/{sample}.freyja_variants.log'
    benchmark:
        '.log/{sample}.freyja_variants.bm'
    run:
        env = str(Path(params.freyja).parents[1])
        shell("""
        set +eu
        source {params.activate} {env}
        {params.freyja} variants {input} --variants {output.var} --depths {output.depth} --ref {params.ref} &> {log}
        """)

rule freyja_demix:
    input:
        var = '6.demix/{sample}.freyja.variants.tsv',
        depth = '6.demix/{sample}.freyja.depths.tsv'
    output:
        '6.demix/{sample}.freyja.demix'
    params:
        freyja = config['path']['freyja'],
        activate = config['path']['activate'],
    log:
        '.log/{sample}.freyja_demix.log'
    benchmark:
        '.log/{sample}.freyja_demix.bm'
    run:
        env = str(Path(params.freyja).parents[1])
        shell("""
        set +eu
        source {params.activate} {env}
        {params.freyja} demix {input.var} {input.depth} --output {output} --confirmedonly &> {log}
        """)

rule abundance:
    input:
        '6.demix/{sample}.freyja.demix'
    output:
        '6.demix/{sample}.pie.png',
        '6.demix/{sample}.abundance.txt'
    benchmark:
        '.log/{sample}.abundance.bm'
    script:
        '../scripts/abundance.py'
