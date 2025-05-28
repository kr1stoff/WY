from pathlib import Path

rule freyja_variants_demix:
    input:
        rules.trim.output
    output:
        var = '6.demix/{sample}.freyja.variants.tsv',
        depth = '6.demix/{sample}.freyja.depths.tsv',
        demix = '6.demix/{sample}.freyja.demix'
    params:
        freyja = config['path']['freyja'],
        activate = config['path']['activate'],
        ref = config['database']['reference']
    log:
        '.log/{sample}.freyja_variants_demix.log'
    benchmark:
        '.log/{sample}.freyja_variants_demix.bm'
    run:
        env = str(Path(params.freyja).parents[1])
        shell("""
        source {params.activate} {env}
        {params.freyja} variants {input} --variants {output.var} --depths {output.depth} --ref {params.ref} &> {log}
        {params.freyja} demix {output.var} {output.depth} --output {output.demix} --confirmedonly &>> {log}
        """)

rule abundance:
    input:
        rules.freyja_variants_demix.output.demix
    output:
        '6.demix/{sample}.pie.png',
        '6.demix/{sample}.abundance.txt'
    benchmark:
        '.log/{sample}.abundance.bm'
    script:
        '../scripts/abundance.py'
