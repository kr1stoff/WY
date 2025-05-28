rule coverage:
    input:
        rules.freyja_variants_demix.output.depth
    output:
        '6.demix/{sample}.coverage.txt',
        '6.demix/{sample}.coverage.png',
    benchmark:
        '.log/{sample}.coverage.bm'
    script:
        '../scripts/genome_coverage.py'
