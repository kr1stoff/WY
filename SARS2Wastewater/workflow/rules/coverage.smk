rule coverage:
    input:
        '6.demix/{sample}.freyja.depths.tsv'
    output:
        '6.demix/{sample}.coverage.txt',
        '6.demix/{sample}.coverage.png',
    benchmark:
        '.log/{sample}.coverage.bm'
    script:
        '../scripts/genome_coverage.py'
