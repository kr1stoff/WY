rule deconvolve:
    input:
        '3.muts/{sample}.ivar.tsv',
        '6.demix/{sample}.abundance.txt'
    output:
        '3.muts/{sample}.deconvolve.tsv',
        '3.muts/{sample}.display.tsv'
    params:
        config['database']['usher_barcodes']
    benchmark:
        '.log/{sample}.deconvolve.bm'
    script:
        '../scripts/deconvolve.py'
