rule deconvolve:
    input:
        rules.snpeff.output.vcf,
        rules.abundance.output[1]
    output:
        '3.muts/{sample}.deconvolve.tsv',
        '3.muts/{sample}.display.tsv'
    params:
        config['database']['usher_barcodes']
    benchmark:
        '.log/{sample}.deconvolve.bm'
    script:
        '../scripts/deconvolve.py'
