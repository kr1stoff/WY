rule pre_nanopore_qc:
    input:
        '.rawdata/{sample}.fastq.gz'
    output:
        '1.qc/{sample}.raw.stats.txt',
        '1.qc/{sample}.raw.stats.png'
    benchmark:
        '.log/{sample}.pre_nanopore_qc.bm'
    script:
        '../scripts/nanopore_fastq_quality.py'

rule pre_nanopore_qc_uncmprs:
    input:
        '.rawdata/{sample}.fastq'
    output:
        '1.qc/{sample}.raw.stats.txt',
        '1.qc/{sample}.raw.stats.png'
    benchmark:
        '.log/{sample}.pre_nanopore_qc.bm'
    script:
        '../scripts/nanopore_fastq_quality.py'


rule post_nanopore_qc:
    input:
        '1.qc/{sample}.cln.fastq'
    output:
        '1.qc/{sample}.cln.stats.txt',
        '1.qc/{sample}.cln.stats.png'
    benchmark:
        '.log/{sample}.post_nanopore_qc.bm'
    script:
        '../scripts/nanopore_fastq_quality.py'
