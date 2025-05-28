from pathlib import Path

rule pangolin:
    input:
        '4.consensus/{sample}.consensus.fa'
    output:
        '5.lineage/{sample}.lineage_report.csv'
    params:
        pangolin = config['path']['pangolin'],
        activate = config['path']['activate']
    threads: 
        4
    log:
        '.log/{sample}.pangolin.log'
    benchmark:
        '.log/{sample}.pangolin.bm'
    run:
        ptho = Path(output[0]).resolve()
        name = ptho.name
        drt = ptho.parent
        env = str(Path(params.pangolin).parents[1])
        shell("""
        set +eu
        source {params.activate} {env}
        pangolin {input} --threads {threads} --outdir {drt} --outfile {name} &> {log}
        """)
