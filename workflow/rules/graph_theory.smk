rule rich_club_coefficient:
    input:
        rules.get_connectome.output.connectome
    output:
        bids_output_dwi(
            atlas="{atlas}",
            desc="{weight}",
            suffix="richclub.csv"""
        )
    log: f"logs/rich_club_coefficient/{'.'.join(wildcards.values())}.{{atlas}}.{{weight}}.log"
    benchmark: f"benchmarks/rich_club_coefficient/{'.'.join(wildcards.values())}.{{atlas}}.{{weight}}.tsv"
    group: "graph_theory"
    threads: 32
    resources:
        mem_mb=64000,
        runtime=75*5,
    params:
        filter_level="1,2,3,4,5",
        normalization=1000,
    shell:
        rich_club_env.script,
        Pyscript(workflow.basedir)(
            script="scripts/rich_club.py",
            params=["filter_level", "normalization"]
        )
