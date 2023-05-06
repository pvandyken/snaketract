# group: nodal_properties
#   group-components: 480

rule rich_club_coefficient:
    input:
        rules.get_connectome.output.connectome
    output:
        bids_output_dwi(
            atlas="{atlas}",
            desc="{weight}",
            suffix="richclub.csv"
        )
    log: f"logs/rich_club_coefficient/{'.'.join(wildcards.values())}.{{atlas}}.{{weight}}.log"
    benchmark: f"benchmarks/rich_club_coefficient/{'.'.join(wildcards.values())}.{{atlas}}.{{weight}}.tsv"
    group: "graph_theory"
    threads: 32
    envmodules:
        'python/3.10'
    resources:
        mem_mb=64000,
        runtime=75*5,
    params:
        filter_level="1,2,3,4,5",
        normalization=1000,
    shell:
        rich_club_env.script,
        pyscript(
            script="scripts/rich_club.py",
            params=["filter_level", "normalization"]
        )


rule nodal_properties:
    input:
        graph=rules.get_connectome.output.connectome,
        metadata=Path(config['bids_dir'], 'participants.tsv')
    output:
        bids_output_dwi(
            rec="{rec}",
            atlas="{atlas}",
            desc="{weight}",
            suffix="graphprops.tsv",
        )
    log: f"logs/nodal_properties/{'.'.join(wildcards.values())}.{{atlas}}.{{weight}}.{{rec}}.log"
    benchmark: f"benchmarks/nodal_properties/{'.'.join(wildcards.values())}.{{atlas}}.{{weight}}.{{rec}}.tsv"
    group: "nodal_properties"
    threads: 2
    resources:
        mem_mb=6000,
        runtime=6,
    envmodules:
        "python/3.8"
    params:
        subject=lambda wcards: wcards['subject'],
        weight=lambda wcards: wcards['weight'],
    shell:
        nodal_props_venv.script,
        pyscript(
            script="scripts/nodal_properties.py",
            input=['graph', 'metadata'],
            params=['subject', 'weight'],
        )

rule merge_nodal_properties:
    input:
        expand(
            rules.nodal_properties.output,
            allow_missing=True,
            weight=config["connectome_weight"],
            **inputs.input_lists['preproc_dwi'],
        )
    output:
        Path(config['output_dir'], 'graph_theory', 'rec-{rec}_atlas-{atlas}_nodes.tsv')
    threads: 1
    run:
        pd.concat(
            pd.read_csv(i, sep="\t", index_col=0) for i in input
        ).reset_index(drop=True).to_csv(str(output), sep="\t")


rule graph_theory_all_atlases:
    input:
        expand(
            rules.merge_nodal_properties.output,
            atlas=config["segmentation"],
            rec=config['tractography']['algorithm']
        )