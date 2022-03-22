LEFT_HEM = "tracts_left_hemisphere"
RIGHT_HEM = "tracts_right_hemisphere"
rule get_hemispheric_tracts:
    input:
        data=rules.separate_clusters_by_hemisphere.output,
        index=config["atlases"]["hemispheric_tracts"]

    output:
        temp(directory(work/"get_hemispheric_tracts"/uid))

    log: f"logs/get_hemispheric_tracts/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/get_hemispheric_tracts/{'.'.join(wildcards.values())}.tsv"
    group: "parcellation"

    threads: 1
    resources:
        mem_mb=1000,
        runtime=30,

    shell: (
        boost(
            datalad.msg(""),
            tar.using(inputs = ["{input.data}"]),
            (
                sh.ShFor(cluster:=sh.ShVar(), _in="`cat {input.index}`").do(
                    f"cp {{input.data}}/{LEFT_HEM}/{cluster} "
                    f"{{output}}/{LEFT_HEM}/{cluster}",

                    f"cp {{input.data}}/{RIGHT_HEM}/{cluster} "
                    f"{{output}}/{RIGHT_HEM}/{cluster}"
                )
            )
        )
    )

rule get_parcellation:
    input:
        tracts=rules.get_hemispheric_tracts.output,
        left_mesh=config["surfaces"]["left"],
        right_mesh=config["surfaces"]["right"],
    output:
        bids_output_dwi(
            atlas="ORG",
            suffix="parcellation.vtk"
        )

    log: f"logs/get_parcellation/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/get_parcellation/{'.'.join(wildcards.values())}.tsv"
    group: "parcellation"

    threads: 8
    resources:
        mem_mb=1000,
        runtime=30,
    shell: (
        boost(
            datalad.msg(""),
            parcellation_env.script,

            "intersection main "
            "{input.left_mesh} {input.right_mesh}"
            f"{{input.tracts}}/{LEFT_HEM} {{input.tracts}}/{RIGHT_HEM} "
            "{output} --threads {threads}"
        )
    )