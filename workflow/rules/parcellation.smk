rule get_hemispheric_tracts:
    input:
        data=rules.separate_clusters_by_hemisphere.output,
        index=config["atlases"]["hemispheric_tracts"]

    output:
        temp(directory(work/"get_hemispheric_tracts"/uid/"{hemi}"))

    log: f"logs/get_hemispheric_tracts/{'.'.join(inputs.input_wildcards['surf'].values())}.log"
    benchmark: f"benchmarks/get_hemispheric_tracts/{'.'.join(inputs.input_wildcards['surf'].values())}.tsv"
    group: "parcellation"
    envmodules:
        "mrtrix/3.0.1",
        "python/3.7",
        "git-annex/8.20200810"

    threads: 1
    resources:
        mem_mb=500,
        runtime=10,

    shell:
        boost(
            datalad.msg(""),
            tar.using(inputs = ["{input.data}"]),
            convert_env.make_venv,
            (
                hemi := sh.ShVar(
                    sh.ShIf("{wildcards.hemi}").eq("L") >> (
                        sh.echo("tracts_left_hemisphere").n()
                    ) >> (
                        sh.echo("tracts_right_hemisphere").n()
                    ),
                    export=True
                ),
                "mkdir -p {output}.tmp",
                sh.ShFor(cluster:=sh.ShVar(), _in="`cat {input.index}`") >> (
                    sh.mv(f"{{input.data}}/{hemi}/{cluster}", f"{{output}}.tmp/{cluster}"),
                ),
                Pyscript(workflow.basedir, python_path=convert_env.python_path)(
                    input={"input": f"{{output}}.tmp/\*.vtp"},
                    output={"output": f"{{output}}/\*.tck"},
                    script="scripts/convert_vtk.py",
                ),
            )
        )


rule get_parcellation:
    input:
        tracts=rules.get_hemispheric_tracts.output,
        mesh=inputs.input_path['surf'],
    output:
        atlas=bids_output_dwi(
            atlas="ORG",
            suffix="parcellation.vtk",
            **inputs.input_wildcards['surf'],
        ),
        connectome=bids_output_dwi(
            atlas="ORG",
            suffix="connectome.pyc",
            **inputs.input_wildcards['surf'],
        ),

    log: f"logs/get_hemisphere_parcellation/{'.'.join(inputs.input_wildcards['surf'].values())}.log"
    benchmark: f"benchmarks/get_hemisphere_parcellation/{'.'.join(inputs.input_wildcards['surf'].values())}.tsv"
    group: "parcellation"
    envmodules:
        "git-annex/8.20200810",
        "python/3.9"


    threads: 8
    resources:
        mem_mb=3000,
        runtime=28,
    shell:
        boost(
            datalad.msg(""),
            parcellation_env.script,

            "intersection main "
            "{input.mesh} {input.tracts} "
            "{output.atlas} {output.connectome} --threads {threads}"
        )
