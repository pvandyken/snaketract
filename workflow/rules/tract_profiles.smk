rule reformat_clusters:
    input:
        rules.separate_clusters_by_hemisphere.output,

    output:
        bids_output_dwi(
            atlas="ORG",
            suffix="clusters.tar.gz"
        )

    log: f"logs/separate_clusters_by_cluster/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/separate_clusters_by_cluster/{'.'.join(wildcards.values())}.tsv"

    envmodules:
        'python/3.7',
        "git-annex/8.20200810"

    group: "cluster_postprocess"
    resources:
        mem_mb=1500,
        runtime=15,
        tmpdir=str(work/"__sn_tmp__"),

    shell:
        boost(
            datalad,
            tar.using(inputs=["{input}"], outputs=["{output}"]),

            "tmpdir={resources.tmpdir}/reformat_clusters; "
            "mkdir -p $tmpdir/vtp-tracts; "
            "mv {input}/tracts_left_hemisphere/* $tmpdir/vtp-tracts; "
            "rename_expr='{{ "
                "number=substr($(NF), match($(NF), /[0-9]{{5}}/), 5); "
                "split($(NF), parts, number); "
                "printf \"%s output/%s%05d%s\n\", parts[1], number+offset, parts[2]"
            "}}'; "

            "find {input}/tracts_right_hemisphere/ -type f | "
            "awk -F'/' -v offset='800' -v output=$tmpdir/vtp-tracts \"$rename_expr\" | "
            "xargs -L 1 mv; "

            "find {input}/tracts_commissural/ -type f | "
            "awk -F'/' -v offset='800' -v output=$tmpdir/vtp-tracts \"$rename_expr\" | "
            "xargs -L 1 mv; "

            + Pyscript(workflow.basedir, wma_env)(
                input={"input": "{resources.tmpdir}/reformat_clusters/vtp-tracts"},
                output={"output": "{resources.tmpdir}/reformat_clusters/vtk-tracts"},
                script="scripts/convert_vtk.py",
            )

            + (
                "find {resources.tmpdir}/reformat_clusters/vtk-tracts -type f | "
                "awk -F'[./]' '{{print $0 \"{output}/\"$(NF-1)\".tck\"}}' | "
                "xargs -L 1 tckconvert"
            )
        )


rule tract_profiles:
    input:
        data=rules.reformat_clusters.output,
        ref=inputs.input_path['t1']
    output:
        shared_work/f"{uid}_tract_profiles.csv"
    log: f"logs/tract_profiles/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/tract_profiles/{'.'.join(wildcards.values())}.tsv"
    group:
    threads: 1
    resources:
        mem_mb=2000,
        runtime=30,
    params:
    shell: (
        boost(
            datalad,
            tar.using(inputs=["{input.data}"]),
            Pyscript(workflow.basedir)(
                "scripts/tract-profiling.py",
                input=["data", "ref"],
                wildcards=["subject"]
            )
        )
    )

rule aggregate_profiles:
    input:
        expand(
            rules.tract_profiles.output,
            **inputs.input_lists['preproc_dwi']
        )
    output:
        config['output_dir'] + "/tract_profiles.csv"
    log: f"logs/aggregate_profiles/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/aggregate_profiles/{'.'.join(wildcards.values())}.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=30,
    params:
    shell: (
        boost(
            Pyscript(workflow.basedir)(
                "scripts/tract-profile-merge.py"
            )
        )
    )