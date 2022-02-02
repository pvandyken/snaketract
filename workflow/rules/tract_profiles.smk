rule reformat_clusters:
    input:
        rules.separate_clusters_by_hemisphere.output,

    output:
        bids_output_dwi(
            atlas="T1w",
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

            "tmpdir={tmpdir}/reformat_clusters; "
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
                input={"input": "{tmpdir}/reformat_clusters/vtp-tracts"},
                output={"output": "{tmpdir}/reformat_clusters/vtk-tracts"}
            )

            + (
                "find {tmpdir}/reformat_clusters/vtk-tracts -type f | "
                "awk -F'[./]' '{{print $0 \"{output}/\"$(NF-1)\".tck\"}}' | "
                "xargs -L 1 tckconvert"
            )
        )
