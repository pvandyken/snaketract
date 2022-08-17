rename_awk_expr = (
    "number=substr($(NF), match($(NF), /[0-9]{{5}}/), 5)",
    "split($(NF), parts, number)",
    'printf "%s "output"/%s%05d%s\\n", $0, parts[1], number+offset, parts[2]'
)


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
        "git-annex/8.20200810",
        "mrtrix"

    group: "profiling"
    resources:
        mem_mb=1000,
        runtime=10,
        tmpdir=str(work/"__sn_tmp__"),

    shell:
        # datalad,
        tar.using(inputs=["{input}"], outputs=["{output}"]),
        wma_env.make_venv,
        (
            tmpdir := sh.ShVar(
                "{resources.tmpdir}/reformat_clusters/{wildcards.subject}"
            ),
            sh.ShTry(
                vtp_dir := sh.ShVar(str(tmpdir)+"/vtp-tracts"),
                sh.mkdir(vtp_dir).p,
                sh.mv("{input}/tracts_left_hemisphere/*", vtp_dir),

                sh.find("{input}/tracts_right_hemisphere/ -type f") |
                sh.awk(*rename_awk_expr).F('/').v(
                    offset='800', output=vtp_dir
                ) |
                "xargs -L 1 mv",

                sh.find("{input}/tracts_commissural/ -type f") |
                sh.awk(*rename_awk_expr).F('/').v(
                    offset='1600', output=vtp_dir
                ) |
                "xargs -L 1 mv",

                Pyscript(workflow.basedir)(
                    input={"input": vtp_dir},
                    output={"output": str(tmpdir)+"/vtk-tracts"},
                    script="scripts/convert_vtk.py",
                    python_path=wma_env.python_path,
                ),

                sh.find(str(tmpdir)+"/vtk-tracts -type f") |
                sh.awk('print $0 " {output}/"$(NF-1)".tck"').F('[./]') |
                "xargs -L 1 tckconvert"
            ).catch(
                f"rm {tmpdir} -rf",
                "false"
            ).els(
                f"rm {tmpdir} -rf"
            )
        )


rule create_r1:
    input:
        data=inputs.input_path["t1_map"],
        mask=inputs.input_path["t1_mask"]
    output:
        bids_output_anat(
            suffix="R1.nii.gz"
        )
    log: f"logs/create_r1/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/create_r1/{'.'.join(wildcards.values())}.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=1,
        tmpdir=str(work/"__sn_tmp__"),
    group: "profiling"
    shell:
        Pyscript(workflow.basedir)(
            "scripts/produce-r1.py",
            input=["data", "mask"]
        )


rule tract_profiles:
    input:
        data=rules.reformat_clusters.output,
        ref=inputs.input_path['t1'],
        r1=rules.create_r1.output,
        fa=bids_output_dwi(model="CSD", suffix="FA.nii.gz")
    output:
        temp(shared_work/f"{uid}_tract_profiles.csv")
    log: f"logs/tract_profiles/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/tract_profiles/{'.'.join(wildcards.values())}.tsv"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=660,
        tmpdir=str(work/"__sn_tmp__"),
    params:
    group: "profiling"
    shell:
        # datalad,
        tar.using(inputs=["{input.data}"]),
        dipy_env.script,
        Pyscript(workflow.basedir)(
            "scripts/tract-profiling.py",
            input=["data", "ref", "r1", "fa"],
            wildcards=["subject"]
        )

rule aggregate_profiles:
    input:
        expand(
            rules.tract_profiles.output,
            **inputs.input_lists['preproc_dwi']
        )
    output:
        config['output_dir'] + "/tract_profiles.csv"
    log: f"logs/aggregate_profiles/all.log"
    benchmark: f"benchmarks/aggregate_profiles/all.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=30,
    shell:
        dipy_env.script,
        Pyscript(workflow.basedir)(
            "scripts/tract-profile-merge.py"
        )
