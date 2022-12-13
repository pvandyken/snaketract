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
    group: "weight"
    shell:
        r1_env.script,
        Pyscript(workflow.basedir)(
            "scripts/produce-r1.py",
            input=["data", "mask"]
        )


rule dtifit_resampled_t1w:
    input:
        dwi=inputs.input_path["preproc_dwi"],
        bvals=inputs.input_path["bval"],
        bvecs=inputs.input_path["bvec"],
        brainmask=inputs.input_path["brainmask"],
    params:
        out_basename=lambda wildcards, output: os.path.join(output.out_folder, "dti"),
    output:
        out_folder=directory(
            bids(
                root=output,
                suffix="dtifit",
                desc="eddy",
                space="T1w",
                res=config["resample_dwi"]["resample_scheme"],
                datatype="dwi",
                **inputs.subj_wildcards
            )
        ),
        out_fa=os.path.join(
            directory(
                bids(
                    root=output,
                    suffix="dtifit",
                    desc="eddy",
                    space="T1w",
                    res=config["resample_dwi"]["resample_scheme"],
                    datatype="dwi",
                    **inputs.subj_wildcards
                )
            ),
            "dti_FA.nii.gz",
        ),
    container:
        'docker://khanlab/neuroglia-core:latest'
    group:
        "weights"
    shell:
        "mkdir -p {output.out_folder} && "
        "dtifit --data={input.dwi} --bvecs={input.bvecs} --bvals={input.bvals} "
        "--mask={input.brainmask} --out={params.out_basename}"