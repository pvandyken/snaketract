# group: weight
#   group-components: 288 

if "t1_map" in inputs:
    rule create_r1:
        input:
            data=inputs.input_path["t1_map"],
            mask=inputs.input_path["t1_mask"]
        output:
            bids_output_anat(
                suffix="R1.nii.gz"
            )
        log: log("create_r1", inputs['preproc_dwi'])
        benchmark: benchmark("create_r1", inputs['preproc_dwi'])
        threads: 1
        envmodules:
            'python/3.10'
        resources:
            mem_mb=1000,
            runtime=1,
        group: "weight"
        shell:
            boost(
                r1_env.script,
                pyscript(
                    "scripts/produce-r1.py",
                    input=["data", "mask"]
                )
            )


rule run_dtifit:
    input:
        dwi=inputs.input_path["preproc_dwi"],
        bvals=inputs.input_path["bval"],
        bvecs=inputs.input_path["bvec"],
        brainmask=inputs.input_path["brainmask"],
    output:
        fa=bids_output_dwi(
            model="tensor",
            desc="FA",
            suffix="mdp.nii.gz",
        ),
        md=bids_output_dwi(
            model="tensor",
            desc="MD",
            suffix="mdp.nii.gz",
        ),
        rd=bids_output_dwi(
            model="tensor",
            desc="RD",
            suffix="mdp.nii.gz",
        ),
        tr=bids_output_dwi(
            model="tensor",
            desc="Tr",
            suffix="mdp.nii.gz",
        ),
        mo=bids_output_dwi(
            model="tensor",
            desc="MO",
            suffix="mdp.nii.gz",
        ),
        l1=bids_output_dwi(
            model="tensor",
            desc="L1",
            suffix="mdp.nii.gz",
        ),
        l2=bids_output_dwi(
            model="tensor",
            desc="L2",
            suffix="mdp.nii.gz",
        ),
        l3=bids_output_dwi(
            model="tensor",
            desc="L3",
            suffix="mdp.nii.gz",
        ),
        v1=bids_output_dwi(
            model="tensor",
            desc="V1",
            suffix="mdp.nii.gz",
        ),
        v2=bids_output_dwi(
            model="tensor",
            desc="V2",
            suffix="mdp.nii.gz",
        ),
        v3=bids_output_dwi(
            model="tensor",
            desc="V3",
            suffix="mdp.nii.gz",
        ),
        s0=bids_output_dwi(
            model="tensor",
            desc="S0",
            suffix="mdp.nii.gz",
        ),

    container:
        'docker://khanlab/neuroglia-core:latest'
    group:
        "weight"
    log: log("run_dtifit", inputs['preproc_dwi'])
    benchmark: benchmark("run_dtifit", inputs['preproc_dwi'])
    resources:
        runtime=10,
    shadow: 'minimal'
    shell: """
        dtifit --data={input.dwi} --bvecs={input.bvecs} --bvals={input.bvals} \
            --mask={input.brainmask} --out=dtifit
        mv dtifit_V1.nii.gz {output.v1}
        mv dtifit_V2.nii.gz {output.v2}
        mv dtifit_V3.nii.gz {output.v3}
        mv dtifit_L1.nii.gz {output.l1}
        mv dtifit_L2.nii.gz {output.l2}
        mv dtifit_L3.nii.gz {output.l3}
        mv dtifit_FA.nii.gz {output.fa}
        mv dtifit_MD.nii.gz {output.md}
        mv dtifit_MO.nii.gz {output.mo}
        mv dtifit_S0.nii.gz {output.s0}
        fslmaths {output.l1} -add {output.l2} -add {output.l3} {output.tr}
        fslmaths {output.l2} -add {output.l3} -div 2 {output.rd}
        """


def _get_image(wildcards):
    if wildcards["weight"][3:] in DIFFUSION_PARAMS:
        return rules.run_dtifit.output.fa
    if wildcards["weight"][3:] == "R1":
        return rules.create_r1.output
    raise ValueError(
        "config key 'connectome_weight' mut be set to '___FA', where ___ is one of "
        f"'avg', 'med', 'min', 'max' currently '{config['connectome_weight']}'"
    )

def _get_stat(wildcards):
    stat = wildcards["weight"][:3]
    mapping = {
        "avg": "mean",
        "med": "median",
        "min": "min",
        "max": "max",
    }
    if stat in mapping:
        return mapping[stat]
    raise ValueError(
        "config key 'connectome_weight' mut be set to '___FA', where ___ is one of "
        f"'avg', 'med', 'min', 'max' currently '{config['connectome_weight']}'"
    )


rule tck_sample:
    input:
        tracks=rules.run_act.output,
        image=_get_image
    output:
        bids_output_dwi(
            rec="{rec}",
            desc="{weight}",
            suffix='tractometry.csv'
        ),
    threads: 4
    resources:
        mem_mb=10000,
        runtime=9,
    envmodules:
        "mrtrix/3.0.1",
        "git-annex/8.20200810"
    log: log("tck_sample", inputs['preproc_dwi'], "weight", "rec")
    benchmark: benchmark("tck_sample", inputs['preproc_dwi'], "weight", "rec")
    group: 'weight'
    params:
        stat=_get_stat
    shell:
        "tcksample {input.tracks} {input.image} {output} "
        "-nthreads {threads} -stat_tck {params.stat} -q"