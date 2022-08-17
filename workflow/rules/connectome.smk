# group: connectome:
#       connected_components: 96
def translate_hemi(path):
    def inner(wcards):
        wcards.hemi = wcards.hemi.lower() + "h"
        return str(path).format(**wcards)
    return inner

rule map_brainnetome_atlas_cortex:
    input:
        atlas=translate_hemi(
            Path(config["atlases"]['brainnetome'], "{hemi}.BN_Atlas.gcs")
        ),
        sphere=translate_hemi(
            Path(config["freesurfer_output"], "sub-{subject}/surf/{hemi}.sphere.reg"),
        ),
        surf=translate_hemi(
            Path(config["freesurfer_output"], "sub-{subject}/surf/{hemi}.white")
        ),
    output:
        bids_output_anat(
            atlas="bn210",
            suffix="dparc.label.gii",
            **inputs.input_wildcards["surf"],
        )
    shadow: "minimal"
    group: "connectome"
    envmodules:
        "freesurfer",
        "connectome-workbench",
        "singularity",
        "git-annex/8.20200810"

    threads: 1
    resources:
        mem_mb=3000,
        runtime=1,
    log: f"logs/map_brainnetome_atlas_cortex/{'.'.join(inputs.input_wildcards['surf'].values())}.log"
    benchmark: f"benchmarks/map_brainnetome_atlas_cortex/{'.'.join(inputs.input_wildcards['surf'].values())}.tsv"
    params:
        hemi=lambda wcards: wcards["hemi"].lower()
    shell:
        env.export(
            SINGULARITYENV_SUBJECTS_DIR = config["freesurfer_output"],
        ),
        (
            "declare -A structure",
            "structure[l]=CORTEX_LEFT",
            "structure[r]=CORTEX_RIGHT",

            "mris_ca_label sub-{wildcards.subject} {params.hemi}h "
            "{input.sphere} $(pwd)/{input.atlas} $(pwd)/tmp.annot",

            "mris_convert --annot $(pwd)/tmp.annot {input.surf} {output}",

            "wb_command -set-structure {output} "
            "${{structure[{params.hemi}]}}"
        )


rule map_brainnetome_atlas_subcortex:
    input:
        atlas=Path(config["atlases"]['brainnetome'], "BN_Atlas_subcortex.gca"),
        txf=Path(
            config["freesurfer_output"], "sub-{subject}/mri/transforms/talairach.m3z"
        ),
        volume=Path(
            config["freesurfer_output"], "sub-{subject}/mri/brain.mgz"
        ),
    output:
        bids_output_anat(
            atlas="bn",
            desc="subcortex",
            suffix="dseg.nii.gz",
            **inputs.input_wildcards["preproc_dwi"],
        )
    shadow: "minimal"
    log: f"logs/map_brainnetome_atlas_subcortex/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/map_brainnetome_atlas_subcortex/{'.'.join(wildcards.values())}.tsv"
    envmodules:
        "freesurfer",
        "mrtrix",
        "git-annex/8.20200810",
        "singularity",
    group: "connectome"
    threads: 1
    resources:
        mem_mb=3000,
        runtime=10,
    shell:
        env.export(
            SINGULARITYENV_SUBJECTS_DIR = config["freesurfer_output"],
        ),
        (
            "mri_ca_label {input.volume} {input.txf} {input.atlas} $(pwd)/out.mgz",
            "mrconvert out.mgz {output} -quiet",
        )


rule map_labels_to_volume_ribbon:
    input:
        label = rules.map_brainnetome_atlas_cortex.output,
        surf = inputs.input_path["surf_mid"],
        anat_ref = inputs.input_path["t1"],
        white_surf = inputs.input_path["surf"],
        pial_surf = inputs.input_path["surf_pial"],
    output:
        temp(work/"map_labels_to_volume_ribbon"/uid/"{hemi}.nii.gz")

    group: "connectome"
    envmodules:
        "connectome-workbench",
        "singularity",
        "git-annex/8.20200810"

    threads: 1
    resources:
        mem_mb=1000,
        runtime=20,
    log: f"logs/map_labels_to_volume_ribbon/{'.'.join(inputs.input_wildcards['surf'].values())}.log"
    benchmark: f"benchmarks/map_labels_to_volume_ribbon/{'.'.join(inputs.input_wildcards['surf'].values())}.tsv"
    shell:
        "wb_command -label-to-volume-mapping "
        "{input.label} {input.surf} {input.anat_ref} {output} "
        "-ribbon-constrained {input.white_surf} {input.pial_surf} -greedy"


rule combine_dseg_hemispheres:
    input:
        left=expand(
            rules.map_labels_to_volume_ribbon.output,
            hemi=["L"],
            allow_missing=True,
        ),
        right=expand(
            rules.map_labels_to_volume_ribbon.output,
            hemi=["R"],
            allow_missing=True,
        ),
    output:
        bids_output_anat(
            atlas="bn210",
            suffix="dseg.nii.gz",
            **inputs.input_wildcards["preproc_dwi"],
        )

    group: "connectome"
    envmodules:
        "git-annex/8.20200810"

    container:
        "docker://khanlab/connectome-workbench:latest"

    threads: 1
    resources:
        mem_mb=500,
        runtime=10,
    log: f"logs/combine_dseg_hemispheres/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/combine_dseg_hemispheres/{'.'.join(wildcards.values())}.tsv"
    shell:
        "wb_command -volume-math "
        "'max(l, 0) + max(r, 0) - (max(l, 0) * (r > 0))' "
        "-var l {input.left} -var r {input.right} {output}"


rule reslice_cortical_seg_to_subcortex_volume_space:
    input:
        cortex=rules.combine_dseg_hemispheres.output,
        subcortex=rules.map_brainnetome_atlas_subcortex.output,
    output:
        temp(work/"reslice_cortical_seg_to_subcortex_volume_space"/uid/".nii.gz")
    log: f"logs/reslice_cortical_seg_to_subcortex_volume_space/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/reslice_cortical_seg_to_subcortex_volume_space/{'.'.join(wildcards.values())}.tsv"
    envmodules:
        "StdEnv/2020",
        "gcc/9.3.0",
        "ants/2.3.5"
    group: "connectome"
    threads: 1
    resources:
        mem_mb=500,
    shell:
        "antsApplyTransforms -i {input.cortex} -r {input.subcortex} -o {output}"


rule combine_dseg_cortex_subcortex:
    input:
        cortex=rules.reslice_cortical_seg_to_subcortex_volume_space.output,
        subcortex=rules.map_brainnetome_atlas_subcortex.output,
    output:
        bids_output_anat(
            atlas="bn246",
            suffix="dseg.nii.gz",
            **inputs.input_wildcards["preproc_dwi"],
        )
    log: f"logs/combine_dseg_cortex_subcortex/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/combine_dseg_cortex_subcortex/{'.'.join(wildcards.values())}.tsv"
    container:
        "docker://khanlab/connectome-workbench:latest"
    group: "connectome"
    threads: 1
    shell:
        "wb_command -volume-math "
        "'round( max(s, 0) + max(c, 0) - (max(s, 0) * (c > 0)) )' "
        "-var s {input.subcortex} -var c {input.cortex} {output} "


def _get_segmentation(wildcards):
    if wildcards["atlas"] == "bn210":
        return rules.combine_dseg_hemispheres.output[0].format(**wildcards)
    if wildcards["atlas"] == "bn246":
        return rules.combine_dseg_cortex_subcortex.output[0].format(**wildcards)
    raise ValueError(
        "config key 'segmentation' mut be set to one of 'bn210' or 'bn246', currently "
        f"'{config['segmentation']}'"
    )

def _get_weights(wcards):
    if wcards["weight"] == "sift2":
        return rules.run_sift2.output.weights.format(**wcards)
    if wcards["weight"][3:] == "FA":
        return rules.tck_sample.output[0].format(**wcards)
    raise ValueError(
        "config key 'connectome_weight' mut be set to 'sift2' or '___FA', where ___ is "
        f"one of 'avg', 'med', 'min', 'max' currently '{config['segmentation']}'"
    )

rule get_connectome:
    input:
        tracks=rules.run_act.output,
        nodes=_get_segmentation,
        tck_weights=_get_weights,
    output:
        connectome=bids_output_dwi(
            atlas="{atlas}",
            desc="{weight}",
            suffix="connectome.csv",
        ),
        assignments=bids_output_dwi(
            atlas="{atlas}",
            desc="{weight}",
            suffix="assignments.csv",
        ),

    group: "connectome"
    envmodules:
        "git-annex/8.20200810",
        "mrtrix"

    threads: 4
    resources:
        mem_mb=5000,
        runtime=2,
    log: f"logs/get_connectome/{'.'.join(wildcards.values())}.{{atlas}}.{{weight}}.log"
    benchmark: f"benchmarks/get_connectome/{'.'.join(wildcards.values())}.{{atlas}}.{{weight}}.tsv"
    shell:
        "tck2connectome {input.tracks} {input.nodes} {output.connectome} "
        "-q -scale_invnodevol -symmetric -keep_unassigned "
        "-tck_weights_in {input.tck_weights} -out_assignments {output.assignments}"
