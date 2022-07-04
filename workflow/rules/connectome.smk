# group: connectome:
#       connected_components: 96
def translate_hemi(path):
    def inner(wcards):
        wcards.hemi = wcards.hemi.lower() + "h"
        return str(path).format(**wcards)
    return inner

rule map_brainnetome_atlas:
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
            atlas="bn",
            suffix="aparc.label.gii",
            **inputs.input_wildcards["surf"],
        )

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
    log: f"logs/map_brainnetome_atlas/{'.'.join(inputs.input_wildcards['surf'].values())}.log"
    benchmark: f"benchmarks/map_brainnetome_atlas/{'.'.join(inputs.input_wildcards['surf'].values())}.tsv"
    shell:
        boost(
            (
                subject_dir := sh.ShVar(
                    config["freesurfer_output"],
                    name="SINGULARITYENV_SUBJECTS_DIR",
                    export=True
                ),
                hemi := sh.ShVar(
                    sh.echo("{wildcards.hemi}") |
                    sh.awk('print tolower($0)')
                ),
                tmp_out := sh.ShVar(
                    work/"map_brainnetome_atlas"/f"{shell_uid('surf')}.annot"
                ),
                "declare -A structure",
                "structure[l]=CORTEX_LEFT",
                "structure[r]=CORTEX_RIGHT",

                f"mkdir -p $(dirname {tmp_out})",

                f"mris_ca_label sub-{{wildcards.subject}} {hemi.escape()}h "
                f"{{input.sphere}} {{input.atlas}} {tmp_out}",

                f"mris_convert --annot {tmp_out} {{input.surf}} {{output}}",

                "wb_command -set-structure {output} "
                f"${{{{structure[{hemi.escape()}]}}}}"
            )
        )


rule map_labels_to_volume_ribbon:
    input:
        label = rules.map_brainnetome_atlas.output,
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
    params:
        tmp_out = str(work/uid/"{hemi}.annot")
    shell:
        # "touch {output}"
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
            atlas="bn",
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
        # "touch {output}"
        "wb_command -volume-math "
        "'max(l, 0) + max(r, 0) - (max(l, 0) * (r > 0))' "
        "-var l {input.left} -var r {input.right} {output}"


rule get_connectome:
    input:
        tracks=rules.run_act.output,
        nodes=rules.combine_dseg_hemispheres.output,
        tck_weights=rules.run_sift2.output.weights,
    output:
        connectome=bids_output_dwi(
            atlas="bn",
            suffix="connectome.csv",
        ),
        assignments=bids_output_dwi(
            atlas="bn",
            desc="tract",
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
    log: f"logs/get_connectome/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/get_connectome/{'.'.join(wildcards.values())}.tsv"
    shell:
        "tck2connectome {input.tracks} {input.nodes} {output.connectome} "
        "-scale_invnodevol -symmetric -keep_unassigned "
        "-tck_weights_in {input.tck_weights} -out_assignments {output.assignments}"
