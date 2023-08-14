"""Rules for getting, transforming, and resampling reference images
"""

# ref_images, path_store = PathStore.from_mapping({
#     "t1w_MNI152": {
#         res: tflow.get(
#             "MNI152NLin6Asym",
#             resolution=res,
#             desc="brain",
#             suffix="T1w",
#             extension=".nii.gz",
#         )
#         for res in range(1,7)
#     },
#     "fa_probseg": resource("FMRIB58_FA_1mm.nii.gz"),
#     "schaefer": {
#         parcels: tflow.get(
#             "MNI152NLin6Asym",
#             resolution=1,
#             atlas="Schaefer2018",
#             desc=f"{parcels}Parcels17Networks",
#             suffix="dseg",
#             extension=".nii.gz",
#         )
#         for parcels in (100, 200, 300, 400, 500, 600, 800, 1000)
#     },
#     "brainnetome246_dseg": ancient(resource(
#         "atlases/atlas-brainnetome246/space-MNI152NLin6ASym_den-1_atlas-brainnetome246_dseg.nii.gz"
#     )),
#     "brainnetome210": ancient(resource(
#         "atlases/atlas-brainnetome210/atlas-brainnetome210_space-fsLR_den-32k_hemi-{hemi}.label.gii"
#     )),

# })

def _get_atlas(wcards):
    return {
        "bn246": ancient(resource(
            "atlases/atlas-brainnetome246/space-MNI152NLin6ASym_den-1_atlas-brainnetome246_dseg.nii.gz"
        )),
    }[wcards["atlas"]]

def _get_dseg(wcards):
    return {
        "bn210": ancient(resource(
            "atlases/atlas-brainnetome210/atlas-brainnetome210_space-fsLR_den-32k_hemi-{hemi}.label.gii"
        )),
    }[wcards["atlas"]]


rule project_atlas_to_cortex:
    input:
        labels=_get_dseg,
        volume=inputs["t1"].path,
        pial=inputs["pial"].path,
        wm=inputs["wm"].path,
    output:
        tempout(
            "project_atlas_to_cortex",
            inputs['wm'],
            ".nii.gz",
            "atlas",
        )
    log:
        log(
            "project_atlas_to_cortex",
            inputs['wm'],
            "atlas",
        )
    benchmark:
        benchmark(
            "project_atlas_to_cortex",
            inputs['wm'],
            "atlas",
        )
    resources:
        runtime=15,
    singularity:
        config["containers"]["workbench"]
    group: "connectome"
    shell:
        """
        wb_command -label-to-volume-mapping \\
            {input.labels} {input.pial} {input.volume} {output} \\
            -ribbon-constrained {input.wm} {input.pial}
        """

rule merge_hemisphere_labels:
    input:
        left=expand(
            rules.project_atlas_to_cortex.output,
            allow_missing=True,
            hemi="L",
        ),  
        right=expand(
            rules.project_atlas_to_cortex.output,
            allow_missing=True,
            hemi="R",
        ),
    output:
        tempout(
            "merge_hemisphere_labels",
            inputs['preproc_dwi'],
            ".nii.gz",
            "atlas",
        )
    log:
        log(
            "merge_hemisphere_labels",
            inputs['preproc_dwi'],
            "atlas",
        )
    benchmark:
        benchmark(
            "merge_hemisphere_labels",
            inputs['preproc_dwi'],
            "atlas",
        )
    resources:
        runtime=10,
    singularity:
        config["containers"]["workbench"]
    group: 'connectome'
    shell:
        """
        wb_command -volume-math \\
            'max(l, 0) + max(r, 0) - (max(l, 0) * (r > 0))' \\
            -var l {input.left} \\
            -var r {input.right} \\
            {output}
        """


rule convert_fnirt_warps_to_itk:
    input:
        ref=_get_atlas,
        warp=inputs["txf_mni_to_t1w"].path
    output:
        tempout(
            "convert_fnirt_warps_to_itk",
            inputs['preproc_dwi'],
            ".nii.gz",
            "atlas",
        )
    log:
        log(
            "convert_fnirt_warps_to_itk",
            inputs['preproc_dwi'],
            "atlas",
        )
    benchmark:
        benchmark(
            "convert_fnirt_warps_to_itk",
            inputs['preproc_dwi'],
            "atlas",
        )
    resources:
        runtime=1,
    singularity:
        config["containers"]["workbench"]
    group: 'connectome'
    shell:
        """
        wb_command -convert-warpfield \\
            -from-fnirt {input.warp} {input.ref} \\
            -to-itk {output}
        """


def _get_t1_to_mni6_transforms(wcards):
    txf_comp = 'txf_mni_to_t1w'
    t1_txf = inputs[txf_comp].path.format(**wcards)
    mni_space = parse_file_entities(t1_txf, config='derivatives')["from"]
    xfm_type = parse_file_entities(t1_txf)['extension']

    if mni_space not in ["MNI152NLin6Asym", "MNI152NLin6ASym"]:
        raise ValueError(
            f"No transform exists to go from '{mni_space}' space to "
            "'MNI152NLin6Asym' space"
        )
        # suppl_txf = tflow.get(
        #     template="MNI152NLin6Asym",
        #     **{"from": mni_space},
        #     mode="image",
        #     suffix="xfm",
        #     extension=".h5"
        # )
        # if not suppl_txf:
        #     raise ValueError(
        #         f"No transform exists to go from '{mni_space}' space to "
        #         "'MNI152NLin6Asym' space"
        #     )
        # if inverse_:
        #     return [t1_txf, suppl_txf]
        # return [suppl_txf, t1_txf]

    # for now, assume all .nii.gz warp files are FSL fnirt warps
    # THIS IS VERY DANGEROUS, and we need better handling
    if xfm_type == ".nii.gz":
        return rules.convert_fnirt_warps_to_itk.output
    return t1_txf
        

rule warp_atlas_to_t1:
    input:
        atlas=_get_atlas,
        txf=_get_t1_to_mni6_transforms,
        t1=inputs["t1"].path,
    output:
        tempout(
            "warp_between_t1_and_mni_via_ants",
            inputs['preproc_dwi'],
            ".nii.gz",
            "atlas",
        )
    log:
        log(
            "warp_between_t1_and_mni_via_ants",
            inputs['preproc_dwi'],
            "atlas",
        )
    benchmark:
        benchmark(
            "warp_between_t1_and_mni_via_ants",
            inputs['preproc_dwi'],
            "atlas",
        )
    resources:
        runtime=2,
    envmodules:
        "StdEnv/2020",
        "gcc/9.3.0",
        "ants/2.3.5",
    group: "connectome"
    shell:
        """
        compose_transforms () {{
            while [[ -n "${{1:-}}" ]]; do
                printf -- '-t %s ' $1
                shift
            done
        }}
        antsApplyTransforms \\
            -d 3 --interpolation NearestNeighbor \\
            -i {input.atlas} -o {output} -r {input.t1} \\
            $(compose_transforms {input.txf})
        """


rule get_bn246_in_subject_space:
    input:
        cortical=expand(
            rules.merge_hemisphere_labels.output,
            allow_missing=True,
            atlas="bn210",
        ),
        subcortical=expand(
            rules.warp_atlas_to_t1.output,
            allow_missing=True,
            atlas="bn246",
        )
    output:
        bids_output_dwi(
            atlas="bn246",
            suffix="dparc.nii.gz",
        ),
    log:
        log(
            "get_bn246_in_subject_space",
            inputs["preproc_dwi"],
        )
    benchmark:
        benchmark(
            "get_bn246_in_subject_space",
            inputs["preproc_dwi"],
        )
    singularity:
        config["containers"]["workbench"]
    group: "connectome"
    shell:
        """
        wb_command -volume-math \\
            'round( s*(s>210) + max(c, 0) - (s*(s>210) * (c>0)) )' \\
            -var s {input.subcortical} \\
            -var c {input.cortical} \\
            {output} 
        """