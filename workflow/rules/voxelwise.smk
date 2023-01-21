rule warp_dmri_to_mni:
    input:
        dmri=inputs.input_path["fa"],
        txf=inputs['txf_t1w_to_mni'].path,
        ref=inputs['mni_t1w_ref'].path,
    output:
        bids_output_dwi(
            space="MNI152NLin2009cAsym",
            suffix="FA.nii.gz",
        )
    log: f"logs/warp_dmri_to_mni/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/warp_dmri_to_mni/{'.'.join(wildcards.values())}.tsv"
    envmodules:
        "StdEnv/2020",
        "gcc/9.3.0",
        "ants/2.3.5",
    group: "wm_mask"
    shell:
        "antsApplyTransforms "
        "-d 3 --interpolation linear "
        "-i {input.dmri} -o {output} -r {input.ref} -t {input.txf}"


rule warp_dmri_to_mni6:
    input:
        dmri=inputs.input_path["fa"],
        mni2009=inputs['txf_t1w_to_mni'].path,
        mni6=tflow.get(
            "MNI152NLin6Asym",
            **{"from":"MNI152NLin2009cAsym"},
            mode="image",
            extension=".h5",
        ),
        ref=tflow.get(
            "MNI152NLin6Asym",
            resolution=2,
            desc="brain",
            suffix="T1w",
            extension=".nii.gz",
        )
    output:
        bids_output_dwi(
            space="MNI152NLin6Asym",
            suffix="FA.nii.gz",
        )
    log: f"logs/warp_dmri_to_mni6/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/warp_dmri_to_mni6/{'.'.join(wildcards.values())}.tsv"
    envmodules:
        "StdEnv/2020",
        "gcc/9.3.0",
        "ants/2.3.5",
    group: "wm_mask"
    shell:
        "antsApplyTransforms "
        "-d 3 --interpolation linear "
        "-i {input.dmri} -o {output} -r {input.ref} "
        "-t {input.mni2009} -t {input.mni6}"



rule transform_dmri_image_to_mask:
    """Finish transforming wm_mask to dmri image

    Calculate and apply the transform needed to make the wm mask exactly
    overlap with the dmri image
    """
    input:
        mask="resources/FMRIB58_FA_1mm.nii.gz",
        dmri=rules.warp_dmri_to_mni.output,
    output:
        temp(
            work/"transform_wm_maks_to_dmri_image"/uid/"txf.nii.gz"
        )
    log: f"logs/transform_dmri_image_to_mask/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/transform_dmri_image_to_mask/{'.'.join(wildcards.values())}.tsv"
    group: "wm_mask"
    envmodules:
        "StdEnv/2020",
        "gcc/9.3.0",
        "ants/2.3.5",
    threads: 4
    resources:
        mem_mb=4000,
        runtime=3,
    shadow: 'minimal'
    shell:
        "antsRegistrationSyNQuick.sh -d 3 -n {threads} "
        "-f {input.mask} -m {input.dmri} -o txf > {log}",

        "mv txfWarped.nii.gz {output}"


rule apply_wm_mask_to_dmri:
    input:
        mask=rules.transform_dmri_image_to_mask.input.mask,
        dmri=rules.transform_dmri_image_to_mask.output,
    output:
        bids_output_dwi(
            space="MNI152NLin2009cAsym",
            label="wm",
            suffix="FA.nii.gz",
        )
    log: f"logs/apply_wm_mask_to_dmri/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/apply_wm_mask_to_dmri/{'.'.join(wildcards.values())}.tsv"
    group: 'wm_mask'
    container:
        "docker://pyushkevich/itksnap:v3.8.2"
    shell:
        "c3d {input.dmri} {input.mask} -thresh 2500 inf 1 0 -times {output}"


def _get_subject_labels_from_group(group):
    filter_list(
        inputs['preproc_dwi'].zip_lists,
        {
            "subject": filter_participants(
                os.path.join(config['bids_dir'], 'participants.tsv'),
                group=group,
            ).map(lambda s: s[4:])
        }
    )



def _get_subjects_from_group(group):
    return expand(
        rules.apply_wm_mask_to_dmri.output,
        zip,
        **filter_list(
            inputs['preproc_dwi'].zip_lists,
            {
                "subject": (
                    filter_participants(
                        Path(config['bids_dir'], 'participants.tsv'),
                        group=group,
                    )['participant_id']
                    .map(lambda s: s[4:])
                )
            }
        )
    )


rule filter_participants_file:
    input:
        Path(config['bids_dir'], 'participants.tsv'),
    output:
        temp(work/"filter_participants_file.tsv")
    log: f"logs/filter_participants_file.log"
    benchmark: f"benchmarks/filter_participants_file.tsv"
    group: "foo"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=30,
    params:
        group=["HC", "FEP"]
    run:
        filter_participants(
            str(input),
            group=params.group,
            participant_id=map(
                lambda id: f"sub-{id}",
                inputs['preproc_dwi'].entities['subject'],
            )
        ).to_csv(str(output), sep="\t")

rule make_voxelwise_design_matrix:
    input:
        rules.filter_participants_file.output
    output:
        mat=temp(work/"make_voxelwise_design_matrix"/"mat.txt"),
        con=temp(work/"make_voxelwise_design_matrix"/"con.txt"),
    log: f"logs/make_voxelwise_design_matrix.log"
    benchmark: f"benchmarks/make_voxelwise_design_matrix.tsv"
    group: "foo"
    threads: 1
    envmodules:
        'python/3.8'
    shell:
        design_matrix_env.script,
        Pyscript(workflow.basedir)(
            script="scripts/design_matrix.py",
            output=["mat", "con"]
        )

rule convert_design_matrices_to_fsl:
    input:
        mat=rules.make_voxelwise_design_matrix.output.mat,
        con=rules.make_voxelwise_design_matrix.output.con,
    output:
        mat=temp(work/"convert_design_matrices_to_fsl.mat"),
        con=temp(work/"convert_design_matrices_to_fsl.con"),
    log: f"logs/convert_design_matrices_to_fsl.log"
    benchmark: f"benchmarks/convert_design_matrices_to_fsl.tsv"
    group: "foo"
    threads: 1
    shell:
        "Text2Vest {input.mat} {output.mat}",
        "Text2Vest {input.con} {output.con}"

rule get_average_fa_image:
    input:
        _get_subjects_from_group(['HC'])
    output:
        Path(config['output_dir'], 'avgFA.nii.gz')
    log: f"logs/get_average_fa_image.log"
    benchmark: f"benchmarks/get_average_fa_image.tsv"
    container:
        "docker://pyushkevich/itksnap:v3.8.2"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=30,
    shell:
        "c3d {input} -mean -o {output}"

# rule get_std_fa_image:
#     input:
#         _get_subjects_from_group(['HC'])
#     output:
#         Path(config['output_dir'], 'avgFA.nii.gz')
#     log: f"logs/get_std_fa_image/{'.'.join(wildcards.values())}.log"
#     benchmark: f"benchmarks/get_std_fa_image/{'.'.join(wildcards.values())}.tsv"
#     container:
#         "docker://pyushkevich/itksnap:v3.8.2"
#     threads: 1
#     resources:
#         mem_mb=1000,
#         runtime=30,
#     shell:
#         "c3d {input} -mean -o {output}"


rule do_voxelwise_comparison:
    input:
        data=_get_subjects_from_group(['HC', 'FEP']),
        design=rules.convert_design_matrices_to_fsl.output.mat,
        contrast=rules.convert_design_matrices_to_fsl.output.con,
        mask=rules.transform_dmri_image_to_mask.input.mask,
    output:
        directory(Path(config['output_dir'], 'voxelwise'))
    log: f"logs/do_voxelwise_comparison.log"
    benchmark: f"benchmarks/do_voxelwise_comparison.tsv"
    params:
        permutations=10000,
    envmodules:
        "StdEnv/2020",
        "gcc/9.3.0",
        "fsl/6.0.4",
    shadow: 'minimal'
    shell:
        "fslmerge -t merged.nii.gz {input.data} ",

        "randomise -i merged.nii.gz -o out "
        "-d {input.design} -t {input.contrast} -m {input.mask} "
        "-n {params.permutations} -T",

        "mkdir {output}",

        """
        for file in out_*; do
          mv "$file" {output}/"${{file#out_}}"
        done
        """
