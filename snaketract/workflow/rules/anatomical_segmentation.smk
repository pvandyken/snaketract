# group segmentation:
#     group_components=288

rule segment_anatomical_image:
    input:
        data=inputs["t1"].path
    output:
        bids_output_anat(
            suffix="dseg.nii.gz"
        )
    log: log("segment_anatomical_image", inputs['t1'])
    benchmark: benchmark("segment_anatomical_image", inputs['t1'])
    group: 'segmentation'
    threads: 4
    container:
        "docker://freesurfer/freesurfer:7.3.1"
    resources:
        mem_mb=14000,
        runtime=14,
    shell: "mri_synthseg --i {input} --o {output} --threads 2"


def _get_segmentation(wcards):
    if "t1_dseg" in inputs:
        return inputs["t1_dseg"].path.format(**wcards)
    return rules.segment_anatomical_image.output[0].format(**wcards)

rule get_5tt_segmentation:
    input:
        segmentation=_get_segmentation,
        lut=ancient(resource("FreeSurferColorLUT.txt")),
    output:
        bids_output_anat(
            suffix="5tt.nii.gz",
        )
    group: "segmentation"
    resources:
        runtime=1,
        mem_mb=4000,
    log: log("get_5tt_segmentation", inputs['t1'])
    benchmark: benchmark("get_5tt_segmentation", inputs['t1'])
    container: 'docker://mrtrix3/mrtrix3:3.0.3'
    shell:
        '5ttgen freesurfer '
        '{input.segmentation} {output} '
        '-scratch {resources.tmpdir} -lut {input.lut} 2> {log}'

rule create_seed_boundary:
    input:
        rules.get_5tt_segmentation.output
    output:
        bids_output_anat(
            suffix="gmwmInterface.nii.gz",
        )
    group: "segmentation"
    log: log("create_seed_boundary", inputs['t1'])
    benchmark: benchmark("create_seed_boundary", inputs['t1'])
    envmodules:
        "mrtrix/3.0.1",
        "git-annex/8.20200810"
    resources:
        runtime=5,
    shell:
        # datalad.msg("Compute boundary between GM/WM for tractography"),
        '5tt2gmwmi {input} {output} 2> {log}'
