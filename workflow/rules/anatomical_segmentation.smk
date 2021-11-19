from snakebids import bids

rule convert_t1_to_mrtrix_format:
    input:
        input_paths['t1']
    output:
        temp(bids_output_anat(
            root=output,
            suffix="t1w.mif"
        ))
    group: "segmentation"
    log: "logs/convert_t1_to_mrtrix_format/{subject}.log"
    envmodules:
        "mrtrix/3.0.1"
    shell:
        'mrconvert {input} {output} 2> {log}'


rule segment_anatomical_image:
    input:
        rules.convert_t1_to_mrtrix_format.output
    output:
        bids_output_anat(
            suffix="5tt.mif",
        )
    group: "segmentation"
    resources:
        mem_mb=2500,
        runtime=20
    log: "logs/segment_anatomical_image/{subject}.log"
    envmodules:
        "mrtrix/3.0.1",
        "StdEnv/2020",
        "gcc/9.3.0",
        "fsl/6.0.4"
    shell:
        '5ttgen fsl {input} {output} -premasked -scratch {resources.tmpdir} 2> {log}'

rule create_seed_boundary:
    input:
        rules.segment_anatomical_image.output
    output:
        bids_output_anat(
            suffix="gmwmInterface.mif",
        )
    group: "segmentation"
    log: "logs/create_seed_boundary/{subject}.log"
    envmodules:
        "mrtrix/3.0.1"
    shell:
        '5tt2gmwmi {input} {output} 2> {log}'