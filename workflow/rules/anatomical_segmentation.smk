from snakebids import bids

wildcards = config['input_wildcards']['preproc_dwi']

work = config['tmpdir'] + '/data'
qc = config['directories']['qc']
output = config['directories']['output']

localrules: convert_t1_to_mrtrix_format

rule convert_t1_to_mrtrix_format:
    input:
        config['input_path']['t1']
    output:
        bids(root=work,
            datatype='anat',
            suffix="t1w.mif",
            **wildcards)
    group: groups.segmentation
    envmodules:
        "mrtrix/3.0.1"
    shell:
        'mrconvert {input} {output}'


rule segment_anatomical_image:
    input:
        bids(root=work,
            datatype='anat',
            suffix="t1w.mif",
            **wildcards)
    output:
        bids(root=work,
            datatype='anat',
            suffix="5tt.mif",
            **wildcards)
    group: groups.segmentation
    resources:
        tmpdir=config['tmpdir']
    benchmark:
        'benchmarks/segment_anatomical_image/{subject}.tsv'
    envmodules:
        "mrtrix/3.0.1",
        "StdEnv/2020",
        "gcc/9.3.0",
        "fsl/6.0.4"
    shell:
        '5ttgen fsl {input} {output} -premasked -scratch {resources.tmpdir}'

rule create_seed_boundary:
    input:
        bids(root=work,
            datatype='anat',
            suffix="5tt.mif",
            **wildcards)
    output:
        bids(root=work,
            datatype='anat',
            suffix="gmwmi.mif",
            **wildcards)
    group: groups.segmentation
    benchmark:
        'benchmarks/create_seed_boundary/{subject}.tsv'
    envmodules:
        "mrtrix/3.0.1"
    shell:
        '5tt2gmwmi {input} {output}'