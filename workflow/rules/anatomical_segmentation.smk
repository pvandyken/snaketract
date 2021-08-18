from snakebids import bids

wildcards = config['input_wildcards']['preproc_dwi']

work = config['directories']['work']
qc = config['directories']['qc']
output = config['directories']['output']

rule convert_t1_to_mrtrix_format:
    input:
        config['input_path']['t1']
    output:
        bids(root=work,
            datatype='anat',
            suffix="t1w.mif",
            **wildcards)
    group: groups.segmentation
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
    shell:
        '5tt2gmwmi {input} {output}'