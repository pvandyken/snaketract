from snakebids import bids

wildcards = config['input_wildcards']['preproc_dwi']

work = config['directories']['output']
qc = config['directories']['qc']
output = config['directories']['output']

localrules: convert_dwi_to_mrtrix_format

rule convert_dwi_to_mrtrix_format:
    input:
        dwi=config['input_path']['preproc_dwi'],
        bvec=config['input_path']['bvec'],
        bval=config['input_path']['bval']
    output:
        temp(bids(root=work,
            datatype='dwi',
            suffix="dwi.mif",
            **wildcards))
    envmodules:
        "mrtrix/3.0.1"
    shell:
        'mrconvert {input.dwi} {output} -fslgrad {input.bvec} {input.bval}'

rule convert_mask_to_mrtrix_format:
    input:
        config['input_path']['brainmask']
    output:
        temp(bids(root=work,
            datatype='dwi',
            suffix="brainmask.mif",
            **wildcards))
    envmodules:
        "mrtrix/3.0.1"
    shell:
        'mrconvert {input} {output}'

rule generate_response_function:
    input:
        dwi=bids(root=work,
            datatype='dwi',
            suffix="dwi.mif",
            **wildcards),
        mask=bids(root=work,
            datatype='dwi',
            suffix="brainmask.mif",
            **wildcards)
    output:
        wm=bids(root=work,
                datatype='dwi',
                suffix='wm.txt',
                **wildcards),
        gm=bids(root=work,
                datatype='dwi',
                suffix='gm.txt',
                **wildcards),
        csf=bids(root=work,
                datatype='dwi',
                suffix='csf.txt',
                **wildcards),
        voxels=bids(root=output,
                datatype='dwi',
                suffix='voxels.mif',
                **wildcards)
    group: groups.response_generation
    resources:
        tmpdir=config['tmpdir']
    benchmark:
        'benchmarks/generate_response_function/{subject}.tsv'
    envmodules:
        "mrtrix/3.0.1"
    shell:
        'dwi2response dhollander {input.dwi} {output.wm} {output.gm} {output.csf} '
        '-voxels {output.voxels} -mask {input.mask} -scratch {resources.tmpdir}'


rule compute_fiber_orientation_densities:
    input:
        dwi=bids(root=work,
                datatype='dwi',
                suffix='dwi.mif',
                **wildcards),
        wm=bids(root=work,
                datatype='dwi',
                suffix='wm.txt',
                **wildcards),
        gm=bids(root=work,
                datatype='dwi',
                suffix='gm.txt',
                **wildcards),
        csf=bids(root=work,
                datatype='dwi',
                suffix='csf.txt',
                **wildcards),
        mask=bids(root=work,
            datatype='dwi',
            suffix="brainmask.mif",
            **wildcards)
    output: 
        wm=bids(root=work,
                datatype='dwi',
                suffix='wmfod.mif',
                **wildcards),
        gm=bids(root=work,
                datatype='dwi',
                suffix='gmfod.mif',
                **wildcards),
        csf=bids(root=work,
                datatype='dwi',
                suffix='csffod.mif',
                **wildcards)
    group: groups.response_generation
    threads: 8
    resources:
        tmpdir=config['tmpdir']
    container:
        'docker://pennbbl/ss3t_beta:0.0.1'
    benchmark:
        'benchmarks/compute_fiber_orientation_densities/{subject}.tsv'
    shell:
        'ss3t_csd_beta1 {input.dwi} '
        '{input.wm} {output.wm} {input.gm} {output.gm} {input.csf} {output.csf} '
        '-mask {input.mask} -nthreads {threads} -scratch {resources.tmpdir}'

rule normalize_fiber_orientation_densities:
    input:
        wm=bids(root=work,
                datatype='dwi',
                suffix='wmfod.mif',
                **wildcards),
        gm=bids(root=work,
                datatype='dwi',
                suffix='gmfod.mif',
                **wildcards),
        csf=bids(root=work,
                datatype='dwi',
                suffix='csffod.mif',
                **wildcards),
        mask=bids(root=work,
            datatype='dwi',
            suffix="brainmask.mif",
            **wildcards)
    output:
        wm=bids(root=work,
                datatype='dwi',
                suffix='wmfod_norm.mif',
                **wildcards),
        gm=bids(root=work,
                datatype='dwi',
                suffix='gmfod_norm.mif',
                **wildcards),
        csf=bids(root=work,
                datatype='dwi',
                suffix='csffod_norm.mif',
                **wildcards)
    group: groups.response_generation
    benchmark:
        'benchmarks/normalize_fiber_orientation_densities/{subject}.tsv'
    envmodules:
        "mrtrix/3.0.1"
    shell:
        'mtnormalize {input.wm} {output.wm} {input.gm} {output.gm} {input.csf} {output.csf} -mask {input.mask}'
