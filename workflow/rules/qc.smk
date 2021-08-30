from snakebids import bids

wildcards = config['input_wildcards']['preproc_dwi']

rule create_qc_reference_image:
    input:
        dwi=config['input_path']['preproc_dwi'],
        bvec=config['input_path']['bvec'],
        bval=config['input_path']['bval']
    output:
        bids(root=qc,
            datatype='dwi',
            suffix="dwi.mif",
            **wildcards)
    envmodules:
        "mrtrix/3.0.1"
    shell:
        'mrconvert {input.dwi} {output} -fslgrad {input.bvec} {input.bval} 2> {log}'

rule generate_qc_FOD_image:
    input:
        wm=bids(root=work,
                datatype='dwi',
                desc="norm",
                suffix='wmfod.mif',
                **wildcards),
        gm=bids(root=work,
                datatype='dwi',
                desc="norm",
                suffix='gmfod.mif',
                **wildcards),
        csf=bids(root=work,
                datatype='dwi',
                desc="norm",
                suffix='csffod.mif',
                **wildcards)
    output:
        bids(root=qc,
            datatype='dwi',
            suffix='vf.mif',
            **wildcards)
    envmodules:
        "mrtrix/3.0.1"
    shell:
        'mrconvert -coord 3 0 {input.wm} - | mrcat {input.csf} {input.gm} - {output}'


rule extract_tractography_subset:
    input:
        bids(root=output,
            datatype='dwi',
            desc="tracks",
            suffix='10M.tck',
            **wildcards)
    output:
        bids(root=qc,
            datatype='dwi',
            desc="tracks",
            suffix='200k.tck',
            **wildcards)
    envmodules:
        "mrtrix/3.0.1"
    shell:
        'tckedit {input} -number 200k {output}'

rule generate_tractography_view_script:
    input:
        dwi=bids(root=qc,
            datatype='dwi',
            suffix="dwi.mif",
            **wildcards),
        tracts=bids(root=qc,
            datatype='dwi',
            desc="tracks",
            suffix='200k.tck',
            **wildcards)
    output:
        bids(root=qc,
            datatype='dwi',
            suffix="viewTractography",
            **wildcards)
    shell:
        'echo "#!/bin/bash\nmrview $(realpath {input.dwi}) -tractography.load $(realpath {input.tracts})" > {output} && chmod +x {output}'

rule generate_5tt_qc_view_script:
    input:
        dwi=bids(root=qc,
            datatype='dwi',
            suffix="dwi.mif",
            **wildcards),
        anat=bids(root=work,
            datatype='anat',
            suffix="gmwmi.mif",
            **wildcards)
    output:
        bids(root=qc,
            datatype='dwi',
            suffix="viewInterface",
            **wildcards)
    shell:
        'echo "#!/bin/bash\nmrview $(realpath {input.dwi}) -overlay.load $(realpath {input.anat})" > {output} && chmod +x {output}'

rule generate_odf_qc_view_script:
    input:
        vf=bids(root=qc,
            datatype='dwi',
            suffix='vf.mif',
            **wildcards),
        odf=bids(root=work,
                datatype='dwi',
                desc="norm",
                suffix='wmfod.mif',
                **wildcards)
    output:
        bids(root=qc,
            datatype='dwi',
            suffix="viewOdf",
            **wildcards)
    shell:
        'echo "#!/bin/bash\nmrview $(realpath {input.vf}) -odf.load_sh $(realpath {input.odf})" > {output} && chmod +x {output}'