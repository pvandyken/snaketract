from snakebids import bids

rule create_qc_reference_image:
    input:
        dwi=input_paths['preproc_dwi'],
        bvec=input_paths['bvec'],
        bval=input_paths['bval']
    output:
        bids(root=qc,
            datatype='dwi',
            suffix="dwi.mif",
            **wildcards)
    envmodules:
        "mrtrix/3.0.1"
    shell:
        'mrconvert {input.dwi} {output} -fslgrad {input.bvec} {input.bval} -quiet'

rule generate_qc_FOD_image:
    input:
        wm=rules.normalize_fiber_orientation_densities.output.wm,
        gm=rules.normalize_fiber_orientation_densities.output.gm,
        csf=rules.normalize_fiber_orientation_densities.output.csf
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
    input: rules.run_act.output
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
        dwi=rules.create_qc_reference_image.output,
        tracts=rules.extract_tractography_subset.output
    output:
        bids(root=qc,
            datatype='dwi',
            suffix="viewTractography",
            **wildcards)
    shell:
        'echo "#!/bin/bash\nmrview $(realpath {input.dwi}) -tractography.load $(realpath {input.tracts})" > {output} && chmod +x {output}'

rule generate_5tt_qc_view_script:
    input:
        dwi=rules.create_qc_reference_image.output,
        anat=rules.create_seed_boundary.output
    output:
        bids(root=qc,
            datatype='dwi',
            suffix="viewInterface",
            **wildcards)
    shell:
        'echo "#!/bin/bash\nmrview $(realpath {input.dwi}) -overlay.load $(realpath {input.anat})" > {output} && chmod +x {output}'

rule generate_odf_qc_view_script:
    input:
        vf=rules.generate_qc_FOD_image.output,
        odf=rules.normalize_fiber_orientation_densities.output.wm
    output:
        bids(root=qc,
            datatype='dwi',
            suffix="viewOdf",
            **wildcards)
    shell:
        'echo "#!/bin/bash\nmrview $(realpath {input.vf}) -odf.load_sh $(realpath {input.odf})" > {output} && chmod +x {output}'

rule create_tractography_png:
    input:
        dwi=input_paths['preproc_dwi'],
        tracts=rules.run_act.output
    output:
        directory(bids(
            root=qc,
            datatype='dwi',
            desc="tracts",
            suffix="imgs",
            **wildcards
        ))
    conda:
        "envs/image-extraction.yaml"
    shell:
        "./scripts/image_extraction.py "
        "--tractography {input.tracts} "
        "{input.dwi} {output}"