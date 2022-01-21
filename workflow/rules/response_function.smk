from snakebids import bids
from lib.shells import FodAlgorithm


rule convert_dwi_to_mrtrix_format:
    input:
        dwi=inputs.input_path['preproc_dwi'],
        bvec=inputs.input_path['bvec'],
        bval=inputs.input_path['bval']

    output:
        temp(work/"response-function"/uid/"dwi.mif")

    envmodules:
        "mrtrix/3.0.1",
        "git-annex/8.20200810"
    log: f"logs/convert_dwi_to_mrtrix_format/{'.'.join(wildcards.values())}.log"
    group: "response_generation"
    shell:
        datalad(
            'mrconvert {input.dwi} {output} -fslgrad {input.bvec} {input.bval} 2> {log}'
        )



rule convert_mask_to_mrtrix_format:
    input:
        inputs.input_path['brainmask']
    output:
        temp(work/"response-function"/uid/"mask.mif")
    envmodules:
        "mrtrix/3.0.1",
        "git-annex/8.20200810"
    log: "logs/convert_mask_to_mrtrix_format/{subject}.log"
    group: "response_generation"
    shell:
        datalad(
            'mrconvert {input} {output} 2> {log}'
        )


rule generate_response_function:
    input:
        dwi=rules.convert_dwi_to_mrtrix_format.output,
        mask=rules.convert_mask_to_mrtrix_format.output
    output:
        wm=bids_output_dwi(
            model="CSD",
            label="wm",
            suffix='fod.txt'
        ),
        gm=bids_output_dwi(
            model="CSD",
            label="gm",
            suffix='fod.txt'
        ),
        csf=bids_output_dwi(
            model="CSD",
            label="csf",
            suffix='fod.txt'
        ),
        voxels=bids_output_dwi(
            model="CSD",
            desc="voxels",
            suffix='fod.mif'
        )
    group: "response_generation"
    log: "logs/generate_response_function/{subject}.log"
    resources:
        runtime=2
    envmodules:
        "mrtrix/3.0.1",
        "git-annex/8.20200810"
    shell:
        datalad.msg("Estimate WM, GM, CSF response functions")(
            'dwi2response dhollander {input.dwi} {output.wm} {output.gm} {output.csf} '
            '-voxels {output.voxels} -mask {input.mask} -scratch {resources.tmpdir} '
            '2> {log}'
        )


rule compute_ss3t_fiber_orientation_densities:
    input:
        dwi=rules.convert_dwi_to_mrtrix_format.output,
        wm=rules.generate_response_function.output.wm,
        gm=rules.generate_response_function.output.gm,
        csf=rules.generate_response_function.output.csf,
        mask=rules.convert_mask_to_mrtrix_format.output
    output:
        wm=bids_output_dwi(
            model="CSD",
            desc="ss3t",
            label="wm",
            suffix='fod.mif'
        ),
        gm=bids_output_dwi(
            model="CSD",
            desc="ss3t",
            label="gm",
            suffix='fod.mif'
        ),
        csf=bids_output_dwi(
            model="CSD",
            desc="ss3t",
            label="csf",
            suffix='fod.mif'
        )
    group: "response_generation"
    threads: 32
    # This still needs to be benchmarked!!
    resources:
        mem_mb=10000,
        runtime=25,
    container:
        'docker://pennbbl/ss3t_beta:0.0.1'
    benchmark:
        'benchmarks/compute_fiber_orientation_densities/{subject}.tsv'
    shell:
        datalad.msg("Compute fod using ss3t algorithm")(
            'ss3t_csd_beta1 {input.dwi} '
            '{input.wm} {output.wm} {input.gm} {output.gm} {input.csf} {output.csf} '
            '-mask {input.mask} -nthreads {threads} -scratch {resources.tmpdir}'
        )



rule compute_ms3t_fiber_orientation_densities:
    input:
        dwi=rules.convert_dwi_to_mrtrix_format.output,
        wm=rules.generate_response_function.output.wm,
        gm=rules.generate_response_function.output.gm,
        csf=rules.generate_response_function.output.csf,
        mask=rules.convert_mask_to_mrtrix_format.output
    output:
        wm=bids_output_dwi(
            model="CSD",
            desc="ms3t",
            label="wm",
            suffix='fod.mif'
        ),
        gm=bids_output_dwi(
            model="CSD",
            desc="ms3t",
            label="gm",
            suffix='fod.mif'
        ),
        csf=bids_output_dwi(
            model="CSD",
            desc="ms3t",
            label="csf",
            suffix='fod.mif'
        )
    group: "response_generation"
    threads: 32
    resources:
        mem_mb=10000,
        runtime=15,
    envmodules:
        "mrtrix/3.0.1",
        "git-annex/8.20200810"
    log: "logs/compute_fiber_orientation_densities/{subject}.log"
    shell:
        datalad.msg("Compute fod using ss3t algorithm")(
            'dwi2fod msmt_csd {input.dwi} '
            '{input.wm} {output.wm} {input.gm} {output.gm} {input.csf} {output.csf} '
            '-mask {input.mask} -nthreads {threads} 2> {log}'
        )


rule normalize_fiber_orientation_densities:
    input:
        wm=FodAlgorithm(
            root=output,
            bval=inputs.input_path['bval'],
            tissue='wm'
        ),
        gm=FodAlgorithm(
            root=output,
            bval=inputs.input_path['bval'],
            tissue='gm'
        ),
        csf=FodAlgorithm(
            root=output,
            bval=inputs.input_path['bval'],
            tissue='csf'
        ),
        mask=rules.convert_mask_to_mrtrix_format.output
    output:
        wm=bids_output_dwi(
            model="CSD",
            desc="norm",
            label="wm",
            suffix='fod.mif'
        ),
        gm=bids_output_dwi(
            model="CSD",
            desc="norm",
            label="gm",
            suffix='fod.mif'
        ),
        csf=bids_output_dwi(
            model="CSD",
            desc="norm",
            label="csf",
            suffix='fod.mif'
        )
    group: "response_generation"
    resources:
        runtime=2
    envmodules:
        "mrtrix/3.0.1",
        "git-annex/8.20200810"
    log: "logs/normalize_fiber_orientation_densities/{subject}.log"
    shell:
        datalad.msg("Normalize fod functions")(
            'mtnormalise '
            '{input.wm} {output.wm} {input.gm} {output.gm} {input.csf} {output.csf} '
            '-mask {input.mask} 2> {log}'
        )
