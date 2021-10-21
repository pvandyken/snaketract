from snakebids import bids
from lib.shells import FodAlgorithm


rule convert_dwi_to_mrtrix_format:
    input:
        dwi=input_paths['preproc_dwi'],
        bvec=input_paths['bvec'],
        bval=input_paths['bval']

    output:
        temp(bids_output_dwi(
            suffix="dwi.mif"
        ))
    
    envmodules:
        "mrtrix/3.0.1"
    log: f"logs/convert_dwi_to_mrtrix_format/{'.'.join(wildcards.values())}.log"
    group: "response_generation"
    shell: 
        'mrconvert {input.dwi} {output} -fslgrad {input.bvec} {input.bval} 2> {log}'


rule convert_mask_to_mrtrix_format:
    input:
        input_paths['brainmask']
    output:
        temp(bids_output_dwi(
            suffix="brainmask.mif"
        ))
    envmodules:
        "mrtrix/3.0.1"
    log: "logs/convert_mask_to_mrtrix_format/{subject}.log"
    group: "response_generation"
    shell:
        'mrconvert {input} {output} 2> {log}'


rule generate_response_function:
    input:
        dwi=rules.convert_dwi_to_mrtrix_format.output,
        mask=rules.convert_mask_to_mrtrix_format.output
    output:
        wm=bids_output_dwi(
            desc="fod",
            suffix='wm.txt'
        ),
        gm=bids_output_dwi(
            desc="fod",
            suffix='gm.txt'
        ),
        csf=bids_output_dwi(
            desc="fod",
            suffix='csf.txt'
        ),
        voxels=bids_output_dwi(
            desc="fod",
            suffix='voxels.mif'
        )
    group: "response_generation"
    log: "logs/generate_response_function/{subject}.log"
    resources:
        runtime=2
    envmodules:
        "mrtrix/3.0.1"
    shell:
        'dwi2response dhollander {input.dwi} {output.wm} {output.gm} {output.csf} '
        '-voxels {output.voxels} -mask {input.mask} -scratch {resources.tmpdir} 2> {log}'


rule compute_ss3t_fiber_orientation_densities:
    input:
        dwi=rules.convert_dwi_to_mrtrix_format.output,
        wm=rules.generate_response_function.output.wm,
        gm=rules.generate_response_function.output.gm,
        csf=rules.generate_response_function.output.csf,
        mask=rules.convert_mask_to_mrtrix_format.output
    output:
        wm=bids_output_dwi(
            desc="fodSs3t",
            suffix='wm.mif'
        ),
        gm=bids_output_dwi(
            desc="fodSs3t",
            suffix='gm.mif'
        ),
        csf=bids_output_dwi(
            desc="fodSs3t",
            suffix='csf.mif'
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
        'ss3t_csd_beta1 {input.dwi} '
        '{input.wm} {output.wm} {input.gm} {output.gm} {input.csf} {output.csf} '
        '-mask {input.mask} -nthreads {threads} -scratch {resources.tmpdir}'


rule compute_ms3t_fiber_orientation_densities:
    input:
        dwi=rules.convert_dwi_to_mrtrix_format.output,
        wm=rules.generate_response_function.output.wm,
        gm=rules.generate_response_function.output.gm,
        csf=rules.generate_response_function.output.csf,
        mask=rules.convert_mask_to_mrtrix_format.output
    output:
        wm=bids_output_dwi(
            desc="fodMs3t",
            suffix='wm.mif'
        ),
        gm=bids_output_dwi(
            desc="fodMs3t",
            suffix='gm.mif'
        ),
        csf=bids_output_dwi(
            desc="fodMs3t",
            suffix='csf.mif'
        )
    group: "response_generation"
    threads: 32
    resources:
        mem_mb=10000,
        runtime=15,
    envmodules:
        "mrtrix/3.0.1"
    log: "logs/compute_fiber_orientation_densities/{subject}.log"
    shell:
        'dwi2fod msmt_csd {input.dwi} '
        '{input.wm} {output.wm} {input.gm} {output.gm} {input.csf} {output.csf} '
        '-mask {input.mask} -nthreads {threads} 2> {log}'


rule normalize_fiber_orientation_densities:
    input:
        wm=FodAlgorithm(
            root=output,
            bval=input_paths['bval'],
            tissue='wm'
        ),
        gm=FodAlgorithm(
            root=output,
            bval=input_paths['bval'],
            tissue='gm'
        ),
        csf=FodAlgorithm(
            root=output,
            bval=input_paths['bval'],
            tissue='csf'
        ),
        mask=rules.convert_mask_to_mrtrix_format.output
    output:
        wm=bids_output_dwi(
            desc='fodNorm',
            suffix='wm.mif'
        ),
        gm=bids_output_dwi(
            desc='fodNorm',
            suffix='gm.mif'
        ),
        csf=bids_output_dwi(
            desc='fodNorm',
            suffix='csf.mif'
        )
    group: "response_generation"
    resources:
        runtime=2
    envmodules:
        "mrtrix/3.0.1"
    log: "logs/normalize_fiber_orientation_densities/{subject}.log"
    shell:
        'mtnormalise {input.wm} {output.wm} {input.gm} {output.gm} {input.csf} {output.csf} '
        '-mask {input.mask} 2> {log}'
