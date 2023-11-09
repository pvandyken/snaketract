from snakebids import bids
from lib.shells import FodAlgorithm

# group response_generation:
#     group_components=14

rule convert_dwi_to_mrtrix_format:
    input:
        dwi=inputs.input_path['preproc_dwi'],
        bvec=inputs.input_path['bvec'],
        bval=inputs.input_path['bval']

    output:
        tempout("convert_dwi_to_mrtrix_format", inputs["preproc_dwi"], ".mif")

    envmodules:
        "mrtrix/3.0.1",
        "git-annex/8.20200810"
    log: log("convert_dwi_to_mrtrix_format", inputs['preproc_dwi'])
    group: "response_generation"
    resources:
        mem_mb=lambda wcards, input: input.size_mb * 3.2
    shell:
        'mrconvert {input.dwi} {output} -fslgrad {input.bvec} {input.bval} 2> {log}'


rule generate_tournier_response:
    input:
        dwi=rules.convert_dwi_to_mrtrix_format.output,
        mask=inputs['brainmask'].path,
    output:
        response=bids_output_dwi(
            model="CSD",
            suffix="fod.txt",
        ),
        voxels=bids_output_dwi(
            model="CSD",
            desc="voxels",
            suffix="fod.nii.gz"
        ),
    log: log("generate_tournier_response", inputs['preproc_dwi'])
    benchmark: benchmark("generate_tournier_response", inputs['preproc_dwi'])
    envmodules:
        "mrtrix/3.0.1",
    group: "response_generation"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
    params:
        lmax=lambda wcards: (
            f"-lmax {lmax}"
            if (lmax := config.get('response_generation_lmax')) else ""
        )
    shell:
        "dwi2response tournier {input.dwi} {output.response} "
        "-mask {input.mask} -voxels {output.voxels} {params.lmax} "
        "-scratch {resources.tmpdir} 2> {log}"


rule generate_dhollander_response:
    input:
        dwi=rules.convert_dwi_to_mrtrix_format.output,
        mask=inputs['brainmask'].path,
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
            suffix='fod.nii.gz'
        )
    group: "response_generation"
    log: log("generate_dhollander_response", inputs['preproc_dwi'])
    resources:
        runtime=2
    envmodules:
        "mrtrix/3.0.1",
        "git-annex/8.20200810"
    params:
        lmax=lambda wcards: (
            f"-lmax {lmax}"
            if (lmax := config.get('response_generation')) else ""
        )
    shell:
        'dwi2response dhollander {input.dwi} '
        '{output.wm} {output.gm} {output.csf} '
        '-voxels {output.voxels} -mask {input.mask} '
        '{params.lmax} '
        '-scratch {resources.tmpdir} '
        '2> {log}'


rule compute_ss3t_fiber_orientation_densities:
    input:
        dwi=rules.convert_dwi_to_mrtrix_format.output,
        wm=rules.generate_dhollander_response.output.wm,
        gm=rules.generate_dhollander_response.output.gm,
        csf=rules.generate_dhollander_response.output.csf,
        mask=inputs["brainmask"].path,
    output:
        wm=bids_output_dwi(
            model="CSD",
            desc="ss3t",
            label="wm",
            suffix='fod.nii.gz'
        ),
        gm=bids_output_dwi(
            model="CSD",
            desc="ss3t",
            label="gm",
            suffix='fod.nii.gz'
        ),
        csf=bids_output_dwi(
            model="CSD",
            desc="ss3t",
            label="csf",
            suffix='fod.nii.gz'
        )
    group: "response_generation"
    threads: 16
    resources:
        mem_mb=1000,
        runtime=96,
    container:
        "docker://vnmd/mrtrix3tissue_5.2.8:latest"
    log: log("compute_ss3t_fiber_orientation_densities", inputs['preproc_dwi'])
    benchmark: benchmark("compute_ss3t_fiber_orientation_densities", inputs['preproc_dwi'])
    shell:
        'ss3t_csd_beta1 {input.dwi} '
        '{input.wm} {output.wm} {input.gm} {output.gm} {input.csf} {output.csf} '
        '-mask {input.mask} -nthreads {threads} -scratch ' + str(work)



rule compute_ms3t_fiber_orientation_densities:
    input:
        dwi=rules.convert_dwi_to_mrtrix_format.output,
        wm=rules.generate_dhollander_response.output.wm,
        gm=rules.generate_dhollander_response.output.gm,
        csf=rules.generate_dhollander_response.output.csf,
        mask=inputs['brainmask'].path
    output:
        wm=bids_output_dwi(
            model="CSD",
            desc="ms3t",
            label="wm",
            suffix='fod.nii.gz'
        ),
        gm=bids_output_dwi(
            model="CSD",
            desc="ms3t",
            label="gm",
            suffix='fod.nii.gz'
        ),
        csf=bids_output_dwi(
            model="CSD",
            desc="ms3t",
            label="csf",
            suffix='fod.nii.gz'
        )
    group: "response_generation"
    threads: 32
    resources:
        mem_mb=10000,
        runtime=15,
    envmodules:
        "mrtrix/3.0.1",
        "git-annex/8.20200810"
    log: log("compute_ms3t_fiber_orientation_densities", inputs['preproc_dwi'])
    benchmark: benchmark("compute_ms3t_fiber_orientation_densities", inputs['preproc_dwi'])
    params:
        ss2t=lambda wcards: (
            not is_multi_shelled(inputs['bval'].path.format(**wcards)) or ""
        )
    shell:
        """
        if [[ -n "{params.ss2t}" ]]; then
            gm_io=''
            touch {output.gm}
        else
            gm_io="{input.gm} {output.gm}"
        fi
        dwi2fod msmt_csd {input.dwi} \
            {input.wm} {output.wm} $gm_io {input.csf} {output.csf} \
            -mask {input.mask} -nthreads {threads} 2> {log}
        """

rule compute_csd_fiber_orientation_densities:
    input:
        response=rules.generate_tournier_response.output['response'],
        mask=inputs['brainmask'].path,
    output:
        bids_output_dwi(
            model="CSD",
            suffix="fod.nii.gz",
        )
    log: log("compute_csd_fiber_orientation_densities", inputs['preproc_dwi'])
    benchmark: benchmark("compute_csd_fiber_orientation_densities", inputs['preproc_dwi'])
    group: "response_generation"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
    shell:
        "dwi2fod csd {input.response} {output} -mask {input.mask} "
        "-nthreads {threads} 2> {log}"


_single_shell_alg = config.get('single_shell_algorithm', 'ss3t')
rule normalize_fiber_orientation_densities:
    input:
        wm=FodAlgorithm(
            root=output,
            bval=inputs.input_path['bval'],
            tissue='wm',
            single_shell_alg=_single_shell_alg,
        ),
        gm=FodAlgorithm(
            root=output,
            bval=inputs.input_path['bval'],
            tissue='gm',
            single_shell_alg=_single_shell_alg,
        ),
        csf=FodAlgorithm(
            root=output,
            bval=inputs.input_path['bval'],
            tissue='csf',
            single_shell_alg=_single_shell_alg,
        ),
        mask=inputs['brainmask'].path
    output:
        wm=bids_output_dwi(
            model="CSD",
            desc="norm",
            label="wm",
            suffix='fod.nii.gz'
        ),
        gm=bids_output_dwi(
            model="CSD",
            desc="norm",
            label="gm",
            suffix='fod.nii.gz'
        ),
        csf=bids_output_dwi(
            model="CSD",
            desc="norm",
            label="csf",
            suffix='fod.nii.gz'
        )
    group: "response_generation"
    resources:
        runtime=2
    envmodules:
        "mrtrix/3.0.1",
        "git-annex/8.20200810"
    log: log("normalize_fiber_orientation_densities", inputs['preproc_dwi'])
    shell:
        # datalad.msg("Normalize fod functions"),
        """
        if [[ -s "{input.gm}" ]]; then
            gm_io="{input.gm} {output.gm}"
        else
            gm_io=''
            touch "{output.gm}"
        fi
        mtnormalise \
            {input.wm} {output.wm} $gm_io {input.csf} {output.csf} \
            -mask {input.mask} 2> {log}
        """
