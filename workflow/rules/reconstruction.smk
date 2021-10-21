from snakebids import bids
from lib.shells import is_multi_shelled

rule run_act:
    input:
        act=rules.segment_anatomical_image.output,
        gmwmi=rules.create_seed_boundary.output,
        fod=rules.normalize_fiber_orientation_densities.output.wm
    output:
        protected(bids_output_dwi(
            desc="tracts",
            suffix='10M.tck'
        ))
    params:
        maxlength=config['tractography']['maxlength'],
        cutoff=config['tractography']['cutoff'],
        num_tracts=config['tractography']['num_tracts']
    threads: 32
    resources:
        runtime=90,
        mem_mb=4000
    log: "logs/run_act/{subject}.log"
    envmodules:
        "mrtrix/3.0.1"
    group: 'act'
    shell:
        'tckgen -act {input.act} -seed_gmwmi {input.gmwmi} '
        '-nthreads {threads} -backtrack '
        '-maxlength {params.maxlength} -cutoff {params.cutoff} -select {params.num_tracts} '
        '{input.fod} {output} 2> {log}'

rule run_sift2:
    input:
        tracks=rules.run_act.output,
        fod=rules.normalize_fiber_orientation_densities.output.wm,
        act=rules.segment_anatomical_image.output
    output:
        weights=bids_output_dwi(
            desc="sift",
            suffix='weights.txt'
        ),
        mu=bids_output_dwi(
            desc="sift",
            suffix='mu.txt'
        ),
        coeffs=bids_output_dwi(
            desc="sift",
            suffix='coeffs.txt'
        )
    threads: 16
    resources:
        mem_mb=10000,
        runtime=9
    envmodules:
        "mrtrix/3.0.1"
    log: "logs/run_sift2/{subject}.log"
    benchmark:
        "benchmarks/run_sift2/sub-{subject}.tsv"
    group: "sift"
    shell: "tcksift2 "
        "-nthreads {threads} "
        "-out_mu {output.mu} -out_coeffs {output.coeffs} "
        "{input.tracks} {input.fod} {output.weights}"