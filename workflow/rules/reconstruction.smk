from snakebids import bids
from lib.shells import is_multi_shelled

rule run_act:
    input:
        act=bids(root=work,
            datatype='anat',
            suffix="5tt.mif",
            **wildcards),
        gmwmi=bids(root=work,
            datatype='anat',
            suffix="gmwmi.mif",
            **wildcards),
        fod=bids(root=work,
                datatype='dwi',
                desc='norm',
                suffix='wmfod.mif',
                **wildcards)
    output:
        protected(bids(root=output,
            datatype='dwi',
            desc="tracks",
            suffix='10M.tck',
            **wildcards))
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
        act=bids(root=work,
            datatype='anat',
            suffix="5tt.mif",
            **wildcards)
    output:
        weights=bids(root=output,
            datatype='dwi',
            desc="sift",
            suffix='weights.txt',
            **wildcards),
        mu=bids(root=output,
            datatype='dwi',
            desc="sift",
            suffix='mu.txt',
            **wildcards),
        coeffs=bids(root=output,
            datatype='dwi',
            desc="sift",
            suffix='coeffs.txt',
            **wildcards)
    threads: 32
    resources:
        mem_mb=125000,
        runtime=15
    log: "logs/run_sift2/{subject}.log"
    benchmark:
        "benchmarks/run_sift2/sub-{subject}.tsv"
    group: "sift"
    shell: "tcksift2 "
        "-nthreads {threads} "
        "-out_mu {output.mu} -out_coeffs {output.coeffs} "
        "{input.tracks} {input.fod} {output.weights}"