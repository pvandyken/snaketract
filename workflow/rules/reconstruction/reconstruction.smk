from snakebids import bids
from lib.shells import is_multi_shelled

# group act:
#     group_components=2

# group sift:
#    group_components=40
rule run_act:
    input:
        act=rules.segment_anatomical_image.output,
        gmwmi=rules.create_seed_boundary.output,
        fod=rules.normalize_fiber_orientation_densities.output.wm
    output:
        bids_output_dwi(
            suffix='tractography.tck'
        )
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
        "mrtrix/3.0.1",
        "git-annex/8.20200810"
    group: 'act'
    shell:
        datalad.msg("Generate anatomically constrained tractography")(
            'tckgen -act {input.act} -seed_gmwmi {input.gmwmi} '
            '-nthreads {threads} -backtrack '
            '-maxlength {params.maxlength} -cutoff {params.cutoff} '
            '-select {params.num_tracts} '
            '{input.fod} {output} 2> {log}'
        )

rule run_sift2:
    input:
        tracks=rules.run_act.output,
        fod=rules.normalize_fiber_orientation_densities.output.wm,
        act=rules.segment_anatomical_image.output
    output:
        weights=bids_output_dwi(
            desc="weights",
            suffix='sift.txt'
        ),
        mu=bids_output_dwi(
            desc="mu",
            suffix='sift.txt'
        ),
        coeffs=bids_output_dwi(
            desc="coeffs",
            suffix='sift.txt'
        )
    threads: 16
    resources:
        mem_mb=10000,
        runtime=9
    envmodules:
        "mrtrix/3.0.1",
        "git-annex/8.20200810"
    log: "logs/run_sift2/{subject}.log"
    benchmark:
        "benchmarks/run_sift2/sub-{subject}.tsv"
    group: "sift"
    shell:
        datalad.msg("Calculate streamline weights")(
            "tcksift2 "
            "-nthreads {threads} "
            "-out_mu {output.mu} -out_coeffs {output.coeffs} "
            "{input.tracks} {input.fod} {output.weights}"
        )

def _get_image(wildcards):
    if wildcards["weight"][3:] == "FA":
        return inputs.input_path["fa"].format(**wildcards)
    raise ValueError(
        "config key 'connectome_weight' mut be set to '___FA', where ___ is one of "
        f"'avg', 'med', 'min', 'max' currently '{config['segmentation']}'"
    )

def _get_stat(wildcards):
    stat = wildcards["weight"][:3]
    mapping = {
        "avg": "mean",
        "med": "median",
        "min": "min",
        "max": "max",
    }
    if stat in mapping:
        return mapping[stat]
    raise ValueError(
        "config key 'connectome_weight' mut be set to '___FA', where ___ is one of "
        f"'avg', 'med', 'min', 'max' currently '{config['segmentation']}'"
    )

rule tck_sample:
    input:
        tracks=rules.run_act.output,
        image=_get_image
    output:
        bids_output_dwi(
            desc="{weight}",
            suffix='tractometry.csv'
        ),
    threads: 1
    resources:
        mem_mb=10000,
        runtime=9
    envmodules:
        "mrtrix/3.0.1",
        "git-annex/8.20200810"
    log: "logs/tck_sample/{subject}.{weight}.log"
    benchmark: "benchmarks/tck_sample/sub-{subject}.{weight}.tsv"
    params:
        stat=_get_stat
    shell:
        "tcksample {input.tracks} {input.image} {output} "
        "-nthreads {threads} -stat_tck {params.stat} -q"
