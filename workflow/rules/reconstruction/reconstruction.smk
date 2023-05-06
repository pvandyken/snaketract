from snakebids import bids
from lib.shells import is_multi_shelled

# group act:
#     group_components=2

# group sift:
#    group_components=40


def _select_fod_algorithm(wcards):
    if (
        config['tractography'].get('single_shell_algorithm', 'ss3t') == "csd"
        and not is_multi_shelled(inputs['bval'].path.format(**wcards))
    ):
        return rules.compute_csd_fiber_orientation_densities.output[0].format(**wcards)
    return rules.normalize_fiber_orientation_densities.output[0].format(**wcards)


rule run_act:
    input:
        act=rules.get_5tt_segmentation.output,
        gmwmi=rules.create_seed_boundary.output,
        fod=_select_fod_algorithm,
    output:
        bids_output_dwi(
            rec="{rec}",
            suffix='tractography.tck'
        )
    params:
        maxlength=config['tractography']['maxlength'],
        cutoff=config['tractography']['cutoff'],
        num_tracts=config['tractography']['num_tracts'],
        rec=lambda wcards: {"iFOD2": "iFOD2"}[wcards['rec']]
    threads: 32
    resources:
        runtime=int(max(config['tractography']['num_tracts'] * 90 / 10_000_000, 10)),
        mem_mb=4000,
    log: "logs/run_act/" + uid + ".{rec}.log"
    envmodules:
        "mrtrix/3.0.1",
        "git-annex/8.20200810"
    group: 'act'
    shell:
        # datalad.msg("Generate anatomically constrained tractography"),
        'tckgen -act {input.act} -seed_gmwmi {input.gmwmi} '
        '-nthreads {threads} -backtrack '
        '-algorithm {params.rec} '
        '-maxlength {params.maxlength} -cutoff {params.cutoff} '
        '-select {params.num_tracts} '
        '{input.fod} {output} 2> {log}'

rule run_sift2:
    input:
        tracks=rules.run_act.output,
        fod=_select_fod_algorithm,
        act=rules.get_5tt_segmentation.output
    output:
        weights=bids_output_dwi(
            rec="{rec}",
            desc="weights",
            suffix='sift.txt'
        ),
        mu=bids_output_dwi(
            rec="{rec}",
            desc="mu",
            suffix='sift.txt'
        ),
        coeffs=bids_output_dwi(
            rec="{rec}",
            desc="coeffs",
            suffix='sift.txt'
        )
    threads: 16
    resources:
        mem_mb=30000,
        runtime=20
    envmodules:
        "mrtrix/3.0.1",
        "git-annex/8.20200810"
    log: "logs/run_sift2/"+uid+".{rec}.log"
    benchmark:
        "benchmarks/run_sift2/"+uid+".{rec}.tsv"
    group: "sift"
    shell:
        # datalad.msg("Calculate streamline weights"),
        "tcksift2 "
        "-nthreads {threads} "
        "-out_mu {output.mu} -out_coeffs {output.coeffs} "
        "{input.tracks} {input.fod} {output.weights}"
