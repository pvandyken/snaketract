# group: connectome:
#       connected_components: 288


def _get_segmentation(wildcards):
    if wildcards["atlas"] == "bn246":
        return expand(
            rules.get_bn246_in_subject_space.output,
            **wildcards,
        )
    if wildcards["atlas"][:3] == "sch":
        return warp_between_t1_and_mni_via_ants(
            group="connectome",
            moving=ref_images['schaefer'][int(wildcards["atlas"][3:])],
            fixed=path_store.register(inputs['t1'].path),
            interpolation="NearestNeighbor",
            **wildcards
        )
    raise ValueError(
        "config key 'segmentation' mut be set to one of 'bn210' or 'bn246', currently "
        f"'{config['segmentation']}'"
    )

def _get_weights(wcards):
    if "log" in wcards["weight"]:
        raise Exception()
    if wcards["weight"] == "sift2":
        return rules.run_sift2.output.weights.format(**wcards)
    if wcards["weight"][3:] in VALID_WEIGHTS:
        return rules.tck_sample.output[0].format(**wcards)
    raise ValueError(
        "config key 'connectome_weight' mut be set to 'sift2' or '___FA', where ___ is "
        f"one of 'avg', 'med', 'min', 'max' currently '{config['segmentation']}'"
    )

def _get_weight_arg(wcards, input):
    if "sift2" in wcards["weight"]:
        scale_invnodevol = (
            " -scale_invnodevol" if wcards["weight"] == "sift2Density" else ""
        )
        return f"-tck_weights_in {input.tck_weights}{scale_invnodevol}"
    if wcards["weight"][3:] in VALID_WEIGHTS:
        return f"-scale_file {input.tck_weights} -stat_edge mean"
    raise ValueError(
        "config key 'connectome_weight' mut be set to 'sift2', 'sift2Density', or "
        "'<stat>FA', where <stat> is one of 'avg', 'med', 'min', 'max' currently "
        f"'{config['segmentation']}'"
    )


rule get_connectome:
    input:
        tracks=rules.run_act.output,
        tck_weights=_get_weights,
        nodes=_get_segmentation,
    output:
        connectome=bids_output_dwi(
            rec="{rec}",
            atlas="{atlas}",
            desc="{weight}",
            suffix="connectome.csv",
        ),
        assignments=bids_output_dwi(
            rec="{rec}",
            atlas="{atlas}",
            desc="{weight}",
            suffix="assignments.csv",
        ),

    group: "connectome"
    envmodules:
        "git-annex/8.20200810",
        "mrtrix"

    threads: 4
    resources:
        mem_mb=5000,
        runtime=int(max(int(config['tractography_numtracts']) * 6 / 10_000_000, 1)),
    log: log("get_connectome", inputs['preproc_dwi'], "atlas", "weight", "rec")
    benchmark: benchmark("get_connectome", inputs['preproc_dwi'], "atlas", "weight", "rec")
    params:
        weight=_get_weight_arg
    shell:
        "tck2connectome {input.tracks} {input.nodes} {output.connectome} "
        "-q -symmetric -keep_unassigned "
        "{params.weight} -out_assignments {output.assignments}"

def _get_log_connectome(wcards):
    if wcards["weight"][:3] != "log":
        return ""
    wcards = {
        **wcards,
        "weight": wcards["weight"][3:],
    }
    result = rules.get_connectome.output.connectome.format(**wcards)
    return result

rule get_log_connectome:
    input:
        connectome=_get_log_connectome
    output:
        connectome=bids_output_dwi(
            rec="{rec}",
            atlas="{atlas}",
            desc="{weight}",
            suffix="connectome.csv",
        ),
    threads: 1
    group: "connectome"
    run:
        data = np.loadtxt(input["connectome"], delimiter=",")
        masked = np.ma.masked_equal(data, 0)
        log = np.ma.log10(masked / masked.min()).filled(0)
        np.savetxt(output["connectome"], log, fmt="%.11f", delimiter=",")

    

    
