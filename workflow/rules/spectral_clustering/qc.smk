from pathlib import Path
from lib.pipenv import PipEnv
from lib.utils import XvfbRun
import itertools as it

wma_env = PipEnv(
    packages = [
        'whitematteranalysis'
    ],
    flags = config["pip-flags"],
    root = Path(work)
)

xvfb_run = XvfbRun(config.get('x11_srv', False))

wma_qc = qc/"whitematteranalysis"

rule qc_tractography_unregistered_overlap:
    input:
        data=rules.collect_registration_output.output.data,
        atlas=config['atlases']['registration_atlas'],

    output:
        directory(wma_qc/"RegTractOverlap/{subject}")

    log: f"logs/qc_tractography_unregistered_overlap/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/qc_tractography_unregistered_overlap/{'.'.join(wildcards.values())}.tsv"

    threads: 1
    resources:
        mem_mb=1000,
        runtime=30,
    shell:
        xvfb_run(
            wma_env.script(
                "wm_quality_control_tract_overlap.py "
                "{input.atlas} {input.data} {output}"
            )
        )


rule qc_tractography_clusters_initial:
    input:
        rules.tractography_spectral_clustering.output
    output:
        directory(wma_qc/"FiberClusters-Initial/{subject}")
    log: f"logs/qc_tractography_clusters/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/qc_tractography_clusters/{'.'.join(wildcards.values())}.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=30,
    shell:
        xvfb_run(
        tar.using(inputs=["{input}"])(
            wma_env.script(
                "wm_quality_control_tractography.py {input} {output}"
            )
        ))



rule qc_tractography_clusters_outliers_removed:
    input:
        rules.remove_cluster_outliers.output
    output:
        directory(wma_qc/"FiberClusters-OutliersRemoved/{subject}")
    log: f"logs/qc_tractography_clusters_outliers_removed/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/qc_tractography_clusters_outliers_removed/{'.'.join(wildcards.values())}.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=30,
    params:
    shell:
        xvfb_run(
        tar.using(inputs=["{input}"])(
            wma_env.script(
                "wm_quality_control_tractography.py {input} {output}"
            )
        ))



rule qc_tractography_anatomical_tracts:
    input:
        rules.assign_to_anatomical_tracts.output
    output:
        directory(wma_qc/"AnatomicalTracts/{subject}")
    log: f"logs/qc_tractography_anatomical_tracts/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/qc_tractography_anatomical_tracts/{'.'.join(wildcards.values())}.tsv"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=30,
    params:
    shell:
        xvfb_run(
        tar.using(inputs=["{input}"])(
            wma_env.script(
                "wm_quality_control_tractography.py {input} {output}"
            )
        ))


def qc_spectral_clustering_collector(*_):
    ret = [*it.chain(
        *(expand(
            r,
            **inputs['input_lists']['preproc_dwi']
        ) for r in [
            rules.qc_tractography_unregistered_overlap.output,
            rules.qc_tractography_clusters_initial.output,
            rules.qc_tractography_clusters_outliers_removed.output,
            rules.qc_tractography_anatomical_tracts.output,
        ])
    )]
    return ret

rule qc_spectral_clustering:
    input: qc_spectral_clustering_collector