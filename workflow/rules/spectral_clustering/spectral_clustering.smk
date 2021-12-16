from pathlib import Path
from functools import partial
import re
from snakebids import bids

from lib.utils import XvfbRun, Tar, display
from lib.pipenv import PipEnv


localrules:
    install_python


# group tract_registration:
#   num_components: 16
#   total_runtime: 3:00
#   total_mem_mb: 480,000
#   cores: 16

# group spectral_clustering:
#   num_components: 24
#   total_runtime: 12:00
#   total_mem_mb: 500,000
#   cores: 32

# group cluster_postprocess:
#   num_components: 160
#   total_runtime: 12:00
#   total_mem_mb: 16,000
#   cores: 32


rule convert_tracts_to_vtk:
    input:
        rules.run_act.output

    output:
        temp(
            work+f"/tractography/{uid}.vtk"
        )

    log: f"logs/convert_tracts_to_vtk/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/convert_tracts_to_vtk/{'.'.join(wildcards.values())}.tsv"

    envmodules: "mrtrix/3.0.1"

    group: "tract_registration"
    resources:
        mem_mb=500,
        runtime=30 # for 10M fibres

    shell: "tckconvert {input} {output}"


# Including {uid} at the end of registration_dir will lead to it appearing twice in the
# path, but this is necessary because Snakemake wants every output to have the wildcards
# at least once. (Main, below, needs wildcards)
registration_dir = work + f"/tractography_registration/{uid}"
registration_files = registration_dir + f"/{uid}/output_tractography"

rule tractography_registration:
    input:
        data=rules.convert_tracts_to_vtk.output[0],
        atlas=config['atlases']['registration_atlas'],

    output:
        main=temp(directory(registration_dir)),
        # These are produced implicitely, but we make them explicit for use in future
        # rules
        data=f"{registration_files}/{uid}_reg.vtk",
        inv_matrix=f"{registration_files}/itk_txform_{uid}.tfm",
        xfm=f"{registration_files}/vtk_txform_{uid}.xfm"

    log: f"logs/tractography_registration/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/tractography_registration/{'.'.join(wildcards.values())}.tsv"

    group: "tract_registration"
    resources:
        mem_mb=60000,
        runtime=44,

    params:
        mode="rigid_affine_fast",
    shell:
        wma_env.script(
            "wm_register_to_atlas_new.py "
            "-mode {params.mode} "
            "{input.data} {input.atlas} {output.main}"
        )


rule collect_registration_output:
    input:
        main=rules.tractography_registration.output.main,
        data=rules.tractography_registration.output.data,
        inv_matrix=rules.tractography_registration.output.inv_matrix,
        xfm=rules.tractography_registration.output.xfm

    output:
        data=bids_output_dwi(
            space="ORG",
            desc="tracts",
            suffix="10M.vtk"
        ),
        inv_matrix=bids_output_dwi(
            space="ORG",
            desc="tracts10M",
            suffix="inverseTransform.tfm"
        ),
        xfm=bids_output_dwi(
            space="ORG",
            desc="tracts10M",
            suffix="vtkTxform.xfm"
        )

    log: f"logs/collect_registration_output/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/collect_registration_output/{'.'.join(wildcards.values())}.tsv"

    group: "tract_registration"
    resources:
        mem_mb=100,
        runtime=1,

    shell:
        (
            "cp {input.data} {output.data} && "
            "cp {input.inv_matrix} {output.inv_matrix} && "
            "cp {input.xfm} {output.xfm}"
        )


rule tractography_spectral_clustering:
    input:
        data=rules.collect_registration_output.output.data,
        atlas=config['atlases']['cluster_atlas'],
    output:
        bids_output_dwi(
            space="ORG",
            desc="clusters800.tar.gz"
        )
    log: f"logs/tractography_spectral_clustering/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/tractography_spectral_clustering/{'.'.join(wildcards.values())}.tsv"

    group: "spectral_clustering"
    threads: 16
    resources:
        mem_mb=250000,
        runtime=60,

    params:
        work_folder=work + "/tractography_clustering",
        results_subfolder=Path(rules.collect_registration_output.output.data).stem
    shell:
        xvfb_run(
        tar.using(outputs=["{output}"])(
            wma_env.script(
                "wm_cluster_from_atlas.py "
                "-j {threads} "
                "{input.data} {input.atlas} {params.work_folder} && "

                "mv {params.work_folder}/{params.results_subfolder} {output}"
            )
        ))


rule remove_cluster_outliers:
    input:
        data=rules.tractography_spectral_clustering.output,
        atlas=config['atlases']['cluster_atlas'],
    output:
        bids_output_dwi(
            space="ORG",
            desc="clusters800",
            suffix="outliersRemoved.tar.gz"
        )
    log: f"logs/remove_cluster_outliers/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/remove_cluster_outliers/{'.'.join(wildcards.values())}.tsv"

    group: "cluster_postprocess"
    threads: 32
    resources:
        mem_mb=10000,
        runtime=4,
    params:
        work_folder=work + "/tractography_outlier_removal",
        results_subfolder=Path(rules.tractography_spectral_clustering.output[0]).name
    shell:
        tar.using(inputs = ["{input.data}"], outputs = ["{output}"])(
            wma_env.script(
                "wm_cluster_remove_outliers.py "
                "-j {threads} "
                "{input.data} {input.atlas} {params.work_folder} && "

                "mv "
                "{params.work_folder}/{params.results_subfolder}_outlier_removed/* {output}/"
            )
        )


rule assess_cluster_location_by_hemisphere:
    input:
        data=rules.remove_cluster_outliers.output,
        atlas=config['atlases']['cluster_atlas'],

    output:
        bids_output_dwi(
            space="ORG",
            desc="clusters800",
            suffix="assignedHemispheres.complete"
        )

    log: f"logs/assess_cluster_location_by_hemisphere/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/assess_cluster_location_by_hemisphere/{'.'.join(wildcards.values())}.tsv"

    group: "cluster_postprocess"
    resources:
        mem_mb=500,
        runtime=13,

    shell:
        tar.using(modify=["{input.data}"])(
            wma_env.script(
                "wm_assess_cluster_location_by_hemisphere.py "
                "{input.data} -clusterLocationFile "
                "{input.atlas}/cluster_hemisphere_location.txt && "

                "touch {output}"
            )
        )


rule transform_clusters_to_subject_space:
    input:
        hemisphereAssignment=rules.assess_cluster_location_by_hemisphere.output,
        data=rules.remove_cluster_outliers.output,
        transform=rules.collect_registration_output.output.inv_matrix,

    output:
        temp(directory(work+"/transformed_clusters/" + uid + "/"))

    log: f"logs/transform_clusters_to_subject_space/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/transform_clusters_to_subject_space/{'.'.join(wildcards.values())}.tsv"

    container: "docker://slicer/slicer-base"
    envmodules: "slicer/4.11.20210226"

    group: "cluster_postprocess"
    resources:
        mem_mb=500,
        runtime=1,

    shell:
        xvfb_run(
        tar.using(inputs=["{input.data}"])(
            wma_env.make_venv(
                f"export PATH={wma_env.bin}:$PATH && "
                f"{wma_env.bin}/wm_harden_transform.py "
                "-i -t {input.transform} "
                "{input.data} {output} $(which Slicer)"
            )
        ))


rule separate_clusters_by_hemisphere:
    input:
        rules.transform_clusters_to_subject_space.output,

    output:
        bids_output_dwi(
            space="individual",
            desc="clusters800",
            suffix="sorted.tar.gz"
        )

    log: f"logs/separate_clusters_by_cluster/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/separate_clusters_by_cluster/{'.'.join(wildcards.values())}.tsv"

    group: "cluster_postprocess"
    resources:
        mem_mb=100,
        runtime=1,

    shell:
        tar.using(outputs=["{output}"])(
            wma_env.script(
                "wm_separate_clusters_by_hemisphere.py {input} {output}"
            )
        )

rule assign_to_anatomical_tracts:
    input:
        data=rules.separate_clusters_by_hemisphere.output,
        atlas=config["atlases"]["cluster_atlas"],

    output:
        bids_output_dwi(
            space="individual",
            desc="tracts",
            suffix="73.tar.gz"
        )

    log: f"logs/assign_to_anatomical_tracts/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/assign_to_anatomical_tracts/{'.'.join(wildcards.values())}.tsv"

    group: "cluster_postprocess"
    resources:
        mem_mb=500,
        runtime=1,

    shell:
        tar.using(inputs=["{input.data}"], outputs=["{output}"])(
            wma_env.script(display(
                "wm_append_clusters_to_anatomical_tracts.py "
                "{input.data} {input.atlas} {output}"
            ))
        )
