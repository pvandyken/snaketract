from pathlib import Path
from functools import partial
import re
from snakebids import bids

from lib.utils import xvfb_run


localrules: 
    transform_clusters_to_subject_space,
    install_whitematteranalysis


rule install_python:
    output: 
        venv=temp(directory(work+"/prepdwi_recon_venv")),
        python=work+"/prepdwi_recon_venv/bin/python"
    envmodules:
        "python/3.7"
    params:
        flags=config["pip-flags"],
        packages="whitematteranalysis"
    shell: 
        (
            "virtualenv --no-download {output.venv}  "
            "{output.python} -m pip install --upgrade pip && "
            "{output.python} -m pip install {params.flags} {params.packages}"
        )
    


rule convert_tracts_to_vtk:
    input: 
        rules.run_act.output
    
    output: 
        temp(
            work+f"/tractography/{uid}.vtk"
        )
    
    log: f"logs/convert_tracts_to_vtk/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/convert_tracts_to_vtk/{'.'.join(wildcards.values())}.tsv"

    group: "spectral_clustering"
    threads: 32
    resources:
        mem_mb=1000,
        runtime=2

    shell: "tckconvert -nthreads {threads} {input} {output}"
    

# Including {uid} at the end of registration_dir will lead to it appearing twice in the
# path, but this is necessary because Snakemake wants every output to have the wildcards
# at least once. (Main, below, needs wildcards)
registration_dir = work + f"/tractography_registration/{uid}"
registration_files = registration_dir + f"/{uid}/output_tractography"

rule tractography_registration:
    input: 
        data=rules.convert_tracts_to_vtk.output[0],
        atlas=config['atlases']['registration_atlas'],
        python=rules.install_python.output.python

    output: 
        main=temp(directory(registration_dir)),
        # These are produced implicitely, but we make them explicit for use in future
        # rules
        data=f"{registration_files}/{uid}_reg.vtk",
        inv_matrix=f"{registration_files}/itk_txform_{uid}.tfm",
        xfm=f"{registration_files}/vtk_txform_{uid}.xfm"

    log: f"logs/tractography_registration/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/tractography_registration/{'.'.join(wildcards.values())}.tsv"

    group: "spectral_clustering"
    resources:
        mem_mb=1000,
        runtime=30,

    params:
        mode="rigid_affine_fast"
    shell: 
        (
            "{input.python} wm_register_to_atlas_new.py "
            "-mode {params.mode} "
            "{input.data} {input.atlas} {output.main}"
        )


rule collect_registration_output:
    input:
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

    group: "spectral_clustering"
    resources:
        mem_mb=500,
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
        python=rules.install_python.output.python
    output: 
        directory(bids_output_dwi(
            space="ORG",
            desc="clusters800"
        ))
    log: f"logs/tractography_spectral_clustering/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/tractography_spectral_clustering/{'.'.join(wildcards.values())}.tsv"
    group: "spectral_clustering"
    threads: 32
    resources:
        mem_mb=1000,
        runtime=30,
        xvfb_run="$([[ -n \"{resources.x11_srv}\"]] && echo xvfb-run)"
    params:
        work_folder="tractography_clustering",
        results_subfolder=Path(rules.collect_registration_output.output.data).stem
    shell:
        (
            f"{xvfb_run(config)}  {{input.python}} wm_cluster_from_atlas.py "
            "-j {threads} "
            "{input.data} {input.atlas} {resources.tmpdir}{params.work_folder} && "

            "mv {resources.tmpdir}/{params.work_folder}/"
                "{params.results_subfolder} {output}"
        )
    

rule remove_cluster_outliers:
    input: 
        data=rules.tractography_spectral_clustering.output,
        atlas=config['atlases']['cluster_atlas'],
        python=rules.install_python.output.python
    output: 
        directory(bids_output_dwi(
            space="ORG",
            desc="clusters800",
            suffix="outliersRemoved"
        ))
    log: f"logs/remove_cluster_outliers/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/remove_cluster_outliers/{'.'.join(wildcards.values())}.tsv"
    group: "spectral_clustering"
    threads: 32
    resources:
        mem_mb=1000,
        runtime=30,
    params:
        work_folder="tractography_outlier_removal",
        results_subfolder=Path(rules.tractography_spectral_clustering.output[0]).stem
    shell: 
        (
            "{input.python} wm_cluster_remove_outliers.py "
            "-j {threads} "
            "{input.data} {input.atlas} {resources.tmpdir}{params.work_folder} && "

            "mv {resources.tmpdir}/{params.work_folder}/"
                "{params.results_subfolder}_outlier_removed {output}"
        )
    

rule assess_cluster_location_by_hemisphere:
    input: 
        data=rules.remove_cluster_outliers.output,
        atlas=config['atlases']['cluster_atlas'],
        python=rules.install_python.output.python

    output: 
        bids_output_dwi(
            space="ORG",
            desc="clusters800",
            suffix="assignedHemispheres.complete"
        )

    log: f"logs/assess_cluster_location_by_hemisphere/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/assess_cluster_location_by_hemisphere/{'.'.join(wildcards.values())}.tsv"

    group: "spectral_clustering"
    resources:
        mem_mb=1000,
        runtime=30,

    shell: 
        (
            "{input.python} wm_assess_cluster_location_by_hemisphere.py "
            "{input.data} -clusterLocationFile "
            "{input.atlas}/cluster_hemisphere_location.txt && "

            "touch {output}"
        )
    
    
rule transform_clusters_to_subject_space:
    input: 
        hemisphereAssignment=rules.assess_cluster_location_by_hemisphere.output,
        data=rules.remove_cluster_outliers.output,
        transform=rules.collect_registration_output.output.inv_matrix,
        python=rules.install_python.output.python

    output: 
        temp(directory(work+"/transformed_clusters/" + uid + "/"))

    log: f"logs/transform_clusters_to_subject_space/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/transform_clusters_to_subject_space/{'.'.join(wildcards.values())}.tsv"

    container: "docker://slicer/slicer-base"
    envmodules: "slicer/4.11.20210226"

    threads: 32
    resources:
        mem_mb=1000,
        runtime=30,

    shell: 
        (
            f"{xvfb_run(config)} {{input.python}} wm_harden_transform.py "
            "-j {threads} -i -t {input.transform} "
            "{input.data} {output} $(which Slicer)"
        )


rule separate_clusters_by_cluster:
    input: 
        data=rules.transform_clusters_to_subject_space.output,
        python=rules.install_python.output.python

    output: 
        directory(bids_output_dwi(
            space="individual",
            desc="clusters800",
            suffix="sorted"
        ))

    log: f"logs/separate_clusters_by_cluster/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/separate_clusters_by_cluster/{'.'.join(wildcards.values())}.tsv"

    group: "spectral_clustering"
    resources:
        mem_mb=1000,
        runtime=30,

    shell: 
        (
            "{input.python} wm_separate_clusters_by_hemisphere.py {input.data} {output}"
        )
    

rule assign_to_anatomical_tracts:
    input: 
        data=rules.separate_clusters_by_cluster.output,
        atlas=config["atlases"]["cluster_atlas"],
        python=rules.install_python.output.python

    output: 
        directory(bids_output_dwi(
            space="individual",
            desc="tracts",
            suffix="73"
        ))

    log: f"logs/assign_to_anatomical_tracts/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/assign_to_anatomical_tracts/{'.'.join(wildcards.values())}.tsv"

    group: "spectral_clustering"
    resources:
        mem_mb=1000,
        runtime=30,

    shell:
        (
            "{input.python} wm_append_clusters_to_anatomical_tracts.py "
            "{input.data} {input.atlas} {output}"
        )
