from pathlib import Path
from snakebids import bids


localrules: transform_clusters_to_subject_space

uid = '.'.join(wildcards.values())

rule convert_tracts_to_vtk:
    input: 
        rules.run_act.output
    
    output: 
        temp(
            work+f"/tractography/{uid}.vtk"
        )
    
    log: f"logs/convert_tracts_to_vtk/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/convert_tracts_to_vtk/{'.'.join(wildcards.values())}.tsv"
    conda:
    container:
    envmodule:
    group:
    threads: 1
    resources:
        mem_mb=1000,
        runtime=2,
    params:
    shell: (
        "tckconvert {input} {output}"
    )


rule tractography_registration:
    input: 
        ata=rules.convert_tracts_to_vtk.output,
        atlas=config['tract-segmentation']['atlas']
    output: 
        temp(directory(
            work+"/tractography_registration"
        ))
    log: f"logs/tractography_registration/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/tractography_registration/{'.'.join(wildcards.values())}.tsv"
    conda:
    container:
    envmodule:
    group:
    threads: 1
    resources:
        mem_mb=1000,
        runtime=30,
    params:
    shell: (
        "wm_register_to_atlas_new.py "
        "-mode {mode} "
        "{input.data} {input.atlas} {output}"
    )

registration_files = rules.tractography_registration.output + f"/{uid}/output_tractography"
rule collect_registration_output:
    input:
        data=f"{registration_files}/{uid}_reg.vtk",
        inv_matrix=f"{registration_files}/itk_txform_{uid}.tfm",
        xfm=f"{registration_files}/vtk_txform_{uid}.tfm"
    output:
        data=bids(
            root=output,
            datatype='dwi',
            desc="registered"
            suffix="10M.vtk",
            **wildcards
        ),
        inv_matrix=bids(
            root=output,
            datatype="dwi",
            desc="inverse",
            suffix="transform.tfm",
            **wildcards
        ),
        xfm=bids(
            root=output,
            datatype="dwi",
            desc
        )
    log: f"logs/collect_registration_output/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/collect_registration_output/{'.'.join(wildcards.values())}.tsv"
    conda:
    container:
    envmodule:
    group:
    threads: 1
    resources:
        mem_mb=500,
        runtime=1,
    shell: (
        "cp {input.data} {ouput.data} && "
        "cp {input.inv_matrix} {output.inv_matrix} && "
        "cp {input.xfm} {output.xfm}"
    )


registration_stem = Path(rules.collect_registration_output.output.data).stem()
rule tractography_spectral_clustering:
    input: 
        data=rules.collect_registration_output.output.data,
        atlas=config['tract-segmentation']['atlas']
    output: 
        bids(
            root=output,
            datatype="dwi",
            desc="clusters",
            suffix="800",
            **wildcards
        )
    log: f"logs/tractography_spectral_clustering/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/tractography_spectral_clustering/{'.'.join(wildcards.values())}.tsv"
    conda:
    container:
    envmodule:
    group:
    threads: 1
    resources:
        mem_mb=1000,
        runtime=30,
    params:
        work_folder="{resources.tmpdir}/tractography_clustering"
    shell: (
        "wm_cluster_from_atlas.py {input.data} {input.atlas} {params.work_folder} && "
        f"mv {{params.work_folder}}/{registration_stem} {{output}}"
    )




rule assess_cluster_location_by_hemisphere:
    input: 
        rules.tractography_spectral_clustering
    output: 
        output
    log: f"logs/assess_cluster_location_by_hemisphere/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/assess_cluster_location_by_hemisphere/{'.'.join(wildcards.values())}.tsv"
    conda:
    container:
    envmodule:
    group:
    threads: 1
    resources:
        mem_mb=1000,
        runtime=30,
    params:
    shell: (
        command
    )


rule transform_clusters_to_subject_space:
    input:
    output:
    log:
    benchmark:
    conda:
    container:
    envmodules:
    group:
    threads:
    resources:
    shell:


rule separate_clusters_by_hemisphere:
    input:
    output:
    log:
    benchmark:
    conda:
    container:
    envmodules:
    group:
    threads:
    resources:
    shell:


rule assign_to_anatomical_tracts:
    input:
    output:
    log:
    benchmark:
    conda:
    container:
    envmodules:
    group:
    threads:
    resources:
    shell: