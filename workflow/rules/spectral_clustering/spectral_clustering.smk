from pathlib import Path


# group convert_tracts_to_vtk:
#   num_components: 192
#   total_runtime: 3:00
#   total_mem_mb: 16_000
#   cores: 32

# group tract_registration:
#   num_components: 64
#   total_runtime: 11:44
#   total_mem_mb: 240_000
#   cores: 4

# group spectral_clustering:
#   num_components: 2
#   total_runtime: 3:00
#   total_mem_mb: 250_000
#   cores: 8

# group cluster_outlier_removal:
#   num_components: 18
#   total_runtime: 3:00
#   total_mem_mb: 3000
#   cores: 32

# group cluster_postprocess:
#   num_components: 180
#   total_runtime: 3:00
#   total_mem_mb: 112_000
#   cores: 32


rule convert_tracts_to_vtk:
    input: rules.run_act.output

    output: temp(shared_work/uid/"tractography.vtk")

    log: f"logs/convert_tracts_to_vtk/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/convert_tracts_to_vtk/{'.'.join(wildcards.values())}.tsv"

    envmodules:
        "mrtrix/3.0.1",
        "git-annex/8.20200810"

    group: "convert_tracts_to_vtk"
    resources:
        mem_mb=500,
        runtime=30 # for 10M fibres

    shell: datalad("tckconvert {input} {output}")


# Including {uid} at the end of registration_dir will lead to it appearing twice in the
# path, but this is necessary because Snakemake wants every output to have the wildcards
# at least once. (Main, below, needs wildcards)
registration_dir = work/"tractography_registration"/uid
registration_output = Path(rules.convert_tracts_to_vtk.output[0]).stem
registration_files = registration_dir/registration_output/"output_tractography"

rule tractography_registration:
    input:
        data=rules.convert_tracts_to_vtk.output,
        atlas=config['atlases']['registration_atlas'],

    output:
        data=temp(shared_work/uid/"registered_tractography.vtk"),
        transform_tfm=bids(
            root=output,
            datatype="dwi",
            _from="T1w",
            to="ORG",
            mode="points",
            suffix="xfm.tfm",
            **wildcards
        ),
        transform_xfm=bids(
            root=output,
            datatype="dwi",
            _from="T1w",
            to="ORG",
            mode="points",
            suffix="xfm.xfm",
            **wildcards
        )

    log: f"logs/tractography_registration/{'.'.join(wildcards.values())}.log"
    benchmark:
        f"benchmarks/tractography_registration/{'.'.join(wildcards.values())}.tsv"

    envmodules:
        'python/3.7',
        "git-annex/8.20200810"

    group: "tract_registration"
    resources:
        mem_mb=60000,
        runtime=44,

    params:
        mode="rigid_affine_fast",

        # Temporary files
        main=str(registration_dir),
        transformed_data=f"{registration_files}/{registration_output}_reg.vtk",
        transform_tfm=f"{registration_files}/itk_txform_{registration_output}.tfm",
        transform_xfm=f"{registration_files}/vtk_txform_{registration_output}.xfm"
    shell:
        boost(
            datalad.msg("Register tractography to ORG atlas"),
            wma_env.script,
            (
                "wm_register_to_atlas_new.py "
                "-mode {params.mode} "
                "{input.data} {input.atlas} {params.main} && "

                "mv {params.transformed_data} {output.data} && "
                "mv {params.transform_tfm} {output.transform_tfm} && "
                "mv {params.transform_xfm} {output.transform_xfm}"
            )
        )


rule tractography_spectral_clustering:
    input:
        data=rules.tractography_registration.output.data,
        atlas=config['atlases']['cluster_atlas'],
    output:
        bids_output_dwi(
            space="ORG",
            atlas="ORG",
            suffix="clusters.tar.gz"
        )
    log: f"logs/tractography_spectral_clustering/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/tractography_spectral_clustering/{'.'.join(wildcards.values())}.tsv"

    envmodules:
        'python/3.7',
        "git-annex/8.20200810"

    group: "spectral_clustering"
    threads: 8
    resources:
        mem_mb=500000,
        runtime=90,

    params:
        work_folder=str(work/"tractography_clustering"/uid),
        results_subfolder=Path(rules.tractography_registration.output.data).stem
    shell:
        boost(
            datalad.msg("Cluster tracts with spectral clustering"),
            xvfb_run,
            tar.using(outputs=["{output}"]),
            wma_env.script,
            (
                "wm_cluster_from_atlas.py "
                "-j {threads} "
                "{input.data} {input.atlas} {params.work_folder} && "

                "mv {params.work_folder}/{params.results_subfolder}/* {output}"
            )
        )


rule remove_cluster_outliers:
    input:
        data=rules.tractography_spectral_clustering.output,
        atlas=config['atlases']['cluster_atlas'],
    output:
        temp(shared_work/uid/"remove-custer-outliers.nii.gz")
    log: f"logs/remove_cluster_outliers/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/remove_cluster_outliers/{'.'.join(wildcards.values())}.tsv"

    envmodules:
        'python/3.7',
        "git-annex/8.20200810"

    group: "cluster_outlier_removal"
    threads: 32
    resources:
        mem_mb=15000,
        runtime=10,
    params:
        work_folder=work/"tractography_outlier_removal",
        results_subfolder=Path(rules.tractography_spectral_clustering.output[0]).name
    shell:
        boost(
            datalad.msg("Remove outliers from clusters"),
            tar.using(inputs = ["{input.data}"], outputs = ["{output}"]),
            wma_env.script,
            (
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
        touch(temp(shared_work/uid/"assigned-clusters-to-hemispheres.complete"))

    log: f"logs/assess_cluster_location_by_hemisphere/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/assess_cluster_location_by_hemisphere/{'.'.join(wildcards.values())}.tsv"

    envmodules:
        'python/3.7',
        "git-annex/8.20200810"

    group: "cluster_postprocess"
    resources:
        mem_mb=1500,
        runtime=10,

    shell:
        boost(
            datalad.msg("Assign cluster locations (left v right hem)"),
            tar.using(modify=["{input.data}"]),
            wma_env.script,
            (
                "wm_assess_cluster_location_by_hemisphere.py "
                "{input.data} -clusterLocationFile "
                "{input.atlas}/cluster_hemisphere_location.txt"
            )
        )

rule transform_clusters_to_subject_space:
    input:
        hemisphereAssignment=rules.assess_cluster_location_by_hemisphere.output,
        data=rules.remove_cluster_outliers.output,
        transform=rules.tractography_registration.output.transform_xfm,

    output:
        temp(directory(work/"transformed_clusters"/uid))

    log: f"logs/transform_clusters_to_subject_space/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/transform_clusters_to_subject_space/{'.'.join(wildcards.values())}.tsv"

    envmodules:
        'python/3.7',
        "git-annex/8.20200810"

    group: "cluster_postprocess"
    resources:
        mem_mb=3500,
        runtime=5,

    shell:
        boost(
            datalad.msg("Convert clusters back to subject T1w space"),
            xvfb_run,
            tar.using(inputs=["{input.data}"]),
            Pyscript(workflow.basedir, wma_env)(
                input=["data", "transform"],
                script="scripts/harden_transform.py"
            )
        )


rule separate_clusters_by_hemisphere:
    input:
        rules.transform_clusters_to_subject_space.output,

    output:
        bids_output_dwi(
            atlas="ORG",
            desc="sorted",
            suffix="clusters.tar.gz"
        )

    log: f"logs/separate_clusters_by_cluster/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/separate_clusters_by_cluster/{'.'.join(wildcards.values())}.tsv"

    envmodules:
        'python/3.7',
        "git-annex/8.20200810"

    group: "cluster_postprocess"
    resources:
        mem_mb=1500,
        runtime=15,

    shell:
        boost(
            datalad.msg("Seperate clusters into folders based on location"),
            tar.using(outputs=["{output}"]),
            wma_env.script,

            "wm_separate_clusters_by_hemisphere.py {input} {output}"
        )


rule assign_to_anatomical_tracts:
    input:
        data=rules.separate_clusters_by_hemisphere.output,
        atlas=config["atlases"]["cluster_atlas"],

    output:
        bids_output_dwi(
            atlas="ORG",
            desc="tracts",
            suffix="clusters.tar.gz"
        )

    log: f"logs/assign_to_anatomical_tracts/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/assign_to_anatomical_tracts/{'.'.join(wildcards.values())}.tsv"

    envmodules:
        'python/3.7',
        "git-annex/8.20200810"

    group: "cluster_postprocess"
    resources:
        mem_mb=2000,
        runtime=6,

    shell:
        boost(
            datalad.msg("Assign clusters to one of 73 anatomical tracts"),
            tar.using(inputs=["{input.data}"], outputs=["{output}"]),
            wma_env.script,
            (
                "wm_append_clusters_to_anatomical_tracts.py "
                "{input.data} {input.atlas} {output}"
            )
        )
