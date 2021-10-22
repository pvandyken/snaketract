rule qc_generate_tractography_images:
    input: 
        expand(
            rules.run_act.output,
            **inputs['input_lists']['preproc_dwi']
        )
    output: 
        output_qc / "tractography"
    log: f"logs/qc_generate_tractography_images/{'.'.join(wildcards.values())}.log"
    benchmark: f"benchmarks/qc_generate_tractography_images/{'.'.join(wildcards.values())}.tsv"
    conda:
    container:
    envmodules:
    group:
    threads: 1
    resources:
        mem_mb=1000,
        runtime=30,
    params:
    shell: (
        command
    )