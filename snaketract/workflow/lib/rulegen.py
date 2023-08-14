
import os
import itertools as it
from workflow.lib.path_store import PathStore

store = PathStore()

def genrule(
    __rulename,
    wildcards=None,
    variable_inputs=None,
    fixed_inputs=None,
    output=None,
    output_root="",
    shell=None,
):
    inputexpr = {
        **{
            key: store(key) for key in variable_inputs
        },
        **fixed_inputs
    }
    outttempl = os.path.join(*it.chain(wildcards, variable_inputs))
    outputexpr = {
        
    }
    f"""
    rule {rulename}:
        input: **{inputexpr}
        output: {outputs}
            tempout(
                "warp_between_t1_and_mni_via_ants",
                inputs['preproc_dwi'],
                ".nii.gz",
                "group",
                "dir",
                "interpolation",
                "moving",
                "fixed",
            )
        log: {log}
            log(
                "warp_between_t1_and_mni_via_ants",
                inputs['preproc_dwi'],
                "moving",
                "fixed",
                "dir",
                "interpolation",
                "group",
            )
        benchmark: {benchmark}
            benchmark(
                "warp_between_t1_and_mni_via_ants",
                inputs['preproc_dwi'],
                "moving",
                "fixed",
                "dir",
                "interpolation",
                "group",
            )
        resources: {resources}
            runtime=2,
        envmodules: {envmodules}
            "StdEnv/2020",
            "gcc/9.3.0",
            "ants/2.3.5",
        group: {group}
        group: lambda wcards: wcards['group']
        shell: {shell}
        shell:
            """
            compose_transforms () {{
                while [[ -n "${{1:-}}" ]]; do
                    printf -- '-t %s ' $1
                    shift
                done
            }}
            antsApplyTransforms \
                -d 3 --interpolation {wildcards.interpolation} \
                -i {input.moving} -o {output} -r {input.fixed} \
                $(compose_transforms {input.txf})
            """
    """
