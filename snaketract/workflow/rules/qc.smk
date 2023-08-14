from snakebids import bids

# rule create_qc_reference_image:
#     input:
#         dwi=inputs.input_path['preproc_dwi'],
#         bvec=inputs.input_path['bvec'],
#         bval=inputs.input_path['bval']
#     output:
#         bids(root=qc,
#             datatype='dwi',
#             suffix="dwi.mif",
#             **wildcards)
#     envmodules:
#         "mrtrix/3.0.1"
#     shell:
#         'mrconvert {input.dwi} {output} -fslgrad {input.bvec} {input.bval} -quiet'

rule generate_qc_FOD_image:
    input:
        wm=rules.normalize_fiber_orientation_densities.output.wm,
        gm=rules.normalize_fiber_orientation_densities.output.gm,
        csf=rules.normalize_fiber_orientation_densities.output.csf
    output:
        bids_output_dwi(
            datatype='qc',
            suffix='vf.nii.gz',
            **wildcards,
        )
    envmodules:
        "mrtrix/3.0.1"
    shell:
        'mrconvert -coord 3 0 {input.wm} - | mrcat {input.csf} {input.gm} - {output}'



# rule generate_tractography_view_script:
#     input:
#         dwi=rules.create_qc_reference_image.output,
#         tracts=rules.extract_tractography_subset.output
#     output:
#         bids(root=qc,
#             datatype='dwi',
#             suffix="viewTractography",
#             **wildcards)
#     shell:
#         'echo "#!/bin/bash\nmrview $(realpath {input.dwi}) -tractography.load $(realpath {input.tracts})" > {output} && chmod +x {output}'

# rule generate_5tt_qc_view_script:
#     input:
#         dwi=rules.create_qc_reference_image.output,
#         anat=rules.create_seed_boundary.output
#     output:
#         bids(root=qc,
#             datatype='dwi',
#             suffix="viewInterface",
#             **wildcards)
#     shell:
#         'echo "#!/bin/bash\nmrview $(realpath {input.dwi}) -overlay.load $(realpath {input.anat})" > {output} && chmod +x {output}'

rule generate_odf_qc_view_script:
    input:
        dwi=rules.generate_qc_FOD_image.output,
        anat=rules.normalize_fiber_orientation_densities.output
    output:
        bids(
            root=qc/"view_scripts",
            datatype='dwi',
            suffix="viewOdf.sh",
            **wildcards
        )
    shell:
        """
        dedent () {{
            xargs -L1 echo
        }}
        echo -e "
            #!/bin/bash
            mrview $(realpath {input.dwi}) -overlay.load $(realpath {input.anat})
        " | dedent > {output}
        chmod +x {output}
        """



rule odf_qc:
    input:
        vf=rules.generate_qc_FOD_image.output,
        odf=rules.normalize_fiber_orientation_densities.output.wm,
        mask=inputs["brainmask"].path,
        planes=script("mrview_planes.sh"),
    output:
        directory(
            tempout(
                "odf_qc",
                inputs["preproc_dwi"],
                "",
            )
        )
    envmodules:
        "StdEnv",
        "mrtrix/3.0.1",
    shadow: 'minimal'
    shell:
        boost(
            xvfb_run,
            """
            mkdir -p {output}
            mrcalc {input.vf} {input.mask} -mult masked.nii.gz
            mrgrid -mask {input.mask} masked.nii.gz crop cropped.nii.gz
            mrview cropped.nii.gz -odf.load_sh {input.odf} \\
                -noannotations -capture.folder {output} \\
                $({input.planes} cropped.nii.gz ::~7 ) \\
                -exit 
            """
        )

rule odf_qc_montage:
    input:
        images=rules.odf_qc.output
    output:
        bids_output_dwi(
            datatype='qc',
            suffix="odf.png",
            **wildcards,
        )
    envmodules:
        "nixpkgs",
        "imagemagick",
    shadow: 'minimal'
    shell:
        """
        montage -font DejaVu-Serif -tile 7x3 -background black -geometry +0+0 \\
            {input}/sagittal* {input}/coronal* {input}/axial* \\
            {output}
        """

rule extract_tractography_subset:
    input: rules.run_act.output
    output:
        bids_output_dwi(
            datatype='qc',
            rec="{rec}",
            desc="200k",
            suffix='tractography.tck',
            **wildcards,
        )
    envmodules:
        "mrtrix/3.0.1"
    shell:
        'tckedit {input} -number 200k {output}'

rule get_tractography_qc:
    input:
        background=inputs['t1'].path,
        tracts=rules.extract_tractography_subset.output,
        planes=script("mrview_planes.sh"),
        mask=inputs['t1_mask'].path,
    output:
        directory(
            tempout(
                "get_tractography_qc",
                inputs["preproc_dwi"],
                "",
                "rec",
            )
        )
    envmodules:
        "StdEnv",
        "mrtrix/3.0.1",
    shadow: 'minimal'
    shell:
        boost(
            xvfb_run,
            """
            mkdir -p {output}
            mrcalc {input.background} {input.mask} -mult masked.nii.gz
            mrgrid -mask {input.mask} masked.nii.gz crop cropped.nii.gz
            mrview cropped.nii.gz -tractography.load {input.tracts} \\
                -noannotations -capture.folder {output} \\
                $({input.planes} cropped.nii.gz ::~7 ) \\
                -exit 
            """
        )

rule get_tractography_qc_montage:
    input:
        images=rules.get_tractography_qc.output
    output:
        bids_output_dwi(
            datatype='qc',
            rec="{rec}",
            desc="200k",
            suffix="tractography.png",
            **wildcards
        )
    envmodules:
        "nixpkgs",
        "imagemagick",
    shadow: 'minimal'
    shell:
        """
        montage -font DejaVu-Serif -tile 7x3 -background black -geometry +0+0 \\
            {input}/sagittal* {input}/coronal* {input}/axial* \\
            {output}
        """

rule qc_scripts:
    input:
        odf_qc=inputs["preproc_dwi"].expand(
            rules.generate_odf_qc_view_script.output[0],
        ),
    output:
        os.path.join(qc, "script_manifest.json"),
    run:
        with open(output[0], "w") as f:
            json.dump(
                {
                    "odf": sorted(input["odf_qc"]),
                },
                f,
            )

rule qc:
    input:
        tract_qc=inputs["preproc_dwi"].expand(
            rules.get_tractography_qc_montage.output[0],
            rec=config["tractography_algorithm"],
        ),
        # odf_qc=inputs["preproc_dwi"].expand(
        #     rules.odf_qc_montage.output[0],
        # ),
    output:
        os.path.join(qc, "data.json"),
    run:
        with open(output[0], "w") as f:
            json.dump(
                {
                    "tracts": {
                        "title": "Tractography",
                        "images": sorted(input["tract_qc"]),
                    },
                    # "odf": {
                    #     "title": "ODFs",
                    #     "images": sorted(input["odf_qc"]),
                    # },
                },
                f,
            )


_qc_app = resource("qc-app.tar.gz")
_qc_scripts = resource("qc.tar.gz")


def _get_tar_contents(file):
    try:
        return [
            p
            for p in sp.check_output(["tar", "-tf", file]).decode().splitlines()
            if p[-1] != "/"
        ]
    except sp.CalledProcessError as err:
        raise Exception("Unable to find qc-app.tar.gz...") from err


rule unpack_qc_app:
    input:
        _qc_app
    output:
        _get_tar_contents(_qc_app),
    shell:
        "tar -xvzf {input}"

rule unpack_qc_scripts:
    input:
        _qc_scripts
    output:
        _get_tar_contents(_qc_scripts),
    shell:
        "tar -xvzf {input}"
