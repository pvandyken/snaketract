from glob import glob
from pathlib import Path
import re
import os

from dipy.io.streamline import load_tractogram, save_tractogram
from dipy.io.vtk import load_vtk_streamlines
from dipy.io.stateful_tractogram import Space, Origin
from dipy.io.stateful_tractogram import StatefulTractogram
import whitematteranalysis as wma
from snakeboost import snakemake_args


def _root_stem(path: Path):
    return os.path.splitext(path)[0]

def load_vtp(
    filename,
    reference,
    to_space=Space.RASMM,
    to_origin=Origin.NIFTI,
    bbox_valid_check=True
):
    """ Load the stateful tractogram from vtp

    Parameters
    ----------
    filename : string
        Filename with valid extension
    reference : Nifti or Trk filename, Nifti1Image or TrkFile, Nifti1Header or
        trk.header (dict), or 'same' if the input is a trk file.
        Reference that provides the spatial attribute.
        Typically a nifti-related object from the native diffusion used for
        streamlines generation

    Returns
    -------
    output : StatefulTractogram
        The tractogram to load (must have been saved properly)
    """
    _, extension = os.path.splitext(filename)
    if extension != '.vtp':
        raise Exception("File must be in .vtp format")


    data_per_point = None
    data_per_streamline = None

    streamlines = load_vtk_streamlines(filename)

    sft = StatefulTractogram(streamlines, reference, Space.RASMM,
                             origin=Origin.NIFTI,
                             data_per_point=data_per_point,
                             data_per_streamline=data_per_streamline)


    sft.to_space(to_space)
    sft.to_origin(to_origin)

    if bbox_valid_check and not sft.is_bbox_in_vox_valid():
        raise ValueError('Bounding box is not valid in voxel space, cannot '
                         'load a valid file if some coordinates are invalid.\n'
                         'Please set bbox_valid_check to False and then use '
                         'the function remove_invalid_streamlines to discard '
                         'invalid streamlines.')

    return sft

def glob_inputs_outputs(in_path: Path, out_path: Path):
    ast_loc = _root_stem(in_path).find('*')

    inputs = [Path(path) for path in glob(str(in_path), recursive=True)]
    var_parts = [
        _root_stem(path)[ast_loc:] for path in inputs
    ]
    out_const_part = _root_stem(out_path)[0:_root_stem(out_path).find("*")]
    outputs = [
        Path(out_const_part + var_part).with_suffix(out_path.suffix)
        for var_part in var_parts
    ]
    return inputs, outputs


def convert_file(in_path: Path, out_path: Path, reference: Path):
    if in_path.suffix == out_path.suffix:
        return 0
    if in_path.suffix == ".vtp":
        tracts = load_vtp(str(in_path), reference)
    else:
        tracts = load_tractogram(str(in_path), reference)
    save_tractogram(tracts, str(out_path))


def main():
    args = snakemake_args(
        input=["in", "--ref"], output=["out"]
    )
    if isinstance(args.input, list):
        data = args.input[0]
        ref = args.input[1]
    else:
        data = args.input.get("input", Path(""))
        ref = Path()

    if isinstance(args.output, list):
        output = args.output[0]
    else:
        output = args.output.get("output", Path(""))

    input_type = data.suffix
    if input_type not in [".vtp", ".vtk"]:
        raise TypeError("Invalid file type as input. Must be one of ['.vtp', 'vtk']")

    output_type = output.suffix
    if output_type not in [".vtp", ".vtk"]:
        raise TypeError("Invalid file type as output. Must be one of ['.vtp', 'vtk']")

    if input_type == output_type:
        raise TypeError(
            f"Input and output types must be different: (both {input_type})"
        )

    data_str = str(data)
    if "*" in data_str:
        if "*" != output.stem[-1]:
            raise TypeError(
                "If input is a glob pattern, the last character of output before the "
                "extension must be an asterisk (*): e.g. output/*.vtp\n"
                f"Got {output}"
            )

        paths, outputs = glob_inputs_outputs(data, output)

        print(f"Converting {len(paths)} files")

        for path, output in zip(paths, outputs):
            output.parent.mkdir(exist_ok=True)
            wma.io.write_polydata(wma.io.read_polydata(str(path)), str(output))
        return 0


    if data.is_dir():
        paths = [*data.glob("*.vtp")]
        out_names = [
            path.with_suffix(".vtk").name
            if path.suffix == ".vtp"
            else path.with_suffix(".vtp")
            for path in paths
        ]
        print(f"Converting {len(paths)} files to vtk")

        if not output.exists():
            output.mkdir()
    else:
        paths = [data]
        out_names = [output.name]
        output = output.parent


    for path, out_name in zip(paths, out_names):
        out_path = output/out_name
        convert_file(path, out_path, ref)

if __name__ == "__main__":
    main()
