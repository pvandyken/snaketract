from glob import glob
import itertools as it
import os
from pathlib import Path
import shutil
import subprocess as sp
import tempfile
from typing import Callable, Dict, Generator, List, Optional, Tuple, Union, cast

from dipy.io.stateful_tractogram import Origin, Space
from dipy.io.stateful_tractogram import StatefulTractogram
from dipy.io.streamline import load_tractogram, save_tractogram
from dipy.io.vtk import load_vtk_streamlines
import fury.io as fio
import more_itertools as itx
import networkx as nx
from snakeboost import snakemake_args


def dipy_convert(src: Path, out: Path):
    def inner(ref: Path):
        save_tractogram(load_tractogram(str(src), str(ref)), str(out))

    return inner


def fury_convert(src: Path, out: Path):
    if out.suffix == ".vtk":
        binary = True
    else:
        binary = False
    fio.save_polydata(fio.load_polydata(str(src)), str(out), binary=binary)


Converter = Union[str, Callable[[Path, Path], Optional[Callable[..., None]]]]
ConversionMapping = Dict[Converter, List[str]]


conversion_map: ConversionMapping = {
    "tckconvert {input} {output}": [".tck", ".vtk", ".bundles"],
    dipy_convert: [".tck", ".trk"],
    fury_convert: [".vtp", ".vtk"],
}


def get_conversion_graph(mapping: ConversionMapping):
    extensions = set(it.chain.from_iterable(mapping.values()))
    G = nx.Graph()
    G.add_nodes_from(extensions)
    for ex1, ex2 in it.combinations(extensions, 2):
        for converter, valid_exts in mapping.items():
            if ex1 in valid_exts and ex2 in valid_exts:
                G.add_edge(ex1, ex2, converter=converter)
    return G


def get_converters(
    G: nx.Graph, source: str, target: str
) -> Generator[Tuple[Converter, str], None, None]:
    path = cast(List[str], nx.shortest_path(G, source, target))
    for ex1, ex2 in itx.windowed(path, 2):
        yield (G[ex1][ex2]["converter"], ex2)  # type: ignore


def _root_stem(path: Path):
    return os.path.splitext(path)[0]


def load_vtp(
    filename,
    reference,
    to_space=Space.RASMM,
    to_origin=Origin.NIFTI,
    bbox_valid_check=True,
):
    """Load the stateful tractogram from vtp

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
    if extension != ".vtp":
        raise Exception("File must be in .vtp format")

    data_per_point = None
    data_per_streamline = None

    streamlines = load_vtk_streamlines(filename)

    sft = StatefulTractogram(
        streamlines,
        reference,
        Space.RASMM,
        origin=Origin.NIFTI,
        data_per_point=data_per_point,
        data_per_streamline=data_per_streamline,
    )

    sft.to_space(to_space)
    sft.to_origin(to_origin)

    if bbox_valid_check and not sft.is_bbox_in_vox_valid():
        raise ValueError(
            "Bounding box is not valid in voxel space, cannot "
            "load a valid file if some coordinates are invalid.\n"
            "Please set bbox_valid_check to False and then use "
            "the function remove_invalid_streamlines to discard "
            "invalid streamlines."
        )

    return sft


def glob_inputs_outputs(in_path: Path, out_path: Path):
    ast_loc = _root_stem(in_path).find("*")

    inputs = [Path(path) for path in glob(str(in_path), recursive=True)]
    var_parts = [_root_stem(path)[ast_loc:] for path in inputs]
    out_const_part = _root_stem(out_path)[0 : _root_stem(out_path).find("*")]
    outputs = [
        Path(out_const_part + var_part).with_suffix(out_path.suffix)
        for var_part in var_parts
    ]
    return inputs, outputs


def convert_file(src: Path, dest: Path, reference: Path = None):
    if src.suffix == dest.suffix:
        return 0
    conversion_graph = get_conversion_graph(conversion_map)
    converters = get_converters(conversion_graph, src.suffix, dest.suffix)
    in_path = src
    in_path_is_temp = False
    for converter, ext in converters:
        if ext == dest.suffix:
            out_path = dest
        else:
            out_path = (
                Path(tempfile.mkdtemp(prefix="convert_tracts."))
                / dest.with_suffix(ext).name
            )
        if callable(converter):
            partial = converter(in_path, out_path)
            if partial:
                if not reference:
                    raise TypeError(
                        f"--ref must be provided when converting from {src.suffix} to "
                        f"{dest.suffix}"
                    )
                if not reference.exists():
                    raise TypeError(f"No file called {reference}")
                partial(reference)
        else:
            cmd = converter.format(input=f"'{in_path}'", output=f"'{out_path}'")
            sp.run(cmd, shell=True)
        if in_path_is_temp:
            shutil.rmtree(in_path.parent)
        in_path = out_path
        in_path_is_temp = True

    return


def main():
    args = snakemake_args(input=["in", "--ref"], output=["out"])
    if isinstance(args.input, list):
        data = args.input[0]
        if len(args.input) == 2:
            ref = args.input[1]
        else:
            ref = None
    else:
        data = args.input.get("input", Path(""))
        ref = None

    if isinstance(args.output, list):
        output = args.output[0]
    else:
        output = args.output.get("output", Path(""))

    input_type = data.suffix
    allowed_types = {*it.chain.from_iterable(conversion_map.values())}
    if input_type not in allowed_types:
        raise TypeError(f"Invalid file type as input. Must be one of {allowed_types}")

    output_type = output.suffix
    if output_type not in allowed_types:
        raise TypeError(f"Invalid file type as output. Must be one of {allowed_types}")

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
            convert_file(path, output, ref)
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
        out_path = output / out_name
        convert_file(path, out_path, ref)


if __name__ == "__main__":
    main()
