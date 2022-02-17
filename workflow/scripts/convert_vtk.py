from glob import glob
from pathlib import Path
import re

import whitematteranalysis as wma
from snakeboost import snakemake_args


def glob_inputs_outputs(in_path: Path, out_path: Path):
    ast_loc = in_path.stem.find('*')
    inputs = [Path(path) for path in glob(str(in_path), recursive=True)]
    var_parts = [
        path.stem[ast_loc:] for path in inputs
    ]
    out_const_part = out_path.stem[0:out_path.stem.find("*")]
    outputs = [
        Path(out_const_part + var_part).with_suffix(out_path.suffix)
        for var_part in var_parts
    ]
    return inputs, outputs


def main():
    args = snakemake_args(
        input=["in"], output=["out"]
    )
    if isinstance(args.input, list):
        data = args.input[0]
    else:
        data = args.input.get("input", Path(""))

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
                "extension must be an asterisk (*): e.g. output/*.vtp"
            )

        paths, outputs = glob_inputs_outputs(data, output)

        for path, output in zip(paths, outputs):
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
        wma.io.write_polydata(wma.io.read_polydata(str(path)), str(out_path))

if __name__ == "__main__":
    main()