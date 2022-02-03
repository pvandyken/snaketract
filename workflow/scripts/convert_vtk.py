from pathlib import Path

import whitematteranalysis as wma
from snakeboost import snakemake_args

if __name__ == "__main__":
    args = snakemake_args()
    if isinstance(args.input, list):
        data = args.input[0]
    else:
        data = Path(args.input.get("input", ""))

    if isinstance(args.output, list):
        output = Path(args.output[0])
    else:
        output = Path(args.output.get("output", ""))

    if data.is_dir():
        paths = [*data.glob("*.vtp")]
        out_names = [
            path.with_suffix(".vtk").name
            if path.suffix == ".vtp"
            else path.with_suffix(".vtp")
            for path in paths
        ]

        if not output.exists():
            output.mkdir()
    else:
        paths = [data]
        out_names = [output.name]
        output = output.parent

    for path, out_name in zip(paths, out_names):
        out_path = output/out_name
        print(out_path)
        wma.io.write_polydata(wma.io.read_polydata(str(path)), str(out_path))