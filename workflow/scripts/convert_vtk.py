from pathlib import Path

import whitematteranalysis as wma
from snakeboost import snakemake_args

if __name__ == "__main__":
    args = snakemake_args()
    assert isinstance(args.input, list)
    assert isinstance(args.output, list)

    data = args.input[0]
    output = Path(args.output[0])

    if data.is_dir():
        paths = [*data.glob("*.vt{p,k}")]
        out_paths = [
            path.with_stem(".vtk") if path.stem == ".vtp" else path.with_stem(".vtp")
            for path in paths
        ]

        if not output.exists():
            output.mkdir()
    else:
        paths = [data]
        out_paths = [output]

    for path, out_path in zip(paths, out_paths):
        wma.io.write_polydata(wma.io.read_polydata(path), out_path)