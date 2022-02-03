from pathlib import Path
from typing import cast
import pandas as pd
from snakeboost import snakemake_args


if __name__ == "__main__":
    args = snakemake_args()
    assert isinstance(args.input, list)
    assert isinstance(args.output, list)

    paths = args.input
    output = Path(args.output[0])

    dataframes = (
        cast(
            pd.DataFrame,
            pd.read_csv(path, index_col=(0, 1))
        ) for path in paths
    )

    merged = pd.concat(dataframes).drop(columns=[0])
    with output.open('w') as f:
        merged.to_csv(f)

