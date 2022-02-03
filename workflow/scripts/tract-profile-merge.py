from pathlib import Path
from typing import cast
import pandas as pd
from snakeboost import snakemake_args




if __name__ == "__main__":
    args = snakemake_args()
    assert isinstance(args.input, list)
    assert isinstance(args.output, list)
    assert isinstance(args.wildcards, dict)

    paths = args.input
    output = Path(args.output[0])

    dataframes = (
        cast(
            pd.DataFrame,
            pd.read_csv(path, iterator=False)
        ).assign(subject=args.wildcards["subject"]) for path in paths
    )

    merged = pd.concat(dataframes)
    merged.set_index(["subject", "cluster"])
    with output.open('w') as f:
        merged.to_csv(f)

