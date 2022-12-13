from __future__ import annotations

import itertools as it
from pathlib import Path

import more_itertools as itx
import numpy as np
import pandas as pd
import sklearn.preprocessing as preproc
import sklearn.impute as impute
import sklearn.compose as compose
import sklearn.pipeline as pipeline
from snakeboost import snakemake_args



def filter_logile(bin: int, num_bins: int = 10):
    def inner(matrix):
        if bin >= num_bins:
            raise ValueError("bin must be less then num_bins")
        masked = np.ma.masked_equal(matrix, 0)
        log = np.ma.log10(masked)
        threshold = 10 ** (((log.max() - log.min()) * bin / num_bins) + log.min())
        return matrix < threshold

    return inner


def main():
    args = snakemake_args(
        input=["participants"],
        output=["out"],
    )
    if isinstance(args.input, dict):
        raise TypeError("Inputs must be specified as a single item")
    if not isinstance(args.output, dict):
        raise TypeError("Outputs must be specified as a dict")

    input = itx.one(args.input)
    partic = pd.read_csv(input, sep='\t')

    groups = ["HC", "FEP"]
    for group in groups:
        partic[group] = (partic['group'] == group).map(int)

    sex_encode = pipeline.make_pipeline(
        preproc.OrdinalEncoder(),
        impute.SimpleImputer(strategy="most_frequent")
    )
    column_preproc = compose.make_column_transformer(
        ("passthrough", groups),
        (sex_encode, ["sex"]),
        (impute.SimpleImputer(), ['age']),
    )

    X = column_preproc.fit_transform(partic)
    np.savetxt(args.output['mat'], X, fmt='%i')
    con = np.array([
        [1, -1, 0, 0],
        [-1, 1, 0, 0]
    ])
    np.savetxt(args.output['con'], con, fmt='%i')
    


if __name__ == "__main__":
    main()

