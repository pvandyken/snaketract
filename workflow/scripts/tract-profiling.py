from pathlib import Path
import numpy as np
from dipy.io.streamline import load_tractogram
import nibabel as nib
import matplotlib.pyplot as plt
import pandas as pd
from dipy.segment.clustering import QuickBundles
from dipy.segment.metric import AveragePointwiseEuclideanMetric, ResampleFeature
import dipy.stats.analysis as dsa
from snakeboost import snakemake_args

import dipy.tracking.streamline as dts

class TractProfile:
    def __init__(self, streamlines, ref):
        cluster = load_tractogram(str(streamlines), str(ref))
        if not cluster:
            raise Exception("Cluster could not be loaded")
        feature = ResampleFeature(nb_points=100)
        metric = AveragePointwiseEuclideanMetric(feature)

        qb = QuickBundles(np.inf, metric=metric)
        cluster_bundle = qb.cluster(cluster.streamlines)

        self.streamlines = dts.orient_by_streamline(
            cluster.streamlines,
            cluster_bundle.centroids[0]
        )

        self.weights = dsa.gaussian_weights(self.streamlines)


    def get_profile(self, img):
        return dsa.afq_profile(
            img.get_fdata(),
            self.streamlines,
            img.affine,
            weights=self.weights
        )


def get_profiles(data, param_maps, ref):
    results = np.empty((len(param_maps), len(data), 100))

    for i, streamlines in enumerate(data):
        profile = TractProfile(streamlines, ref)
        for j, param_map in enumerate(param_maps.values()):
            results[j][i] = profile.get_profile(param_map)
    return results.mean(axis=2)


if __name__ == "__main__":
    args = snakemake_args()
    assert isinstance(args.input, dict)
    assert isinstance(args.output, list)

    data = args.input["data"]
    ref_img = args.input["ref_img"]
    output = Path(args.output[0])

    parameter_maps = {
        "FA": nib.load(args.input["FA"]),
        "R1": nib.load(args.input["R1"])
    }

    if data.is_dir():
        paths = [*data.glob("*.tck")]
    else:
        raise FileNotFoundError("Input must be a directory")

    profiles = get_profiles(paths, parameter_maps, ref_img)
    profile_table = pd.DataFrame(
        {
            key: data for key, data in zip(parameter_maps, profiles)
        }
    )
    with output.open('w') as f:
        profile_table.to_csv(f)
