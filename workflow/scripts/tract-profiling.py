from pathlib import Path
import re
import numpy as np
from dipy.io.streamline import load_tractogram
import nibabel as nib
import pandas as pd
from dipy.segment.clustering import QuickBundles
from dipy.segment.metric import AveragePointwiseEuclideanMetric, ResampleFeature
import dipy.stats.analysis as dsa
from snakeboost import snakemake_args

import dipy.tracking.streamline as dts

class TractProfile:
    def __init__(self, streamlines, ref):
        cluster = load_tractogram(str(streamlines), str(ref))
        if cluster is False:
            raise Exception(f"Cluster {streamlines} could not be loaded")

        if not cluster:
            self.streamlines = None
            return
        feature = ResampleFeature(nb_points=100)
        metric = AveragePointwiseEuclideanMetric(feature)

        qb = QuickBundles(np.inf, metric=metric)
        cluster_bundle = qb.cluster(cluster.streamlines)

        self.cluster_bundle = cluster_bundle
        self.streamlines = cluster.streamlines

        self.weights = dsa.gaussian_weights(self.streamlines)


    def get_profile(self, img):
        return dsa.afq_profile(
            img.get_fdata(),
            dts.orient_by_streamline(
                self.streamlines,
                self.cluster_bundle.centroids[0]
            ),
            img.affine,
            weights=self.weights
        )


def get_profiles_and_streamlines(paths, param_maps, ref):
    profiles = np.empty((len(param_maps), len(paths), 100))
    streamlines = np.empty(len(paths))

    for i, path in enumerate(paths):
        profile = TractProfile(path, ref)
        if profile.streamlines:
            streamlines[i] = len([*profile.streamlines])
            for j, param_map in enumerate(param_maps.values()):
                profiles[j, i] = profile.get_profile(param_map)
        else:
            streamlines[i] = 0
            profiles[:, i, :] = 0
    return profiles.mean(axis=2), streamlines



if __name__ == "__main__":
    args = snakemake_args()
    assert isinstance(args.input, dict)
    assert isinstance(args.output, list)
    assert isinstance(args.wildcards, dict)

    data = args.input["data"]
    ref_img = args.input["ref"]
    output = Path(args.output[0])

    parameter_maps = {
        "FA": nib.load(args.input["FA"]),
        "R1": nib.load(args.input["R1"])
    }

    if data.is_dir():
        paths = [*data.glob("*.tck")]
    else:
        raise FileNotFoundError("Input must be a directory")

    profiles, streamlines = get_profiles_and_streamlines(paths, parameter_maps, ref_img)
    cluster_numbers = (
        match[0] for match in (
            re.search(r'(?<=cluster_)\d{5}(?=\.tck$)', str(path)) for path in paths
        ) if match
    )
    profile_table = pd.DataFrame(
        {
            key: data for key, data in zip(parameter_maps, profiles)
        }
    ).assign(
        cluster=pd.Series(cluster_numbers),
        subject=args.wildcards["subject"],
        streamlines=pd.Series(streamlines)
    ).set_index(["subject", "cluster"])
    with output.open('w') as f:
        profile_table.to_csv(f)
