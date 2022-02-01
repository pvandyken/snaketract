from pathlib import Path
import numpy as np
from dipy.viz import window, actor, colormap as cmap
from dipy.io.streamline import load_tractogram
import nibabel as nib
import matplotlib.pyplot as plt
from dipy.segment.clustering import QuickBundles
from dipy.segment.metric import AveragePointwiseEuclideanMetric, ResampleFeature
import dipy.stats.analysis as dsa

import dipy.tracking.streamline as dts


rand_cluster = next(clusters.iterdir())

class TractProfile:
    def __init__(self, cluster, ref)
        cluster = load_tractogram(str(streamlines), str(ref))
        feature = ResampleFeature(nb_points=100)
        metric = AveragePointwiseEuclideanMetric(feature)

        qb = QuickBundles(np.inf, metric=metric)
        cluster_bundle = qb.cluster(cluster.streamlines)

        self.streamlines = dts.orient_by_streamline(
            cluster.streamlines,
            cluster_bundle.centroids[0]
        )

        self.weights = dsa.gaussian_weights(self.streamlines)


    def get_profile(self, img_path):
        img = nib.load(img_path)

        return dsa.afq_profile(
            img.get_fdata(),
            self.streamlines,
            img.affine,
            weights=self.weights
        )




