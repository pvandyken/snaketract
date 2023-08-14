from itertools import groupby
from tokenize import group
from snakebids import bids

class FodAlgorithm:
    def __init__(self, tissue, bval, root, single_shell_alg="ss3t"):
        self.tissue = tissue
        self.bval = bval
        self.root = root
        self.single_shell_alg = single_shell_alg

    def __call__(self, wildcards):
        bval = self.bval.format(**wildcards)
        if is_multi_shelled(bval):
            algorithm = 'ms3t'
        else:
            algorithm = self.single_shell_alg
        return bids(
            root=self.root,
            datatype='dwi',
            model="CSD",
            space="T1w",
            desc=algorithm,
            label=self.tissue,
            suffix="fod.nii.gz",
            **wildcards
        )

def is_multi_shelled(bval_file: str):
    with open(bval_file) as bval_s:
        bvals = filter(None, next(bval_s).split())
        ranks = set()
        ranks = { rank for rank, _ in groupby(bvals, key=_shell_rank) }
        for rank, _ in groupby(bvals, key=_shell_rank):
            ranks.add(rank)

        return len(ranks) > 2


def _shell_rank(bval):
    return round(int(float(bval)) / 500) * 500
