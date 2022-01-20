from itertools import groupby
from tokenize import group
from snakebids import bids

class FodAlgorithm:
    def __init__(self, tissue, bval, root):
        self.tissue = tissue
        self.bval = bval
        self.root = root

    def __call__(self, wildcards):
        sub = wildcards['subject']
        bval = self.bval.format(subject=sub)
        if is_multi_shelled(bval):
            algorithm = 'ms3t'
        else:
            algorithm = 'ss3t'
        return bids(
            root=self.root,
            datatype='dwi',
            model="CSD",
            space="T1w",
            desc=algorithm,
            label=self.tissue,
            suffix="fod.mif",
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
    return round(int(bval) / 500) * 500
