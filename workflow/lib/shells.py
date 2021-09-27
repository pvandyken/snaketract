import toolz as tz
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
        return bids(root=self.root,
                datatype='dwi',
                desc=algorithm,
                suffix=f"{self.tissue}fod.mif",
                **wildcards)

def is_multi_shelled(bval_file: str):
    with open(bval_file) as bval_s:
        bvals = filter(None, tz.first(bval_s).split(" "))
        ranks = tz.groupby(_shell_rank, bvals)
        return len(ranks) > 2

        
def _shell_rank(bval):
    return round(int(bval) / 500) * 500
    