import toolz as tz

def is_multi_shelled(bval_file: str):
    with open(bval_file) as bval_s:
        bvals = filter(None, tz.first(bval_s).split(" "))
        ranks = tz.groupby(_shell_rank, bvals)
        print(ranks)
        return len(ranks) > 2

        
def _shell_rank(bval):
    return round(int(bval) / 500) * 500