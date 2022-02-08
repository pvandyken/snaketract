from pathlib import Path
import numpy as np
import nibabel as nib
from snakeboost import snakemake_args


if __name__ == "__main__":
    args = snakemake_args()
    assert isinstance(args.input, dict)
    assert isinstance(args.output, list)
    assert isinstance(args.wildcards, dict)

    data = args.input["data"]
    mask = args.input["mask"]
    output = Path(args.output[0])

    img = nib.load(data)
    mask_data = nib.load(t1_mask_path).get_fdata().astype(bool)
    masked = np.ma.masked_where(~mask_data, data)
    r1_data = (1/masked).data
    r1_data[~mask_data] = 0

    r1 = nib.Nifti1Image(r1, t1map.affine)
    nib.save(r1, output)
