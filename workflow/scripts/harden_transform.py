from typing import TYPE_CHECKING, Union, Any
from pathlib import Path
from glob import glob

import whitematteranalysis as wma
from snakeboost import snakemake_args
import vtk
if TYPE_CHECKING:
    vtk: Any

def open_transform(filename: Union[Path, str], inverse: bool = False):
    filename = Path(filename)
    if filename.suffix == '.xfm':
        reader = vtk.vtkMNITransformReader()
        reader.SetFileName(str(filename))
        reader.Update()
        transform = reader.GetTransform()
    elif filename.suffix == '.img':
        reader = vtk.vtkImageReader()
        reader.SetFileName(str(filename))
        reader.Update()
        coeffs = reader.GetOutput()
        transform = vtk.vtkBSplineTransform()
        transform.SetCoefficients(coeffs)
        print(coeffs)
        print(transform)
    else:
        f = open(filename, 'r')
        transform = vtk.vtkTransform()
        matrix = vtk.vtkMatrix4x4()
        for i in range(0,4):
            for j in range(0,4):
                matrix_val = float(f.readline())
                matrix.SetElement(i,j, matrix_val)
        transform.SetMatrix(matrix)
        del matrix
    if inverse:
        transform.Inverse()
    return transform


if __name__ == "__main__":
    args = snakemake_args()
    assert isinstance(args.input, dict) and len(args.input) == 2, (
        "Incorrect number of inputs provided"
    )

    data = args.input["data"]

    transform_path = args.input["transform"]
    inverse = True
    assert isinstance(args.output, list) and len(args.output) == 1, (
        "Incorrect number of outputs provided"
    )
    output = Path(args.output[0])
    transform = open_transform(transform_path, inverse=inverse)

    if data.is_dir():
        paths = data.glob("*.vtp")

        if not output.exists():
            output.mkdir()


        for path in paths:
            wma.io.transform_polydata_from_disk_using_transform_object(
                str(path), transform, str(output/path.name)
            )
    else:
        wma.io.transform_polydata_from_disk_using_transform_object(
            str(data), transform, str(output)
        )