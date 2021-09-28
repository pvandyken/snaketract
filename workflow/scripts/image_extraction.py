import argparse
import sys
from typing import Literal, Tuple, Optional, Union, List, TYPE_CHECKING
import nibabel as nib
from dipy.io.image import load_nifti
from dipy.io.streamline import load_tractogram
from dipy.viz import window, actor, colormap as cmap
from numpy import ndarray
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

if TYPE_CHECKING:
    from snakemake.jobs import Job
    snakemake: Job

Direction = Union[
    Literal["axial"],
    Literal["coronal"],
    Literal["sagittal"]
]

def streamline_actor(tractogram):
    streamlines = tractogram.streamlines
    color = cmap.line_colors(streamlines)
    return actor.line(streamlines, color)

def movie_calc(rotation, length):
    fps = 26
    num_frames = length * fps
    rotation_per_frame = rotation / num_frames
    return num_frames, rotation_per_frame

px = 1/plt.rcParams['figure.dpi'] # type: ignore

class SliceView:
    def __init__(self,
            nifti: Union[Path, str],
            output_folder: Path,
            file_prefix: str,
            tractography: Optional[Union[Path, str]] = None,
            overlays: List[Union[Path, str]] = [],
        ):
        self.output_folder = output_folder
        self.file_prefix = file_prefix

        if self._validate_nifti(Path(nifti)):
            self.img = nib.load(str(nifti))
            self.data = self.img.get_fdata()
            self.mids = [mid(self.data, d) for d in [
                "sagittal",
                "coronal",
                "axial"
            ]]
            if not self._validate_rgb(self.data):
                raise Exception("rgb data must have 3 dimensions")

        else:
            raise Exception("Must provide a nifti file")

        if tractography:
            if self._validate_tck(Path(tractography)):
                if x := load_tractogram(str(tractography), self.img):
                    self.tractography = x
                else:
                    raise Exception("Tractography could not be loaded")
            else:
                raise Exception("Tractography must be a .tck file")
        else:
            self.tractography = None

        if all([self._validate_nifti(Path(overlay)) for overlay in overlays]):
            self.overlays = [nib.load(str(overlay)) for overlay in overlays]
        else:
            raise Exception("Overlay must be a .nii[.gz] file")

        self._brightness = None

    def _validate_nifti(self, file: Path):
        nii = file.suffix == ".nii"
        niigz = file.suffixes[-1] == ".gz" and file.suffixes[-2] == ".nii"
        if nii or niigz:
            return True
        else:
            return False

    def _validate_rgb(self, data: ndarray):
        if len(data.shape) == 4 and not data.shape[3] == 3:
            return False
        else:
            return True

    def _validate_tck(self, file: Path):
        if file.suffix == ".tck":
            return True
        return False

    def _open_nifti(self, path: Path):
        return load_nifti(str(path), return_img=True)

    @property
    def brightness(self):
        return self._brightness

    @brightness.setter
    def brightness(self, values):
        self._brightness = values

    def get_brightness(self, data):
        if self.brightness:
            mean, std = data[data>0].mean(), data[data>0].std()
            return mean + self.brightness[0] * std, mean + self.brightness[1] * std
        else:
            diff = np.max(data) - np.min(data[data>0])
            return np.min(data[data>0]), np.max(data)

    def tract_actor(self):
        if self.tractography:
            streamlines = self.tractography.streamlines
            color = cmap.line_colors(streamlines)
            return actor.line(streamlines, color)
        else:
            raise Exception("No tractography was defined")

    def render_3d(self, view: Direction, size: Tuple[int, int]=(1000,1000)) -> ndarray:
        scene = window.Scene()
        a = actor.slicer(self.data, self.img.affine, self.brightness)
        scene.add(self.tract_actor())
        scene.reset_camera_tight()
        if view == "sagittal" or view == "s":
            a.display(mid(a, "sagittal"), None, None)
            scene.azimuth(-90)
            scene.roll(90)
        if view == "coronal" or view == "c":
            a.display(None, mid(a, "coronal"), None)
            scene.elevation(90)
            scene.roll(180)
        if view == "axial" or view == "a":
            a.display(None, None, mid(a, "axial"))
        scene.add(a)
        #n_frames, az_angle = movie_calc(90, 5)
        #return window.record(self.scene, out_path="mov/mov", n_frames=n_frames, size=size, az_ang=az_angle, path_numbering=True)
        return window.snapshot(scene, size=size)

    def render_2d(self, img, view: Direction, size: Tuple[int, int]=(1000,1000)):
        scene = window.Scene()
        a = self.get_section(img.get_fdata(), img.affine, view)
        print(self.data.shape)
        scene.add(a)
        for overlay, hue in zip(self.overlays, (0.1, 0.5)):
            overlay_arr = self.get_section(overlay.get_fdata(), overlay.affine, view, hue=hue)
            print(overlay_arr.shape)
            #o_data = np.ma.masked_equal(overlay_arr, 0)
            scene.add(overlay_arr)
        scene.reset_camera_tight()
        window.show(scene)
        return window.snapshot(scene, size=size)

    def get_section(self, data, affine, view: Direction, hue: Optional[float] = None):
        if hue:
            m = cmap.colormap_lookup_table(hue_range=(hue, hue), value_range=(0,1), scale_range=(np.min(data), np.max(data)))
        else:
            m = None
        a = actor.slicer(data, affine, self.get_brightness(data), lookup_colormap=m)
        if view == "sagittal":
            i = self.mids[0]
            a.display(i, None, None)
        if view == "coronal":
            i = self.mids[1]
            a.display(None, i, None)
        if view == "axial":
            i = self.mids[2]
            a.display(None, None, i)
        return a

    def _get_section(self, data, view: Direction):
        i = mid(data, view)
        if view == "sagittal":
            return data[i, :, :, ...].T
        if view == "coronal":
            return data[:, i, :, ...].T
        if view == "axial":
            return data[:, :, i, ...].T



    def save_fig(self, view: Direction, size: Tuple[int, int]):
        if self.tractography:
            img_arr = self.render_3d(view)
        else:
            img_arr = self._get_section(self.data, view)

        plt.subplots(figsize=(size[0]*px, size[1]*px))[1].set_axis_off()
        plt.imshow(img_arr, cmap='gray', origin='lower')

        if self.file_prefix:
            name = f"{view}.png"
        else:
            name = self.file_prefix + f"_{view}.png"
        p = self.output_folder / name
        plt.savefig(p, bbox_inches="tight")

    def show_fig(self, view: Direction, size: Tuple[int, int]):
        if self.tractography:
            img_arr = self.render_3d(view)
        else:
            img_arr = self._get_section(self.data, view)
        plt.subplots(figsize=(size[0]*px, size[1]*px))[1].set_axis_off()
            #ma = np.ma.masked_less(img, 0.01)
        plt.imshow(img_arr, cmap="hot", origin='lower')

        plt.show()

def mid(data, direction: Direction):
    i = [
        "sagittal",
        "coronal",
        "axial"
    ].index(direction)

    return data.shape[i] // 2


def main():

    if snakemake:
        dwi = snakemake.input["dwi"]
        tracts = snakemake.input["tracts"]
        output = snakemake.output
        prefix = ""

    else:
        parser = argparse.ArgumentParser(description="Extract png images from MRI")
        parser.add_argument('input', type=Path)
        parser.add_argument('outputFolder', type=Path)
        parser.add_argument('--tractography', type=Path)
        parser.add_argument('--file-prefix', type=str, help="String to prepend on file outputs")

        parsed = parser.parse_args(sys.argv[1:])
        dwi = parsed.input
        tracts = parsed.tractography
        output = parsed.outputFolder
        prefix = parsed.filePrefix

    view = SliceView(
        dwi,
        tractography=tracts,
        output_folder=output,
        file_prefix=prefix
    )

    for v in ["axial", "sagittal", "coronal"]:
        view.save_fig(v, size=(1000, 1000))

if __name__ == "__main__":
    main()