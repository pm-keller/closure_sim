"""

beam.py

Authors: Pascal M. Keller
Affiliation: Cavendish Astrophysics
Created on: 04/02/2022

Description:
    Tools for adding a beam response to data models

"""

import os
import sys
import h5py
import copy
import numpy as np
from pyuvdata import UVBeam
import astropy.units as u

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

sys.path.append(os.getcwd())
from util import dircos_to_altaz


class BeamModel:
    def __init__(self, uvbeam=None, beam_image=None):
        self.uvbeam = uvbeam
        self.beam_image = beam_image

        if isinstance(beam_image, np.ndarray):
            self.shape = beam_image.shape

    def copy(self):
        """
        Return a copy of self
        """
        return copy.deepcopy(self)

    def read_uvbeam(self, path, **kwargs):
        """
        Read a beam as UVBeam object
        """
        beam = UVBeam(**kwargs)
        beam.read_beamfits(path)
        beam.efield_to_power()
        beam.interpolation_function = "az_za_simple"
        beam.freq_interp_kind = "cubic"
        beam.peak_normalize()

        self.uvbeam = beam

    def read_image(self, path):
        """
        Read beam image data from file
        """
        with h5py.File(path, "r") as f:
            self.beam_image = f["Data"][()]
            self.half_width = f["Half-Width"][()]
            self.freq_array = f["Frequency Array"][()] * u.MHz
            self.shape = self.beam_image.shape

    def write_image(self, path, fname, **kwargs):
        """
        Write beam image data to file
        """
        assert isinstance(self.beam_image, np.ndarray), "Beam image not specified!"

        f = h5py.File(path, "a")
        f.create_dataset("Data", data=self.beam_image)
        f.create_dataset("Frequency Array", data=self.freq_array.to(u.MHz).value)
        f.create_dataset("Half-Width", data=self.half_width)
        f.close()

    def make_image(self, Naxis, half_width, freq, inplace=True):
        """
        Make a beam image at given frequencies
        """
        assert self.uvbeam != None, "Beam not specified!"

        l_array = np.linspace(-half_width, half_width, Naxis)
        freq_array = np.atleast_1d(freq.to(u.Hz).value)

        # get direction cosine meshgrid
        ll, mm = np.meshgrid(l_array, l_array)

        # convert direction cosines to alt/az coordinates
        alt, az = dircos_to_altaz(ll, mm)
        alt, az = alt.reshape(-1), az.reshape(-1)

        # zenith angle
        za = np.pi / 2 - alt

        # interpolate beam
        beam_im = self.uvbeam.interp(az, za, freq_array=freq_array, reuse_spline=True)
        beam_im = beam_im[0][0, 0, 0].reshape((-1, Naxis, Naxis)).real

        if inplace:
            BeamObj = self
        else:
            self.copy()

        BeamObj.beam_image = np.moveaxis(beam_im, 0, -1)
        BeamObj.half_width = np.max(np.abs(l_array))
        BeamObj.freq_array = np.atleast_1d(freq)
        BeamObj.shape = (Naxis, Naxis, len(freq))

        return BeamObj

    def plot_image(self, ax=None, freq_idx=0, path=None, **kwargs):
        """
        Plot a 2D projected beam image at a given frequency
        """
        assert isinstance(self.beam_image, np.ndarray), "Beam image not specified!"

        hw = self.half_width

        if ax == None:
            fig, ax = plt.subplots()

        im = ax.imshow(
            self.beam_image[:, :, freq_idx],
            extent=(-hw, hw, -hw, hw),
            cmap="rainbow",
            norm=LogNorm(),
            **kwargs,
        )

        freq = self.freq_array[freq_idx].to(u.MHz)
        ax.set_title("Frequency: {:.2f}".format(freq))
        ax.set_xlabel("l")
        ax.set_ylabel("m")
        plt.colorbar(im, label="Normalised Beam Power")

        # save plot
        if isinstance(path, str):
            plt.savefig(path)

        return ax