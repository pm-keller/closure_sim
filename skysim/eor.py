"""

eor.py

Authors: Pascal M. Keller
Affiliation: Cavendish Astrophysics
Created on: 04/02/2022

Description:
    Tools for simulating the response to the HI signal.
    Use with astropy units!

"""

import os
import sys
import h5py
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.cosmology import Planck15

sys.path.append(os.getcwd())
from util import c, kb, f_HI, ztof, ftoz, extract_from_str


class EORModel:
    def __init__(
        self,
        eor=None,
        dimensions=None,
        redshift_start=None,
        units="mK",
        cosmology=Planck15,
    ):
        self.eor = eor
        self.z0 = redshift_start
        self._set_shape(np.shape(eor))
        self._set_dimensions(dimensions)
        self.cosmo = cosmology
        self.units = units

    def _set_shape(self, shape):

        if np.isscalar(shape):
            self.shape = 3 * (shape,)
        else:
            self.shape = shape

    def _set_dimensions(self, dimensions):
        if not isinstance(dimensions, list):
            self.dimensions = 3 * [dimensions]
        else:
            self.dimensions = dimensions

    def read(
        self, path, shape=None, dimensions=None, redshift_start=None, dtype=np.float32
    ):
        """
        Read EoR model from file
        """

        fext = os.path.splitext(path)[-1]

        try:
            if fext == ".h5":
                with h5py.File(path, "r") as f:
                    self.eor = f["Data"][()]
                    self.dimensions = list(f["Dimensions"][()] * u.Gpc)
                    self.z0 = f["Redshift Start"][()]
                    self.shape = self.eor.shape
            else:
                self.shape = shape
                self._set_shape(shape)
                self._set_dimensions(dimensions)
                self.z0 = redshift_start
                self.eor = np.fromfile(path, dtype=dtype).reshape(self.shape)
        except:
            raise Exception("Read file failed! Try different path, shape or dtype")

    def write(self, path, fname="-", **kwargs):
        """
        Write EoR model to file
        """
        assert isinstance(self.eor, np.ndarray), "EoR cuboid not specified!"

        f = h5py.File(path, "a")
        f.create_dataset("Data", data=self.eor)
        f.create_dataset(
            "Dimensions", data=[d.to(u.Gpc).value for d in self.dimensions]
        )
        f.create_dataset("Redshift Start", data=self.z0)
        f.close()

    def convert_units(self, units="Jy"):
        """
        Convert between mK and Jy
        """

        if units == self.units:
            print(f"EoR is already in units of {units}!")
        else:
            freq = self.get_freq_array()
            lam = c / freq

            # steradian per pixel
            rad1 = self.get_angular_resolution(axis=0)
            rad2 = self.get_angular_resolution(axis=1)
            sterad = rad1 * rad2

            # conversion factor
            Jy_to_mK = (lam ** 2 / (2 * kb)).to(u.mK / u.Jy) / sterad

            if units == "Jy":
                self.eor = self.eor / Jy_to_mK.value
                self.units = units
            elif units == "mK":
                self.eor = self.eor * Jy_to_mK.value
                self.units = units
            else:
                print("Units must be one of ('mK', 'Jy')")

    def get_z_array(self):
        """
        Get redshifts
        """
        assert (
            np.isscalar(self.z0) and self.dimensions != None
        ), "Start redshift or dimensions not specified!"

        nz = self.shape[-1]
        z_array = np.zeros(nz)
        z_array[0] = self.z0

        # comoving size of cell
        dR = self.dimensions[-1] / nz

        for i in range(1, nz):
            # redshift of previous cell
            z_prev = z_array[i - 1]

            # comoving distance per redshift interval
            dr_dz = c / self.cosmo.H(z_prev)

            # redshift of current cell
            z_array[i] = z_prev + dR / dr_dz

        return z_array

    def get_freq_array(self, z_array=None):
        """
        Compute frequencies corresponding to redshift axis
        """
        if z_array == None:
            z_array = self.get_z_array()
        return ztof(z_array, f_HI)

    def select_freq_range(self, fmin, fmax):
        """
        Select a frequency range to keep
        """
        freq = self.get_freq_array()
        idx = np.where((freq > fmin) & (freq < fmax))[0]

        freq = freq[idx]
        self.eor = self.eor[:, :, idx]
        self.shape = self.eor.shape

        z = ftoz(freq[[0, -1]])
        Dc_start, Dc_stop = self.cosmo.comoving_distance(z)
        self.dimensions[-1] = Dc_stop - Dc_start
        self.z0 = z[0].value

    def smooth(self, kernel, **kwargs):
        """
        Smooth EoR with the given kernel
        """
        self.eor = ndimage.convolve(self.eor, kernel, **kwargs)

    def downsample(self, n=(1, 1, 1)):
        """
        Downsample EoR cuboid
        """
        assert isinstance(self.eor, np.ndarray), "EoR model not specified!"

        if np.isscalar(n):
            n = 3 * (n,)

        self.eor = self.eor[:: n[0], :: n[1], :: n[2]]
        self._set_shape(self.eor.shape)

    def get_comoving_distance(self):
        """
        Get comoving distances corresponding to the redshift axis
        """
        return self.cosmo.comoving_distance(self.get_z_array())

    def get_comoving_coordinates(self):
        """
        Get comoving coordinate axes of EoR cuboid
        """
        assert isinstance(
            self.dimensions, list
        ), "Dimensions not specified! Must be list."

        Lx, Ly, Lz = self.dimensions
        Dx = np.linspace(-Lx, Lx, self.shape[0]) / 2
        Dy = np.linspace(-Ly, Ly, self.shape[1]) / 2
        Dz = self.get_comoving_distance()

        return Dx, Dy, Dz

    def get_comoving_mesh(self):
        """
        Get comoving coordinate meshgrid of EoR cuboid
        """
        Dx, Dy, Dz = self.get_comoving_coordinates()
        return np.meshgrid(Dx, Dy, Dz)

    def get_angular_resolution(self, axis=0):
        """
        Get angular resoluion of EoR cuboid
        Note that this is redshift dependent, so an array is returned.
        """
        Dz = self.get_comoving_distance()
        dx = self.dimensions[axis] / self.shape[axis]
        return np.arcsin(dx / np.sqrt(Dz ** 2 + dx ** 2))

    def get_dircos_array(self, axis=0):
        """
        Get array with direction cosines of a given axis of the EoR cuboid
        """
        dl = np.sin(self.get_angular_resolution(axis))
        l_array = np.arange(self.shape[axis]) * np.mean(dl)

        return l_array - np.mean(l_array)

    def get_angular_radius(self, axis=0):
        """
        Get the angular radius that the EoR subtends on the sky along a given axis
        """
        l_array = self.get_dircos_array(axis)
        return np.max(np.arcsin(np.abs(l_array)))

    def get_sky_coodinates(self):
        """
        Get sky coorinates axes: l, m and f,
        where l and m are the direction cosines and f is the frequency.
        The flat sky/narrow beam approximations are applied!
        """

        f_array = self.get_freq_array()
        l_array = self.get_dircos_array(axis=0)
        m_array = self.get_dircos_array(axis=1)

        return l_array, m_array, f_array

    def _tile(self, n=1, axis=0):
        """
        Tile the EoR cuboid along a given axis
        """
        self.eor = np.concatenate(n * (self.eor,), axis=axis)
        self.shape = self.eor.shape
        self.dimensions[axis] = n * self.dimensions[axis]

    def tile(self, n=(1, 1, 1)):
        """
        Tile the EoR cuboid along various axes
        """
        if np.isscalar(n):
            n = 3 * (n,)
        if len(n) == 2:
            n = n + (1,)

        for i in range(3):
            self._tile(n[i], axis=i)

    def plot(self, ax=None, z_index=0, path=None, **kwargs):
        """
        Plot a slice of the EoR cuboid
        """
        Dx, Dy, Dz = self.dimensions
        z = self.get_z_array()[z_index]
        freq = self.get_freq_array(z)

        if ax == None:
            fig, ax = plt.subplots()

        im = ax.imshow(
            self.eor[:, :, z_index],
            extent=(0, Dx.value, 0, Dy.value),
            cmap="afmhot",
            **kwargs,
        )

        ax.set_xlabel(f"X ({str(Dx.unit)})")
        ax.set_ylabel(f"Y ({str(Dy.unit)})")
        ax.set_title("Redshift: {:.2f}, Frequency: {:.2f}".format(z, freq.to(u.MHz)))

        if self.units == "mK":
            plt.colorbar(im, label=r"$\delta T_\mathrm{B}$ (mK)")
        else:
            plt.colorbar(im, label=r"Flux Density (Jy)")

        # save plot
        if isinstance(path, str):
            plt.savefig(path)

        return ax
