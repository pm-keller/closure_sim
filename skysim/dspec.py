"""

dspec.py

Authors: Pascal M. Keller
Affiliation: Cavendish Astrophysics
Created on: 14/02/2022

Description:
    Tools for working with delay transformed data
    Intended for use with astropy units

"""

import numpy as np
import copy
import sys
import os

from scipy.signal.windows import get_window

import astropy.units as u
from astropy.cosmology import Planck18

import matplotlib.pyplot as plt

sys.path.append(os.getcwd())

from ft import fft, xps
from util import c, kb, f_HI, ftoz, comoving_depth, NoneType


class DSpec:
    def __init__(
        self,
        data,
        noise=None,
        freq_array=None,
        window="blackmanharris",
        axis=-1,
        units="Jy",
        cosmology=Planck18,
    ):
        self.data = data
        self.noise = noise
        self.N_axis = data.shape[axis]
        self.axis = axis
        self.freq_array = freq_array
        self.freq = np.mean(freq_array)
        self.bandwidth = np.abs(freq_array[-1] - freq_array[0])
        self.fs = self.N_axis / self.bandwidth
        self.units = units
        self.cosmo = cosmology

        if isinstance(window, str):
            self.window = get_window(window, self.N_axis)
        else:
            self.window = window

        dspec = fft(data, window, inverse=True, shift=True, axis=axis)
        self.dspec = dspec * self.N_axis / self.fs

    def __mul__(self, other):
        """
        Returns cross-power spectrum object
        """
        return XPSpec(self.dspec, other.dspec, self.freq_array)

    def get_delays(self):
        """
        Get delay array
        """
        return (np.arange(self.N_axis) - self.N_axis / 2) * self.fs / self.N_axis

    def get_k_parallel(self):
        """
        Get parallel k-modes
        """
        z = ftoz(self.freq)
        delay = self.get_delays()

        return 2 * np.pi * delay * f_HI * self.cosmo.H(z) / c / (1 + z) ** 2

    def get_xps(self, axis=None):
        """
        Get cross power spectra
        """
        return XPSpec(self.dspec, self.dspec, self.freq_array, self.cosmo)

    def get_xps_noise(self, axis=None):
        """
        Get noise cross power spectra
        """
        assert isinstance(self.noise, np.ndarray), "Noise must be ndarray!"

        return XPSpec(self.noise, self.noise, self.freq_array, axis)

    def downsample(self, n):
        """
        Downsample delay spectrum
        """
        dspec = np.moveaxis(self.dspec, self.axis, 0)
        self.dspec = np.moveaxis(dspec[::n], 0, self.axis)
        self.N_axis = self.dspec.shape[self.axis]


class XPSpec(DSpec):
    def __init__(
        self,
        dspec1=None,
        dspec2=None,
        freq_array=None,
        cosmology=Planck18,
    ):
        self.freq_array = freq_array
        self.freq = np.mean(freq_array)
        self.bandwidth = np.abs(freq_array[-1] - freq_array[0]) / 2
        self.N_axis = np.shape(freq_array)[-1]
        self.fs = self.N_axis / self.bandwidth
        self.cosmo = cosmology

        if not isinstance(dspec1, NoneType) and not isinstance(dspec2, NoneType):
            assert dspec1.shape == dspec2.shape, "Input data must be equal in shape!"
            self.xps = self._make_xps(dspec1, dspec2, axis=None)
        else:
            self.xps = None

    def copy(self):
        """
        Return a copy of self
        """
        return copy.deepcopy(self)

    def _make_xps(self, dspec1, dspec2, axis=None):
        """
        Generate a matrix with all possible cross-power spectra along a given axis.
        """

        if axis == None:
            return dspec1 * dspec2.conjugate()

        elif isinstance(axis, int):
            dspec1 = np.moveaxis(dspec1, axis, 0)
            dspec2 = np.moveaxis(dspec1, axis, 0)

            N = dspec1.shape[0]

            xps = [
                [dspec1[i] * dspec2[j].conjugate() for i in range(N)] for j in range(N)
            ]
            return np.moveaxis(xps, (0, 1), (axis, axis + 1))
        else:
            raise Exception("Invalid axis!")
    
    def get_cosmo_scaling(self, Ae, closure_scaling=False):
        """ 
        Get power spectrum scaling for a given cosmology and beam area
        """
        lam = c / self.freq
        Jy_to_mK = (lam ** 2 / (2 * kb)).to(u.mK / u.Jy)
        z = ftoz(self.freq, f_HI)
        Dc = self.cosmo.comoving_distance(z)
        dDc = comoving_depth(self.bandwidth, self.freq, f_HI, self.cosmo)

        j1 = Ae / (lam ** 2 * self.bandwidth)
        j2 = Dc ** 2 * dDc / self.bandwidth

        scaling = Jy_to_mK ** 2 * j1 * j2
        
        if closure_scaling:
            scaling *= 2 / 3
        
        return scaling

    def apply_cosmo_units(self, Ae, **kwargs):
        """
        Apply cosmological units to delay power spectrum
        """
        scaling = self.get_cosmo_scaling(Ae, **kwargs)
        self.xps = (self.xps * scaling).to(u.mK ** 2 * u.Mpc ** 3)

        return self

    def plot_power(self, ax=None, path=None, **kwargs):
        """
        Plot cross-power spectrum
        """

        delay = self.get_delays().to(u.us).value
        k_prll = self.get_k_parallel().to(1 / u.Mpc).value

        if ax == None:
            fig, ax = plt.subplots()

        ax.plot(k_prll, 2/3 * self.xps.real, **kwargs)
        ax.set_xlabel(r"$\kappa_{||}$ (pseudo h Mpc$^{-1}$)")
        ax.set_xlim([-k_prll[-1], k_prll[-1]])

        # delay axis
        dax = ax.twiny()
        dax.set_xlabel(r"Delay $(\mu \mathrm{s})$")
        dax.set_xlim([-delay[-1], delay[-1]])

        ax.set_yscale("log")
        ax.set_ylabel(
            r"$\frac{2}{3}P_\bigtriangledown\left(\kappa_{||}\right)$ (pseudo $\mathrm{mK}^2\mathrm{h}^{-3}\mathrm{Mpc}^3$)"
        )

        ax.minorticks_on()
        dax.minorticks_on()

        # save plot
        if isinstance(path, str):
            plt.savefig(path)

        return ax