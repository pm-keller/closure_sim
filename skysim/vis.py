"""
vis.py

Authors: Pascal M. Keller
Affiliation: Cavendish Astrophysics
Created on: 07/02/2022

"""

import numpy as np
import h5py
import copy
import sys
import os

from scipy.signal.windows import get_window

import astropy.units as u

import matplotlib.pyplot as plt

sys.path.append(os.getcwd())
from antennas import Array

from util import c, NoneType


class VisData(Array):
    def __init__(self, ant_pos=None, ant_num=None, freq_array=None, **kwargs):
        self.vis = np.atleast_2d(None)
        self.freq_array = freq_array
        Array.__init__(self, ant_pos, ant_num, **kwargs)

    def copy(self):
        """
        Return a copy of self
        """
        return copy.deepcopy(self)

    def isempty(self):
        """
        Check if visibility data is empty
        """
        return None in self.vis

    def add_axis(self, axis=1):
        """
        Add an axis to visibility data
        """
        axis = np.atleast_1d(axis)

        for ax in axis:
            assert (
                ax > 0 and ax < self.vis.ndim + ax and ax != -1
            ), "Axes 0 and -1 are reserved for different baselines and frequencies respectively. Choose a different axis."

            vis = self.vis[np.newaxis]
            self.vis = np.moveaxis(vis, 0, ax)

    def append(self, VisObj, axis):
        """
        Append a visibility data set along axis
        """
        if self.isempty():
            self.vis = VisObj.vis
        else:
            try:
                self.vis = np.concatenate([self.vis, VisObj.vis], axis)
            except:
                raise Exception("Shapes are incompatible!")

    def _make_vis(self, bls_vectors, sky, coords):
        """
        Compute vsibilities given some baseline vectors, a sky model with image coordinates (l, m).
        """

        # image coordinates
        l_array, m_array = coords

        # frequency in Hz
        freq_array = self.freq_array.to(u.Hz)

        # if baseline vectors are 3D, append zeros to coordinate axis
        if bls_vectors.shape[-1] == 3:
            zeros = np.zeros_like(l_array)
            s_vectors = np.array([l_array, m_array, zeros])
        else:
            s_vectors = np.array([l_array, m_array])

        # compute visibility phases of sky components
        b_dot_s = np.dot(bls_vectors, s_vectors.swapaxes(0, -2))
        phases = np.array([2 * np.pi * b_dot_s * freq / c for freq in freq_array])
        phases = np.moveaxis(phases, 0, -1)

        # sum over all sky components
        sum_axes = tuple(range(phases.ndim - s_vectors.ndim, s_vectors.ndim))
        vis = (np.exp(1j * phases) * sky).sum(sum_axes)

        return vis

    def vis_from_image(self, sky, coords, inplace=True):
        """
        Compute visibilities from sky image
        """
        assert not isinstance(self.bls_vectors, NoneType), "Baselines not specified!"
        assert not isinstance(self.freq_array, NoneType), "Frequencies not specified!"

        if inplace:
            VisObj = self
        else:
            VisObj = self.copy()

        VisObj.vis = self._make_vis(self.bls_vectors, sky, coords)

        return VisObj

    def write_vis(
        self, path, data_list=None, names=None, axis_order=None, overwrite=False
    ):
        """
        Write visibility data to HDF5 file
        """

        if not isinstance(axis_order, NoneType):
            default_order = tuple(np.arange(0, vis.ndim))
            vis = np.moveaxis(self.vis, default_order, axis_order)
        else:
            vis = self.vis

        if overwrite and os.path.exists(path):
            os.remove(path)

        f = h5py.File(path, "w")

        f.create_dataset("Antenna Positions", data=self.ant_pos)
        f.create_dataset("Antenna Numbers", data=self.ant_num)
        f.create_dataset("Baseline Vectors", data=self.bls_vectors)
        f.create_dataset("Baseline Pairs", data=self.bls_pairs)
        f.create_dataset("Frequencies", data=self.freq_array)
        f.create_dataset("Visibilities", data=vis)

        if not isinstance(data_list, NoneType):
            assert len(data_list) == len(names), "Data and names must have same length!"

            for data, name in zip(data_list, names):
                f.create_dataset(str(name), data=data)
        f.close()

    def read_vis(self, path, names=None):
        """
        Read visibility data from HDF5 file
        """
        with h5py.File(path, "r") as f:
            self.ant_pos = f["Antenna Positions"][()]
            self.ant_num = f["Antenna Numbers"][()]
            self.bls_vectors = f["Baseline Vectors"][()]
            self.bls_pairs = f["Baseline Pairs"][()]
            self.freq_array = f["Frequencies"][()] * u.MHz
            self.vis = f["Visibilities"][()]

            if not isinstance(names, NoneType):
                names = np.atleast_1d(names)
                self.header_data = dict()

                for name in names:
                    self.header_data.update({name: f[name][()]})


    def cut_nindex(self, lower=0, upper=0, axis=-1):
        """ 
        Trim visibility data along a given axis
        """

        vis = np.moveaxis(self.vis, axis, 0)

        if upper == 0:
            upper = -vis.shape[0]

        vis = vis[lower:-upper]
        self.vis = np.moveaxis(vis, 0, axis)

        if axis == -1:
            self.freq_array = self.freq_array[lower:-upper]


    def _get_bispec(self, triad):
        """
        Get the bispectrum for a given triad
        """

        triad_pairs = self.get_triad_pairs(triad)
        bispec = np.ones(self.vis.shape[1:], dtype=complex)

        for bls_pair in triad_pairs:
            idx, signs = self.get_bls_idx(bls_pair, return_signs=True)

            if signs[0] == 1:
                bispec *= self.vis[idx][0]
            elif signs[0] == -1:
                bispec *= self.vis[idx][0].conjugate()
        return bispec

    def get_bispec(self, triads=None):
        """
        Get bispectra
        """
        if triads == None:
            triads = self.get_triads()

        triads = np.atleast_2d(triads)
        bispec = np.ones((len(triads),) + self.vis.shape[1:], dtype=complex)

        for i, triad in enumerate(triads):
            bispec[i] = self._get_bispec(triad)

        return np.squeeze(bispec)

    def get_closure_phases(self, triads=None):
        """
        Get closure phases
        """
        return np.angle(self.get_bispec(triads))

    def get_eicp(self, triads=None, veff=False, veff_cut=0):
        """
        Get complex exponential of closure phase
        """
        if veff:
            amp = self.get_veff(triads)
            amp[np.where(amp < veff_cut)] = np.nan
        else:
            amp = 1

        cp = np.moveaxis(self.get_closure_phases(triads), -1, 0)
        return np.moveaxis(amp * np.exp(1j * cp), 0, -1)

    def _get_veff(self, triad, window):
        """
        Get the effective visibility for a triad
        """

        window_sum = window.sum(axis=-1)

        triad_pairs = self.get_triad_pairs(triad)
        idx = self.get_bls_idx(triad_pairs)
        v1, v2, v3 = self.vis[idx]
        v1 = np.sum(np.abs(v1) * window, axis=-1) / window_sum
        v2 = np.sum(np.abs(v2) * window, axis=-1) / window_sum
        v3 = np.sum(np.abs(v3) * window, axis=-1) / window_sum

        return np.sqrt(1 / (1 / v1**2 + 1 / v2**2 + 1 / v3**2))

    def get_veff(self, triads=None, window="blackmanharris"):
        """
        Get the effective visibilities for given triads
        """
        if triads == None:
            triads = self.get_triads()

        triads = np.atleast_2d(triads)
        veff = []

        if isinstance(window, str):
            window = get_window(window, len(self.freq_array))

        for triad in triads:
            veff += [self._get_veff(triad, window)]

        return np.squeeze(veff)

    def plot_vis_spec(
        self, bls_pair, index=None, ax=None, entity="abs", path=None, **kwargs
    ):
        """
        Plot visibility spectrum of baseline pairs
        """
        entities = ["abs", "real", "imag"]
        assert entity in entities, f"entity must be one of {entities}"

        if ax == None:
            fig, ax = plt.subplots()

        idx, sign = self.get_bls_idx(bls_pair, return_signs=True)

        if index != None:
            idx = (idx, *tuple(index))

        vis = self.vis[idx][0]

        assert vis.ndim == 1, "Visibility spectrum must be 1D. Set index."

        if sign[0] == -1:
            vis = vis.conjugate()

        freq_array = self.freq_array.to(u.MHz).value

        if entity == "abs":
            ax.plot(freq_array, np.abs(vis), label=f"{bls_pair}", **kwargs)
            ax.set_ylim(bottom=0)
            ax.set_ylabel(r"$|V|$ (Jy)")
        if entity == "real":
            ax.plot(freq_array, vis.real, label=f"{bls_pair}", **kwargs)
            ax.set_ylabel(r"$\mathrm{Re}(V)$ (Jy)")
        if entity == "imag":
            ax.plot(freq_array, vis.imag, label=f"{bls_pair}", **kwargs)
            ax.set_ylabel(r"$\mathrm{Im}(V)$ (Jy)")

        ax.set_xlim([np.min(freq_array), np.max(freq_array)])
        ax.set_xlabel("Frequency (MHz)")

        # save plot
        if isinstance(path, str):
            plt.savefig(path)

        return ax

    def plot_closure_spec(
        self, triad, index=None, ax=None, entity="abs", path=None, **kwargs
    ):
        """
        Plot visibility spectrum of baseline pairs
        """
        if ax == None:
            fig, ax = plt.subplots()

        freq_array = self.freq_array.to(u.MHz).value

        if entity == "abs":
            ax.plot(freq_array, self.get_closure_phases(triad)[index], label=f"{triad}", **kwargs)
            ax.set_ylim(bottom=0)
            ax.set_ylabel(r"Closure Phase (rad)")

        ax.set_xlim([np.min(freq_array), np.max(freq_array)])
        ax.set_ylim([-np.pi, np.pi])
        ax.set_xlabel("Frequency (MHz)")

        # save plot
        if isinstance(path, str):
            plt.savefig(path)

        return ax
