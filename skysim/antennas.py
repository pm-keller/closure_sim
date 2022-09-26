"""

array.py

Authors: Pascal M. Keller
Affiliation: Cavendish Astrophysics
Created on: 04/02/2022

Description:

Notes:
    Use with astropy units!

"""

import os
import sys
import numpy as np
from itertools import product
import matplotlib.pyplot as plt

import astropy.units as u

sys.path.append(os.getcwd())
from util import NoneType


class Array:
    def __init__(self, ant_pos=None, ant_num=None, ant_diameter=None, ant_flags=None, include_autos=False):
        self.ant_pos = np.atleast_2d(ant_pos)
        self.ant_num = ant_num
        self.Nants = len(self.ant_pos)
        self.ant_diameter = ant_diameter
        self.include_autos = include_autos

        if isinstance(ant_num, NoneType):
            self.ant_num = np.arange(self.ant_pos.shape[0])

        if isinstance(ant_flags, NoneType):
            ant_flags = np.zeros_like(ant_num).astype(bool)
        self.ant_flags = ant_flags

        if not isinstance(ant_pos, NoneType):
            self._get_bls()
            

    def isempty(self):
        """
        Check if array data is empty
        """
        return None in self.ant_pos

    def _get_bls(self):
        """
        Get baselines associated with antenna positions
        """
        bls_vectors = self.ant_pos[np.newaxis] - self.ant_pos[:, np.newaxis]

        if self.include_autos:
            idx = np.triu_indices(self.Nants, k=0)
        else:
            idx = np.triu_indices(self.Nants, k=1)

        self.bls_vectors = bls_vectors[idx]
        self.bls_pairs = np.array(
            list(zip(self.ant_num[idx[0]], self.ant_num[idx[1]]))
        )
        self.Nbls = len(self.bls_pairs)


    def read_ant_file(self, path, antnum=True, **kwargs):
        """
        Read antenna positions from file
        """
        data = np.loadtxt(path)

        # get antenna numbers from file
        if antnum:
            ant_num = data[:, 0]
            ant_pos = data[:, 1:] * u.m
        else:
            ant_pos = data[1:] * u.m 
            ant_num = np.arange(ant_pos.shape[0])

        self.__init__(ant_pos=ant_pos, ant_num=ant_num, **kwargs)

    def select(self, idx):
        """
        Select a set of antennas
        """

        self.ant_pos = self.ant_pos[idx]
        self.ant_num = self.ant_num[idx]
        self.Nants = len(self.ant_pos)
        self._get_bls()

    def get_uvw_array(self):
        """
        Returns an uvw-array
        """
        if self.bls_vectors.shape[-1] == 3:
            return self.bls_vectors
        else:
            zeros = np.zeros(self.bls_vectors.shape[:-1])
            zeros = zeros[..., np.newaxis]
            return np.append(self.bls_vectors, zeros, axis=-1)

    def remove_redundancies(self, ndec=2):
        """
        Remove all redundant baselines
        """
        idx = np.unique(
            self.bls_vectors.round(decimals=ndec), axis=0, return_index=True
        )[1]
        self.bls_vectors = self.bls_vectors[idx]
        self.bls_pairs = self.bls_pairs[idx]

        ant_num = np.unique(self.bls_pairs.reshape(-1))
        self.ant_num, idx1, idx2 = np.intersect1d(
            ant_num, self.ant_num, return_indices=True
        )
        self.ant_pos = self.ant_pos[idx2]
        self.Nants = len(self.ant_pos)
        self.Nbls = len(self.bls_pairs)

    def get_bls_idx(self, bls_pairs, return_signs=False):
        """
        Get the index values corresponding to a list of given baseline pairs
        """
        bls_pairs = np.atleast_2d(bls_pairs)
        idx, signs = [], []

        for bls_pair in bls_pairs.tolist():
            ant1, ant2 = bls_pair

            if ant1 > ant2:
                signs += [-1]
                bls_pair = [ant2, ant1]
            else:
                signs += [1]

            idx += [self.bls_pairs.tolist().index(bls_pair)]

        if return_signs:
            return np.array(idx), np.array(signs)
        else:
            return np.array(idx)

    def get_bls_len(self):
        """
        Get baseline lengths
        """
        return np.sqrt((self.bls_vectors ** 2).sum(axis=-1))

    def get_triads(self):
        """
        Get all possible triads of array
        """
        anum = self.ant_num
        triads = np.array(np.meshgrid(anum, anum, anum)).reshape(-1, 3)
        triads = np.unique(triads, axis=0)

        idx = [np.unique(triad).shape[0] == 3 for triad in triads]
        return triads[idx]

    def get_triad_pairs(self, triad):
        """
        Get baseline pairs associated with a list of antenna numbers
        """

        return np.array([[triad[i], triad[(i + 1) % 3]] for i in range(3)])

    def get_ant_pos(self, ant_num):
        """
        Get all antenna positions for a list of given antennas
        """

        ant_num, _idx, idx = np.intersect1d(ant_num, self.ant_num, return_indices=True)
        return self.ant_pos[idx]

    def plot_antennas(self, ax=None, color="k"):
        """
        Plot antenna positions
        """

        if ax == None:
            fig, ax = plt.subplots()

        x, y = self.ant_pos.value.T[[0, 1]]

        for ant_num, flag, xi, yi in zip(self.ant_num, self.ant_flags, x, y):
            if flag:
                alpha = 0.0
            else:
                alpha = 1.0
            cc = plt.Circle((xi , yi), self.ant_diameter / 2 , alpha=alpha, color=color, fc="None") 
            ax.set_aspect(1) 
            ax.add_artist(cc)
            ax.text(xi, yi, ant_num, ha="center", va="center", alpha=alpha, clip_on=True, fontsize=7)

        ax.set_xlabel("East-West position (m)")
        ax.set_ylabel("North-South position (m)")

        return ax

    def plot_triads(self, triads, ax=None):
        """
        Plot antenna positions
        """

        if ax == None:
            fig, ax = plt.subplots()

        x, y = self.ant_pos.value.T[[0, 1]]

        for i, triad in enumerate(triads):
            pos = self.get_ant_pos(triad)[:, :2].value
            cc = plt.Polygon(pos, color="k", fc="None") 
            ax.set_aspect(1) 
            ax.add_artist(cc)
            ax.text(np.mean(pos, axis=0)[0], np.mean(pos, axis=0)[1], i, ha="center", va="center", clip_on=True)

        ax.set_xlabel("East-West position (m)")
        ax.set_ylabel("North-South position (m)")

        return ax