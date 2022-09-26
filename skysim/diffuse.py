"""

diffuse.py

Authors: Pascal M. Keller
Affiliation: Cavendish Astrophysics
Created on: 04/02/2022

Description:
    Tools for working with diffuse radio emission.
    Use with astropy units!

"""

import os
import sys
import h5py
import numpy as np
import healpy as hp
from pygsm import GlobalSkyModel2016

import astropy.units as u 
from astropy.coordinates import SkyCoord

import matplotlib.pyplot as plt

sys.path.append(os.getcwd())

from util import NoneType, altaz, dircos, HERA_LAT


gstr = ["galactic", "G", "g"]
estr = ["equatorial", "E", "e"]

class DiffuseModel():
    def __init__(self):
        self.data = None

    def read_gsm(self, freq_array, freq_unit="MHz", unit="MJysr", resolution="low"):
        """
        Read global sky model using PyGSM 
        """
        gsm = GlobalSkyModel2016(freq_unit=freq_unit, unit=unit, resolution=resolution)
        self.data = gsm.generate(freq_array)
        self.freq = freq_array
        self.nest = False
        self.frame = "G"

        return gsm

    def read_healpix(self, path, freq, frame="G", nest=True):
        """
        Read diffuse model in healapix format.
        """
        ext = os.path.splitext(path)[-1]
        
        assert frame in gstr + estr, "Frame must be 'G', 'C' or 'E'."

        self.frame = frame
        self.nest = nest
        self.freq = np.array([freq])

        if ext in [".txt", ".dat"]:
            self.data = np.loadtxt(path)
        else:
            raise Exception("File extension not recognised!")

    def get_equatorial(self):
        """
        Get equatorial coordinates Dec/RA. 
        """

        assert not isinstance(self.data, NoneType), "Data array is empty."

        npix = self.data.shape[-1]
        nside = hp.npix2nside(npix)
        lat, lon = hp.pix2ang(nside, np.arange(npix), nest=self.nest)
        lat = np.pi / 2 - lat
        print(np.min(lat), np.max(lat), np.min(lon), np.max(lon))
        c = SkyCoord(l=lon*u.rad, b=lat*u.rad, frame='galactic')
        c = c.transform_to("icrs")

        return c.ra.deg, c.dec.deg
    
    def get_alt_az(self, lst, lat=HERA_LAT):
        """
        Get Alt/Az coordinates
        """
        ra, dec = self.get_equatorial()
        aaz = altaz(ra * u.deg, dec * u.deg, lat, lst.value * 15 * u.deg, units=u.deg)

        return aaz

    def get_za_az(self, lst, lat=HERA_LAT):
        """
        Get zenith angle and azimuth coordinates
        """
        alt, az = self.get_alt_az(lst, lat)
        za = np.pi / 2 * u.rad - alt

        return za, az

    def get_image_coordinates(self, lst, lat=HERA_LAT):
        """
        Get image plane coordinates
        """
        aaz = self.get_alt_az(lst, lat)

        return dircos(*aaz)

    def get_spec_data(self, freq_array, alpha=0.0):
        """
        Get diffuse data for a set of frequencies given a spectral index alpha.
        """
        Nfreq = len(freq_array)
        data = [self.data * (freq_array[i] / self.freq[0]) ** alpha for i in range(Nfreq)]
        return np.moveaxis(data, 0, -1)

    def get_sky_data(self, lst, lat=HERA_LAT, coordinates="image"):
        """
        Get sky at given LST and frequencies. Sky includes data above horizon only. 
        """
        za, az = self.get_za_az(lst, lat)
        idx = za < np.pi / 2 * u.rad
        sky = np.atleast_2d(self.data)[:, idx]

        if coordinates == "za/az":
            coord = (za, az)
        elif coordinates == "alt/az":
            coord = self.get_alt_az(lst, lat)
        elif coordinates == "image":
            coord = self.get_image_coordinates(lst, lat)
        else:
            raise Exception("Coordinates must be 'za/az', 'alt/az' ord 'image'.")

        coord = (coord[0][idx], coord[1][idx])

        return sky, coord
