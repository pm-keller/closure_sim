"""

gleam.py

Authors: Pascal M. Keller
Affiliation: Cavendish Astrophysics
Created on: 28/01/2022

Description:
    Tools for working with the GLEAM catalogue
    Use astropy units!

"""


import os
import sys
import copy
import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.table import Table
from astroquery.vizier import Vizier
import astropy.coordinates as coord
import astropy.units as u

sys.path.append(os.getcwd())

from util import altaz, dircos, make_image

class GLEAMData:
    def __init__(self, columns=["*"], filters={}):
        """
        Parameters
        ----------
        columns : list, optional
            Columns to keep, by default ["*"]
        filters : dict, optional
            VizieR filters, by default {}
        """

        self.columns = columns
        self.filters = filters
        self.sources = Table(names=('RAJ2000', 'DEJ2000', 'Fintfit200', 'alpha'))
        self.ra = None
        self.dec = None
        self.radius = None
        self.catalog_names = {"EGC": "VIII/100/gleamegc", "GAL": "VIII/102/gleamgal"}

    def copy(self):
        """
        Return a copy of self
        """
        return copy.deepcopy(self)

    def _merge(self, tables):
        """ 
        Merge rows of two FITS tables
        """

        table_1 = tables[0]
        
        if len(tables) > 1:
            for table in tables[1:]:
                for row in table:
                    table_1.add_row(row)
        
        return table_1

    def _set_region(self, ra, dec, radius):
        """
        Set region variables
        """
        if ra.unit == u.h:
            ra = ra.value * 15 * u.deg

        self.ra = ra
        self.dec = dec
        self.radius = radius

    def query_all(self):
        """
        Query whole GLEAM catalogue
        """
        res = []

        for catalog_name in self.catalog_names.values():
            vizier = Vizier(self.columns, self.filters, catalog_name, row_limit=-1)
            res.append(vizier.get_catalogs(catalog_name)[catalog_name])

        self.sources = self._merge(res)

    def query_region(self, ra, dec, radius):
        """
        Query sky region
        """
        res = []
        coordinates = coord.SkyCoord(ra=ra, dec=dec, frame="icrs")

        for catalog_name in self.catalog_names.values():
            vizier = Vizier(self.columns, self.filters, catalog_name, row_limit=-1)
            table = vizier.query_region(coordinates, radius=radius)
            
            if len(table) > 0:
                res.append(table[catalog_name])

        self.sources = self._merge(res)
        self._set_region(ra, dec, radius)

    def select_region(self, ra, dec, radius, inplace=False):
        """
        Select sources within a region of the sky
        """
        assert "RAJ2000" in self.sources.keys() and "DEJ2000" in self.sources.keys()

        srcs_coord = np.array([self.sources["RAJ2000"], self.sources["DEJ2000"]])
        alt = altaz(*srcs_coord * u.deg, dec, ra, units=u.deg)[0]

        if inplace:
            gleam_obj = self
        else:
            gleam_obj = self.copy()

        gleam_obj.sources = gleam_obj.sources[alt > 90 * u.deg - radius.to(u.deg)]
        gleam_obj._set_region(ra, dec, radius)

        return gleam_obj

    def add_source(self, ra, dec, f0, F0, alpha):
        """
        Add a source the souce list
        """

        F = F0 * (200 * u.MHz / f0) ** alpha
        row = {
            "RAJ2000": ra.to(u.deg).value,
            "DEJ2000": dec.to(u.deg).value,
            "Fintfit200": F.to(u.Jy).value,
            "alpha": alpha,
        }
        self.sources.add_row(row)

    def apply_filters(self):
        """
        Apply filters to source list
        """

        for column in self.filters:
            symbol = self.filters[column][0]
            value = float(self.filters[column][1:])
            col = self.sources[column]

            if symbol == "<":
                self.sources = self.sources[np.where(col.data < value)]
            elif symbol == ">":
                self.sources = self.sources[np.where(col.data > value)]

            if self.columns != ["*"]:
                self.sources = self.sources[self.columns]

    def read(self, path="./gleam.fits"):
        """
        Read GLEAM sources from FITS file
        """

        hdu_list = fits.open(path)
        self.sources = Table(hdu_list[1].data)
        self.apply_filters()

    def write(self, path="./gleam.fits", **kwargs):
        """
        Write GLEAM sources to FITS file
        """
        self.sources.write(path, **kwargs)


    def get_freq(self, return_names=False):
        """
        Get frequencies of the integrated flux densities
        """

        keys = self.sources.keys()
        column_names = [name for name in keys if "Fint" in name and name[4:].isdigit()]
        freq = np.array([int(name[4:]) for name in column_names])

        if return_names:
            return freq, column_names
        else:
            return freq

    def _spec(self, freq, F0, alpha):
        """
        Generate spectrum characterised by spectral index
        """
        return F0 * (freq / 200) ** alpha

    def _spec_residuals(self, coeff, freq, F):
        """
        Residuals of spectral model
        """
        return F - self._spec(freq, *coeff)

    def _fit_spectral_index(self, freq, F):
        """
        Fit spectral index
        """

        popt = leastsq(self._spec_residuals, x0=(np.min(F), 0), args=(freq, F))[0]
        return popt

    def fit_spectral_index(self):
        """
        Fit a spectral index to GLEAM fluxes where masked
        """
        assert "alpha" in self.sources.keys() and "Fintfit200" in self.sources.keys()

        freq, column_names = self.get_freq(return_names=True)

        for source in self.sources:
            F = np.array([source[name] for name in column_names])

            if np.ma.is_masked(source["alpha"]):
                source[["Fintfit200", "alpha"]] = self._fit_spectral_index(freq, F)

    def get_flux_array(self, freq_array):
        """
        Compute flux densities for an array of frequencies using spectral indices
        """

        f0 = 200 * u.MHz
        F0 = self.sources["Fintfit200"] * u.Jy
        alpha = self.sources["alpha"]

        flux_array = np.array(
            [F0 * (freq / f0) ** alpha for freq in freq_array.to(u.MHz)]
        )

        return flux_array

    def get_alt_az(self):
        """
        Get Alt/Az coordinates
        """
        assert "RAJ2000" in self.sources.keys() and "DEJ2000" in self.sources.keys()

        srcs_coord = np.array([self.sources["RAJ2000"], self.sources["DEJ2000"]])
        aaz = altaz(*srcs_coord * u.deg, self.dec, self.ra, units=u.deg)

        return aaz

    def get_za_az(self):
        """
        Get zenith angle and azimuth coordinates
        """
        alt, az = self.get_alt_az()
        za = np.pi / 2 * u.rad - alt

        return za, az

    def get_image_coordinates(self):
        """
        Get image plane coordinates
        """
        aaz = self.get_alt_az()

        return dircos(*aaz)

    def get_beam_sky(self, uvbeam, freq_array, pol=0):
        """
        Return a beam weighted GELAM sky
        """
        za, az = self.get_za_az()
        za, az = za.to(u.rad).value, az.to(u.rad).value

        flux_array = self.get_flux_array(freq_array)
        freq_array = freq_array.to(u.Hz).value

        beam = uvbeam.interp(az, za, freq_array=freq_array, reuse_spline=True)[0][
            0, 0, pol
        ]

        return np.moveaxis(beam * flux_array, 0, -1)

    def make_sky_image(self, freq_array, n=256, **kwargs):
        """
        Make a spectral image of GLEAM source
        """
        assert isinstance(
            self.radius, u.Quantity
        ), "Cant make all-sky image. Use select_region or query_region."
        width = np.sin(self.radius)
        im_range = (-width, width)

        l_array, m_array = self.get_image_coordinates()
        flux_array = self.get_flux_array(freq_array)

        l_grid, m_grid = 2 * (np.linspace(-width, width, n),)
        im = np.zeros((n, n, len(freq_array)))

        for i, (l, m) in enumerate(zip(l_array, m_array)):
            li, mi = np.argmin(np.abs(l_grid - l)), np.argmin(np.abs(m_grid - m))
            im[li, mi] += flux_array[..., i]

        return im

    def plot(self, freq, n=256, ax=None, **kwargs):
        """
        Plot GLEAM sources
        """

        if ax == None:
            fig, ax = plt.subplots()

        im = self.make_sky_image(np.atleast_1d(freq), n)
        hw = np.sin(self.radius)

        implt = ax.imshow(
            im[:, :, 0],
            cmap="bone",
            vmin=0.0,
            extent=[-hw, hw, -hw, hw],
            interpolation="gaussian",
            **kwargs
        )
        title = "RA/DEC: {:.2f}/{:.2f}, Radius: {:.2f}, Frequency: {:.2f}"
        ax.set_title(
            title.format(
                self.ra.value / 15 * u.h,
                self.dec,
                self.radius.to(u.deg),
                freq.to(u.MHz),
            )
        )
        ax.set_xlabel("l")
        ax.set_ylabel("m")

        plt.colorbar(implt, label="Flux Density (Jy)")

        return ax