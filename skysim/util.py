"""

util.py

Authors: Pascal M. Keller
Affiliation: Cavendish Astrophysics
Created on: Fr 28th January 2022

Description:
    Tools and Utilities

"""

import re
import numpy as np
import astropy.units as u
from astropy.cosmology import Planck18

f_HI = 1.4204057517667e9 * u.Hz  # H I rest frequency [Hz]
c = 2.99792458e8 * u.m / u.s  # speed of light [m/s]
kb = 1.3806503e-23 * u.J / u.K  # Boltzmann [J/K]
HERA_LAT = -30.72138329631366 * u.deg # Latitude of HERA array [deg]

NoneType = type(None)


def omit_units(quants, units):
    """
    Omit units of a list of quantities
    """

    for i, quant in enumerate(quants):
        if isinstance(quant, u.Quantity):
            quants[i] = quant.to(units[i]).value

    return quants


def extract_from_str(string, dt):
    """
    Extract floats from string
    """
    str_list = re.findall(r"\d*\.\d+|\d+", string)
    return [dt(s) for s in str_list]


def ftoz(fobs, fem=f_HI):
    """
    Frequency to redshift
    """
    return (fem - fobs) / fobs


def ztof(z, fem=f_HI):
    """
    Redshift to frequency
    """
    return fem / (z + 1)


def comoving_depth(B, fc, f0=f_HI, cosmology=Planck18):
    """
    Convert a bandwidth to a comoving depth
    """
    z = ftoz(fc, f0)
    dD = B * c * (1 + z) ** 2 / f0 / cosmology.H(z)

    return dD.to(u.Mpc)


def comoving_depth_to_bandwidth(dD, fc, f0=f_HI, cosmology=Planck18):
    """
    Convert comoving depth to bandwidth
    """
    z = ftoz(fc, f0)
    B = dD * f0 * cosmology.H(z) / c / (1 + z) ** 2

    return B


def altaz(ra, dec, lat, lst, units=u.deg):
    """
    Transform J2000 RA/DEC to Alt/Az given Latitude and LST.
    """

    # local hour angle
    h = lst - ra

    # altitude
    alt = np.arcsin(np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(h))

    # azimuth
    x = -np.sin(lat) * np.cos(dec) * np.cos(h) + np.cos(lat) * np.sin(dec)
    y = np.cos(dec) * np.sin(h)
    az = np.reshape(-np.arctan2(y, x), -1)

    # ensure azimuth takes values in [0, 2*pi]
    idx = np.where(az < 0)
    az[idx] = az[idx] + 2 * np.pi * u.rad

    if len(az) == 1:
        az = az[0]
    if type(units) is not tuple:
        units = 2 * (units,)

    return alt.to(units[0]), az.to(units[1])


def dircos(alt, az):
    """
    Compute direction cosines
    """

    l = np.cos(alt) * np.sin(az)
    m = np.cos(alt) * np.cos(az)

    return l, m


def dircos_to_altaz(l, m):
    """
    Convert direction cosines to alt/az coordinates
    """

    alt = np.arccos(np.sqrt(l ** 2 + m ** 2))
    az = np.arctan(l / m)

    # ensure azimuth takes values in [0, 2*pi]
    idx = np.where(az < 0)
    az[idx] = az[idx] + 2 * np.pi

    return alt, az


def make_image(x, y, z, xrange=None, yrange=None, L=(1024, 1024)):
    """
    Make an image of z-data with (x, y) coordinates
    """

    if type(L) is not tuple:
        L = 2 * (L,)
    if xrange == None:
        xrange = (np.min(x), np.max(x))
    if yrange == None:
        yrange = (np.min(y), np.max(y))

    x_grid = np.linspace(xrange[0], xrange[1], L[0])
    y_grid = np.linspace(yrange[0], yrange[1], L[1])
    xx, yy = np.meshgrid(x_grid, y_grid)
    zz = np.zeros(xx.shape)

    for i in range(len(z)):
        l = np.argmin(np.abs(x_grid - x[i]))
        k = np.argmin(np.abs(y_grid - y[i]))

        zz[k, l] = zz[k, l] + z[i]

    return xx, yy, zz