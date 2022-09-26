"""

make_scaling.py

Authors: Pascal M. Keller
Affiliation: Cavendish Astrophysics
Created on: 15/02/2022

Description:
    Compute power spectrum scaling

    Example:
    python3 make_scaling.py -p "${vispath}"

"""

import argparse
import numpy as np
import os
import sys

from scipy import interpolate

import astropy.units as u
from astropy.cosmology import Planck18

sys.path.append(os.path.join(os.getcwd(), "skysim/"))
from dspec import XPSpec
from util import c


parser = argparse.ArgumentParser()
parser.add_argument("-v", "--veffpath", help="Path to effective visibility file.", type=str)
parser.add_argument("-a", "--omegapath", help="Path to integrated squared beam area file.", type=str)
parser.add_argument("-o", "--outfile", help="Path to write scaling coefficients to.", type=str)
parser.add_argument("--fmin", help="Minimum frequency in MHz", type=float)
parser.add_argument("--fmax", help="Maximum frequency in MHz", type=float)
parser.add_argument("--nf", help="Number of frequencies", type=int)

args = parser.parse_args()

# load data
veff = np.loadtxt(args.veffpath) * u.Jy
#veff[np.where(veff < 5 * u.Jy)] *= np.nan 
Omega2 = np.loadtxt(args.omegapath)

# compute beam area
BW = (args.fmax - args.fmin) * 1e6
fB = np.linspace(100, 200, 1024)
fb = np.linspace(args.fmin, args.fmax, args.nf)
Ae = np.mean((c.value / fb * 1e-6)**2 / interpolate.interp1d(fB, Omega2, kind="cubic")(fb))

# compute scaling
ps = XPSpec(freq_array=fb * u.MHz, cosmology=Planck18)
scaling = veff ** 2 * ps.get_cosmo_scaling(Ae * u.m**2, closure_scaling=True) / u.s**2
scaling = scaling.to(u.mK ** 2 * u.Mpc ** 3).value * Planck18.h**3

# write to file
np.savetxt(args.outfile, scaling, header="A [mK^2 * Mpc^3]")

print("All done \n")