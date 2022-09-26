"""

make_gleam_catalog.py

Authors: Pascal M. Keller
Affiliation: Cavendish Astrophysics
Created on: 15/02/2022

Description:
    Query GLEAM data and write to disc.
    
    Example: 

    python3 make_gleam_catalog.py -o "${gleam_file}" --Fmin 0.05

"""

import argparse
import os
import sys

import astropy.units as u

sys.path.append(os.path.join(os.getcwd(), "skysim/"))
from gleam import GLEAMData

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfile", help="File to write GLEAM data to.", type=str)
parser.add_argument("--Fmin", help="Minimum flux density at 151 MHz in Jy", type=float)

args = parser.parse_args()

columns = [
    "RAJ2000",
    "DEJ2000",
    "Fint143",
    "Fint151",
    "Fint158",
    "Fint166",
    "Fint174",
    "Fintfit200",
    "alpha",
]
filters = {"Fint151": f">{args.Fmin}"}

gleam = GLEAMData(columns, filters)

print("\nQuerying GLEAM data")
gleam.query_all()

print("Fitting spectral indices")
gleam.fit_spectral_index()

print(f"Writing GLEAM data to {args.outfile}\n")
gleam.write(args.outfile, overwrite=True)

print("All done.")