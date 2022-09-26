"""

make_beam_image.py

Authors: Pascal M. Keller
Affiliation: Cavendish Astrophysics
Created on: 08/02/2022

Description:
    Make a spectral image of an antenna beam on the image plane (l, m)
    and save to disc. 
    
    Example: 
    
    python3 make_beam_image.py -i "${beam_file}" -o "${outfile}" -w 0.2625456 --np 256 --nf 167 --fmin 152.3 --fmax 167.94

"""

import argparse
import numpy as np
import os
import sys

import astropy.units as u

sys.path.append(os.path.join(os.getcwd(), "skysim/"))
from beam import BeamModel

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", help="Path to beam FITS-file.", type=str)
parser.add_argument("-o", "--outfile", help="File to write beam image to.", type=str)
parser.add_argument("-w", help="Half-width of the beam image.", type=float)
parser.add_argument("--np", help="Number of pixels per axis.", type=int)
parser.add_argument("--nf", help="Number of frequencies.", type=int)
parser.add_argument("--fmin", help="Minimum frequency in MHz.", type=float)
parser.add_argument("--fmax", help="Maximum frequency in MHz.", type=float)

args = parser.parse_args()

# make data header
path, fname = os.path.split(args.infile)

# spectral cuboid coordinates
freq_array = np.linspace(args.fmin, args.fmax, args.nf) * u.MHz

# generate beam image and write to output directory
model = BeamModel()

print(f"Reading beam FITS {args.infile}")
model.read_uvbeam(args.infile)

print("Generating beam image")
model.make_image(args.np, args.w, freq_array)

print(f"Writing beam FITS to {args.outfile} \n")
model.write_image(args.outfile, fname)

print("All done. \n")