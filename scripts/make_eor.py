"""

make_eor.py

Authors: Pascal M. Keller
Affiliation: Cavendish Astrophysics
Created on: 08/02/2022

Description:
    Tile EoR model and downsize by smoothing and downsampling. Write to disc.
    
    Example: 
    
    python3 make_eor.py -i "${eor_file}" -o "${outfile}" -d 1.6 -z 5.0 --np 1024 --fmin 152.3 --fmax 167.94 --ntile 3 -r 256

"""

import argparse
import numpy as np
import os
import sys

import astropy.units as u

sys.path.append(os.path.join(os.getcwd(), "skysim/"))
from eor import EORModel

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", help="Path to EoR model.", type=str)
parser.add_argument("-o", "--outfile", help="File to write new EoR model to.", type=str)
parser.add_argument(
    "-d", "--dimension", help="Dimension of EoR cube in Gpc.", default=1.6, type=float
)
parser.add_argument(
    "-z", "--redshift", help="Start redshift of EoR cube", default=5, type=float
)
parser.add_argument("--np", help="Number of pixels per axis.", default=1024, type=int)
parser.add_argument("--fmin", help="Minimum frequency in MHz.", type=float)
parser.add_argument("--fmax", help="Maximum frequency in MHz.", type=float)
parser.add_argument(
    "--ntile",
    help="Number of times to tile EoR along first two axes",
    default=3,
    type=int,
)
parser.add_argument(
    "-r", "--resolution", help="Resolution of EoR cube", default=256, type=int
)


args = parser.parse_args()

# Generate an EoR model
model = EORModel()

print(f"Reading EoR: {args.infile}")
model.read(args.infile, args.np, args.dimension * u.Gpc, args.redshift)

print(f"Selecting frequency range {args.fmin}-{args.fmax} MHz")
model.select_freq_range(args.fmin * u.MHz, args.fmax * u.MHz)

print(f"Tiling EoR cuboid with n={args.ntile}")
model.tile(n=(args.ntile, args.ntile, 1))

nbox = int(model.shape[0] / args.resolution)
print(f"Smoothing EoR with box kernel of shape {(nbox, nbox, 1)}")
kernel = np.ones((nbox, nbox, 1)) / nbox ** 2
model.smooth(kernel)

print(f"Downsample EoR to shape {(args.resolution, args.resolution, model.shape[-1])}")
model.downsample((nbox, nbox, 1))

print(f"Write EoR to file {args.outfile} \n")

# write to file
path, fname = os.path.split(args.infile)
model.write(args.outfile, fname)

print(f"All done. \n")
