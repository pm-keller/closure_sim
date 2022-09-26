"""

make_veff.py

Authors: Pascal M. Keller
Affiliation: Cavendish Astrophysics
Created on: 15/02/2022

Description:
    Compute power effective visibilities

    Example:
    python3 make_veff.py -g ${gleamfile} -b ${brightfile} -o ${outfile} -n 16

"""

import argparse
import numpy as np
import os
import sys

import astropy.units as u

sys.path.append(os.path.join(os.getcwd(), "skysim/"))
from vis import VisData


parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gleampath", help="Path to GLEAM visibilities.", type=str)
parser.add_argument("-b", "--brightpath", help="Path to visibilities of bright sources.", type=str)
parser.add_argument("-o", "--outfile", help="Name of file to write effective visibilities to.", type=str)
parser.add_argument("-n", help="Number of neighbouring Veff to average", type=int)
parser.add_argument("-l", help="Number of lower channels to cut", type=int, default=0)
parser.add_argument("-u", help="Number of upper channels to cut", type=int, default=0)

args = parser.parse_args()

# get gleam visibilities
gleam = VisData()
gleam.read_vis(args.gleampath)

# get visibilities of bright sources
bright = VisData()
bright.read_vis(args.brightpath)

# add visibilities
vis = gleam.copy()
vis.vis += bright.vis

# trim visibility data along frequency axis
vis.cut_nindex(lower=args.l, upper=args.u, axis=-1)

# get effective visibilities
veff = vis.get_veff()

# average effective visibilities
#veff = np.array([np.mean(veff[:, i*args.n:(i+1)*args.n], axis=1) for i in range(veff.shape[1] // args.n)])

# write to file
np.savetxt(args.outfile, veff, header="V_eff [Jy]")

print("All done \n")