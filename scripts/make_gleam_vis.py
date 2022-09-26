"""

make_gleam_vis.py

Authors: Pascal M. Keller
Affiliation: Cavendish Astrophysics
Created on: 15/02/2022

Description:
    Simulate visibilitites of a sky populated with GLEAM point sources.

    Example:
    python3 make_gleam_vis.py -g "${gleamfile}" -o "${outfile}" -b "${beamfile}" -a "${antfile}" --nf 161 --fmin 152.3 --fmax 167.94 --nt 752 --tmin 4.0 --tmax 6.25

"""

import argparse
import numpy as np
import h5py
import os
import sys

import astropy.units as u

sys.path.append(os.path.join(os.getcwd(), "skysim/"))
from gleam import GLEAMData
from beam import BeamModel
from vis import VisData

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gleamfile", help="Path to GLEAM file.", type=str)
parser.add_argument("-o", "--outfile", help="File to write visibilities to.", type=str)
parser.add_argument("-b", "--beamfile", help="Path to beam FITS-file.", type=str)
parser.add_argument(
    "-a", "--antfile", help="Path to file containing antenna positions.", type=str
)
parser.add_argument(
    "-r", "-radius", help="Sky radius to query in deg", default=15.0, type=float
)
parser.add_argument("--nf", help="Number of frequencies.", default=None, type=int)
parser.add_argument(
    "--fmin", help="Minimum frequency in MHz.", default=None, type=float
)
parser.add_argument(
    "--fmax", help="Maximum frequency in MHz.", default=None, type=float
)
parser.add_argument("--nt", help="Number of times.", type=int)
parser.add_argument("--tmin", help="Minimum LST in h.", type=float)
parser.add_argument("--tmax", help="Maximum LST in h.", type=float)
parser.add_argument(
    "--lat",
    help="Latitude of telescope in deg.",
    default=-30.72138329631366,
    type=float,
)

args = parser.parse_args()

print(f"\nReading GLEAM data {args.gleamfile}")
gleam = GLEAMData()
gleam.read(args.gleamfile)

print(f"Reading Beam FITS {args.beamfile}\n")
beam = BeamModel()
beam.read_uvbeam(args.beamfile)

if args.nf == None:
    print("Frequency parameters not specified. Using beam frequencies. \n")
    freq_array = (beam.uvbeam.freq_array[0] * u.Hz).to(u.MHz)

    if args.fmin != None and args.fmax != None:
        idx = np.where(
            (freq_array >= args.fmin * u.MHz) & (freq_array <= args.fmax * u.MHz)
        )[0]
        freq_array = freq_array[idx]
else:
    freq_array = np.linspace(args.fmin, args.fmax, args.nf) * u.MHz

# generate LST array
lst_array = np.linspace(args.tmin, args.tmax, args.nt)

# generate visibility data object
VisObj = VisData()
VisObj.read_ant_file(args.antfile)
VisObj.freq_array = freq_array
VisObj.remove_redundancies()
ds = VisObj.copy()

# iterate through LST data
for i, lst in enumerate(lst_array):
    print("LST={:.3f} h, {:d}/{:d}".format(lst, i, args.nt))
    gleam_region = gleam.select_region(15 * lst * u.deg, args.lat * u.deg, args.r * u.deg)
    gleam_coords = gleam_region.get_image_coordinates()
    gleam_vis = VisObj.copy()

    # polarisations
    for pol in range(2):
        gleam_sky = gleam_region.get_beam_sky(beam.uvbeam, freq_array, pol)
        gleam_vis_pol = VisObj.copy()
        gleam_vis_pol.vis_from_image(gleam_sky, gleam_coords)
        gleam_vis_pol.add_axis(axis=1)
        gleam_vis.append(gleam_vis_pol, axis=1)

    gleam_vis.add_axis(axis=2)
    ds.append(gleam_vis, axis=2)

# get closure phase and effective visibility
eicp = ds.get_eicp()
veff = ds.get_veff()

print(f"\nWriting visibility data to file {args.outfile} \n")
data_list = [lst_array, eicp, veff]
names = ["LST", "eicp", "Veff"]
ds.write_vis(args.outfile, data_list, names, overwrite=True)

print("All done \n")