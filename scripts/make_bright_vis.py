"""

make_bright_vis.py

Authors: Pascal M. Keller
Affiliation: Cavendish Astrophysics
Created on: 15/02/2022

Description:
    Simulate visibilitites of a sky populated with bright point sources which are excluded from the GLEAM catalog (e.g. Fornax A.)

    Example:
    python3 make_bright_vis.py -o "${outfile}" -b "${beamfile}" -a "${antfile}" --nf 167 --fmin 152.3 --fmax 167.94 --nt 752 --tmin 4.0 --tmax 6.25

"""

import argparse
import numpy as np
import os
import sys

import astropy.units as u
import astropy.coordinates as coord

sys.path.append(os.path.join(os.getcwd(), "skysim/"))
from gleam import GLEAMData
from beam import BeamModel
from vis import VisData


parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfile", help="File to write visibilities to.", type=str)
parser.add_argument("-b", "--beamfile", help="Path to beam FITS-file.", type=str)
parser.add_argument(
    "-a", "--antfile", help="Path to file containing antenna positions.", type=str
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

# sources to include in model (c.f. Table 2 in Hurley-Walker et al. 2017)
sources = {
    "Fornax A Core": {
        "Coord": coord.SkyCoord("03h22m43s", "-37d12m02s"),
        "Flux": 12 * u.Jy,
        "Frequency": 154 * u.MHz,
        "alpha": -0.88,
    },
    "Fornax A West": {
        "Coord": coord.SkyCoord("03h21m17s", "-37d9m10s"),
        "Flux": 478 * u.Jy,
        "Frequency": 154 * u.MHz,
        "alpha": -0.77,
    },
    "Fornax A East": {
        "Coord": coord.SkyCoord("03h24m12s", "-37d14m30s"),
        "Flux": 260 * u.Jy,
        "Frequency": 154 * u.MHz,
        "alpha": -0.77,
    },
    "Pictor A": {
        "Coord": coord.SkyCoord("05h19m50s", "-45d46m44s"),
        "Flux": 390 * u.Jy,
        "Frequency": 200 * u.MHz,
        "alpha": -0.99,
    },
    "Hydra A": {
        "Coord": coord.SkyCoord("09h18m06s", "-12d05m44s"),
        "Flux": 280 * u.Jy,
        "Frequency": 200 * u.MHz,
        "alpha": -0.96,
    },
    "Hercules A": {
        "Coord": coord.SkyCoord("16h51m08s", "04d59m33s"),
        "Flux": 377 * u.Jy,
        "Frequency": 200 * u.MHz,
        "alpha": -1.07,
    },
    "Virgo A": {
        "Coord": coord.SkyCoord("12h30m49s", "12d23m28s"),
        "Flux": 861 * u.Jy,
        "Frequency": 200 * u.MHz,
        "alpha": -0.86,
    },
    "Crab": {
        "Coord": coord.SkyCoord("05h34m32s", "22d00m52s"),
        "Flux": 1340 * u.Jy,
        "Frequency": 200 * u.MHz,
        "alpha": -0.22,
    },
    "Cygnus A": {
        "Coord": coord.SkyCoord("19h59m28s", "40d44m02s"),
        "Flux": 7920 * u.Jy,
        "Frequency": 200 * u.MHz,
        "alpha": -0.78,
    },
    "Cassiopeia A": {
        "Coord": coord.SkyCoord("23h23m28s", "58d48m42s"),
        "Flux": 11900 * u.Jy,
        "Frequency": 200 * u.MHz,
        "alpha": -0.41,
    },
}


# insert sources in gleam catalog
srcs = GLEAMData()
srcs.dec = args.lat * u.deg

for source in sources:
    source_coord, F0, f0, alpha = sources[source].values()
    ra = source_coord.ra.deg * u.deg
    dec = source_coord.dec.deg * u.deg
    srcs.add_source(ra, dec, f0, F0, alpha)

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
    srcs.ra = lst * 15 * u.deg
    alt = srcs.get_alt_az()[0]
    idx = np.where(alt > 0)[0]
    coords = srcs.get_image_coordinates()
    coords = (coords[0][idx], coords[1][idx])
    vis = VisObj.copy()

    # polarisations
    for pol in range(2):
        sky = srcs.get_beam_sky(beam.uvbeam, freq_array, pol)[idx]
        vis_pol = VisObj.copy()
        vis_pol.vis_from_image(sky, coords)
        vis_pol.add_axis(axis=1)
        vis.append(vis_pol, axis=1)

    vis.add_axis(axis=2)
    ds.append(vis, axis=2)

# get closure phase and effective visibility
eicp = ds.get_eicp()
veff = ds.get_veff()

print(f"\nWriting visibility data to file {args.outfile} \n")
data_list = [lst_array, eicp, veff]
names = ["LST", "eicp", "Veff"]
ds.write_vis(args.outfile, data_list, names, overwrite=True)

print("All done \n")