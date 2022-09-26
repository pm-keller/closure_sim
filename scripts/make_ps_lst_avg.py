import os
import sys
import h5py
import numpy as np
from scipy.signal.windows import get_window
import astropy.units as u
import astropy.coordinates as coord
from multiprocessing import Pool
import argparse

sys.path.append(os.path.join(os.getcwd(), "skysim/"))

from beam import BeamModel
from gleam import GLEAMData
from eor import EORModel
from vis import VisData
import dspec


parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gleamfile", help="Path to GLEAM file.", type=str)
parser.add_argument("-e", "--eorfile", help="Path to EoR file.", type=str)
parser.add_argument("-b", "--beamfile", help="Path to beam FITS-file.", type=str)
parser.add_argument("-o", "--outfile", help="File to write visibilities to.", type=str)
parser.add_argument(
    "-a", "--antfile", help="Path to file containing antenna positions.", type=str
)
parser.add_argument("--nt", help="Number of times.", type=int)
parser.add_argument("--tmin", help="Minimum LST in h.", type=float)
parser.add_argument("--tmax", help="Maximum LST in h.", type=float)
parser.add_argument("--Nmax", help="Maximum number of neighbouring LST's to integrate.", type=int)
parser.add_argument(
    "--lat",
    help="Latitude of telescope in deg.",
    default=-30.72138329631366,
    type=float,
)
args = parser.parse_args()


def sliding_average(data, n, axis):
    """
    Sliding average along axis 
    """
    data = np.moveaxis(data, axis, 0)
    data = np.array([np.nanmean(data[i:i+n], axis=0) for i in range(data.shape[0]-n)])
    data = np.exp(1j * np.angle(data)) * np.nanmean(np.abs(data))
    return np.moveaxis(data, 0, axis)

# read EoR
eor = EORModel()
eor.read(args.eorfile)
eor.convert_units("Jy")

# get coordinates
l_array, m_array, freq_array = eor.get_sky_coodinates()
coord_grid = np.meshgrid(l_array, m_array[args.Nmax//2:-args.Nmax//2])

# read beam
beam = BeamModel()
beam.read_image(args.beamfile)

# read GLEAM data
radius = eor.get_angular_radius()
freq_array = eor.get_freq_array()
gleam = GLEAMData(filters = {"Fint151":">0.1"})
gleam.read(args.gleamfile)

# bright sources that are not included in GLEAM
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
for source in sources.values():
    source_coord, F0, f0, alpha = source.values()
    ra = source_coord.ra.deg * u.deg
    dec = source_coord.dec.deg * u.deg
    gleam.add_source(ra, dec, f0, F0, alpha)

# prepare visibility data
VisObj = VisData()
VisObj.read_ant_file(args.antfile)
VisObj.remove_redundancies()
VisObj.freq_array = freq_array

# get window function (convolved Blackman-Harris window)
bh = get_window("blackmanharris", len(freq_array))
bh = bh / np.sqrt(np.mean(bh**2))
bh2 = np.convolve(bh, bh)[::2]
bh2 = bh2 / np.sqrt(np.mean(bh2**2))

ps_list = []
ps_units = u.mK ** 2 * u.Mpc ** 3

lst_array = np.linspace(args.tmin, args.tmax, args.nt)

def psavg(lst):
    with open('printfile.txt', 'a') as f:
        f.write("{:.4f}\n".format(lst))

    # make gleam images
    gleam_region = gleam.select_region(lst * 15 * u.deg, args.lat * u.deg, radius)    
    gleam_im = gleam_region.make_sky_image(freq_array, n=256)

    # make skies
    sky = eor.eor + gleam_im
    skies = np.array([sky[i:-(args.Nmax-i)] * beam.beam_image[args.Nmax//2:-args.Nmax//2] for i in range(args.Nmax)])

    # compute sky visibilities
    ds = VisObj.copy()

    for i, sky in enumerate(skies):
        sky_vis = VisObj.copy()
        sky_vis.vis_from_image(sky, coord_grid)
        sky_vis.add_axis(axis=1)
        ds.append(sky_vis, axis=1)
    
    eicp = ds.get_eicp(veff=True, veff_cut=5)

    # compute closure phase averages
    eicp_averages = [sliding_average(eicp, n, 0) for n in range(1, args.Nmax)]

    # compute power spectra
    ps_avg = []
    for eicp_avg in eicp_averages:
        ps = np.zeros_like(eicp_avg) * ps_units

        for i, eicp in enumerate(eicp_avg):
            DSpec = dspec.DSpec(eicp*u.Jy, freq_array=freq_array, window=bh2)
            ps[i] = DSpec.get_xps().apply_cosmo_units(212.550498068445 * u.m**2).xps

        ps_avg.append(np.nanmean(ps, axis=0))
    
    return ps_avg

with Pool(processes=16) as pool:
    ps_list = pool.map(psavg, lst_array)

ps_list = np.array(ps_list)
DSpec = dspec.DSpec(np.zeros_like(freq_array), freq_array=freq_array, window=bh2)

if os.path.exists(args.outfile):
    os.remove(args.outfile)

f = h5py.File(args.outfile, "a")
f.create_dataset("Data", data=ps_list)
f.create_dataset("LST Array", data=lst_array)
f.create_dataset("Delay", data=DSpec.get_delays().to(u.us).value)
f.close()
