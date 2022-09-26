#!/bin/bash
#SBATCH -p hera
#SBATCH -J gleam_vis
#SBATCH -o gleam_vis.out
#SBATCH -t 12:00:00
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user pmk46@cam.ac.uk

beamfile="/lustre/aoc/projects/hera/pkeller/HERA-Beams/NicolasFagnoniBeams/NF_HERA_Dipole_efield_beam.fits"
gleamfile="/users/pkeller/code/ClosureSim/data/gleam.fits"
pythonfile="/users/pkeller/code/ClosureSim/scripts/make_gleam_vis.py"

antfile="/users/pkeller/code/ClosureSim/data/antenna_positions_EQ28.dat"
outfile="/users/pkeller/code/ClosureSim/data/vis_EQ28_FornaxA.h5"
/lustre/aoc/projects/hera/pkeller/anaconda3/bin/python3 "${pythonfile}" -g "${gleamfile}" -o "${outfile}" -b "${beamfile}" -a "${antfile}" --nf 1024 --fmin 100 --fmax 200 --nt 1 --tmin 3.366 --tmax 3.367
