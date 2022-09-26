#!/bin/bash
#SBATCH -p hera
#SBATCH -J bright_vis
#SBATCH -o bright_vis.out
#SBATCH -t 24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user pmk46@cam.ac.uk

beamfile="/lustre/aoc/projects/hera/pkeller/HERA-Beams/NicolasFagnoniBeams/NF_HERA_Dipole_efield_beam.fits"
pythonfile="/users/pkeller/code/ClosureSim/scripts/make_bright_vis.py"

antfile="/users/pkeller/code/ClosureSim/data/antenna_positions_EQ28.dat"

outfile="/users/pkeller/code/ClosureSim/data/vis_bright_EQ28_J0136.h5"
/lustre/aoc/projects/hera/pkeller/anaconda3/bin/python3 "${pythonfile}" -o "${outfile}" -b "${beamfile}" -a "${antfile}" --nf 1024 --fmin 100 --fmax 200 --nt 1 --tmin 1.600 --tmax 1.601
