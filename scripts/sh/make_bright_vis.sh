#!/bin/bash
#SBATCH -p hera
#SBATCH -J bright_vis
#SBATCH -o bright_vis.out
#SBATCH -t 12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user pmk46@cam.ac.uk

beamfile="/lustre/aoc/projects/hera/pkeller/HERA-Beams/NicolasFagnoniBeams/NF_HERA_Dipole_efield_beam.fits"
pythonfile="/users/pkeller/code/ClosureSim/scripts/make_bright_vis.py"

antfile="/users/pkeller/code/ClosureSim/data/antenna_positions_EQ14.dat"

outfile="/users/pkeller/code/ClosureSim/data/vis_bright_EQ14_FAB2.h5"
#/lustre/aoc/projects/hera/pkeller/anaconda3/bin/python3 "${pythonfile}" -o "${outfile}" -b "${beamfile}" -a "${antfile}" --nf 161 --fmin 152.25 --fmax 167.94 --nt 835 --tmin 21.5 --tmax 24.0

outfile="/users/pkeller/code/ClosureSim/data/vis_bright_EQ14_FBB2.h5"
/lustre/aoc/projects/hera/pkeller/anaconda3/bin/python3 "${pythonfile}" -o "${outfile}" -b "${beamfile}" -a "${antfile}" --nf 161 --fmin 152.25 --fmax 167.94 --nt 669 --tmin 0.75 --tmax 2.75

outfile="/users/pkeller/code/ClosureSim/data/vis_bright_EQ14_FCB2.h5"
#/lustre/aoc/projects/hera/pkeller/anaconda3/bin/python3 "${pythonfile}" -o "${outfile}" -b "${beamfile}" -a "${antfile}" --nf 161 --fmin 152.25 --fmax 167.94 --nt 752 --tmin 4.0 --tmax 6.25

#outfile="/users/pkeller/code/ClosureSim/data/vis_bright_EQ14_FDB2.h5"
#/lustre/aoc/projects/hera/pkeller/anaconda3/bin/python3 "${pythonfile}" -o "${outfile}" -b "${beamfile}" -a "${antfile}" --nf 161 --fmin 152.25 --fmax 167.94 --nt 1003 --tmin 6.25 --tmax 9.25

#outfile="/users/pkeller/code/ClosureSim/data/vis_bright_EQ14_FEB2.h5"
#/lustre/aoc/projects/hera/pkeller/anaconda3/bin/python3 "${pythonfile}" -o "${outfile}" -b "${beamfile}" -a "${antfile}" --nf 161 --fmin 152.25 --fmax 167.94 --nt 1838 --tmin 9.25 --tmax 14.75

antfile="/users/pkeller/code/ClosureSim/data/antenna_positions_EQ28.dat"

outfile="/users/pkeller/code/ClosureSim/data/vis_bright_EQ28_FAB2.h5"
#/lustre/aoc/projects/hera/pkeller/anaconda3/bin/python3 "${pythonfile}" -o "${outfile}" -b "${beamfile}" -a "${antfile}" --nf 161 --fmin 152.25 --fmax 167.94 --nt 835 --tmin 21.5 --tmax 24.0

outfile="/users/pkeller/code/ClosureSim/data/vis_bright_EQ28_FBB2.h5"
/lustre/aoc/projects/hera/pkeller/anaconda3/bin/python3 "${pythonfile}" -o "${outfile}" -b "${beamfile}" -a "${antfile}" --nf 161 --fmin 152.25 --fmax 167.94 --nt 669 --tmin 0.75 --tmax 2.75

outfile="/users/pkeller/code/ClosureSim/data/vis_bright_EQ28_FCB2.h5"
#/lustre/aoc/projects/hera/pkeller/anaconda3/bin/python3 "${pythonfile}" -o "${outfile}" -b "${beamfile}" -a "${antfile}" --nf 161 --fmin 152.25 --fmax 167.94 --nt 752 --tmin 4.0 --tmax 6.25

#outfile="/users/pkeller/code/ClosureSim/data/vis_bright_EQ28_FDB2.h5"
#/lustre/aoc/projects/hera/pkeller/anaconda3/bin/python3 "${pythonfile}" -o "${outfile}" -b "${beamfile}" -a "${antfile}" --nf 161 --fmin 152.25 --fmax 167.94 --nt 1003 --tmin 6.25 --tmax 9.25

#outfile="/users/pkeller/code/ClosureSim/data/vis_bright_EQ28_FEB2.h5"
#/lustre/aoc/projects/hera/pkeller/anaconda3/bin/python3 "${pythonfile}" -o "${outfile}" -b "${beamfile}" -a "${antfile}" --nf 161 --fmin 152.25 --fmax 167.94 --nt 1838 --tmin 9.25 --tmax 14.75