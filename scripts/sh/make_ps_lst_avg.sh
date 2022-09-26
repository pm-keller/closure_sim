#!/bin/bash
#SBATCH -p hera
#SBATCH -J ps_avg
#SBATCH -o ps_avg.out
#SBATCH -t 1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user pmk46@cam.ac.uk

eorfile="/users/pkeller/code/ClosureSim/data/eor_FG_B2_256.h5"
beamfile="/users/pkeller/code/ClosureSim/data/beam_B2_256.h5"
gleamfile="/users/pkeller/code/ClosureSim/data/gleam.fits"
pythonfile="/users/pkeller/code/ClosureSim/scripts/make_ps_lst_avg.py"
antfile="/users/pkeller/code/ClosureSim/data/antenna_positions_EQ14.dat"
outfile="/users/pkeller/code/ClosureSim/data/ps_lst_avg_EQ14_FCB2.h5"

/lustre/aoc/projects/hera/pkeller/anaconda3/bin/python3 "${pythonfile}" -g "${gleamfile}" -e "${eorfile}" -b "${beamfile}" -a "${antfile}" -o "${outfile}" --nt 1000 --tmin 4.0 --tmax 6.25 --Nmax 10
