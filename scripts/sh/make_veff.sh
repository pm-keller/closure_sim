#!/bin/bash
#SBATCH -p hera
#SBATCH -J veff
#SBATCH -o veff.out
#SBATCH -t 1:00:00
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user pmk46@cam.ac.uk

pythonfile="/users/pkeller/code/ClosureSim/scripts/make_veff.py"

for trclass in EQ14 EQ28
do 
    for field in A B C
    do
        echo $trclass $field
        gleamfile="/users/pkeller/code/ClosureSim/data/vis_${trclass}_F${field}B2.h5"
        brightfile="/users/pkeller/code/ClosureSim/data/vis_bright_${trclass}_F${field}B2.h5"
        outfile="/users/pkeller/code/H1C_IDR3.2/data/veff_${trclass}_F${field}B2_model.dat"
        /lustre/aoc/projects/hera/pkeller/anaconda3/bin/python3 ${pythonfile} -g ${gleamfile} -b ${brightfile} -o ${outfile} -n 16 -l 85
    done
done
