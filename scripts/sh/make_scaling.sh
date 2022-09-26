#!/bin/bash
#SBATCH -p hera
#SBATCH -J scaling
#SBATCH -o scaling.out
#SBATCH -t 1:00:00
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user pmk46@cam.ac.uk

pythonfile="/users/pkeller/code/ClosureSim/scripts/make_scaling.py"

for trclass in EQ14 EQ28
do 
    for field in A B C D E
    do
        echo $trclass $field
        veffpath="/users/pkeller/code/H1C_IDR3.2/data/veff_cal_${trclass}_F${field}B2.dat"
        omegapath="/lustre/aoc/projects/hera/pkeller/HERA-Beams/NicolasFagnoniBeams/Omega2_Dipole.dat"
        outfile="/users/pkeller/code/H1C_IDR3.2/data/scaling_${trclass}_F${field}B2.dat"
        /lustre/aoc/projects/hera/pkeller/anaconda3/bin/python3 ${pythonfile} -v ${veffpath} -a ${omegapath} -o ${outfile} --nf 76 --fmin 160.59108434 --fmax 167.94
    done
done
