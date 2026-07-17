#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --threads-per-core=1
#SBATCH --output=dftd4-4000-%j.out
#SBATCH --error=dftd4-4000-%j.out
#SBATCH --mem=0
module load gcc/12.2.0
export OMP_NUM_THREADS=24
./dftd4/build/app/dftd4 _MOLECULES/4000.xyz --func pbe -s --mbdscale 0 --noedisp
