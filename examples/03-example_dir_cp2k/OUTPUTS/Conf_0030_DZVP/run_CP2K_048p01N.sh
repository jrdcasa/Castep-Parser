#!/bin/bash
#SBATCH --partition=generic
#SBATCH --nodes=1
#SBATCH --tasks-per-node=48
#SBATCH --job-name=01-SiO2-1_01Al_P01b_noCl_screen_Conf_0030
#SBATCH --reservation=curso

WD=`pwd`

# ============== LOAD MODULES ================
module purge
module load GCC/11.2.0  OpenMPI/4.1.1 CP2K/2023.1

# ============== CP2K ================
mpirun -n 48 cp2k.popt input.dat >output_opt.dat
