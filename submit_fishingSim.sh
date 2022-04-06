#!/bin/bash
#SBATCH -J Fishing
#SBATCH -n 36
#SBATCH --time=00:02:00
module load gnu openmpi
mpirun /home/vvo2/HPC-Course/Viikko2/FishingSimulation/fishingSim

