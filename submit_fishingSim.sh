#!/bin/bash
#SBATCH -J Fishing
#SBATCH -n 36
module load gnu openmpi
mpirun /home/vvo2/HPC-Course/Viikko2/FishingSimulation/fishingSim

