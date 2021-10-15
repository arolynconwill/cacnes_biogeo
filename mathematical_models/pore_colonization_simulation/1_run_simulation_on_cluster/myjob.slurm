#!/bin/bash
#SBATCH --job-name PoreSimSimple.SM.main
#SBATCH -p sched_mem1TB,defq
#SBATCH -n 1
#SBATCH --time=09:00:00
#SBATCH -o main_out.txt
#SBATCH -e main_err.txt
#SBATCH --mem=2000
#SBATCH --mail-user=aconwill@mit.edu
#SBATCH --mail-type=ALL


module add mit/matlab/2018a; 
matlab -singleCompThread -r "simple_pore_diffusion"


echo Done!!!
