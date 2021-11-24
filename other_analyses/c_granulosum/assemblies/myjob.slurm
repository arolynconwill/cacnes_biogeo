#!/bin/bash
#SBATCH -p sched_mem1TB,defq
#SBATCH -n 1
#SBATCH --time=1-00:00:00
#SBATCH -o masterout.txt
#SBATCH -e mastererr.txt
#SBATCH --mem=2000
#SBATCH --mail-user=PLACEHOLDER@mit.edu
#SBATCH --mail-type=ALL

bash snakemakeslurm.sh

echo Done!!!
