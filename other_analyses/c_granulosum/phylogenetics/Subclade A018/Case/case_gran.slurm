#!/bin/bash
#SBATCH -p sched_mem1TB,defq
#SBATCH -n 1
#SBATCH --time=2-00:00:00
#SBATCH -o masterout.txt
#SBATCH -e mastererr.txt
#SBATCH --mem=8000


bash snakemakeslurm.sh

echo Done!!!
