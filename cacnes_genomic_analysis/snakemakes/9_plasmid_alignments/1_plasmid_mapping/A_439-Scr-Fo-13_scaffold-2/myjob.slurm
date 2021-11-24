#!/bin/bash
#SBATCH --job-name MapScaff1.SM.main
#SBATCH -p sched_mem1TB,defq
#SBATCH -n 1
#SBATCH --time=1-00:00:00
#SBATCH -o mainout.txt
#SBATCH -e mainerr.txt
#SBATCH --mem=2000
#SBATCH --mail-user=aconwill@mit.edu
#SBATCH --mail-type=ALL
#SBATCH --exclude=node325

bash snakemakeslurm.sh

echo Done!!!
