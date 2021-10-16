#!/bin/bash
#SBATCH --job-name=LinAM.main
#SBATCH -p sched_mem1TB,defq
#SBATCH -n 1
#SBATCH --time=2-00:00:00
#SBATCH -o mainout.txt
#SBATCH -e mainerr.txt
#SBATCH --mem=2000
#SBATCH --mail-user=aconwill@mit.edu
#SBATCH --mail-type=ALL

bash snakemakeslurm.sh

echo Done!!!
