#!/bin/bash
#SBATCH -p sched_mem1TB
#SBATCH -n 1
#SBATCH --time=70:00:00
#SBATCH -o mainout.txt
#SBATCH -e mainerr.txt
#SBATCH --mem=500G
#SBATCH -x node325
#SBATCH --mail-user=aconwill@mit.edu
#SBATCH --mail-type=ALL

module add mit/matlab/2015b

matlab -r make_tree_clade_cluster
echo Done!!!
