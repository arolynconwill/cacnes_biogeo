#!/bin/bash
#SBATCH -p sched_mem1TB,defq
#SBATCH -n 1
#SBATCH --time=30:00:00
#SBATCH -o log_2_out.txt
#SBATCH -e log_2_err.txt
#SBATCH --mem=128000


bash scripts/build_classifer_template.sh

echo Done!!!
