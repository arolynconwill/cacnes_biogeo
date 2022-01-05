#!/bin/bash
#SBATCH -p sched_mem1TB,defq
#SBATCH -n 1
#SBATCH --time=30:00:00
#SBATCH -o log_3_out.txt
#SBATCH -e log_3_err.txt
#SBATCH --mem=128000


bash scripts/run_classifer.sh

echo Done!!!
