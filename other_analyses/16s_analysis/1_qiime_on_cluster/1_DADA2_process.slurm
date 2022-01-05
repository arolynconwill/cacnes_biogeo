#!/bin/bash
#SBATCH -p sched_mem1TB,defq
#SBATCH -n 8
#SBATCH --time=30:00:00
#SBATCH -o log_1_out.txt
#SBATCH -e log_1_err.txt
#SBATCH --mem=256000


bash scripts/data_import_and_dada2.sh

echo Done!!!
