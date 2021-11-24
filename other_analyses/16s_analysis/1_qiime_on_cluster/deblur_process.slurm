#!/bin/bash
#SBATCH -p sched_mem1TB,defq
#SBATCH -n 8
#SBATCH --time=30:00:00
#SBATCH -o masterout.txt
#SBATCH -e mastererr.txt
#SBATCH --mem=256000


bash data_import_and_deblur.sh

echo Done!!!
