#!/bin/bash
#SBATCH -p sched_mem1TB,defq
#SBATCH -n 1
#SBATCH --time=30:00:00
#SBATCH -o masterout.txt
#SBATCH -e mastererr.txt
#SBATCH --mem=128000


bash build_classifer_template.sh

echo Done!!!
