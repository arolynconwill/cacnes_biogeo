#!/bin/bash 
#SBATCH -p sched_mem1TB,defq 
#SBATCH -n 1 
#SBATCH --time=1-00:00:00 
#SBATCH -o mainout.txt 
#SBATCH -e mainerr.txt 
#SBATCH --mem=2000 
#SBATCH --mail-user=aconwill@mit.edu 
#SBATCH --mail-type=ALL 
#SBATCH -D /scratch/mit_lieberman/projects/aro_cacnes_biogeo/6_lineage_assembly_case/clade_16 

bash snakemakeslurm.sh 

echo Done!!! 

