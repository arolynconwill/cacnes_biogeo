Overview:

A separate case step is necessary for each lineage.


Instructions:

1. Run make_clade_dirs.m 
2. Run case step snakemake in each directory "6_lineage_assembly_case/snakemake_directories/clade_*"


How to run all at once:

for i in `seq 1 53`; do sbatch /scratch/mit_lieberman/projects/aro_cacnes_biogeo/6_lineage_assembly_case/clade_$i/myjob_$i.slurm; done 