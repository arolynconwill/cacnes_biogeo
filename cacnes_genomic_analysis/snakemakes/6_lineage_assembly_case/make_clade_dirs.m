%% OBJECTIVE

% Make snakemake directories for each clade
% Additional versions at the bottom to specify the directory (for batching
% via loop)


%% Path

path(path,'matlab_functions')


%% Directory setup (local)

% Directory with csv files
dir_csvs = 'samples_case_csvs';

% Directory with base files:
dir_skeleton = 'snakemake_skeleton';

% Directory for putting all the clade directories:
dir_clades = 'snakemake_directories';


%% Directories (cluster)

dir_cluster_case_step = '/scratch/mit_lieberman/projects/aro_cacnes_biogeo/6_lineage_assembly_case/';


%% Generate samples_case.csv's

% Call matlab function
make_samples_case_csv


%% Make directory and fill it
% Alternative that specifies working directory

% NEW VERSION: allows for batching all clades at once since it specifies
% the working directory in the sbatch script


% Loop through clades
for c=1:53

    this_clade = c;
    dir_this_clade = [ 'clade_' num2str(this_clade) ];

    % Make directory
    if ~exist( [ dir_clades '/' dir_this_clade ], 'dir')
        mkdir( [ dir_clades '/' dir_this_clade ] )
    end

    % Snakefile
    % % Update csv file name
    copyfile( [ dir_skeleton '/' 'Snakefile' ], ...
        [ dir_clades '/' dir_this_clade '/' 'Snakefile' ] );

    % samples_case_clade_N.csv
    copyfile( [ dir_csvs '/' 'samples_case_clade_' num2str(this_clade) '.csv' ], ...
        [ dir_clades '/' dir_this_clade '/' 'samples_case.csv' ] );

    % myjob_N.slurm
    fileID = fopen([ dir_clades '/' dir_this_clade '/' 'myjob_' num2str(this_clade) '.slurm' ],'w');
    fprintf(fileID,'#!/bin/bash \n');
    fprintf(fileID,'#SBATCH -p sched_mem1TB,defq \n');
    fprintf(fileID,'#SBATCH -n 1 \n');
    fprintf(fileID,'#SBATCH --time=1-00:00:00 \n');
    fprintf(fileID,'#SBATCH -o mainout.txt \n');
    fprintf(fileID,'#SBATCH -e mainerr.txt \n');
    fprintf(fileID,'#SBATCH --mem=2000 \n');
    fprintf(fileID,'#SBATCH --mail-user=aconwill@mit.edu \n');
    fprintf(fileID,'#SBATCH --mail-type=ALL \n');
    fprintf(fileID,['#SBATCH -D ' dir_cluster_case_step 'clade_' num2str(this_clade) ' \n']);
    fprintf(fileID,'\n');
    fprintf(fileID,'bash snakemakeslurm.sh \n');
    fprintf(fileID,'\n');
    fprintf(fileID,'echo Done!!! \n');
    fprintf(fileID,'\n');
    fclose(fileID);

    % snakemakeslurm.sh
    copyfile( [ dir_skeleton '/' 'snakemakeslurm.sh' ], ...
        [ dir_clades '/' dir_this_clade '/' 'snakemakeslurm.sh' ] );

    % cluster.slurm.json
    if c<=4 % big clades need more compute time
        copyfile( [ dir_skeleton '/' 'cluster_bigclade.slurm.json' ], ...
            [ dir_clades '/' dir_this_clade '/' 'cluster.slurm.json' ] );
    else
        copyfile( [ dir_skeleton '/' 'cluster.slurm.json' ], ...
            [ dir_clades '/' dir_this_clade '/' 'cluster.slurm.json' ] );
    end

end



